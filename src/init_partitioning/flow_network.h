/*******************************************************************************
 * MIT License
 *
 * This file is part of HeiProMap.
 *
 * Copyright (C) 2025 Henning Woydt <henning.woydt@informatik.uni-heidelberg.de>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#ifndef HEIPROMAP_FLOW_NETWORK_H
#define HEIPROMAP_FLOW_NETWORK_H

#include <climits>
#include <vector>
#include <cstring>

#include "../definitions.h"
#include "../utility/aligned_array.h"
#include "../datastructures/csr_graph.h"

namespace HeiProMap {
    /**
     * Flow network of a sub graph with undirected edges (two paired arcs with the edge weight as capacity) on the
     * nodes 0..nl-1, and an incremental Boykov-Kolmogorov max flow between two growing terminal sets.
     *
     * s and t are implicit: a source terminal (root of the source tree) has an infinite arc s -> v, a sink
     * terminal one v -> t. A simple augmenting path never returns to s or leaves t, so no capacity is stored.
     *
     * Search trees: label 0 free, 1 in the source tree, 2 in the sink tree. pa[v] is the arc (in v's list) to the
     * parent, -1 for a root (terminal), -3 for an orphan. At termination of bk_maxflow the source tree is exactly the
     * set reachable from s and the sink tree the set that reaches t in the residual graph.
     *
     * All buffers are AlignedArray members that only grow (AlignedArray::initialize reallocates only if the request
     * exceeds the capacity). Their size() is the capacity, the logical sizes are nl, N and na.
     */
    class FlowNetwork {
    protected:
        static constexpr s64 INF = (s64) 1 << 60;

        struct Arc {
            int to;
            int rev;
            s64 cap;
        };

        int nl = 0; // number of vertices
        int N = 0; // = nl, size of the per vertex arrays
        size_t na = 0; // number of arcs
        AlignedArray<Arc> arcs; // residual network
        AlignedArray<Arc> base_arcs; // network without flow
        AlignedArray<int> first; // arcs of v: [first[v], first[v + 1])
        AlignedArray<weight_t> vw;
        weight_t total_w = 0;

        // Boykov-Kolmogorov state
        AlignedArray<char> label, isact;
        struct Verified {
            int time; // bk_time at which the distance was verified
            int dist; // distance to the root of the tree
        };
        AlignedArray<int> pa, act, orph;
        AlignedArray<Verified> verified;
        int bk_time = 0;
        // IBFS (first solve from scratch)
        AlignedArray<int> ib_d, ib_cur;
        std::vector<int> ib_fs, ib_ft, ib_nx, ib_q;
        std::vector<std::vector<int>> ib_bk;
        std::vector<int> *ib_next = nullptr, *ib_otherf = nullptr;
        int act_h = 0, act_t = 0, orph_n = 0;

        // scratch of build_network
        AlignedArray<int> deg_buf, pos_buf;

        /// builds the network of the vertices verts[0, cnt); loc[v] is the local index of a vertex of the graph or -1
        void build_network(const CSRGraph &g, const vertex_t *verts, const int cnt, const int *loc) {
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "build_network");
            nl = cnt;
            N = nl;

            vw.initialize(nl);
            deg_buf.initialize(N);
            pos_buf.initialize(N);
            first.initialize(N + 1);
            total_w = 0;
            for (int i = 0; i < N; ++i) { deg_buf[i] = 0; }
            for (int u = 0; u < nl; ++u) {
                const vertex_t gu = verts[u];
                vw[u] = g.v_weight(gu);
                total_w += vw[u];
                for (size_t i = g.neighborhoods[gu]; i < g.neighborhoods[gu + 1]; ++i) {
                    if (loc[g.edges_v[i]] >= 0) { deg_buf[u]++; }
                }
            }

            first[0] = 0;
            for (int i = 0; i < N; ++i) { first[i + 1] = first[i] + deg_buf[i]; }
            na = (size_t) first[N];
            arcs.initialize(na);
            base_arcs.initialize(na);
            for (int i = 0; i < N; ++i) { pos_buf[i] = first[i]; }

            for (int u = 0; u < nl; ++u) {
                const vertex_t gu = verts[u];
                for (size_t i = g.neighborhoods[gu]; i < g.neighborhoods[gu + 1]; ++i) {
                    const int v = loc[g.edges_v[i]];
                    if (v < 0 || u >= v) { continue; } // each undirected edge once
                    const s64 w = g.edge_w(i);
                    const int a = pos_buf[u]++;
                    const int b = pos_buf[v]++;
                    arcs[a] = {v, b, w};
                    arcs[b] = {u, a, w};
                }
            }
            std::memcpy(base_arcs.get_ptr(), arcs.get_ptr(), na * sizeof(Arc));

            label.initialize(N, 0);
            isact.initialize(N, 0);
            pa.initialize(N, -1);
            act.initialize(N + 1);
            orph.initialize(N);
            verified.initialize(N, Verified{0, 0});
            ib_d.initialize(N);
            ib_cur.initialize(N);
        }

        // ---------------------------------------------------------------- Boykov-Kolmogorov
        // (all functions below keep the search trees between calls; terminals may be added at any time)
        void bk_activate(const int v) {
            if (isact[v]) { return; }
            isact[v] = 1;
            act[act_t] = v;
            act_t = (act_t + 1 == N + 1) ? 0 : act_t + 1;
        }

        void bk_orphan(const int v) {
            if (pa[v] == -3) { return; }
            pa[v] = -3;
            orph[orph_n++] = v;
        }

        // A terminal is a root of a tree. Free vertices have label 0 (and pa -1).
        bool is_terminal(const int v) const { return label[v] != 0 && pa[v] == -1; }

        // empty trees (no terminals), nothing active
        void reset_trees() {
            std::fill_n(label.get_ptr(), nl, 0);
            std::fill_n(pa.get_ptr(), nl, -1);
            std::fill_n(isact.get_ptr(), nl, 0);
            act_h = act_t = 0;
            orph_n = 0;
        }

        // v becomes a terminal of tree `lab` (1 source, 2 sink) while the trees are still empty
        void set_terminal(const int v, const char lab) {
            label[v] = lab;
            pa[v] = -1;
        }

        // trees from scratch: the terminals set since reset_trees are the roots
        void bk_start() {
            for (int v = 0; v < nl; ++v) { if (label[v]) { bk_activate(v); } }
        }

        // v becomes a source terminal while the trees are kept
        void bk_add_source(const int v) {
            if (label[v] == 2) {
                for (int i = first[v]; i < first[v + 1]; ++i) {
                    const int q = arcs[i].to;
                    if (label[q] == 2 && pa[q] == arcs[i].rev) { bk_orphan(q); }
                }
            }
            label[v] = 1;
            pa[v] = -1;
            bk_activate(v);
        }

        void bk_add_sink(const int v) {
            if (label[v] == 1) {
                for (int i = first[v]; i < first[v + 1]; ++i) {
                    const int q = arcs[i].to;
                    if (label[q] == 1 && pa[q] == arcs[i].rev) { bk_orphan(q); }
                }
            }
            label[v] = 2;
            pa[v] = -1;
            bk_activate(v);
        }

        // tries to find a new parent for the orphan p, otherwise frees it
        void bk_adopt(const int p) {
            if (pa[p] != -3) { return; } // became a terminal (root) while waiting
            const char lab = label[p];
            int best = -1, best_d = INT_MAX;
            for (int i = first[p]; i < first[p + 1]; ++i) {
                const int q = arcs[i].to;
                if (label[q] != lab) { continue; }
                const bool res = lab == 1 ? arcs[arcs[i].rev].cap > 0 : arcs[i].cap > 0;
                if (!res) { continue; }
                // distance of q to its root, valid only if the chain reaches a root without orphans
                int d = 0, x = q;
                while (true) {
                    if (verified[x].time == bk_time) { d += verified[x].dist; break; }
                    const int a = pa[x];
                    if (a == -1) { d += 1; break; }
                    if (a == -3) { d = INT_MAX; break; }
                    x = arcs[a].to;
                    ++d;
                }
                if (d == INT_MAX) { continue; }
                // mark the chain as verified for this time stamp
                int y = q, dd = d;
                while (verified[y].time != bk_time) {
                    verified[y] = Verified{bk_time, dd--};
                    const int a = pa[y];
                    if (a == -1) { break; }
                    y = arcs[a].to;
                }
                if (d < best_d) { best_d = d; best = i; }
            }
            if (best >= 0) {
                pa[p] = best;
                verified[p] = Verified{bk_time, best_d + 1};
                return;
            }
            label[p] = 0;
            for (int i = first[p]; i < first[p + 1]; ++i) {
                const int q = arcs[i].to;
                if (label[q] != lab) { continue; }
                const bool can_grow = lab == 1 ? arcs[arcs[i].rev].cap > 0 : arcs[i].cap > 0;
                if (can_grow) { bk_activate(q); }
                if (pa[q] == arcs[i].rev) { bk_orphan(q); }
            }
        }

        void bk_process_orphans() {
            ++bk_time;
            while (orph_n > 0) { bk_adopt(orph[--orph_n]); }
        }

        // pushes flow along s-tree path to a, arc ba (a -> b), t-tree path from b
        s64 bk_augment(const int ba, const int a_node) {
            const int b = arcs[ba].to;
            s64 bn = arcs[ba].cap;
            for (int v = a_node; pa[v] != -1; v = arcs[pa[v]].to) { bn = std::min(bn, arcs[arcs[pa[v]].rev].cap); }
            for (int v = b; pa[v] != -1; v = arcs[pa[v]].to) { bn = std::min(bn, arcs[pa[v]].cap); }
            arcs[ba].cap -= bn;
            arcs[arcs[ba].rev].cap += bn;
            for (int v = a_node; pa[v] != -1;) {
                const int up = pa[v];
                const int par = arcs[up].to;
                Arc &pr = arcs[arcs[up].rev];
                pr.cap -= bn;
                arcs[up].cap += bn;
                if (pr.cap == 0) { bk_orphan(v); }
                v = par;
            }
            for (int v = b; pa[v] != -1;) {
                const int up = pa[v];
                const int par = arcs[up].to;
                arcs[up].cap -= bn;
                arcs[arcs[up].rev].cap += bn;
                if (arcs[up].cap == 0) { bk_orphan(v); }
                v = par;
            }
            return bn;
        }


        // ---------------------------------------------------------------- IBFS
        // Incremental breadth-first search max flow (Goldberg, Hed, Kaplan, Tarjan, Werneck 2011) from scratch for the
        // terminals set since reset_trees (no bk_start). Uses the same tree representation as BK (label, pa), and at
        // the end completes both trees, so bk_maxflow can continue from the result (incremental piercing).
        // Returns the flow; at `limit` it returns early (the state is then not usable).

        // v left its tree or was relabeled: its children become orphans
        void ibfs_orphan_children(const int v, const char lab) {
            for (int j = first[v]; j < first[v + 1]; ++j) {
                const int w = arcs[j].to;
                if (label[w] == lab && pa[w] == arcs[j].rev) { bk_orphan(w); }
            }
        }

        // adoption of all orphans; tr = tree grown in this pass (levels up to D + 1), the other one has levels up to Do
        void ibfs_adopt(const char tr, const int D, const int Do) {
            auto proc = [&](const int v) {
                if (pa[v] != -3) { return; }
                const char lab = label[v];
                const int dv = ib_d[v];
                const int maxlev = lab == tr ? D + 1 : Do;
                // scan from the current arc on: stop at an admissible arc (new parent, same label); otherwise remember
                // the first arc to the lowest neighbour. The arcs before the current arc come first in the list order.
                bool found = false;
                int best = -1, bd = INT_MAX;
                const int c0 = ib_cur[v];
                for (int i = c0; i < first[v + 1] && !found; ++i) {
                    const Arc &a = arcs[i];
                    const int u = a.to;
                    if (label[u] != lab || ib_d[u] >= bd) { continue; }
                    if (lab == 1 ? arcs[a.rev].cap > 0 : a.cap > 0) {
                        best = i;
                        bd = ib_d[u];
                        if (bd == dv - 1) { found = true; }
                    }
                }
                if (found) {
                    pa[v] = best;
                    ib_cur[v] = best;
                    return;
                }
                {
                    int bestB = -1, bdB = INT_MAX;
                    for (int i = first[v]; i < c0; ++i) {
                        const Arc &a = arcs[i];
                        const int u = a.to;
                        if (label[u] != lab || ib_d[u] >= bdB) { continue; }
                        if (lab == 1 ? arcs[a.rev].cap > 0 : a.cap > 0) { bestB = i; bdB = ib_d[u]; }
                    }
                    if (bdB <= bd) { best = bestB; bd = bdB; }
                }
                if (best < 0 || bd + 1 > maxlev) {
                    label[v] = 0;
                    pa[v] = -1;
                    ibfs_orphan_children(v, lab);
                    return;
                }
                pa[v] = best;
                ib_cur[v] = best;
                if (bd + 1 != dv) {
                    ib_d[v] = bd + 1;
                    ibfs_orphan_children(v, lab);
                    if (lab == tr) { if (bd + 1 == D + 1) { ib_next->push_back(v); } }
                    else if (bd + 1 == Do) { ib_otherf->push_back(v); }
                }
            };
            size_t lvl = ib_bk.size();
            auto drain = [&]() {
                while (orph_n > 0) {
                    const int w = orph[--orph_n];
                    const size_t l = ib_d[w];
                    if (l >= ib_bk.size()) { ib_bk.resize(l + 1); }
                    ib_bk[l].push_back(w);
                    if (l < lvl) { lvl = l; }
                }
            };
            drain();
            while (lvl < ib_bk.size()) {
                if (ib_bk[lvl].empty()) { ++lvl; continue; }
                const int v = ib_bk[lvl].back();
                ib_bk[lvl].pop_back();
                if (pa[v] == -3) { proc(v); }
                drain();
            }
        }

        // one pass: grows tree tr by one level. 0 = limit reached, 1 = continue, 2 = no vertex at the next level
        int ibfs_pass(const char tr, int &D, const int Do, std::vector<int> &f, std::vector<int> &other_f,
                      s64 &flow, const s64 limit) {
            const char other = 3 - tr;
            // frontier: the vertices of the tree at level D (entries may be stale or duplicated)
            size_t m = 0;
            for (size_t x = 0; x < f.size(); ++x) {
                const int v = f[x];
                if (label[v] == tr && ib_d[v] == D && pa[v] != -3 && !isact[v]) { isact[v] = 1; f[m++] = v; }
            }
            f.resize(m);
            for (size_t x = 0; x < m; ++x) { isact[f[x]] = 0; }
            ib_nx.clear();
            ib_next = &ib_nx;
            ib_otherf = &other_f;
            for (size_t x = 0; x < f.size(); ++x) {
                const int v = f[x];
                if (label[v] != tr || ib_d[v] != D) { continue; }
                int i = first[v];
                while (i < first[v + 1]) {
                    if (label[v] != tr || ib_d[v] != D) { break; }
                    const Arc &a = arcs[i];
                    const bool res = tr == 1 ? a.cap > 0 : arcs[a.rev].cap > 0;
                    if (!res) { ++i; continue; }
                    const int w = a.to;
                    if (label[w] == 0) {
                        label[w] = tr;
                        pa[w] = a.rev;
                        ib_d[w] = D + 1;
                        ib_cur[w] = first[w];
                        ib_nx.push_back(w);
                        ++i;
                    } else if (label[w] == other) {
                        flow += bk_augment(tr == 1 ? i : a.rev, tr == 1 ? v : w);
                        if (flow >= limit) { return 0; }
                        ibfs_adopt(tr, D, Do);
                    } else {
                        ++i;
                    }
                }
            }
            // next frontier
            m = 0;
            for (size_t x = 0; x < ib_nx.size(); ++x) {
                const int v = ib_nx[x];
                if (label[v] == tr && ib_d[v] == D + 1 && pa[v] != -3 && !isact[v]) { isact[v] = 1; ib_nx[m++] = v; }
            }
            ib_nx.resize(m);
            for (size_t x = 0; x < m; ++x) { isact[ib_nx[x]] = 0; }
            if (m == 0) { return 2; }
            f.swap(ib_nx);
            ++D;
            return 1;
        }

        // grows the tree of `lab` over all residual arcs (BFS) so that it is exactly the reachable (S) / reaching (T) set
        void ibfs_complete(const char lab) {
            ib_q.clear();
            for (int v = 0; v < nl; ++v) { if (label[v] == lab) { ib_q.push_back(v); } }
            for (size_t h = 0; h < ib_q.size(); ++h) {
                const int x = ib_q[h];
                for (int i = first[x]; i < first[x + 1]; ++i) {
                    const int w = arcs[i].to;
                    if (label[w] != 0) { continue; }
                    if (lab == 1 ? arcs[i].cap > 0 : arcs[arcs[i].rev].cap > 0) {
                        label[w] = lab;
                        pa[w] = arcs[i].rev;
                        ib_q.push_back(w);
                    }
                }
            }
        }

        s64 ibfs_maxflow(const s64 limit = INF) {
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "ibfs");
            s64 flow = 0;
            ib_fs.clear();
            ib_ft.clear();
            for (int v = 0; v < nl; ++v) {
                if (label[v] == 1) { ib_d[v] = 0; ib_cur[v] = first[v]; ib_fs.push_back(v); }
                else if (label[v] == 2) { ib_d[v] = 0; ib_cur[v] = first[v]; ib_ft.push_back(v); }
            }
            orph_n = 0;
            act_h = act_t = 0;
            int Ds = 0, Dt = 0;
            while (!ib_fs.empty() && !ib_ft.empty()) {
                int r;
                if (ib_fs.size() <= ib_ft.size()) { r = ibfs_pass(1, Ds, Dt, ib_fs, ib_ft, flow, limit); }
                else { r = ibfs_pass(2, Dt, Ds, ib_ft, ib_fs, flow, limit); }
                if (r == 0) { return flow; }
                if (r == 2) { break; }
            }
            ibfs_complete(1);
            ibfs_complete(2);
            return flow;
        }

        // grows the trees until no active vertex is left; returns the added flow (early return at `limit`)
        s64 bk_maxflow(const s64 limit = INF) {
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "maxflow");
            s64 flow = 0;
            bk_process_orphans();
            while (act_h != act_t) {
                const int p = act[act_h];
                if (label[p] == 0) {
                    isact[p] = 0;
                    act_h = (act_h + 1 == N + 1) ? 0 : act_h + 1;
                    continue;
                }
                int ba = -1, an = -1;
                if (label[p] == 1) {
                    for (int i = first[p]; i < first[p + 1]; ++i) {
                        if (arcs[i].cap <= 0) { continue; }
                        const int q = arcs[i].to;
                        if (label[q] == 0) {
                            label[q] = 1;
                            pa[q] = arcs[i].rev;
                            bk_activate(q);
                        } else if (label[q] == 2) {
                            ba = i;
                            an = p;
                            break;
                        }
                    }
                } else {
                    for (int i = first[p]; i < first[p + 1]; ++i) {
                        if (arcs[arcs[i].rev].cap <= 0) { continue; }
                        const int q = arcs[i].to;
                        if (label[q] == 0) {
                            label[q] = 2;
                            pa[q] = arcs[i].rev;
                            bk_activate(q);
                        } else if (label[q] == 1) {
                            ba = arcs[i].rev;
                            an = q;
                            break;
                        }
                    }
                }
                if (ba >= 0) {
                    flow += bk_augment(ba, an);
                    if (flow >= limit) { return flow; }
                    bk_process_orphans();
                } else {
                    isact[p] = 0;
                    act_h = (act_h + 1 == N + 1) ? 0 : act_h + 1;
                }
            }
            return flow;
        }
    };
} // namespace HeiProMap

#endif // HEIPROMAP_FLOW_NETWORK_H
