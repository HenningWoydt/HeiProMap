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
#include <cstring>

#include "../definitions.h"
#include "../utility/aligned_array.h"
#include "../datastructures/csr_graph.h"

namespace HeiProMap {
    /**
     * Flow network of a sub graph with undirected edges (two paired arcs with the edge weight as capacity) on the
     * nodes 0..nl-1, and an incremental Boykov-Kolmogorov max flow between two growing terminal sets.
     *
     * s and t are implicit: a vertex with in_s[v] has an infinite arc s -> v, a vertex with in_t[v] one v -> t. A
     * simple augmenting path never returns to s or leaves t, so no capacity is stored for these arcs.
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

        AlignedArray<char> in_s, in_t;
        AlignedArray<int> src_list; // vertices with in_s
        int src_n = 0;

        // Boykov-Kolmogorov state
        AlignedArray<char> label, isact;
        AlignedArray<int> pa, act, orph, ts, dst;
        int bk_time = 0;
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

            src_list.initialize(N);
            label.initialize(N, 0);
            isact.initialize(N, 0);
            pa.initialize(N, -1);
            act.initialize(N + 1);
            orph.initialize(N);
            ts.initialize(N, 0);
            dst.initialize(N, 0);
            in_s.initialize(N);
            in_t.initialize(N);
        }

        void add_source(const int v) {
            in_s[v] = 1;
            src_list[src_n++] = v;
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

        // trees from scratch: terminals are the roots
        void bk_init() {
            for (int v = 0; v < nl; ++v) {
                label[v] = in_s[v] ? 1 : (in_t[v] ? 2 : 0);
                pa[v] = -1;
                isact[v] = 0;
            }
            act_h = act_t = 0;
            orph_n = 0;
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
            if (in_s[p] || in_t[p]) { return; } // became a terminal (root) while waiting
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
                    if (ts[x] == bk_time) { d += dst[x]; break; }
                    const int a = pa[x];
                    if (a == -1) { d += 1; break; }
                    if (a == -3) { d = INT_MAX; break; }
                    x = arcs[a].to;
                    ++d;
                }
                if (d == INT_MAX) { continue; }
                // mark the chain as verified for this time stamp
                int y = q, dd = d;
                while (ts[y] != bk_time) {
                    ts[y] = bk_time;
                    dst[y] = dd--;
                    const int a = pa[y];
                    if (a == -1) { break; }
                    y = arcs[a].to;
                }
                if (d < best_d) { best_d = d; best = i; }
            }
            if (best >= 0) {
                pa[p] = best;
                ts[p] = bk_time;
                dst[p] = best_d + 1;
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
