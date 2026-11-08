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

#ifndef HEIPROMAP_FLOW_BISECTION_H
#define HEIPROMAP_FLOW_BISECTION_H

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <random>
#include <vector>
#include <cstdlib>
#include <iostream>

#include "../definitions.h"
#include "../utility/aligned_array.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"

namespace HeiProMap {
    struct FlowBisectionConfiguration {
        // number of far-apart seed vertices; all pairs of seeds are tried
        u64 num_seeds = 6;
        // stop a pair as soon as its flow value reaches the best feasible cut found so far (exact)
        bool prune = true;
        // both terminals are first grown (alternating BFS) to this fraction of the sub graph weight
        f64 init_fraction = 0.1;
        // A pair is dropped as soon as its flow reaches prune_factor * (best feasible cut). 1.0 is exact (a
        // pair can never end below its flow), smaller values also drop pairs that might win by a few percent.
        f64 prune_factor = 0.95;
        // Only the best keep_fraction of all (pair, orientation) jobs, ranked by the flow after their first
        // solve (a lower bound of their cut), are continued. 1.0 continues all of them.
        f64 keep_fraction = 0.5;
        // Boykov-Kolmogorov max flow (search trees are kept when terminals are added); false: Dinic
        bool use_bk = true;
    };

    /**
     * Initial partitioning by flow based recursive bisection.
     *
     * Multisection (highest hierarchy level first, then the next, ...) where each level is a
     * k'-partition computed by recursive bisection of the block range into two halves. Each bisection:
     *  1. picks `num_seeds` far-apart vertices (farthest point sampling on hop distance),
     *  2. for every pair (a, b) puts a into the source set and b into the sink set and
     *     solves max-flow. While the min cut is not balanced (the side weight has to be in
     *     [W - K2*lmax, K1*lmax]), the lighter side is grown by piercing vertices next to
     *     its cut (preferring vertices that do not create an augmenting path and that are
     *     far from the other terminal) and the flow is augmented incrementally,
     *  3. keeps the balanced cut with minimum weight over all pairs (or the least violated one)
     *  and recurses on both sides.
     * The min cut can only grow when terminals grow, so the first balanced cut of a pair is
     * the best one of that pair.
     *
     * All work buffers are AlignedArray members that only ever grow (AlignedArray::initialize
     * reallocates only if the request exceeds the capacity), so after the first (largest) sub
     * problem nothing is allocated or freed anymore. Their size() is the capacity, the logical
     * sizes are tracked explicitly (nl, N, na, ...).
     */
    class FlowBisectionPartitioner {
        static constexpr s64 INF = (s64) 1 << 60;

        struct Arc {
            int to;
            int rev;
            s64 cap;
        };

        struct Cand {
            int far_flag; // 1 if the vertex is on the other side's region
            int neg_dist; // - hop distance to the other terminal
            u64 rnd;
            int v;

            bool operator<(const Cand &o) const {
                if (far_flag != o.far_flag) { return far_flag < o.far_flag; }
                if (neg_dist != o.neg_dist) { return neg_dist < o.neg_dist; }
                if (rnd != o.rnd) { return rnd < o.rnd; }
                return v < o.v;
            }
        };

        // fixed capacity stack of trivially destructible items
        template<typename T>
        struct Stack {
            AlignedArray<T> a;
            size_t n = 0;

            void reset(const size_t capacity) {
                a.initialize(capacity);
                n = 0;
            }

            void push(const T &x) { a[n++] = x; }
            T pop() { return a[--n]; }
            bool empty() const { return n == 0; }
        };

        // one (seed pair, orientation) of a bisection
        struct Job {
            int a, b; // terminals: a -> source, b -> sink
            weight_t lo, hi; // window of the source side weight
            int ia, ib; // seed indices (rows of dist_rows)
            bool mirrored;
            s64 lower; // flow after the first solve: a lower bound on the cut of this job
        };

        struct Task {
            size_t b, e; // range in the order array
            partition_t lo, hi; // block range
        };

        struct PairResult {
            bool valid = false;
            weight_t violation = std::numeric_limits<weight_t>::max();
            s64 cut = std::numeric_limits<s64>::max();
            AlignedArray<char> in_source;

            void reset() {
                valid = false;
                violation = std::numeric_limits<weight_t>::max();
                cut = std::numeric_limits<s64>::max();
            }

            // O(1), no allocation
            void swap(PairResult &o) {
                std::swap(valid, o.valid);
                std::swap(violation, o.violation);
                std::swap(cut, o.cut);
                std::swap(in_source, o.in_source);
            }
        };

        // Flow network of the current sub problem on the nodes 0..nl-1. s and t are implicit: a vertex with
        // in_s[v] has an infinite arc s -> v, a vertex with in_t[v] one v -> t. A simple augmenting path never
        // returns to s or leaves t, so no capacity has to be stored for these arcs.
        int nl = 0;
        int N = 0; // = nl, size of the per vertex arrays
        AlignedArray<int> src_list; // vertices with in_s
        // Boykov-Kolmogorov search trees. label: 0 free, 1 in the source tree, 2 in the sink tree. pa[v]: arc
        // (in v's list) to the parent, -1 root (terminal), -3 orphan. At termination the source tree is the set
        // reachable from s and the sink tree the set that reaches t in the residual graph.
        bool use_bk = true;
        AlignedArray<char> label, isact;
        AlignedArray<int> pa, act, orph, ts, dst;
        int bk_time = 0;
        int act_h = 0, act_t = 0, orph_n = 0;
        int src_n = 0;
        int dlevel = -1; // level of the sink terminals of the current Dinic phase
        size_t na = 0; // number of arcs
        AlignedArray<Arc> arcs;
        AlignedArray<Arc> base_arcs;
        AlignedArray<int> first; // N + 1
        AlignedArray<int> ptr;
        AlignedArray<int> level;
        AlignedArray<int> queue; // generic BFS queue
        AlignedArray<int> bq; // queue of the Dinic BFS (levels are reset through it)
        int last_qt = 0;
        AlignedArray<int> path;
        size_t pn = 0;
        AlignedArray<weight_t> vw;
        weight_t total_w = 0;

        AlignedArray<char> in_s, in_t, reach_s, reach_t;
        AlignedArray<int> stamp;
        int stamp_id = 0;

        // scratch
        AlignedArray<int> deg_buf, pos_buf; // build_network
        AlignedArray<char> owner; // grow_initial_terminals
        AlignedArray<int> gq_s, gq_t;
        AlignedArray<Cand> cand;
        size_t cand_n = 0;
        AlignedArray<int> seeds, tmp_dist, min_d, dist_rows; // bisect: dist_rows = n_seeds x nl
        AlignedArray<char> chosen;
        AlignedArray<char> first_part;
        PairResult best_res, mir_res, scr_res;
        AlignedArray<vertex_t> scratch_order; // stable partition of a range
        AlignedArray<weight_t> prefix, mw;
        AlignedArray<size_t> blk_b, blk_e; // result of kway_partition: range of every block
        AlignedArray<int> loc;
        Stack<Task> task_stack;
        // snapshots of the networks after the first solve of every job of a bisection
        AlignedArray<Arc> slot_arcs;
        AlignedArray<char> slot_in_s, slot_in_t, slot_label;
        AlignedArray<int> slot_pa;
        AlignedArray<Job> jobs;
        AlignedArray<int> job_order;

        std::mt19937_64 rng;
        bool debug = false;
        u64 pair_salt = 0;
        bool prune = true;
        f64 init_fraction = 0.1;
        f64 prune_factor = 1.0;
        f64 keep_fraction = 1.0;

        int tail(const int a) const { return arcs[arcs[a].rev].to; }

        void build_network(const CSRGraph &g, const vertex_t *verts, const int cnt) {
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

            ptr.initialize(N);
            level.initialize(N, -1);
            last_qt = 0;
            queue.initialize(N);
            bq.initialize(N);
            path.initialize(N);
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
            reach_s.initialize(N);
            reach_t.initialize(N);
            stamp.initialize(N, 0);
            stamp_id = 0;
            owner.initialize(N);
            gq_s.initialize(N);
            gq_t.initialize(N);
            cand.initialize(N);
            tmp_dist.initialize(N);
            min_d.initialize(N);
            chosen.initialize(N);
            first_part.initialize(N);
        }

        void add_source(const int v) {
            in_s[v] = 1;
            src_list[src_n++] = v;
        }

        // Dinic level graph from all source terminals (level 0). Stops at the first level that contains a
        // sink terminal (dlevel); only sink terminals of that level end an augmenting path.
        bool bfs() {
            // reset only what the previous search touched
            for (int i = 0; i < last_qt; ++i) { level[bq[i]] = -1; }
            int qh = 0, qt = 0;
            dlevel = -1;
            for (int i = 0; i < src_n; ++i) {
                const int v = src_list[i];
                level[v] = 0;
                ptr[v] = first[v];
                bq[qt++] = v;
            }
            while (qh < qt) {
                const int u = bq[qh++];
                if (dlevel != -1 && level[u] >= dlevel) { break; } // deeper levels are not used by Dinic
                for (int i = first[u]; i < first[u + 1]; ++i) {
                    const Arc &a = arcs[i];
                    if (a.cap > 0 && level[a.to] == -1) {
                        level[a.to] = level[u] + 1;
                        ptr[a.to] = first[a.to];
                        bq[qt++] = a.to;
                        if (dlevel == -1 && in_t[a.to]) { dlevel = level[a.to]; }
                    }
                }
            }
            last_qt = qt;
            return dlevel != -1;
        }

        // augments the current flow to a maximum flow, returns the added flow.
        // afterwards level[v] != -1 iff v is reachable from s in the residual graph.
        // If the added flow reaches `limit` it returns early (the caller prunes the pair anyway).
        s64 maxflow(const s64 limit = INF) {
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "maxflow");
            s64 flow = 0;
            while (bfs()) {
                for (int si = 0; si < src_n; ++si) {
                    int u = src_list[si];
                    if (level[u] == -1) { continue; } // dead end
                    pn = 0;
                    while (true) {
                        if (in_t[u] && level[u] == dlevel) {
                            s64 bn = INF;
                            for (size_t i = 0; i < pn; ++i) { bn = std::min(bn, arcs[path[i]].cap); }
                            size_t cut = pn;
                            for (size_t i = 0; i < pn; ++i) {
                                Arc &a = arcs[path[i]];
                                a.cap -= bn;
                                arcs[a.rev].cap += bn;
                                if (a.cap == 0 && cut == pn) { cut = i; }
                            }
                            flow += bn;
                            if (flow >= limit) { return flow; }
                            u = tail(path[cut]);
                            pn = cut;
                            continue;
                        }
                        bool advanced = false;
                        for (int &i = ptr[u]; i < first[u + 1]; ++i) {
                            const Arc &a = arcs[i];
                            if (a.cap > 0 && level[a.to] == level[u] + 1) {
                                advanced = true;
                                break;
                            }
                        }
                        if (advanced) {
                            path[pn++] = ptr[u];
                            u = arcs[ptr[u]].to;
                        } else {
                            level[u] = -1; // dead end
                            if (pn == 0) { break; }
                            const int a = path[--pn];
                            u = tail(a);
                            ++ptr[u];
                        }
                    }
                }
            }
            return flow;
        }

        // ---------------------------------------------------------------- Boykov-Kolmogorov
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

        // vertices that can still reach t in the residual graph
        void compute_reach_t() {
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "reach_t");
            std::fill_n(reach_t.get_ptr(), N, 0);
            int qh = 0, qt = 0;
            for (int v = 0; v < nl; ++v) {
                if (in_t[v]) {
                    reach_t[v] = 1;
                    queue[qt++] = v;
                }
            }
            while (qh < qt) {
                const int x = queue[qh++];
                for (int i = first[x]; i < first[x + 1]; ++i) {
                    const Arc &a = arcs[i];
                    if (!reach_t[a.to] && arcs[a.rev].cap > 0) {
                        reach_t[a.to] = 1;
                        queue[qt++] = a.to;
                    }
                }
            }
        }

        // hop distance inside the sub graph, unreachable = large; dist has nl entries
        void hop_distance(const int src, int *dist) {
            constexpr int FAR = 1 << 29;
            std::fill_n(dist, nl, FAR);
            int qh = 0, qt = 0;
            dist[src] = 0;
            queue[qt++] = src;
            while (qh < qt) {
                const int u = queue[qh++];
                for (int i = first[u]; i < first[u + 1]; ++i) {
                    const int v = arcs[i].to;
                    if (v < nl && dist[v] == FAR) {
                        dist[v] = dist[u] + 1;
                        queue[qt++] = v;
                    }
                }
            }
        }

        static weight_t violation(const weight_t a, const weight_t lo, const weight_t hi) {
            if (a < lo) { return lo - a; }
            if (a > hi) { return a - hi; }
            return 0;
        }

        // Candidates to pierce into one side: vertices next to `reach` (the side's current region)
        // that are not terminals yet. Sorted so that the best are first.
        void collect_candidates(const char *reach, const char *other_reach, const int *dist_other) {
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "collect_candidates");
            cand_n = 0;
            ++stamp_id;
            for (int v = 0; v < nl; ++v) {
                if (!reach[v]) { continue; }
                for (int i = first[v]; i < first[v + 1]; ++i) {
                    const int w = arcs[i].to;
                    if (w >= nl || reach[w] || in_s[w] || in_t[w] || stamp[w] == stamp_id) { continue; }
                    stamp[w] = stamp_id;
                    cand[cand_n++] = Cand{other_reach[w] ? 1 : 0, -dist_other[w], rng(), w};
                }
            }
            if (cand_n == 0) {
                // disconnected: any free vertex
                for (int w = 0; w < nl; ++w) {
                    if (reach[w] || in_s[w] || in_t[w]) { continue; }
                    cand[cand_n++] = Cand{other_reach[w] ? 1 : 0, -dist_other[w], rng(), w};
                }
            }
            std::sort(cand.get_ptr(), cand.get_ptr() + cand_n);
        }

        // pierces candidates into the side until `batch` weight is added without exceeding max_side
        bool pierce(const bool source_side, const weight_t side_w, const weight_t max_side, const weight_t batch) {
            weight_t added = 0;
            bool any = false;
            for (size_t ci = 0; ci < cand_n; ++ci) {
                const int w = cand[ci].v;
                if (side_w + added + vw[w] > max_side) { continue; }
                if (source_side) {
                    add_source(w);
                    if (use_bk) { bk_add_source(w); }
                } else {
                    in_t[w] = 1;
                    if (use_bk) { bk_add_sink(w); }
                }
                added += vw[w];
                any = true;
                if (added >= batch) { break; }
            }
            return any;
        }

        // Grows the terminals from a and b by alternating BFS (lighter side first) until both reach
        // `target` weight (bounded by what the window allows) or no vertex is left to claim.
        void grow_initial_terminals(const int a, const int b, const weight_t target_s, const weight_t target_t) {
            std::fill_n(owner.get_ptr(), nl, 0);
            size_t ns = 0, nt = 0;
            gq_s[ns++] = a;
            gq_t[nt++] = b;
            size_t hs = 0, ht = 0;
            owner[a] = 1;
            owner[b] = 2;
            weight_t ws = vw[a], wt = vw[b];
            while (true) {
                const bool s_open = ws < target_s && hs < ns;
                const bool t_open = wt < target_t && ht < nt;
                if (!s_open && !t_open) { break; }
                const bool take_s = s_open && (!t_open || ws <= wt);
                AlignedArray<int> &q = take_s ? gq_s : gq_t;
                size_t &qn = take_s ? ns : nt;
                size_t &h = take_s ? hs : ht;
                const int u = q[h++];
                for (int i = first[u]; i < first[u + 1]; ++i) {
                    const int v = arcs[i].to;
                    if (v >= nl || owner[v]) { continue; }
                    if (take_s ? (ws + vw[v] > target_s) : (wt + vw[v] > target_t)) { continue; }
                    owner[v] = take_s ? 1 : 2;
                    q[qn++] = v;
                    (take_s ? ws : wt) += vw[v];
                }
            }
            for (int v = 0; v < nl; ++v) {
                if (owner[v] == 1 && !in_s[v]) {
                    add_source(v);
                } else if (owner[v] == 2 && !in_t[v]) {
                    in_t[v] = 1;
                }
            }
        }

        // Sets up the network of a job (terminals a, b pre-grown), solves the first max flow and saves the
        // network in `slot`. The flow value is a lower bound for the cut the job can end with.
        void prepare_job(Job &j, const int slot) {
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "prepare_job");
            std::memcpy(arcs.get_ptr(), base_arcs.get_ptr(), na * sizeof(Arc));
            std::fill_n(in_s.get_ptr(), N, 0);
            std::fill_n(in_t.get_ptr(), N, 0);
            src_n = 0;
            add_source(j.a);
            in_t[j.b] = 1;
            if (init_fraction > 0.0) {
                const weight_t target = (weight_t) (init_fraction * (f64) total_w);
                grow_initial_terminals(j.a, j.b, std::min(target, j.hi), std::min(target, total_w - j.lo));
            }
            if (use_bk) {
                bk_init();
                j.lower = bk_maxflow();
                std::memcpy(slot_label.get_ptr() + (size_t) slot * N, label.get_ptr(), N);
                std::memcpy(slot_pa.get_ptr() + (size_t) slot * N, pa.get_ptr(), N * sizeof(int));
            } else {
                j.lower = maxflow();
            }
            std::memcpy(slot_arcs.get_ptr() + (size_t) slot * na, arcs.get_ptr(), na * sizeof(Arc));
            std::memcpy(slot_in_s.get_ptr() + (size_t) slot * N, in_s.get_ptr(), N);
            std::memcpy(slot_in_t.get_ptr() + (size_t) slot * N, in_t.get_ptr(), N);
        }

        // Continues a job from its snapshot: pierce the lighter side until the cut is balanced.
        // The result is merged into `best` (scr_res is used as scratch).
        void run_job(const Job &j, const int slot, const int *dist_a, const int *dist_b, PairResult &best,
                     const s64 cut_bound = INF) {
            const int a = j.a, b = j.b;
            const weight_t lo = j.lo, hi = j.hi;
            rng.seed(pair_salt ^ ((u64) a * 1000003ULL + (u64) b * 7919ULL + 17));
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "run_pair");
            std::memcpy(arcs.get_ptr(), slot_arcs.get_ptr() + (size_t) slot * na, na * sizeof(Arc));
            std::memcpy(in_s.get_ptr(), slot_in_s.get_ptr() + (size_t) slot * N, N);
            std::memcpy(in_t.get_ptr(), slot_in_t.get_ptr() + (size_t) slot * N, N);
            src_n = 0;
            for (int v = 0; v < nl; ++v) { if (in_s[v]) { src_list[src_n++] = v; } }
            if (use_bk) {
                std::memcpy(label.get_ptr(), slot_label.get_ptr() + (size_t) slot * N, N);
                std::memcpy(pa.get_ptr(), slot_pa.get_ptr() + (size_t) slot * N, N * sizeof(int));
                std::fill_n(isact.get_ptr(), N, 0);
                act_h = act_t = 0;
                orph_n = 0;
            }

            const weight_t avg_w = std::max<weight_t>(1, total_w / std::max(1, nl));
            s64 flow = j.lower;
            PairResult &cur = scr_res;
            cur.reset();

            for (int iter = 0; iter < nl + 2; ++iter) {
                if (use_bk) {
                    flow += bk_maxflow(cut_bound - flow);
                    if (flow >= cut_bound) { return; }
                    for (int v = 0; v < nl; ++v) {
                        reach_s[v] = label[v] == 1;
                        reach_t[v] = label[v] == 2;
                    }
                } else {
                    flow += maxflow(cut_bound - flow);
                    // the cut never shrinks when terminals grow: this pair cannot beat a known feasible cut
                    if (flow >= cut_bound) { return; }
                    for (int v = 0; v < nl; ++v) { reach_s[v] = level[v] != -1; }
                    compute_reach_t();
                }

                weight_t w_rs = 0, w_ct = 0;
                for (int v = 0; v < nl; ++v) {
                    if (reach_s[v]) { w_rs += vw[v]; }
                    if (reach_t[v]) { w_ct += vw[v]; }
                }
                const weight_t a1 = w_rs; // source closest cut
                const weight_t a2 = total_w - w_ct; // sink closest cut
                const weight_t v1 = violation(a1, lo, hi);
                const weight_t v2 = violation(a2, lo, hi);
                const weight_t v = std::min(v1, v2);

                if (!cur.valid || v < cur.violation) {
                    cur.valid = true;
                    cur.violation = v;
                    cur.cut = flow;
                    cur.in_source.initialize(nl);
                    if (v1 <= v2) {
                        for (int x = 0; x < nl; ++x) { cur.in_source[x] = reach_s[x]; }
                    } else {
                        for (int x = 0; x < nl; ++x) { cur.in_source[x] = !reach_t[x]; }
                    }
                }
                if (v == 0) { break; }

                // decide which side to grow
                bool grow_source;
                if (a1 > hi) {
                    grow_source = false; // source region alone is too heavy
                } else if (a2 < lo) {
                    grow_source = true; // sink region alone is too heavy
                } else {
                    grow_source = w_rs <= w_ct;
                }

                bool done = false;
                for (int attempt = 0; attempt < 2 && !done; ++attempt) {
                    if (grow_source) {
                        const weight_t deficit = std::max<weight_t>(1, lo - w_rs);
                        collect_candidates(reach_s.get_ptr(), reach_t.get_ptr(), dist_b);
                        done = pierce(true, w_rs, hi, std::max(avg_w, deficit / 8));
                    } else {
                        const weight_t deficit = std::max<weight_t>(1, (total_w - hi) - w_ct);
                        collect_candidates(reach_t.get_ptr(), reach_s.get_ptr(), dist_a);
                        done = pierce(false, w_ct, total_w - lo, std::max(avg_w, deficit / 8));
                    }
                    grow_source = !grow_source;
                }
                if (!done) { break; }
            }

            if (!cur.valid) { return; }
            const bool better = !best.valid
                                || cur.violation < best.violation
                                || (cur.violation == best.violation && cur.cut < best.cut);
            if (better) { best.swap(cur); }
        }

        s64 bound(const PairResult &best) const {
            return (prune && best.valid && best.violation == 0) ? (s64) ((f64) best.cut * prune_factor) : INF;
        }

        // splits the sub problem; first_part[v] = 1 if local vertex v belongs to the first part
        void bisect(const weight_t lo, const weight_t hi, const bool symmetric, const u64 num_seeds) {
            std::fill_n(first_part.get_ptr(), nl, 1);
            if (nl < 2) { return; }

            // farthest point sampling of seeds
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "seed_selection");
            const int n_seeds = (int) std::min<u64>(std::max<u64>(2, num_seeds), (u64) nl);
            seeds.initialize(n_seeds);
            dist_rows.initialize((size_t) n_seeds * (size_t) nl);
            int ns = 0;
            hop_distance((int) (rng() % (u64) nl), tmp_dist.get_ptr());
            int s0 = 0;
            for (int v = 0; v < nl; ++v) { if (tmp_dist[v] > tmp_dist[s0]) { s0 = v; } }
            std::fill_n(min_d.get_ptr(), nl, std::numeric_limits<int>::max());
            std::fill_n(chosen.get_ptr(), nl, 0);
            int next = s0;
            while (ns < n_seeds) {
                int *row = dist_rows.get_ptr() + (size_t) ns * (size_t) nl;
                seeds[ns++] = next;
                chosen[next] = 1;
                hop_distance(next, row);
                int arg = -1;
                for (int v = 0; v < nl; ++v) {
                    min_d[v] = std::min(min_d[v], row[v]);
                    if (!chosen[v] && (arg < 0 || min_d[v] > min_d[arg])) { arg = v; }
                }
                if (arg < 0) { break; }
                next = arg;
            }

            best_res.reset();
            pair_salt = rng();
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "all_pairs");

            // all jobs: every seed pair, and the other orientation if the sides differ (seeds[j] must end up
            // in the first part: source side weight window mirrored)
            const size_t max_jobs = (size_t) ns * (size_t) (ns - 1);
            jobs.initialize(max_jobs);
            job_order.initialize(max_jobs);
            size_t nj = 0;
            for (int i = 0; i < ns; ++i) {
                for (int j = i + 1; j < ns; ++j) {
                    jobs[nj++] = Job{seeds[i], seeds[j], lo, hi, i, j, false, 0};
                    if (!symmetric) { jobs[nj++] = Job{seeds[j], seeds[i], total_w - hi, total_w - lo, j, i, true, 0}; }
                }
            }
            slot_arcs.initialize(nj * na);
            slot_in_s.initialize(nj * (size_t) N);
            slot_in_t.initialize(nj * (size_t) N);
            if (use_bk) {
                slot_label.initialize(nj * (size_t) N);
                slot_pa.initialize(nj * (size_t) N);
            }
            for (size_t q = 0; q < nj; ++q) {
                prepare_job(jobs[q], (int) q);
                job_order[q] = (int) q;
            }
            // most promising (smallest lower bound) first: finds a good feasible cut early, which then
            // prunes the remaining jobs
            std::stable_sort(job_order.get_ptr(), job_order.get_ptr() + nj,
                             [&](const int x, const int y) { return jobs[x].lower < jobs[y].lower; });

            const size_t run_jobs = std::max<size_t>(1, (size_t) std::ceil((double) nj * keep_fraction));
            for (size_t q = 0; q < run_jobs; ++q) {
                const int slot = job_order[q];
                const Job &job = jobs[slot];
                const s64 cb = bound(best_res);
                if (job.lower >= cb) { continue; } // cannot beat the best feasible cut
                const int *da = dist_rows.get_ptr() + (size_t) job.ia * (size_t) nl;
                const int *db = dist_rows.get_ptr() + (size_t) job.ib * (size_t) nl;
                if (!job.mirrored) {
                    run_job(job, slot, da, db, best_res, cb);
                } else {
                    mir_res.reset();
                    run_job(job, slot, da, db, mir_res, cb);
                    if (mir_res.valid) {
                        for (int x = 0; x < nl; ++x) { mir_res.in_source[x] = !mir_res.in_source[x]; }
                        const bool better = !best_res.valid
                                            || mir_res.violation < best_res.violation
                                            || (mir_res.violation == best_res.violation && mir_res.cut < best_res.cut);
                        if (better) { best_res.swap(mir_res); }
                    }
                }
            }

            // the number of random draws inside the pairs depends on pruning: restore a fixed stream
            rng.seed(pair_salt + 1);

            if (best_res.valid) { std::memcpy(first_part.get_ptr(), best_res.in_source.get_ptr(), nl); }
        }

        /**
         * Computes a k'-partition of the vertices order[begin, end) (induced subgraph of g) by recursive
         * flow based bisection, in place: afterwards block i owns order[blk_b[i], blk_e[i]) (the blocks
         * appear in increasing order). Block i may weigh at most max_weights[i]. `loc` is -1 outside
         * of this function.
         */
        void kway_partition(const CSRGraph &g, vertex_t *order, const size_t begin, const size_t end,
                            const partition_t k_prime, const weight_t *max_weights,
                            const FlowBisectionConfiguration &config) {
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "kway_partition");
            blk_b.initialize(k_prime);
            blk_e.initialize(k_prime);
            for (partition_t i = 0; i < k_prime; ++i) { blk_b[i] = begin; blk_e[i] = begin; }
            prefix.initialize(k_prime + 1);
            prefix[0] = 0;
            for (partition_t i = 0; i < k_prime; ++i) { prefix[i + 1] = prefix[i] + max_weights[i]; }
            scratch_order.initialize(end - begin);

            task_stack.reset(k_prime + 2);
            task_stack.push({begin, end, 0, k_prime});

            while (!task_stack.empty()) {
                const Task task = task_stack.pop();
                const partition_t K = task.hi - task.lo;
                const size_t cnt = task.e - task.b;

                if (K == 1 || cnt < 2) {
                    // fewer vertices than needed: everything goes to the first block of the range
                    blk_b[task.lo] = task.b;
                    blk_e[task.lo] = task.e;
                    continue;
                }

                const partition_t mid = task.lo + (K + 1) / 2;
                const weight_t max1 = prefix[mid] - prefix[task.lo];
                const weight_t max2 = prefix[task.hi] - prefix[mid];

                vertex_t *verts = order + task.b;
                for (size_t i = 0; i < cnt; ++i) { loc[verts[i]] = (int) i; }
                build_network(g, verts, (int) cnt);
                for (size_t i = 0; i < cnt; ++i) { loc[verts[i]] = -1; }

                const weight_t lo_w = std::max<weight_t>(0, total_w - max2);
                bisect(lo_w, max1, max1 == max2, config.num_seeds);
                if (debug) {
                    weight_t w1 = 0;
                    size_t n1 = 0;
                    for (int v = 0; v < nl; ++v) { if (first_part[v]) { w1 += vw[v]; ++n1; } }
                    const weight_t w2 = total_w - w1;
                    std::cerr << "[FB]     bisect blocks[" << task.lo << "," << task.hi << ") n=" << nl
                              << " W=" << total_w << " max=" << max1 << "|" << max2
                              << " window=[" << lo_w << "," << max1 << "] -> " << w1 << "|" << w2
                              << " (n " << n1 << "|" << (nl - n1) << ")"
                              << ((w1 > max1 || w2 > max2) ? "  VIOLATED" : "") << std::endl;
                }

                // stable partition of the range: first part to the front
                size_t n1 = 0, n2 = 0;
                for (size_t i = 0; i < cnt; ++i) {
                    if (first_part[i]) { verts[n1++] = verts[i]; } else { scratch_order[n2++] = verts[i]; }
                }
                std::memcpy(verts + n1, scratch_order.get_ptr(), n2 * sizeof(vertex_t));

                task_stack.push({task.b, task.b + n1, task.lo, mid});
                task_stack.push({task.b + n1, task.e, mid, task.hi});
            }
        }

    public:
        /**
         * Multisection driven by flow based recursive bisection: split along the highest hierarchy
         * level first, then the next one, down to the lowest. Every part of a level may weigh at most
         * (#PEs below it) * lmax.
         */
        static void partition(const CSRGraph &g,
                              PartitionManager &p_manager,
                              const std::vector<partition_t> &hierarchy,
                              const f64 imbalance,
                              const FlowBisectionConfiguration &config,
                              const u64 seed) {
            FlowBisectionPartitioner p;
            p.debug = std::getenv("FLOW_BISECTION_DEBUG") != nullptr;
            p.prune = config.prune;
            p.init_fraction = config.init_fraction;
            p.prune_factor = config.prune_factor;
            p.keep_fraction = config.keep_fraction;
            p.use_bk = config.use_bk;
            p.rng.seed(seed * 0x9E3779B97F4A7C15ULL + 12345);

            const size_t l = hierarchy.size();
            AlignedArray<partition_t> stride; // PEs below one part of level j
            stride.initialize(l);
            partition_t k = 1;
            size_t stack_cap = 2;
            for (size_t j = 0; j < l; ++j) {
                stride[j] = k;
                k *= hierarchy[j];
                stack_cap += hierarchy[j];
            }
            const weight_t lmax = (weight_t) std::ceil((1.0 + imbalance) * (f64) g.g_weight / (f64) k);

            p.loc.initialize(g.n, -1);
            AlignedArray<vertex_t> order;
            order.initialize(g.n);
            for (vertex_t u = 0; u < g.n; ++u) { order[u] = u; }

            struct Item {
                size_t b, e; // range in order
                size_t level; // index into hierarchy, split into hierarchy[level] parts
                partition_t offset;
            };
            Stack<Item> stack;
            stack.reset(stack_cap);
            stack.push({0, (size_t) g.n, l - 1, 0});

            while (!stack.empty()) {
                const Item item = stack.pop();

                const partition_t kk = hierarchy[item.level];
                p.mw.initialize(kk, (weight_t) stride[item.level] * lmax);
                p.kway_partition(g, order.get_ptr(), item.b, item.e, kk, p.mw.get_ptr(), config);

                if (p.debug) {
                    weight_t parent = 0, mx = 0, mn = std::numeric_limits<weight_t>::max();
                    size_t over = 0;
                    for (partition_t i = 0; i < kk; ++i) {
                        weight_t w = 0;
                        for (size_t x = p.blk_b[i]; x < p.blk_e[i]; ++x) { w += g.v_weight(order[x]); }
                        parent += w;
                        mx = std::max(mx, w);
                        mn = std::min(mn, w);
                        over += w > p.mw[i];
                    }
                    std::cerr << "[FB] level " << item.level << " offset " << item.offset << ": split weight "
                              << parent << " (n=" << (item.e - item.b) << ") into " << kk << " parts, max allowed "
                              << p.mw[0] << ", part weights min=" << mn << " max=" << mx
                              << " overloaded=" << over << std::endl;
                }
                for (partition_t i = 0; i < kk; ++i) {
                    if (item.level == 0) {
                        for (size_t x = p.blk_b[i]; x < p.blk_e[i]; ++x) {
                            const vertex_t u = order[x];
                            p_manager.set(u, g.v_weight(u), item.offset + i);
                        }
                    } else {
                        stack.push({p.blk_b[i], p.blk_e[i], item.level - 1, item.offset + i * stride[item.level]});
                    }
                }
            }
        }
    };
} // namespace HeiProMap

#endif // HEIPROMAP_FLOW_BISECTION_H
