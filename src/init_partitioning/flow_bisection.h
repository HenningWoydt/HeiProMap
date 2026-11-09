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
#include <cmath>
#include <cstdlib>
#include <limits>
#include <random>
#include <vector>

#include "flow_network.h"
#include "../datastructures/partition_manager.h"

namespace HeiProMap {
    struct FlowBisectionConfiguration {
        // number of far-apart seed vertices; all pairs of seeds are tried
        u64 num_seeds = 6;        // top hierarchy level
        u64 num_seeds_mid = 2;    // levels between top and bottom
        u64 num_seeds_bottom = 2; // bottom level
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
        // The first `bound_jobs` jobs are solved completely and the best of them is run; its cut bounds the
        // first solves of all other jobs (stopped early). 0 = solve all first solves completely.
        u64 bound_jobs = 4;
        // Large sub problems get fewer seeds: a bisection with `arcs` arcs uses at most 1 + seed_budget / arcs seeds
        // (the star around the first seed costs about one first solve per seed, and a solve scales with the arcs).
        // 0 = no limit.
        u64 seed_budget = 300000;
        // Every hierarchy level gets its own imbalance (same formula as the multisection): the upper levels are
        // tight and leave slack to the lower ones. false: every part may weigh stride * lmax (global imbalance).
        bool adaptive_imbalance = true;
        // Seed diversity: every seed is drawn uniformly from the `seed_pool` vertices that are farthest from the
        // seeds so far (1 = always the farthest one, deterministic up to the start vertex).
        u64 seed_pool = 1;
        // Flow refinement of every bisection: a corridor around the best cut is freed (the `refine_band` fraction of
        // the weight closest to the cut), everything outside is terminal, and the flow is solved again (with the
        // same balanced piercing). A better cut is accepted and the corridor tried again, otherwise the corridor is
        // widened (0.2, 0.5, 1, 2, 4 times `refine_band`); after an improvement it starts again with the narrowest. At most
        // `refine_rounds` improvements (it normally stops earlier, when no corridor improves the cut); 0 = off.
        u64 refine_rounds = 100;
        f64 refine_band = 0.1;
        // After the recursive bisection of a level: passes of flow refinement between all pairs of adjacent blocks
        // (largest shared cut first), each pair as a bisection of the union of both blocks. 0 = off.
        u64 kway_rounds = 3;
        // Bisections of fewer than `small_n` vertices are cheap: they use `small_seeds` seeds (no seed budget) and
        // continue every job (no keep_fraction, exact pruning). small_n = 0 = off.
        u64 small_n = 300;
        u64 small_seeds = 8;
        // levels below the top: -1 = same as above, otherwise these values replace kway_rounds / small_seeds /
        // refine_rounds there (the top level decides most of the QAP, the lower levels are cheap to settle)
        s64 low_kway_rounds = -1, low_small_seeds = 2, low_refine_rounds = -1;
        bool small_exact = true; // small bisections: keep_fraction 1 and prune_factor 1
        // seed_budget of the top hierarchy level (its few, large bisections decide most of the QAP); 0 = no limit
        u64 top_seed_budget = 0;
        // seed pairs that do not contain the first seed: only this fraction (the farthest apart) is tried; 1 = all
        f64 pair_fraction = 1.0;
        int pair_proxy = 0; // 0 = hop distance, 1 = sum of the seed ball cuts
        // move vertices out of overloaded blocks after every k'-partition (smallest cut increase first)
        bool repair = true;
        // first max flow of every seed pair with IBFS instead of BK (identical result, currently slower)
        bool ibfs = false;
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
     *     its cut (preferring vertices that do not create an augmenting path and that increase
     *     the cut the least) and the flow is augmented incrementally,
     *  3. keeps the balanced cut with minimum weight over all pairs (or the least violated one)
     *  and recurses on both sides.
     * The min cut can only grow when terminals grow, so the first balanced cut of a pair is
     * the best one of that pair.
     *
     * The flow network and the max flow are in FlowNetwork. All work buffers are AlignedArray members that only
     * ever grow, so after the first (largest) sub problem nothing is allocated or freed anymore.
     */
    class FlowBisectionPartitioner : FlowNetwork {
        // A piercing step adds about 1/PIERCE_DIVISOR of the weight that the side still lacks to reach the window
        // (at least one average vertex). Larger steps need fewer flow computations, smaller ones overshoot less.
        static constexpr weight_t PIERCE_DIVISOR = 2;

        struct Cand {
            int far_flag; // 1 if the vertex is on the other side's region
            int neg_dist; // - hop distance to the other terminal
            s64 delta; // increase of the cut when the vertex becomes a terminal of the side
            u64 rnd;
            int v;

            // best first: vertices outside the other side's tree, then the smallest cut increase, then the
            // largest distance to the other terminal
            bool operator<(const Cand &o) const {
                if (far_flag != o.far_flag) { return far_flag < o.far_flag; }
                if (delta != o.delta) { return delta < o.delta; }
                if (neg_dist != o.neg_dist) { return neg_dist < o.neg_dist; }
                if (rnd != o.rnd) { return rnd < o.rnd; }
                return v < o.v;
            }
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

        // scratch
        AlignedArray<char> owner; // marks of grow_initial_terminals, pick_seeds and collect_candidates
        AlignedArray<int> gq_s, gq_t; // BFS queues (also used by hop_distance)
        AlignedArray<Cand> cand;
        size_t cand_n = 0;
        std::vector<int> seeds;
        AlignedArray<int> tmp_dist, dist_rows; // bisect: dist_rows = n_seeds x nl
        PairResult best_res, mir_res, scr_res;
        AlignedArray<vertex_t> scratch_order; // stable partition of a range
        std::vector<weight_t> prefix, mw;
        std::vector<size_t> blk_b, blk_e; // result of kway_partition: range of every block
        AlignedArray<int> loc;
        std::vector<int> blk_of; // refine_kway: block of a vertex, -1 outside
        std::vector<Task> task_stack;
        // snapshots of the networks after the first solve of every job of a bisection
        AlignedArray<Arc> slot_arcs;
        AlignedArray<char> slot_label;
        AlignedArray<int> slot_pa;
        std::vector<Job> jobs;
        std::vector<int> job_order;
        std::vector<s64> lam; // lam[i * ns + k]: known lower bound of the min cut of seeds i and k

        std::vector<Cand> heap;
        AlignedArray<s64> cdelta;
        AlignedArray<u64> crnd;
        char pq_side = 0, pq_other = 0;
        const int *pq_dist = nullptr;
        std::mt19937_64 rng;
        u64 pair_salt = 0;
        bool prune = true;
        f64 init_fraction = 0.1;
        f64 prune_factor = 1.0;
        f64 keep_fraction = 1.0;
        size_t bound_jobs = 0;
        size_t seed_budget = 0;
        size_t seed_pool = 1;
        size_t refine_rounds = 0;
        f64 refine_band = 0.1;
        size_t kway_rounds = 0;
        bool use_ibfs = false;
        bool repair = true;
        f64 pair_fraction = 1.0;
        int pair_proxy = 0; // 0 = hop distance, 1 = sum of the seed ball cuts
        size_t small_n = 0, small_seeds = 0;
        bool small_exact = true;
        AlignedArray<int> ref_q, ref_d, ref_d1, ref_d2;
        AlignedArray<char> ref_core;
        std::vector<std::pair<int, int> > far_top; // (distance, vertex), farthest first

        /// network of the vertices verts[0, cnt) (loc must be set for them) and the buffers of the partitioner
        void setup_network(const CSRGraph &g, const vertex_t *verts, const int cnt) {
            build_network(g, verts, cnt, loc.get_ptr());
            owner.initialize(N);
            gq_s.initialize(N);
            gq_t.initialize(N);
            cand.initialize(N);
            cdelta.initialize(N);
            crnd.initialize(N);
            tmp_dist.initialize(N);
        }

        // hop distance inside the sub graph, unreachable = large; dist has nl entries
        void hop_distance(const int src, int *dist) {
            constexpr int FAR = 1 << 29;
            std::fill_n(dist, nl, FAR);
            int qh = 0, qt = 0;
            dist[src] = 0;
            gq_s[qt++] = src;
            while (qh < qt) {
                const int u = gq_s[qh++];
                for (int i = first[u]; i < first[u + 1]; ++i) {
                    const int v = arcs[i].to;
                    if (v < nl && dist[v] == FAR) {
                        dist[v] = dist[u] + 1;
                        gq_s[qt++] = v;
                    }
                }
            }
        }

        static weight_t violation(const weight_t a, const weight_t lo, const weight_t hi) {
            if (a < lo) { return lo - a; }
            if (a > hi) { return a - hi; }
            return 0;
        }

        // Initial candidates to pierce into one side: vertices next to the side's current region that are not
        // terminals yet, with their score (see Cand).
        // how much the cut grows when w (outside of the side) becomes a terminal of `side`: its edges to the rest
        // become cut edges, its edges into the side stop being cut
        s64 cut_increase(const int w, const char side) const {
            s64 d = 0;
            for (int i = first[w]; i < first[w + 1]; ++i) {
                d += label[arcs[i].to] == side ? -base_arcs[i].cap : base_arcs[i].cap;
            }
            return d;
        }

        void collect_candidates(const char side, const char other_side, const int *dist_other) {
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "collect_candidates");
            cand_n = 0;
            std::fill_n(owner.get_ptr(), nl, 0); // candidates already collected
            for (int v = 0; v < nl; ++v) {
                if (label[v] != side) { continue; }
                for (int i = first[v]; i < first[v + 1]; ++i) {
                    const int w = arcs[i].to;
                    if (w >= nl || label[w] == side || is_terminal(w) || owner[w]) { continue; }
                    owner[w] = 1;
                    cand[cand_n++] = Cand{label[w] == other_side ? 1 : 0, -dist_other[w], cut_increase(w, side), rng(), w};
                }
            }
            if (cand_n == 0) {
                // disconnected: any free vertex
                for (int w = 0; w < nl; ++w) {
                    if (label[w] == side || is_terminal(w)) { continue; }
                    cand[cand_n++] = Cand{label[w] == other_side ? 1 : 0, -dist_other[w], cut_increase(w, side), rng(), w};
                }
            }
            pq_side = side;
            pq_other = other_side;
            pq_dist = dist_other;
            for (size_t ci = 0; ci < cand_n; ++ci) {
                cdelta[cand[ci].v] = cand[ci].delta;
                crnd[cand[ci].v] = cand[ci].rnd;
            }
        }

        // Pierces vertices into the side until `batch` weight is added without exceeding max_side. The frontier is
        // dynamic: when a vertex joins the side, the cut increase of its neighbours changes and new vertices become
        // adjacent. Min-heap with lazy deletion (an entry is stale if its score differs from the current score).
        bool pierce(const bool source_side, const weight_t side_w, const weight_t max_side, const weight_t batch) {
            const char side = pq_side;
            weight_t added = 0;
            bool any = false;
            auto greater = [](const Cand &a, const Cand &b) { return b < a; };
            heap.assign(cand.get_ptr(), cand.get_ptr() + cand_n);
            std::make_heap(heap.begin(), heap.end(), greater);
            while (!heap.empty()) {
                std::pop_heap(heap.begin(), heap.end(), greater);
                const Cand c = heap.back();
                heap.pop_back();
                const int w = c.v;
                if (label[w] == side || is_terminal(w) || c.delta != cdelta[w]) { continue; }
                if (side_w + added + vw[w] > max_side) { continue; }
                if (source_side) {
                    bk_add_source(w);
                } else {
                    bk_add_sink(w);
                }
                added += vw[w];
                any = true;
                if (added >= batch) { break; }
                for (int i = first[w]; i < first[w + 1]; ++i) {
                    const int x = arcs[i].to;
                    if (label[x] == side || is_terminal(x)) { continue; }
                    if (!owner[x]) {
                        owner[x] = 1;
                        crnd[x] = rng();
                    }
                    const s64 d = cut_increase(x, side);
                    cdelta[x] = d;
                    heap.push_back(Cand{label[x] == pq_other ? 1 : 0, -pq_dist[x], d, crnd[x], x});
                    std::push_heap(heap.begin(), heap.end(), greater);
                }
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
                if (label[v] == 0 && owner[v] != 0) { set_terminal(v, (char) owner[v]); }
            }
        }

        // Sets up the network of a job (terminals a, b pre-grown), solves the first max flow and saves the
        // network in `slot`. The flow value is a lower bound for the cut the job can end with.
        // `limit`: the first solve stops once its flow reaches it (the job is pruned then: lower >= limit)
        void prepare_job(Job &j, const int slot, const s64 limit = INF) {
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "prepare_job");
            std::memcpy(arcs.get_ptr(), base_arcs.get_ptr(), na * sizeof(Arc));
            reset_trees();
            set_terminal(j.a, 1);
            set_terminal(j.b, 2);
            if (init_fraction > 0.0) {
                const weight_t target = (weight_t) (init_fraction * (f64) total_w);
                grow_initial_terminals(j.a, j.b, std::min(target, j.hi), std::min(target, total_w - j.lo));
            }
            if (use_ibfs) {
                j.lower = ibfs_maxflow(limit);
            } else {
                bk_start();
                j.lower = bk_maxflow(limit);
            }
            if (j.lower >= limit) { return; }
            std::memcpy(slot_label.get_ptr() + (size_t) slot * N, label.get_ptr(), N);
            std::memcpy(slot_pa.get_ptr() + (size_t) slot * N, pa.get_ptr(), N * sizeof(int));
            std::memcpy(slot_arcs.get_ptr() + (size_t) slot * na, arcs.get_ptr(), na * sizeof(Arc));
        }

        // snapshot of the mirrored job (b, a) from the state of the solved job (a, b) that is still in the
        // working arrays
        void make_mirror_snapshot(const Job &j, Job &jm, const int slot) {
            jm.lower = j.lower;
            Arc *da = slot_arcs.get_ptr() + (size_t) slot * na;
            for (size_t e = 0; e < na; ++e) {
                da[e].to = arcs[e].to;
                da[e].rev = arcs[e].rev;
                da[e].cap = arcs[arcs[e].rev].cap;
            }
            char *dl = slot_label.get_ptr() + (size_t) slot * N;
            for (int v = 0; v < N; ++v) { dl[v] = label[v] ? (char) (3 - label[v]) : 0; }
            std::memcpy(slot_pa.get_ptr() + (size_t) slot * N, pa.get_ptr(), N * sizeof(int));
        }

        // Continues a job from its snapshot: pierce the lighter side until the cut is balanced.
        // The result is merged into `best` (scr_res is used as scratch).
        void run_job(const Job &j, const int slot, const int *dist_a, const int *dist_b, PairResult &best,
                     const s64 cut_bound = INF) {
            const int a = j.a, b = j.b;
            rng.seed(pair_salt ^ ((u64) a * 1000003ULL + (u64) b * 7919ULL + 17));
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "run_pair");
            std::memcpy(arcs.get_ptr(), slot_arcs.get_ptr() + (size_t) slot * na, na * sizeof(Arc));
            std::memcpy(label.get_ptr(), slot_label.get_ptr() + (size_t) slot * N, N);
            std::memcpy(pa.get_ptr(), slot_pa.get_ptr() + (size_t) slot * N, N * sizeof(int));
            std::fill_n(isact.get_ptr(), N, 0);
            act_h = act_t = 0;
            orph_n = 0;
            continue_job(j.lo, j.hi, j.lower, dist_a, dist_b, best, cut_bound);
        }

        // Continues from a solved max flow (value flow0, trees and residual network in the working arrays): pierces
        // the lighter side until the cut is balanced; the result is merged into `best`.
        void continue_job(const weight_t lo, const weight_t hi, const s64 flow0, const int *dist_a, const int *dist_b,
                          PairResult &best, const s64 cut_bound) {
            const weight_t avg_w = std::max<weight_t>(1, total_w / std::max(1, nl));
            s64 flow = flow0;
            PairResult &cur = scr_res;
            cur.reset();

            for (int iter = 0; iter < nl + 2; ++iter) {
                flow += bk_maxflow(cut_bound - flow);
                // the cut never shrinks when terminals grow: this pair cannot beat a known feasible cut
                if (flow >= cut_bound) { return; }

                weight_t w_rs = 0, w_ct = 0;
                for (int v = 0; v < nl; ++v) {
                    if (label[v] == 1) { w_rs += vw[v]; }
                    if (label[v] == 2) { w_ct += vw[v]; }
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
                        for (int x = 0; x < nl; ++x) { cur.in_source[x] = label[x] == 1; }
                    } else {
                        for (int x = 0; x < nl; ++x) { cur.in_source[x] = label[x] != 2; }
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
                        collect_candidates(1, 2, dist_b);
                        done = pierce(true, w_rs, hi, std::max(avg_w, deficit / PIERCE_DIVISOR));
                    } else {
                        const weight_t deficit = std::max<weight_t>(1, (total_w - hi) - w_ct);
                        collect_candidates(2, 1, dist_a);
                        done = pierce(false, w_ct, total_w - lo, std::max(avg_w, deficit / PIERCE_DIVISOR));
                    }
                    grow_source = !grow_source;
                }
                if (!done) { break; }
            }

            if (cur.valid) { merge(best, cur); }
        }

        // keeps the better of the two results in best (less violation, then smaller cut)
        static void merge(PairResult &best, PairResult &cur) {
            const bool better = !best.valid
                                || cur.violation < best.violation
                                || (cur.violation == best.violation && cur.cut < best.cut);
            if (better) { best.swap(cur); }
        }

        s64 bound(const PairResult &best) const {
            return (prune && best.valid && best.violation == 0) ? (s64) ((f64) best.cut * prune_factor) : INF;
        }

        // Continues a prepared job and merges its cut into best_res. Jobs whose first flow already reaches the
        // prune bound are dropped.
        void run_slot(const int slot) {
            const Job &job = jobs[slot];
            const s64 cb = bound(best_res);
            if (job.lower >= cb) { return; } // cannot beat the best feasible cut
            const int *da = dist_rows.get_ptr() + (size_t) job.ia * (size_t) nl;
            const int *db = dist_rows.get_ptr() + (size_t) job.ib * (size_t) nl;
            if (!job.mirrored) {
                run_job(job, slot, da, db, best_res, cb);
                return;
            }
            mir_res.reset();
            run_job(job, slot, da, db, mir_res, cb);
            if (mir_res.valid) {
                for (int x = 0; x < nl; ++x) { mir_res.in_source[x] = !mir_res.in_source[x]; }
                merge(best_res, mir_res);
            }
        }

        // a random one of the seed_pool vertices (not yet a seed) with the largest distance d; -1 if none is left
        int pick_far(const int *d) {
            far_top.clear();
            for (int v = 0; v < nl; ++v) {
                if (owner[v]) { continue; }
                if (far_top.size() == seed_pool && d[v] <= far_top.back().first) { continue; }
                auto it = far_top.begin();
                while (it != far_top.end() && it->first >= d[v]) { ++it; }
                far_top.insert(it, {d[v], v});
                if (far_top.size() > seed_pool) { far_top.pop_back(); }
            }
            if (far_top.empty()) { return -1; }
            return far_top.size() == 1 ? far_top[0].second : far_top[rng() % far_top.size()].second;
        }

        /// farthest point sampling on the hop distance; fills seeds and dist_rows, returns the number of seeds
        int pick_seeds(const u64 num_seeds) {
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "seed_selection");
            const int n_seeds = (int) std::min<u64>(std::max<u64>(2, num_seeds), (u64) nl);
            seeds.resize(n_seeds);
            dist_rows.initialize((size_t) n_seeds * (size_t) nl);
            int ns = 0;
            hop_distance((int) (rng() % (u64) nl), tmp_dist.get_ptr());
            std::fill_n(owner.get_ptr(), nl, 0);
            int s0 = pick_far(tmp_dist.get_ptr());
            // tmp_dist now holds the distance to the nearest seed, owner marks the seeds
            std::fill_n(tmp_dist.get_ptr(), nl, std::numeric_limits<int>::max());
            int next = s0;
            while (ns < n_seeds) {
                int *row = dist_rows.get_ptr() + (size_t) ns * (size_t) nl;
                seeds[ns++] = next;
                owner[next] = 1;
                hop_distance(next, row);
                for (int v = 0; v < nl; ++v) { tmp_dist[v] = std::min(tmp_dist[v], row[v]); }
                const int arg = pick_far(tmp_dist.get_ptr());
                if (arg < 0) { break; }
                next = arg;
            }
            return ns;
        }

        // cut of the BFS ball of weight init_fraction * total_w around v (proxy for how cheap it is to separate v)
        s64 ball_cut(const int v0) {
            const weight_t target = std::max<weight_t>(1, (weight_t) (init_fraction * (f64) total_w));
            std::vector<char> in(nl, 0);
            std::vector<int> q;
            q.push_back(v0);
            in[v0] = 1;
            weight_t w = 0;
            for (size_t h = 0; h < q.size() && w < target; ++h) {
                const int u = q[h];
                w += vw[u];
                for (int i = first[u]; i < first[u + 1]; ++i) {
                    const int x = arcs[i].to;
                    if (x < nl && !in[x]) { in[x] = 1; q.push_back(x); }
                }
            }
            // vertices queued but beyond the target weight are not part of the ball
            const size_t ball_n = [&]() { weight_t a = 0; size_t c = 0; while (c < q.size() && a < target) { a += vw[q[c]]; ++c; } return c; }();
            for (size_t h = ball_n; h < q.size(); ++h) { in[q[h]] = 0; }
            s64 cut = 0;
            for (size_t h = 0; h < ball_n; ++h) {
                const int u = q[h];
                for (int i = first[u]; i < first[u + 1]; ++i) { if (arcs[i].to < nl && !in[arcs[i].to]) { cut += base_arcs[i].cap; } }
            }
            return cut;
        }

        /// all jobs: every seed pair, and the other orientation if the sides differ (seeds[j] must end up in the
        /// first part: source side weight window mirrored). Allocates the snapshot slots.
        size_t make_jobs(const weight_t lo, const weight_t hi, const bool symmetric, const int ns) {
            const size_t max_jobs = (size_t) ns * (size_t) (ns - 1);
            jobs.resize(max_jobs);
            job_order.resize(max_jobs);
            size_t nj = 0;
            // pairs without the first seed: only the pair_fraction best ones, ranked by a cheap proxy of their cut:
            // pair_proxy 0: hop distance (farther is better); 1: sum of the cuts of the BFS balls around both seeds
            std::vector<s64> pscore;
            s64 thr = 0;
            bool has_thr = false;
            auto score = [&](const int i, const int j) -> s64 {
                if (pair_proxy == 0) { return -(s64) dist_rows[(size_t) i * (size_t) nl + (size_t) seeds[j]]; }
                return pscore[i] + pscore[j];
            };
            if (pair_fraction < 1.0 && ns > 2) {
                if (pair_proxy != 0) { pscore.resize(ns); for (int i = 0; i < ns; ++i) { pscore[i] = ball_cut(seeds[i]); } }
                std::vector<s64> ds;
                for (int i = 1; i < ns; ++i) {
                    for (int j = i + 1; j < ns; ++j) { ds.push_back(score(i, j)); }
                }
                std::sort(ds.begin(), ds.end());
                const size_t keep_pairs = std::max<size_t>(1, (size_t) std::ceil(pair_fraction * (f64) ds.size()));
                thr = ds[std::min(keep_pairs, ds.size()) - 1];
                has_thr = true;
            }
            for (int i = 0; i < ns; ++i) {
                for (int j = i + 1; j < ns; ++j) {
                    if (i > 0 && has_thr && score(i, j) > thr) { continue; }
                    jobs[nj++] = Job{seeds[i], seeds[j], lo, hi, i, j, false, 0};
                    if (!symmetric) { jobs[nj++] = Job{seeds[j], seeds[i], total_w - hi, total_w - lo, j, i, true, 0}; }
                }
            }
            slot_arcs.initialize(nj * na);
            slot_label.initialize(nj * (size_t) N);
            slot_pa.initialize(nj * (size_t) N);
            return nj;
        }

        // index in job_order[0, m) of the job with the smallest lower bound
        size_t best_ranked(const size_t m) const {
            size_t b = 0;
            for (size_t x = 1; x < m; ++x) { if (jobs[job_order[x]].lower < jobs[job_order[b]].lower) { b = x; } }
            return b;
        }

        /**
         * Solves the first max flow of the jobs; job_order[0, m) are the jobs that still have to be ranked.
         *  - The jobs of the star around seed 0 (all pairs with seed 0) are solved completely, the best of them is
         *    run at once (returned in ran_slot): its cut bounds all later first solves, which stop as soon as their
         *    flow reaches the prune bound.
         *  - Transitivity of min cuts: lambda(i,k) >= min(lambda(i,j), lambda(j,k)) for disjoint terminal sets. A
         *    pair whose bound already reaches the prune bound is not solved at all (it would be dropped anyway).
         *  - The mirrored job of a pair has the same pre-grown balls (if the ball fits both windows), and the
         *    reversed flow of a max flow is a max flow with source and sink swapped: its snapshot is derived from
         *    the first one (residual capacities exchanged with the reverse arcs, tree labels swapped, the parent
         *    arcs stay valid).
         */
        size_t prepare_jobs(const size_t nj, const int ns, const bool symmetric, const bool can_mirror,
                            int &ran_slot) {
            const size_t star = bound_jobs > 0 ? (size_t) (ns - 1) * (symmetric ? 1 : 2) : 0;
            const bool pair_sym = symmetric || can_mirror;
            lam.assign((size_t) ns * (size_t) ns, (s64) -1);
            auto lam_at = [&](const int i, const int k) -> s64 & { return lam[(size_t) i * ns + k]; };

            size_t m = 0;
            ran_slot = -1;
            for (size_t q = 0; q < nj; ++q) {
                const bool bounded = star > 0 && q >= star;
                if (star > 0 && q == star) {
                    ran_slot = job_order[best_ranked(m)];
                    run_slot(ran_slot);
                }
                Job &job = jobs[q];
                const s64 lim = bounded ? bound(best_res) : INF;
                if (bounded && pair_sym && lim != INF) {
                    s64 lb = -1;
                    for (int x = 0; x < ns; ++x) {
                        const s64 u = lam_at(job.ia, x), w = lam_at(x, job.ib);
                        if (u >= 0 && w >= 0) { lb = std::max(lb, std::min(u, w)); }
                    }
                    if (lb >= lim) {
                        job.lower = lb;
                        lam_at(job.ia, job.ib) = lam_at(job.ib, job.ia) = lb;
                        job_order[m++] = (int) q;
                        if (can_mirror && !job.mirrored && q + 1 < nj) {
                            jobs[q + 1].lower = lb;
                            job_order[m++] = (int) (q + 1);
                            ++q;
                        }
                        continue;
                    }
                }
                prepare_job(job, (int) q, lim);
                if (pair_sym) { lam_at(job.ia, job.ib) = lam_at(job.ib, job.ia) = job.lower; }
                job_order[m++] = (int) q;
                if (can_mirror && !job.mirrored && q + 1 < nj) {
                    if (job.lower < lim) {
                        make_mirror_snapshot(job, jobs[q + 1], (int) q + 1);
                    } else {
                        jobs[q + 1].lower = job.lower;
                    }
                    job_order[m++] = (int) (q + 1);
                    ++q;
                }
            }
            if (star > 0 && ran_slot < 0) { // fewer jobs than the star
                ran_slot = job_order[best_ranked(m)];
                run_slot(ran_slot);
            }
            if (ran_slot >= 0) { // already run: leave it out of the ranking
                for (size_t x = 0; x < m; ++x) { if (job_order[x] == ran_slot) { job_order[x] = job_order[--m]; break; } }
            }
            return m;
        }

        // hop distance from the vertices with mark[v] == m (several sources); order (optional) gets the visit order
        void multi_hop_distance(const char *mark, const char m, int *dist, int *order) {
            constexpr int FAR = 1 << 29;
            std::fill_n(dist, nl, FAR);
            int qh = 0, qt = 0;
            for (int v = 0; v < nl; ++v) { if (mark[v] == m) { dist[v] = 0; order[qt++] = v; } }
            while (qh < qt) {
                const int u = order[qh++];
                for (int i = first[u]; i < first[u + 1]; ++i) {
                    const int v = arcs[i].to;
                    if (v < nl && dist[v] == FAR) {
                        dist[v] = dist[u] + 1;
                        order[qt++] = v;
                    }
                }
            }
        }

        // Flow refinement of best_res (feasible): free a corridor around the cut, make the rest of both sides
        // terminal and solve again. See FlowBisectionConfiguration::refine_rounds.
        void refine_bisection(const weight_t lo, const weight_t hi) {
            if (refine_rounds == 0 || !best_res.valid || best_res.violation != 0 || best_res.cut <= 0 || nl < 4) { return; }
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "refine_bisection");
            ref_q.initialize(N);
            ref_d.initialize(N);
            ref_d1.initialize(N);
            ref_d2.initialize(N);
            ref_core.initialize(N);
            AlignedArray<char> &side = ref_core; // 1 source side, 2 sink side (before the corridor is cut out)
            // corridor widths: tried from narrow to wide; after an improvement start again with the narrowest one
            static constexpr f64 BAND_FACTORS[] = {0.2, 0.5, 1.0, 2.0, 4.0};
            constexpr size_t N_BANDS = sizeof(BAND_FACTORS) / sizeof(BAND_FACTORS[0]);
            size_t band_idx = 0;
            f64 f = refine_band * BAND_FACTORS[0];
            size_t improvements = 0;
            auto next_band = [&]() { f = ++band_idx < N_BANDS ? refine_band * BAND_FACTORS[band_idx] : 2.0; };
            while (improvements < refine_rounds && f < 0.95 && best_res.cut > 0) {
                // boundary vertices: sources of the BFS that orders the vertices by distance to the cut
                for (int v = 0; v < nl; ++v) { side[v] = best_res.in_source[v] ? 1 : 2; }
                int qt = 0;
                for (int v = 0; v < nl; ++v) {
                    for (int i = first[v]; i < first[v + 1]; ++i) {
                        if (side[arcs[i].to] != side[v]) { ref_q[qt++] = v; break; }
                    }
                }
                if (qt == 0) { return; }
                std::fill_n(ref_d.get_ptr(), nl, 1 << 29);
                for (int x = 0; x < qt; ++x) { ref_d[ref_q[x]] = 0; }
                for (int qh = 0; qh < qt; ++qh) {
                    const int u = ref_q[qh];
                    for (int i = first[u]; i < first[u + 1]; ++i) {
                        const int v = arcs[i].to;
                        if (v < nl && ref_d[v] == (1 << 29)) { ref_d[v] = ref_d[u] + 1; ref_q[qt++] = v; }
                    }
                }
                // corridor = the closest vertices up to f * total weight; terminals = the rest of each side
                const f64 corridor_w = f * (f64) total_w;
                weight_t acc = 0;
                for (int x = 0; x < qt && (f64) acc < corridor_w; ++x) {
                    acc += vw[ref_q[x]];
                    side[ref_q[x]] = 0;
                }
                int n1 = 0, n2 = 0;
                for (int v = 0; v < nl; ++v) { n1 += side[v] == 1; n2 += side[v] == 2; }
                if (n1 == 0 || n2 == 0) { next_band(); continue; }

                multi_hop_distance(side.get_ptr(), 1, ref_d1.get_ptr(), ref_q.get_ptr());
                multi_hop_distance(side.get_ptr(), 2, ref_d2.get_ptr(), ref_q.get_ptr());

                std::memcpy(arcs.get_ptr(), base_arcs.get_ptr(), na * sizeof(Arc));
                reset_trees();
                for (int v = 0; v < nl; ++v) { if (side[v]) { set_terminal(v, side[v]); } }
                bk_start();
                const s64 before = best_res.cut;
                const s64 flow = bk_maxflow(before);
                if (flow >= before) { next_band(); continue; }
                rng.seed(pair_salt ^ (0x9E3779B97F4A7C15ULL * (improvements + 1)));
                continue_job(lo, hi, flow, ref_d1.get_ptr(), ref_d2.get_ptr(), best_res, before);
                if (best_res.valid && best_res.violation == 0 && best_res.cut < before) {
                    ++improvements;
                    band_idx = 0;
                    f = refine_band * BAND_FACTORS[0];
                } else {
                    next_band();
                }
            }
        }

        // splits the sub problem; afterwards best_res.in_source[v] = 1 if local vertex v belongs to the first part
        void bisect(const weight_t lo, const weight_t hi, const bool symmetric, const u64 num_seeds) {
            best_res.reset();
            best_res.in_source.initialize(nl);
            std::fill_n(best_res.in_source.get_ptr(), nl, 1);
            if (nl < 2) { return; }

            u64 seeds_cap = num_seeds;
            const bool small = small_n > 0 && (size_t) nl < small_n;
            const f64 keep_saved = keep_fraction, pf_saved = prune_factor;
            if (small) {
                seeds_cap = std::max<u64>(num_seeds, small_seeds);
                if (small_exact) {
                    keep_fraction = 1.0;
                    prune_factor = 1.0;
                }
            } else if (seed_budget > 0) {
                seeds_cap = std::min<u64>(num_seeds, std::max<u64>(2, 1 + seed_budget / std::max<size_t>(1, na)));
            }
            const int ns = pick_seeds(seeds_cap);
            pair_salt = rng();
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "all_pairs");

            const size_t nj = make_jobs(lo, hi, symmetric, ns);
            const weight_t ball = (weight_t) (init_fraction * (f64) total_w);
            const bool can_mirror = !symmetric && ball <= std::min(hi, total_w - lo);
            int ran_slot = -1;
            const size_t m = prepare_jobs(nj, ns, symmetric, can_mirror, ran_slot);

            // most promising (smallest lower bound) first: finds a good feasible cut early, which then prunes the
            // remaining jobs. Only the best keep_fraction of all jobs is continued.
            std::stable_sort(job_order.begin(), job_order.begin() + m,
                             [&](const int x, const int y) { return jobs[x].lower < jobs[y].lower; });
            const size_t run_jobs = std::max<size_t>(1, (size_t) std::ceil((double) nj * keep_fraction));
            const size_t to_run = ran_slot >= 0 ? (run_jobs > 1 ? run_jobs - 1 : 0) : run_jobs;
            for (size_t r = 0; r < std::min(to_run, m); ++r) { run_slot(job_order[r]); }

            keep_fraction = keep_saved;
            prune_factor = pf_saved;
            refine_bisection(lo, hi);

            // the number of random draws inside the pairs depends on pruning: restore a fixed stream
            rng.seed(pair_salt + 1);

        }

        /**
         * Computes a k'-partition of the vertices order[begin, end) (induced subgraph of g) by recursive
         * flow based bisection, in place: afterwards block i owns order[blk_b[i], blk_e[i]) (the blocks
         * appear in increasing order). Block i may weigh at most max_weights[i]. `loc` is -1 outside
         * of this function.
         */
        void kway_partition(const CSRGraph &g, vertex_t *order, const size_t begin, const size_t end,
                            const partition_t k_prime, const weight_t *max_weights,
                            const u64 num_seeds) {
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "kway_partition");
            blk_b.assign(k_prime, begin);
            blk_e.assign(k_prime, begin);
            prefix.resize(k_prime + 1);
            prefix[0] = 0;
            for (partition_t i = 0; i < k_prime; ++i) { prefix[i + 1] = prefix[i] + max_weights[i]; }
            scratch_order.initialize(end - begin);

            task_stack.clear();
            task_stack.push_back({begin, end, 0, k_prime});

            while (!task_stack.empty()) {
                const Task task = task_stack.back();
                task_stack.pop_back();
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
                setup_network(g, verts, (int) cnt);
                for (size_t i = 0; i < cnt; ++i) { loc[verts[i]] = -1; }

                const weight_t lo_w = std::max<weight_t>(0, total_w - max2);
                bisect(lo_w, max1, max1 == max2, num_seeds);

                // stable partition of the range: first part to the front
                size_t n1 = 0, n2 = 0;
                for (size_t i = 0; i < cnt; ++i) {
                    if (best_res.in_source[i]) { verts[n1++] = verts[i]; } else { scratch_order[n2++] = verts[i]; }
                }
                std::memcpy(verts + n1, scratch_order.get_ptr(), n2 * sizeof(vertex_t));

                task_stack.push_back({task.b, task.b + n1, task.lo, mid});
                task_stack.push_back({task.b + n1, task.e, mid, task.hi});
            }
            if (repair && k_prime > 1) { repair_kway(g, order, begin, end, k_prime, max_weights); }
            if (kway_rounds > 0 && k_prime > 2) { refine_kway(g, order, begin, end, k_prime, max_weights); }
        }

        // Repairs overloaded blocks of the finished k'-partition (a bisection could not hit its weight window):
        // repeatedly moves the vertex with the smallest cut increase from an overloaded block to a neighbouring block
        // where it fits (or at least lowers the overload). Returns true if anything moved.
        bool repair_kway(const CSRGraph &g, vertex_t *order, const size_t begin, const size_t end,
                         const partition_t k_prime, const weight_t *max_weights) {
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "repair_kway");
            std::vector<weight_t> bw(k_prime, 0);
            bool any_over = false;
            for (partition_t b = 0; b < k_prime; ++b) {
                for (size_t x = blk_b[b]; x < blk_e[b]; ++x) { bw[b] += g.v_weight(order[x]); }
                any_over |= bw[b] > max_weights[b];
            }
            if (!any_over) { return false; }
            if (blk_of.size() < (size_t) g.n) { blk_of.assign(g.n, -1); }
            for (partition_t b = 0; b < k_prime; ++b) {
                for (size_t x = blk_b[b]; x < blk_e[b]; ++x) { blk_of[order[x]] = (int) b; }
            }
            std::vector<s64> conn(k_prime, 0);
            std::vector<partition_t> touched;
            bool moved_any = false;
            while (true) {
                // most overloaded block
                partition_t ob = k_prime;
                weight_t worst = 0;
                for (partition_t b = 0; b < k_prime; ++b) {
                    if (bw[b] - max_weights[b] > worst) { worst = bw[b] - max_weights[b]; ob = b; }
                }
                if (ob == k_prime) { break; }
                // best move out of ob: (fits, cut increase, larger weight)
                vertex_t best_u = 0;
                int best_c = -1;
                bool best_fits = false;
                s64 best_delta = 0;
                weight_t best_w = 0;
                for (size_t x = begin; x < end; ++x) {
                    const vertex_t u = order[x];
                    if (blk_of[u] != (int) ob) { continue; }
                    const weight_t wu = g.v_weight(u);
                    touched.clear();
                    for (size_t e = g.neighborhoods[u]; e < g.neighborhoods[u + 1]; ++e) {
                        const int c = blk_of[g.edges_v[e]];
                        if (c < 0) { continue; }
                        if (conn[c] == 0) { touched.push_back((partition_t) c); }
                        conn[c] += g.edge_w(e);
                    }
                    for (const partition_t c: touched) {
                        if (c != ob) {
                            const weight_t after = bw[c] + wu;
                            const bool fits = after <= max_weights[c];
                            // otherwise it must at least lower the overload of the pair
                            const bool helps = std::max<weight_t>(0, after - max_weights[c]) < std::max<weight_t>(0, bw[ob] - max_weights[ob]) && wu <= bw[ob] - max_weights[ob] + (max_weights[c] - std::min(bw[c], max_weights[c]));
                            if (fits || helps) {
                                const s64 delta = conn[ob] - conn[c];
                                if (best_c < 0 || (fits && !best_fits) || (fits == best_fits && (delta < best_delta || (delta == best_delta && wu > best_w)))) {
                                    best_u = u; best_c = (int) c; best_fits = fits; best_delta = delta; best_w = wu;
                                }
                            }
                        }
                    }
                    for (const partition_t c: touched) { conn[c] = 0; }
                }
                if (best_c < 0) { break; }
                bw[ob] -= best_w;
                bw[best_c] += best_w;
                blk_of[best_u] = best_c;
                moved_any = true;
                // the order array is rebuilt at the end; keep the scan consistent by tracking blk_of only
            }
            if (moved_any) {
                std::vector<size_t> pos(k_prime + 1, 0);
                for (size_t x = begin; x < end; ++x) { ++pos[blk_of[order[x]] + 1]; }
                for (partition_t b = 0; b < k_prime; ++b) { pos[b + 1] += pos[b]; }
                scratch_order.initialize(end - begin);
                std::vector<size_t> fill(pos.begin(), pos.end() - 1);
                for (size_t x = begin; x < end; ++x) { const vertex_t u = order[x]; scratch_order[fill[blk_of[u]]++] = u; }
                for (size_t x = begin; x < end; ++x) { order[x] = scratch_order[x - begin]; }
                for (partition_t b = 0; b < k_prime; ++b) { blk_b[b] = begin + pos[b]; blk_e[b] = begin + pos[b + 1]; }
            }
            for (size_t x = begin; x < end; ++x) { blk_of[order[x]] = -1; }
            return moved_any;
        }

        // Pairwise flow refinement of the finished k'-partition (blocks = ranges blk_b/blk_e of order). Every pair of
        // adjacent blocks is refined as a bisection of its union (the edges to other blocks are cut either way).
        void refine_kway(const CSRGraph &g, vertex_t *order, const size_t begin, const size_t end,
                         const partition_t k_prime, const weight_t *max_weights) {
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "refine_kway");
            if (blk_of.size() < (size_t) g.n) { blk_of.assign(g.n, -1); }
            for (partition_t b = 0; b < k_prime; ++b) {
                for (size_t x = blk_b[b]; x < blk_e[b]; ++x) { blk_of[order[x]] = (int) b; }
            }
            std::vector<std::vector<vertex_t> > members(k_prime);
            for (partition_t b = 0; b < k_prime; ++b) { members[b].assign(order + blk_b[b], order + blk_e[b]); }
            std::vector<weight_t> shared;
            std::vector<vertex_t> verts;
            // a pair is only refined again if one of its blocks changed since its last refinement
            std::vector<u64> ver(k_prime, 1), seen((size_t) k_prime * k_prime, 0);
            bool changed_any = false;
            for (size_t round = 0; round < kway_rounds; ++round) {
                shared.assign((size_t) k_prime * k_prime, 0);
                for (size_t x = begin; x < end; ++x) {
                    const vertex_t u = order[x];
                    for (size_t e = g.neighborhoods[u]; e < g.neighborhoods[u + 1]; ++e) {
                        const vertex_t v = g.edges_v[e];
                        const int bu = blk_of[u], bv = blk_of[v];
                        if (bv >= 0 && bu < bv) { shared[(size_t) bu * k_prime + bv] += g.edge_w(e); }
                    }
                }
                std::vector<std::pair<weight_t, std::pair<int, int> > > pairs;
                for (partition_t a = 0; a < k_prime; ++a) {
                    for (partition_t b = a + 1; b < k_prime; ++b) {
                        if (shared[(size_t) a * k_prime + b] > 0) { pairs.push_back({shared[(size_t) a * k_prime + b], {(int) a, (int) b}}); }
                    }
                }
                std::sort(pairs.begin(), pairs.end(), [](const auto &x, const auto &y) { return x.first > y.first; });
                bool improved = false;
                for (const auto &pr: pairs) {
                    const int a = pr.second.first, b = pr.second.second;
                    const u64 stamp = ver[a] * 1000003ULL + ver[b];
                    if (seen[(size_t) a * k_prime + b] == stamp) { continue; }
                    seen[(size_t) a * k_prime + b] = stamp;
                    verts.clear();
                    verts.insert(verts.end(), members[a].begin(), members[a].end());
                    verts.insert(verts.end(), members[b].begin(), members[b].end());
                    const int cnt = (int) verts.size();
                    if (cnt < 4) { continue; }
                    for (int i = 0; i < cnt; ++i) { loc[verts[i]] = i; }
                    setup_network(g, verts.data(), cnt);
                    for (int i = 0; i < cnt; ++i) { loc[verts[i]] = -1; }
                    best_res.reset();
                    best_res.in_source.initialize(nl);
                    weight_t w1 = 0;
                    s64 cut = 0;
                    for (int v = 0; v < nl; ++v) { best_res.in_source[v] = v < (int) members[a].size() ? 1 : 0; }
                    for (int v = 0; v < nl; ++v) {
                        if (best_res.in_source[v]) { w1 += vw[v]; }
                        for (int i = first[v]; i < first[v + 1]; ++i) {
                            if (arcs[i].to < nl && best_res.in_source[arcs[i].to] != best_res.in_source[v] && best_res.in_source[v]) {
                                cut += base_arcs[i].cap;
                            }
                        }
                    }
                    const weight_t lo_w = std::max<weight_t>(0, total_w - max_weights[b]);
                    const weight_t hi_w = max_weights[a];
                    if (w1 < lo_w || w1 > hi_w || cut <= 0) { continue; }
                    best_res.valid = true;
                    best_res.violation = 0;
                    best_res.cut = cut;
                    pair_salt = rng();
                    refine_bisection(lo_w, hi_w);
                    if (!(best_res.valid && best_res.violation == 0 && best_res.cut < cut)) { continue; }
                    improved = changed_any = true;
                    ++ver[a];
                    ++ver[b];
                    members[a].clear();
                    members[b].clear();
                    for (int v = 0; v < nl; ++v) {
                        const int blk = best_res.in_source[v] ? a : b;
                        members[blk].push_back(verts[v]);
                        blk_of[verts[v]] = blk;
                    }
                }
                if (!improved) { break; }
            }
            if (changed_any) {
                size_t pos = begin;
                for (partition_t b = 0; b < k_prime; ++b) {
                    blk_b[b] = pos;
                    for (const vertex_t u: members[b]) { order[pos++] = u; }
                    blk_e[b] = pos;
                }
            }
            for (partition_t b = 0; b < k_prime; ++b) { for (const vertex_t u: members[b]) { blk_of[u] = -1; } }
        }

    public:
        /// copies the parameters of the configuration (the partitioner is reusable and not thread safe)
        void configure(const FlowBisectionConfiguration &config) {
            prune = config.prune;
            init_fraction = config.init_fraction;
            prune_factor = config.prune_factor;
            keep_fraction = config.keep_fraction;
            bound_jobs = (size_t) config.bound_jobs;
            seed_budget = (size_t) config.seed_budget;
            seed_pool = (size_t) std::max<u64>(1, config.seed_pool);
            refine_rounds = (size_t) config.refine_rounds;
            refine_band = config.refine_band;
            kway_rounds = (size_t) config.kway_rounds;
            pair_fraction = config.pair_fraction;
            pair_proxy = config.pair_proxy;
            repair = config.repair;
            use_ibfs = config.ibfs;
            small_n = (size_t) config.small_n;
            small_exact = config.small_exact;
            small_seeds = (size_t) config.small_seeds;
        }

        /// deterministic random stream of a part
        void reseed(const u64 seed, const u64 part) {
            rng.seed((seed + 1) * 0x9E3779B97F4A7C15ULL ^ (part + 1) * 0xC2B2AE3D27D4EB4FULL);
        }

        /// loc must have g.n entries (all -1)
        void init_loc(const vertex_t n) { loc.initialize(n, -1); }

        /// result of the last split_item: block i owns order[block_begin(i), block_end(i))
        size_t block_begin(const partition_t i) const { return blk_b[i]; }
        size_t block_end(const partition_t i) const { return blk_e[i]; }

        /**
         * Splits the vertices order[b, e) (one part of hierarchy level `level`) into hierarchy[level] parts, each
         * of which may weigh at most (#PEs below it) * lmax (or less with the adaptive imbalance). `stride[j]` is
         * the number of PEs below one part of level j, k the total number of PEs.
         */
        void split_item(const CSRGraph &g, vertex_t *order, const size_t b, const size_t e, const size_t level,
                        const std::vector<partition_t> &hierarchy, const partition_t *stride, const partition_t k,
                        const weight_t lmax, const f64 imbalance, const FlowBisectionConfiguration &config,
                        const u64 top_budget) {
            const partition_t kk = hierarchy[level];
            weight_t part_max = (weight_t) stride[level] * lmax;
            if (config.adaptive_imbalance) {
                // imbalance of this item such that the remaining (level + 1) levels can share the global slack
                weight_t item_w = 0;
                for (size_t x = b; x < e; ++x) { item_w += g.v_weight(order[x]); }
                const f64 local_k = (f64) stride[level] * (f64) kk; // PEs below this item
                f64 imb = ((1.0 + imbalance) * local_k * (f64) g.g_weight) / ((f64) k * (f64) std::max<weight_t>(1, item_w));
                imb = std::max(0.0, std::pow(imb, 1.0 / (f64) (level + 1)) - 1.0);
                const weight_t adaptive = (weight_t) std::ceil((1.0 + imb) * (f64) item_w / (f64) kk);
                part_max = std::min(part_max, adaptive);
            }
            mw.assign(kk, part_max);
            const size_t top = hierarchy.size() - 1;
            const u64 ns = level == top ? config.num_seeds : (level == 0 ? config.num_seeds_bottom : config.num_seeds_mid);
            seed_budget = (size_t) (level == top ? top_budget : config.seed_budget);
            const size_t kr = kway_rounds, ss = small_seeds, rr = refine_rounds;
            if (level != top) {
                if (config.low_kway_rounds >= 0) { kway_rounds = (size_t) config.low_kway_rounds; }
                if (config.low_small_seeds >= 0) { small_seeds = (size_t) config.low_small_seeds; }
                if (config.low_refine_rounds >= 0) { refine_rounds = (size_t) config.low_refine_rounds; }
            }
            kway_partition(g, order, b, e, kk, mw.data(), ns);
            kway_rounds = kr;
            small_seeds = ss;
            refine_rounds = rr;
        }

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
            p.configure(config);
            p.rng.seed(seed * 0x9E3779B97F4A7C15ULL + 12345);

            const size_t l = hierarchy.size();
            std::vector<partition_t> stride(l); // PEs below one part of level j
            partition_t k = 1;
            for (size_t j = 0; j < l; ++j) {
                stride[j] = k;
                k *= hierarchy[j];
            }
            const weight_t lmax = (weight_t) std::ceil((1.0 + imbalance) * (f64) g.g_weight / (f64) k);

            p.init_loc(g.n);
            AlignedArray<vertex_t> order;
            order.initialize(g.n);
            for (vertex_t u = 0; u < g.n; ++u) { order[u] = u; }

            struct Item {
                size_t b, e; // range in order
                size_t level; // index into hierarchy, split into hierarchy[level] parts
                partition_t offset;
            };
            std::vector<Item> stack;
            stack.push_back({0, (size_t) g.n, l - 1, 0});

            while (!stack.empty()) {
                const Item item = stack.back();
                stack.pop_back();
                p.split_item(g, order.get_ptr(), item.b, item.e, item.level, hierarchy, stride.data(), k, lmax,
                             imbalance, config, config.top_seed_budget);
                const partition_t kk = hierarchy[item.level];
                for (partition_t i = 0; i < kk; ++i) {
                    if (item.level == 0) {
                        for (size_t x = p.blk_b[i]; x < p.blk_e[i]; ++x) {
                            const vertex_t u = order[x];
                            p_manager.set(u, g.v_weight(u), item.offset + i);
                        }
                    } else {
                        stack.push_back({p.blk_b[i], p.blk_e[i], item.level - 1, item.offset + i * stride[item.level]});
                    }
                }
            }
        }
    };
} // namespace HeiProMap

#endif // HEIPROMAP_FLOW_BISECTION_H
