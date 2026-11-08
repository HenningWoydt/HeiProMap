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
        u64 num_seeds_mid = 3;    // levels between top and bottom
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
     * The flow network and the max flow are in FlowNetwork. All work buffers are AlignedArray members that only
     * ever grow, so after the first (largest) sub problem nothing is allocated or freed anymore.
     */
    class FlowBisectionPartitioner : FlowNetwork {
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

        AlignedArray<int> queue, stamp;
        int stamp_id = 0;

        // scratch
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
        AlignedArray<char> slot_label;
        AlignedArray<int> slot_pa;
        AlignedArray<Job> jobs;
        AlignedArray<int> job_order;
        AlignedArray<s64> lam; // lam[i * ns + k]: known lower bound of the min cut of seeds i and k

        std::mt19937_64 rng;
        u64 pair_salt = 0;
        bool prune = true;
        f64 init_fraction = 0.1;
        f64 prune_factor = 1.0;
        f64 keep_fraction = 1.0;
        size_t bound_jobs = 0;

        /// network of the vertices verts[0, cnt) (loc must be set for them) and the buffers of the partitioner
        void setup_network(const CSRGraph &g, const vertex_t *verts, const int cnt) {
            build_network(g, verts, cnt, loc.get_ptr());
            queue.initialize(N);
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
        void collect_candidates(const char side, const char other_side, const int *dist_other) {
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "collect_candidates");
            cand_n = 0;
            ++stamp_id;
            for (int v = 0; v < nl; ++v) {
                if (label[v] != side) { continue; }
                for (int i = first[v]; i < first[v + 1]; ++i) {
                    const int w = arcs[i].to;
                    if (w >= nl || label[w] == side || is_terminal(w) || stamp[w] == stamp_id) { continue; }
                    stamp[w] = stamp_id;
                    cand[cand_n++] = Cand{label[w] == other_side ? 1 : 0, -dist_other[w], rng(), w};
                }
            }
            if (cand_n == 0) {
                // disconnected: any free vertex
                for (int w = 0; w < nl; ++w) {
                    if (label[w] == side || is_terminal(w)) { continue; }
                    cand[cand_n++] = Cand{label[w] == other_side ? 1 : 0, -dist_other[w], rng(), w};
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
                    bk_add_source(w);
                } else {
                    bk_add_sink(w);
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
            bk_start();
            j.lower = bk_maxflow(limit);
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
            const weight_t lo = j.lo, hi = j.hi;
            rng.seed(pair_salt ^ ((u64) a * 1000003ULL + (u64) b * 7919ULL + 17));
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "run_pair");
            std::memcpy(arcs.get_ptr(), slot_arcs.get_ptr() + (size_t) slot * na, na * sizeof(Arc));
            std::memcpy(label.get_ptr(), slot_label.get_ptr() + (size_t) slot * N, N);
            std::memcpy(pa.get_ptr(), slot_pa.get_ptr() + (size_t) slot * N, N * sizeof(int));
            std::fill_n(isact.get_ptr(), N, 0);
            act_h = act_t = 0;
            orph_n = 0;

            const weight_t avg_w = std::max<weight_t>(1, total_w / std::max(1, nl));
            s64 flow = j.lower;
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
                        done = pierce(true, w_rs, hi, std::max(avg_w, deficit / 8));
                    } else {
                        const weight_t deficit = std::max<weight_t>(1, (total_w - hi) - w_ct);
                        collect_candidates(2, 1, dist_a);
                        done = pierce(false, w_ct, total_w - lo, std::max(avg_w, deficit / 8));
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

        /// farthest point sampling on the hop distance; fills seeds and dist_rows, returns the number of seeds
        int pick_seeds(const u64 num_seeds) {
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
            return ns;
        }

        /// all jobs: every seed pair, and the other orientation if the sides differ (seeds[j] must end up in the
        /// first part: source side weight window mirrored). Allocates the snapshot slots.
        size_t make_jobs(const weight_t lo, const weight_t hi, const bool symmetric, const int ns) {
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
            lam.initialize((size_t) ns * (size_t) ns);
            std::fill_n(lam.get_ptr(), (size_t) ns * (size_t) ns, (s64) -1);
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

        // splits the sub problem; first_part[v] = 1 if local vertex v belongs to the first part
        void bisect(const weight_t lo, const weight_t hi, const bool symmetric, const u64 num_seeds) {
            std::fill_n(first_part.get_ptr(), nl, 1);
            if (nl < 2) { return; }

            const int ns = pick_seeds(num_seeds);
            best_res.reset();
            pair_salt = rng();
            HEIPROMAP_PROFILE_SCOPE("partition", "FlowBisection", "all_pairs");

            const size_t nj = make_jobs(lo, hi, symmetric, ns);
            const weight_t ball = (weight_t) (init_fraction * (f64) total_w);
            const bool can_mirror = !symmetric && ball <= std::min(hi, total_w - lo);
            int ran_slot = -1;
            const size_t m = prepare_jobs(nj, ns, symmetric, can_mirror, ran_slot);

            // most promising (smallest lower bound) first: finds a good feasible cut early, which then prunes the
            // remaining jobs. Only the best keep_fraction of all jobs is continued.
            std::stable_sort(job_order.get_ptr(), job_order.get_ptr() + m,
                             [&](const int x, const int y) { return jobs[x].lower < jobs[y].lower; });
            const size_t run_jobs = std::max<size_t>(1, (size_t) std::ceil((double) nj * keep_fraction));
            const size_t to_run = ran_slot >= 0 ? (run_jobs > 1 ? run_jobs - 1 : 0) : run_jobs;
            for (size_t r = 0; r < std::min(to_run, m); ++r) { run_slot(job_order[r]); }

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
                            const u64 num_seeds) {
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
                setup_network(g, verts, (int) cnt);
                for (size_t i = 0; i < cnt; ++i) { loc[verts[i]] = -1; }

                const weight_t lo_w = std::max<weight_t>(0, total_w - max2);
                bisect(lo_w, max1, max1 == max2, num_seeds);

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
            p.prune = config.prune;
            p.init_fraction = config.init_fraction;
            p.prune_factor = config.prune_factor;
            p.keep_fraction = config.keep_fraction;
            p.bound_jobs = (size_t) config.bound_jobs;
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
                const size_t top = hierarchy.size() - 1;
                const u64 ns = item.level == top ? config.num_seeds
                               : (item.level == 0 ? config.num_seeds_bottom : config.num_seeds_mid);
                p.kway_partition(g, order.get_ptr(), item.b, item.e, kk, p.mw.get_ptr(), ns);

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
