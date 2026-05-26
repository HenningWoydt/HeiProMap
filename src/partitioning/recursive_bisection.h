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

#ifndef HEIPROMAP_RECURSIVE_BISECTION_H
#define HEIPROMAP_RECURSIVE_BISECTION_H

#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <vector>

#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../definitions.h"
#include "../utility/random_engine.h"

namespace HeiProMap {

    /**
     * Partitions a graph into k blocks using recursive bisection.
     *
     * Each bisection step uses greedy interleaved BFS region growing:
     *   1. Two seeds are chosen via furthest-point seeding (the pair of vertices
     *      that are maximally far apart inside the current vertex subset).
     *   2. Two BFS frontiers grow simultaneously from the seeds. At each step
     *      the lighter frontier is expanded, naturally producing balanced halves.
     *   3. Any vertices unreachable from both seeds (disconnected components)
     *      are greedily assigned to whichever side is lighter.
     *
     * The public `partition()` method supports `kappa` independent trials with
     * different random seeds, keeping the result with the lowest edge-cut.
     */
    class RecursiveBisectionPartitioner {
    private:
        // ---------------------------------------------------------------
        // Working arrays (allocated once, reused across recursive calls)
        // ---------------------------------------------------------------

        /**
         * Side marker for every vertex in the graph.
         *   0  = assigned to left region
         *   1  = assigned to right region
         *   2  = unassigned, belongs to the current bisection subset
         *   3  = not part of the current bisection subset (sibling subproblem)
         *
         * The invariant maintained by recurse() is:
         *   - on entry  : side[v] == 2  for every v in `vertices`
         *                 side[v] == 3  for every other vertex
         *   - on exit   : contents inside `vertices` are undefined
         *                 (parent restores them as needed)
         */
        std::vector<u8> side;

        /**
         * BFS distance array, always reset to -1 after each bfs_furthest() call.
         * Uses -1 as the "unvisited" sentinel so we never need a full reset.
         */
        std::vector<s32> dist;

        RandomEngine *rng = nullptr;

        // ---------------------------------------------------------------
        // Seed selection
        // ---------------------------------------------------------------

        /**
         * BFS from `start` within the current subset (side[v] != 3).
         * Returns the vertex with the greatest BFS depth — the "furthest point"
         * in the subset reachable from `start`.
         *
         * The dist[] array is reset to -1 for all visited vertices before returning,
         * so no global reset is needed between calls.
         */
        vertex_t bfs_furthest(const CSRGraph &g, vertex_t start) {
            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "bfs_furthest");

            std::vector<vertex_t> visited;
            visited.reserve(64);

            std::queue<vertex_t> q;
            dist[start] = 0;
            visited.push_back(start);
            q.push(start);

            vertex_t furthest = start;

            while (!q.empty()) {
                vertex_t u = q.front();
                q.pop();

                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    vertex_t v = g.edges_v[i];

                    // Only traverse within the current subset
                    if (dist[v] == -1 && side[v] != 3) {
                        dist[v] = dist[u] + 1;
                        visited.push_back(v);
                        q.push(v);

                        if (dist[v] > dist[furthest]) {
                            furthest = v;
                        }
                    }
                }
            }

            // Reset dist for every visited vertex (O(visited) instead of O(g.n))
            for (vertex_t v : visited) {
                dist[v] = -1;
            }

            return furthest;
        }

        // ---------------------------------------------------------------
        // Bisection
        // ---------------------------------------------------------------

        /**
         * Greedy interleaved BFS bisection.
         *
         * Precondition : side[v] == 2 for all v in `vertices`
         *                side[v] == 3 for all other vertices
         * Postcondition: side[v] == 0 (left) or 1 (right) for all v in `vertices`
         *
         * @param g               The graph.
         * @param vertices        Vertices belonging to this bisection subproblem.
         * @param target_weight_left  Target total vertex weight for the left side.
         * @param left_out        Output: vertices assigned to the left side.
         * @param right_out       Output: vertices assigned to the right side.
         */
        void bisect(const CSRGraph &g,
                    const std::vector<vertex_t> &vertices,
                    weight_t target_weight_left,
                    std::vector<vertex_t> &left_out,
                    std::vector<vertex_t> &right_out) {
            // ---- Furthest-point seed selection ----
            vertex_t seed_0;
            {
                HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "bisect_overhead");
                seed_0 = vertices[rng->get_u64() % vertices.size()];
            }

            vertex_t seed_1 = bfs_furthest(g, seed_0);

            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "bisect_grow");
            // Fallback: if the graph is a single isolated vertex or all vertices
            // are in one connected component of size 1, seed_1 == seed_0.
            // Pick the last vertex in the list as an alternative.
            if (seed_1 == seed_0 && vertices.size() > 1) {
                seed_1 = (vertices.front() != seed_0) ? vertices.front() : vertices.back();
            }

            // Total weight of this subset (precomputed by caller, but we need it
            // here to derive the right-side target)
            weight_t total_weight   = 0;
            for (vertex_t v : vertices) { total_weight += g.v_weights[v]; }
            weight_t target_weight_right = total_weight - target_weight_left;

            // ---- Interleaved BFS growing ----
            side[seed_0] = 0;
            side[seed_1] = 1;
            weight_t w0 = g.v_weights[seed_0];
            weight_t w1 = g.v_weights[seed_1];

            // Handle the degenerate case where both seeds are the same vertex
            if (seed_0 == seed_1) {
                side[seed_0] = 0;
                w0 = g.v_weights[seed_0];
                w1 = 0;
            }

            std::queue<vertex_t> q0, q1;
            q0.push(seed_0);
            if (seed_1 != seed_0) { q1.push(seed_1); }

            while (!q0.empty() || !q1.empty()) {
                // Determine which side can still grow
                const bool can_grow_0 = !q0.empty() && w0 < target_weight_left;
                const bool can_grow_1 = !q1.empty() && w1 < target_weight_right;

                if (!can_grow_0 && !can_grow_1) { break; }

                // Expand the lighter eligible side; ties broken towards left
                const bool grow_left = (can_grow_0 && can_grow_1)
                                       ? (w0 <= w1)
                                       : can_grow_0;

                std::queue<vertex_t> &q = grow_left ? q0 : q1;
                const u8              s = grow_left ? u8(0) : u8(1);
                weight_t             &w = grow_left ? w0 : w1;

                vertex_t u = q.front();
                q.pop();

                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    vertex_t v = g.edges_v[i];

                    if (side[v] == 2) { // unassigned vertex inside this subset
                        side[v] = s;
                        w += g.v_weights[v];
                        q.push(v);
                    }
                }
            }

            // ---- Assign remaining unassigned vertices (disconnected components) ----
            // Any vertex that neither BFS frontier reached goes to the lighter side.
            for (vertex_t v : vertices) {
                if (side[v] == 2) {
                    const u8 s = (w0 <= w1) ? u8(0) : u8(1);
                    side[v]    = s;
                    if (s == 0) { w0 += g.v_weights[v]; }
                    else        { w1 += g.v_weights[v]; }
                }
            }

            // ---- Collect results ----
            left_out.clear();
            right_out.clear();
            left_out.reserve(vertices.size() / 2 + 1);
            right_out.reserve(vertices.size() / 2 + 1);

            for (vertex_t v : vertices) {
                if (side[v] == 0) { left_out.push_back(v); }
                else              { right_out.push_back(v); }
            }
        }

        // ---------------------------------------------------------------
        // Recursion
        // ---------------------------------------------------------------

        /**
         * Recursively partitions `vertices` into `k` blocks starting at `block_start`.
         *
         * @param g             The graph.
         * @param pm            PartitionManager to write assignments into (via set()).
         * @param vertices      Vertices in this subproblem (modified during recursion).
         * @param block_start   First block ID allocated to this subproblem.
         * @param k             Number of blocks to partition into.
         * @param total_weight  Sum of vertex weights in `vertices`.
         */
        void recurse(const CSRGraph &g,
                     PartitionManager &pm,
                     std::vector<vertex_t> &vertices,
                     partition_t block_start,
                     partition_t k,
                     weight_t total_weight) {
            // Base case: assign all vertices in this subset to block_start
            if (k == 1) {
                HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "recurse_base_case");
                for (vertex_t v : vertices) {
                    pm.set(v, g.v_weights[v], block_start);
                }
                return;
            }

            partition_t k_left, k_right;
            weight_t target_left;
            {
                HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "recurse_overhead");
                // Split k into two halves; left gets the larger slice when k is odd
                k_left  = k / 2;
                k_right = k - k_left;

                // Target weight for the left half: proportional to its share of blocks
                target_left = (weight_t) std::llround((f64) total_weight * (f64) k_left / (f64) k);
            }

            // Bisect the current vertex set
            std::vector<vertex_t> left_v, right_v;
            bisect(g, vertices, target_left, left_v, right_v);

            weight_t w_left = 0, w_right = 0;
            {
                HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "recurse_overhead");
                // Compute actual weights of each half for the recursive calls
                for (vertex_t v : left_v)  { w_left  += g.v_weights[v]; }
                for (vertex_t v : right_v) { w_right += g.v_weights[v]; }

                // ---- Recurse left ----
                // Hide the right subset from the left recursion's BFS calls
                for (vertex_t v : right_v) { side[v] = 3; }
                for (vertex_t v : left_v)  { side[v] = 2; } // reset to unassigned
            }

            recurse(g, pm, left_v, block_start, k_left, w_left);

            {
                HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "recurse_overhead");
                // ---- Recurse right ----
                // Hide the (now processed) left subset from the right recursion
                for (vertex_t v : left_v)  { side[v] = 3; }
                for (vertex_t v : right_v) { side[v] = 2; } // reset to unassigned
            }

            recurse(g, pm, right_v, block_start + k_left, k_right, w_right);
        }

        // ---------------------------------------------------------------
        // Edge-cut evaluation (for multi-trial selection)
        // ---------------------------------------------------------------

        /**
         * Computes the total edge-cut of a partition (each cut edge counted once).
         * O(n + m).
         */
        static weight_t compute_edge_cut(const CSRGraph &g, const PartitionManager &pm) {
            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "compute_edge_cut");

            weight_t cut = 0;
            for (vertex_t u = 0; u < g.n; ++u) {
                const partition_t u_id = pm[u];
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    const vertex_t v = g.edges_v[i];
                    if (pm[v] != u_id) {
                        cut += g.edges_w[i];
                    }
                }
            }
            return cut / 2; // each edge is counted from both endpoints
        }

    public:
        // ---------------------------------------------------------------
        // Public API
        // ---------------------------------------------------------------

        /**
         * Partition `g` into `k` balanced blocks using recursive bisection with
         * greedy BFS region growing. The result is written into `out_pm`.
         *
         * If `kappa > 1`, runs `kappa` independent trials with different random
         * seeds derived from `seed`, and keeps the trial with the lowest edge-cut.
         *
         * @param g        Input graph.
         * @param out_pm   Output PartitionManager. Must already be initialized
         *                 with `out_pm.initialize(g.n, k, 0)` before calling.
         * @param k        Number of blocks.
         * @param seed     Base random seed.
         * @param kappa    Number of independent trials (default 1).
         */
        void partition(const CSRGraph &g,
                       PartitionManager &out_pm,
                       partition_t k,
                       u64 seed,
                       u64 kappa = 1) {
            RandomEngine rand_engine(seed);
            rng = &rand_engine;

            std::vector<vertex_t> all_vertices(g.n);
            weight_t best_cut = std::numeric_limits<weight_t>::max();

            // Reusable PartitionManager for each trial
            PartitionManager trial_pm;

            {
                HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "allocate");
                // Allocate/resize working arrays
                side.assign(g.n, u8(2));
                dist.assign(g.n, s32(-1));
                std::iota(all_vertices.begin(), all_vertices.end(), vertex_t(0));
                trial_pm.initialize(g.n, k, 0); // start with all block weights at 0
            }

            for (u64 trial = 0; trial < kappa; ++trial) {
                {
                    HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "trial_reset");
                    // Reset side array and block weights for this trial
                    std::fill(side.begin(), side.end(), u8(2));
                    trial_pm.reset_weights(); // zeroes bweights and n_vertices
                }

                recurse(g, trial_pm, all_vertices, 0, k, g.g_weight);

                const weight_t cut = compute_edge_cut(g, trial_pm);

                HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "copy_best");
                if (cut < best_cut) {
                    best_cut = cut;
                    out_pm.copy_from(trial_pm);
                }
            }
        }
    };

} // namespace HeiProMap

#endif // HEIPROMAP_RECURSIVE_BISECTION_H
