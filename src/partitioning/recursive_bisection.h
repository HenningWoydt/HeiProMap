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

#include "../utility/indexed_max_heap.h"

namespace HeiProMap {

    enum struct BisectionMethod {
        BFS,
        GGG,
        HYBRID
    };

    /**
     * Partitions a graph into k blocks using recursive bisection.
     *
     * Each bisection step uses either BFS region growing, Greedy Graph Growing (GGG),
     * or a Hybrid approach that tries both and picks the best local result.
     *
     * BFS region growing:
     *   1. Two seeds are chosen via furthest-point seeding (the pair of vertices
     *      that are maximally far apart inside the current vertex subset).
     *   2. Two BFS frontiers grow simultaneously from the seeds. At each step
     *      the lighter frontier is expanded, naturally producing balanced halves.
     *   3. Any vertices unreachable from both seeds (disconnected components)
     *      are greedily assigned to whichever side is lighter.
     *
     * Greedy Graph Growing (GGG):
     *   1. A random seed is chosen.
     *   2. A priority queue stores all boundary vertices of the growing region,
     *      ordered by their "gain" (connectivity to the growing region).
     *   3. Vertices are added to the region until the target weight is reached.
     *
     * Hybrid:
     *   1. Runs `kappa` trials of BFS and `kappa` trials of GGG.
     *   2. Picks the split with the lowest local edge-cut and recurses.
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
         */
        std::vector<u8> side;

        /**
         * BFS distance array.
         */
        std::vector<s32> dist;

        /**
         * Gains for greedy refinement.
         */
        std::vector<s64> gains;

        /**
         * Lock markers for greedy refinement.
         */
        std::vector<u8> locked;

        /**
         * IndexedMaxHeap for GGG and refinement.
         */
        IndexedMaxHeap<weight_t> pq;

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
        void bisect_bfs(const CSRGraph &g,
                        const std::vector<vertex_t> &vertices,
                        weight_t target_weight_left,
                        std::vector<vertex_t> &left_out,
                        std::vector<vertex_t> &right_out) {
            if (vertices.empty()) {
                left_out.clear();
                right_out.clear();
                return;
            }

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

        /**
         * Greedy Graph Growing (GGG) bisection.
         *
         * @param g               The graph.
         * @param vertices        Vertices belonging to this bisection subproblem.
         * @param target_weight_left  Target total vertex weight for the left side.
         * @param left_out        Output: vertices assigned to the left side.
         * @param right_out       Output: vertices assigned to the right side.
         */
        void bisect_ggg(const CSRGraph &g,
                        const std::vector<vertex_t> &vertices,
                        weight_t target_weight_left,
                        std::vector<vertex_t> &left_out,
                        std::vector<vertex_t> &right_out) {
            if (vertices.empty()) {
                left_out.clear();
                right_out.clear();
                return;
            }

            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "bisect_ggg");

            pq.clear();

            // Pick a random seed vertex
            vertex_t seed = vertices[rng->get_u64() % vertices.size()];
            side[seed] = 0;
            weight_t w0 = g.v_weights[seed];

            // Initialize PQ with neighbors of seed
            for (size_t i = g.neighborhoods[seed]; i < g.neighborhoods[seed + 1]; ++i) {
                vertex_t v = g.edges_v[i];
                if (side[v] == 2) {
                    pq.push_increment(v, g.edges_w[i]);
                }
            }

            // Grow the region
            while (!pq.empty() && w0 < target_weight_left) {
                vertex_t u = pq.top_key();
                pq.pop();

                side[u] = 0;
                w0 += g.v_weights[u];

                // Update neighbors
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    vertex_t v = g.edges_v[i];
                    if (side[v] == 2) {
                        pq.push_increment(v, g.edges_w[i]);
                    }
                }
            }

            // Assign remaining vertices to right
            for (vertex_t v : vertices) {
                if (side[v] == 2) {
                    side[v] = 1;
                }
            }

            // Collect results
            left_out.clear();
            right_out.clear();
            left_out.reserve(vertices.size() / 2 + 1);
            right_out.reserve(vertices.size() / 2 + 1);

            for (vertex_t v : vertices) {
                if (side[v] == 0) { left_out.push_back(v); }
                else              { right_out.push_back(v); }
            }
        }

        /**
         * Refines the current bisection using a greedy local search (FM-like).
         *
         * @param g               The graph.
         * @param vertices        Vertices belonging to this bisection subproblem.
         * @param target_weight_left  Target total vertex weight for the left side.
         * @param imbalance       Allowed imbalance (slack).
         */
        void refine_bisection(const CSRGraph &g,
                              const std::vector<vertex_t> &vertices,
                              weight_t target_weight_left,
                              f64 imbalance) {
            if (vertices.empty()) return;

            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "refine_bisection");

            weight_t total_weight = 0;
            weight_t w_left = 0;
            for (vertex_t u : vertices) {
                total_weight += g.v_weights[u];
                if (side[u] == 0) { w_left += g.v_weights[u]; }
            }
            weight_t w_right = total_weight - w_left;

            // Simple greedy refinement: move vertices that improve the cut and keep balance
            bool improved = true;
            int max_passes = 10; // Limit passes to prevent excessive work

            while (improved && max_passes-- > 0) {
                improved = false;
                pq.clear();
                std::fill(locked.begin(), locked.end(), 0);

                // Calculate initial gains for all vertices in this subset
                for (vertex_t u : vertices) {
                    s64 internal_deg = 0;
                    s64 external_deg = 0;
                    u8 u_side = side[u];

                    for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                        vertex_t v = g.edges_v[i];
                        if (side[v] == 3) continue; // outside subset

                        if (side[v] == u_side) {
                            internal_deg += g.edges_w[i];
                        } else {
                            external_deg += g.edges_w[i];
                        }
                    }

                    gains[u] = external_deg - internal_deg;
                    if (gains[u] > 0) {
                        pq.push(u, (weight_t)(gains[u] + 1000000000)); 
                    }
                }

                while (!pq.empty()) {
                    vertex_t u = pq.top_key();
                    s64 u_gain = (s64)pq.top() - 1000000000;
                    pq.pop();

                    if (u_gain <= 0) break;
                    if (locked[u]) continue;

                    u8 from = side[u];
                    u8 to = 1 - from;
                    weight_t u_w = g.v_weights[u];

                    // Check balance constraint (with provided imbalance slack)
                    if (from == 0) { // moving 0 -> 1
                        if (w_right + u_w > (total_weight - target_weight_left) * (1.0 + imbalance)) continue; 
                    } else { // moving 1 -> 0
                        if (w_left + u_w > target_weight_left * (1.0 + imbalance)) continue;
                    }

                    // Move vertex
                    side[u] = to;
                    if (from == 0) { w_left -= u_w; w_right += u_w; }
                    else           { w_left += u_w; w_right -= u_w; }
                    
                    locked[u] = 1;
                    improved = true;

                    // Update gains of neighbors incrementally
                    for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                        vertex_t v = g.edges_v[i];
                        if (side[v] == 3 || locked[v]) continue;

                        // Correct gain update:
                        // If v is on side 'from' (u's old side), u just moved to 'to' (external side for v).
                        // So v's external degree increases, internal decreases -> Gain(v) increases.
                        if (side[v] == from) {
                            gains[v] += 2 * (s64)g.edges_w[i];
                        } else {
                            gains[v] -= 2 * (s64)g.edges_w[i];
                        }

                        if (gains[v] > 0) {
                            pq.push_update(v, (weight_t)(gains[v] + 1000000000));
                        } else if (pq.entry_exists(v)) {
                            pq.update(v, 1000000000); 
                        }
                    }
                }
            }
        }

        /**
         * Computes the edge-cut of the current bisection within the current subset.
         */
        weight_t compute_bisection_cut(const CSRGraph &g, const std::vector<vertex_t> &vertices) {
            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "compute_bisection_cut");

            weight_t cut = 0;
            for (vertex_t u : vertices) {
                u8 u_side = side[u];
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    vertex_t v = g.edges_v[i];
                    // Only count edges within the current subset that span the bisection
                    if (side[v] != 3 && side[v] != u_side) {
                        cut += g.edges_w[i];
                    }
                }
            }
            return cut / 2; // each edge is counted twice
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
         * @param method        Bisection method to use.
         * @param kappa         Number of trials per bisection step.
         * @param imbalance     Allowed imbalance (slack).
         */
        void recurse(const CSRGraph &g,
                     PartitionManager &pm,
                     std::vector<vertex_t> &vertices,
                     partition_t block_start,
                     partition_t k,
                     weight_t total_weight,
                     BisectionMethod method,
                     u64 kappa,
                     f64 imbalance) {
            // Safety: handle empty vertex set
            if (vertices.empty()) return;

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

            // Bisect the current vertex set using multiple trials
            weight_t best_bisect_cut = std::numeric_limits<weight_t>::max();
            std::vector<vertex_t> best_left_v, best_right_v;

            auto perform_trials = [&](BisectionMethod m) {
                for (u64 trial = 0; trial < kappa; ++trial) {
                    for (vertex_t v : vertices) { side[v] = 2; }
                    std::vector<vertex_t> left_v_trial, right_v_trial;
                    if (m == BisectionMethod::BFS) {
                        bisect_bfs(g, vertices, target_left, left_v_trial, right_v_trial);
                    } else {
                        bisect_ggg(g, vertices, target_left, left_v_trial, right_v_trial);
                    }

                    // Perform greedy FM refinement with slack
                    refine_bisection(g, vertices, target_left, imbalance);

                    weight_t cut = compute_bisection_cut(g, vertices);
                    if (cut < best_bisect_cut) {
                        best_bisect_cut = cut;
                        best_left_v.clear();
                        best_right_v.clear();
                        for (vertex_t u : vertices) {
                            if (side[u] == 0) best_left_v.push_back(u);
                            else              best_right_v.push_back(u);
                        }
                    }
                }
            };

            if (method == BisectionMethod::HYBRID) {
                perform_trials(BisectionMethod::BFS);
                perform_trials(BisectionMethod::GGG);
            } else {
                perform_trials(method);
            }

            // Ensure side[] reflects the best split for recursion
            for (vertex_t v : best_left_v)  { side[v] = 0; }
            for (vertex_t v : best_right_v) { side[v] = 1; }

            weight_t w_left = 0, w_right = 0;
            {
                HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "recurse_overhead");
                // Compute actual weights of each half for the recursive calls
                for (vertex_t v : best_left_v)  { w_left  += g.v_weights[v]; }
                for (vertex_t v : best_right_v) { w_right += g.v_weights[v]; }

                // ---- Recurse left ----
                // Hide the right subset from the left recursion
                for (vertex_t v : best_right_v) { side[v] = 3; }
                for (vertex_t v : best_left_v)  { side[v] = 2; } // reset to unassigned
            }

            recurse(g, pm, best_left_v, block_start, k_left, w_left, method, kappa, imbalance);

            {
                HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "recurse_overhead");
                // ---- Recurse right ----
                // Hide the (now processed) left subset from the right recursion
                for (vertex_t v : best_left_v)  { side[v] = 3; }
                for (vertex_t v : best_right_v) { side[v] = 2; } // reset to unassigned
            }

            recurse(g, pm, best_right_v, block_start + k_left, k_right, w_right, method, kappa, imbalance);
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
         * greedy BFS region growing or GGG. The result is written into `out_pm`.
         *
         * Performs `kappa` trials per bisection step and picks the best local cut.
         *
         * @param g        Input graph.
         * @param out_pm   Output PartitionManager. Must already be initialized
         *                 with `out_pm.initialize(g.n, k, 0)` before calling.
         * @param k        Number of blocks.
         * @param seed     Base random seed.
         * @param kappa    Number of independent trials per step (default 1).
         * @param imbalance Allowed imbalance slack (default 0.03).
         * @param method   Bisection method to use.
         */
        void partition(const CSRGraph &g,
                       PartitionManager &out_pm,
                       partition_t k,
                       u64 seed,
                       u64 kappa = 1,
                       f64 imbalance = 0.03,
                       BisectionMethod method = BisectionMethod::BFS) {
            RandomEngine rand_engine(seed);
            rng = &rand_engine;

            std::vector<vertex_t> all_vertices(g.n);
            {
                HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "allocate");
                // Allocate/resize working arrays
                side.assign(g.n, u8(2));
                dist.assign(g.n, s32(-1));
                pq.initialize(g.n);
                std::iota(all_vertices.begin(), all_vertices.end(), vertex_t(0));
                out_pm.reset_weights(); // zeroes bweights and n_vertices
            }

            recurse(g, out_pm, all_vertices, 0, k, g.g_weight, method, kappa, imbalance);
        }
    };

} // namespace HeiProMap

#endif // HEIPROMAP_RECURSIVE_BISECTION_H
