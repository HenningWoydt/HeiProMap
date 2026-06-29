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
#include <unordered_set>

#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../datastructures/subgraph_extractor.h"
#include "../definitions.h"
#include "../utility/random_engine.h"
#include "../utility/indexed_max_heap.h"
#include "../utility/translation_table.h"
#include "../utility/utils.h"

#include "../refinement/label_propagation_refinement.h"
#include "../refinement/quotient_graph_refinement.h"
#include "../refinement/flow_based_refinement.h"
#include "../refinement/swap_refinement.h"

namespace HeiProMap {
    enum struct BisectionMethod {
        BFS,
        GGG,
        HYBRID
    };

    struct RecursiveBisectionConfiguration {
        u64 kappa = 1;
        bool use_full_refine = false;
        BisectionMethod method = BisectionMethod::HYBRID;
        SwapRefinementConfiguration swap_config = SwapRefinementConfiguration("Swap");

        RecursiveBisectionConfiguration() {
            swap_config.enabled = true;
            swap_config.max_iteration = 5;
        }
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
        std::vector<s32> dist;
        IndexedMaxHeap<weight_t> pq;

        RandomEngine rng;

        // Refinement instances
        SwapRefinement swap_refine;

        partition_t original_k = 0;
        std::vector<partition_t> k_vec;
        std::vector<weight_t> dist_vec;
        DistanceOracle d_oracle;

        void refine_pm(const CSRGraph &g,
                       PartitionManager &pm,
                       const AlignedArray<weight_t> &lmax_constraints,
                       const RecursiveBisectionConfiguration &config) {
            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "refine_pm");
            if (pm.k <= 1) return;

            if (k_vec.empty()) {
                k_vec = std::vector<partition_t>(1, original_k);
                dist_vec = std::vector<weight_t>(1, 1);
                d_oracle.initialize(k_vec, dist_vec);

                swap_refine.initialize(0, config.swap_config);
            }

            if (config.swap_config.enabled) {
                swap_refine.refine(const_cast<CSRGraph &>(g), d_oracle, pm, lmax_constraints, g.uniform_v_weights, g.uniform_e_weights);
            }
        }


        void perform_full_refinement(const CSRGraph &g,
                                     std::vector<u8> &side,
                                     weight_t lmax_left,
                                     weight_t lmax_right,
                                     const RecursiveBisectionConfiguration &config) {
            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "full_refine");

            // Setup 2-way environment for the current subgraph
            p_manager_t sub_pm;
            sub_pm.initialize(g.n, 2, g.g_weight);
            sub_pm.reset_weights();
            for (vertex_t u = 0; u < g.n; ++u) {
                sub_pm.set(u, g.v_weights[u], side[u]);
            }

            AlignedArray<weight_t> lmax_constraints;
            lmax_constraints.initialize(2);
            lmax_constraints[0] = lmax_left;
            lmax_constraints[1] = lmax_right;

            refine_pm(g, sub_pm, lmax_constraints, config);

            // Update side array from refined partition
            for (vertex_t u = 0; u < g.n; ++u) {
                side[u] = (u8) sub_pm[u];
            }
        }

        void evaluate_trial(const CSRGraph &g,
                            weight_t lmax_left,
                            weight_t lmax_right,
                            const RecursiveBisectionConfiguration &config,
                            std::vector<u8> &current_side,
                            std::vector<u8> &best_side,
                            weight_t &best_bisect_cut,
                            bool &best_is_balanced) {
            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "evaluate_trial");

            if (config.use_full_refine) {
                perform_full_refinement(g, current_side, lmax_left, lmax_right, config);
            }

            weight_t cut = compute_bisection_cut(g, current_side);

            weight_t w0 = 0, w1 = 0;
            for (vertex_t u = 0; u < g.n; ++u) {
                if (current_side[u] == 0) w0 += g.v_weights[u];
                else w1 += g.v_weights[u];
            }
            bool is_balanced = (w0 <= lmax_left && w1 <= lmax_right);

            bool update_best = false;
            if (best_side.empty()) {
                update_best = true;
            } else if (is_balanced && !best_is_balanced) {
                update_best = true;
            } else if (is_balanced == best_is_balanced && cut < best_bisect_cut) {
                update_best = true;
            }

            if (update_best) {
                best_bisect_cut = cut;
                best_is_balanced = is_balanced;
                best_side = current_side;
            }
        }

        vertex_t bfs_furthest(const CSRGraph &g, vertex_t start) {
            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "bfs_furthest");

            std::vector<vertex_t> visited;
            visited.reserve(std::min((vertex_t) 64, g.n));

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

                    if (dist[v] == -1) {
                        dist[v] = dist[u] + 1;
                        visited.push_back(v);
                        q.push(v);

                        if (dist[v] > dist[furthest]) {
                            furthest = v;
                        }
                    }
                }
            }

            // Reset dist for every visited vertex
            for (vertex_t v: visited) {
                dist[v] = -1;
            }

            return furthest;
        }

        void run_bfs_trials(const CSRGraph &g, u64 kappa, weight_t lmax_left, weight_t lmax_right, const RecursiveBisectionConfiguration &config, std::vector<u8> &best_side, weight_t &best_cut, bool &best_is_balanced) {
            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "run_bfs_trials");
            if (g.n == 0 || kappa == 0) return;

            std::vector<u8> side(g.n, 0);
            std::vector<u8> visited(g.n, 0);
            weight_t target_w1 = (weight_t) ((f64) g.g_weight * (f64) lmax_right / (f64) (lmax_left + lmax_right));

            for (u64 trial = 0; trial < kappa; ++trial) {
                vertex_t seed = rng.get_u64() % g.n;
                vertex_t seed_1 = bfs_furthest(g, seed);

                std::fill(side.begin(), side.end(), 0);
                std::fill(visited.begin(), visited.end(), 0);
                visited[seed_1] = 1;
                side[seed_1] = 1;
                weight_t w1 = g.v_weights[seed_1];

                std::queue<vertex_t> q;
                q.push(seed_1);

                while (!q.empty() && w1 < target_w1) {
                    vertex_t u = q.front();
                    q.pop();

                    for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                        vertex_t v = g.edges_v[i];
                        if (visited[v] == 0) {
                            visited[v] = 1;
                            if (w1 + g.v_weights[v] <= lmax_right) {
                                side[v] = 1;
                                w1 += g.v_weights[v];
                            }
                            q.push(v);
                        }
                    }
                }

                evaluate_trial(g, lmax_left, lmax_right, config, side, best_side, best_cut, best_is_balanced);
            }
        }

        void run_ggg_trials(const CSRGraph &g, u64 kappa, weight_t lmax_left, weight_t lmax_right, const RecursiveBisectionConfiguration &config,
                            std::vector<u8> &best_side, weight_t &best_cut, bool &best_is_balanced) {
            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "run_ggg_trials");
            if (g.n == 0 || kappa == 0) return;

            std::vector<weight_t> weighted_degrees(g.n, 0);
            for (vertex_t u = 0; u < g.n; ++u) {
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    weighted_degrees[u] += g.edges_w[i];
                }
            }

            u64 attempt_threshold = 10;
            std::vector<vertex_t> used_seeds;
            for (u64 trial = 0; trial < kappa; ++trial) {
                u64 attempts = 0;
                while (attempts < attempt_threshold) {
                    vertex_t temp_s = rng.get_u64() % g.n;
                    if (std::find(used_seeds.begin(), used_seeds.end(), temp_s) == used_seeds.end()) {
                        used_seeds.push_back(temp_s);
                        break;
                    }
                    attempts++;
                }
                if (attempts == attempt_threshold) break;
            }

            std::vector<u8> side(g.n, 0);
            std::vector<u8> visited(g.n, 0);
            weight_t target_w1 = (weight_t) ((f64) g.g_weight * (f64) lmax_right / (f64) (lmax_left + lmax_right));

            for (vertex_t seed: used_seeds) {
                pq.clear();
                std::fill(side.begin(), side.end(), 0);
                std::fill(visited.begin(), visited.end(), 0);
                visited[seed] = 1;
                side[seed] = 1;
                weight_t w1 = g.v_weights[seed];

                auto add_or_update = [&](vertex_t v, weight_t w) {
                    if (pq.entry_exists(v)) pq.increment(v, 2 * w);
                    else {
                        pq.push(v, 2 * w - weighted_degrees[v]);
                    }
                };

                for (size_t i = g.neighborhoods[seed]; i < g.neighborhoods[seed + 1]; ++i) {
                    vertex_t v = g.edges_v[i];
                    if (visited[v] == 0) add_or_update(v, g.edges_w[i]);
                }

                while (!pq.empty() && w1 < target_w1) {
                    vertex_t u = pq.top_key();
                    pq.pop();
                    if (visited[u] == 1) continue;
                    visited[u] = 1;

                    if (w1 + g.v_weights[u] <= lmax_right) {
                        side[u] = 1;
                        w1 += g.v_weights[u];
                    }

                    for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                        vertex_t v = g.edges_v[i];
                        if (visited[v] == 0) add_or_update(v, g.edges_w[i]);
                    }
                }

                evaluate_trial(g, lmax_left, lmax_right, config, side, best_side, best_cut, best_is_balanced);
            }
        }

        void perform_single_bisection(const CSRGraph &g,
                                      weight_t lmax_left,
                                      weight_t lmax_right,
                                      const RecursiveBisectionConfiguration &config,
                                      std::vector<u8> &best_side,
                                      weight_t &best_bisect_cut,
                                      bool &best_is_balanced) {
            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "perform_single_bisection");
            if (g.n == 0) {
                best_side.clear();
                best_bisect_cut = 0;
                best_is_balanced = true;
                return;
            }

            best_side.clear();
            best_bisect_cut = std::numeric_limits<weight_t>::max();
            best_is_balanced = false;

            if (config.method == BisectionMethod::HYBRID) {
                u64 kappa_bfs = config.kappa / 2;
                u64 kappa_ggg = config.kappa - kappa_bfs;
                run_bfs_trials(g, kappa_bfs, lmax_left, lmax_right, config, best_side, best_bisect_cut, best_is_balanced);
                run_ggg_trials(g, kappa_ggg, lmax_left, lmax_right, config, best_side, best_bisect_cut, best_is_balanced);
            } else if (config.method == BisectionMethod::BFS) {
                run_bfs_trials(g, config.kappa, lmax_left, lmax_right, config, best_side, best_bisect_cut, best_is_balanced);
            } else {
                run_ggg_trials(g, config.kappa, lmax_left, lmax_right, config, best_side, best_bisect_cut, best_is_balanced);
            }
        }

        weight_t compute_bisection_cut(const CSRGraph &g, const std::vector<u8> &side) {
            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "compute_bisection_cut");
            weight_t cut = 0;
            for (vertex_t u = 0; u < g.n; ++u) {
                u8 u_side = side[u];
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    vertex_t v = g.edges_v[i];
                    cut += (side[v] != u_side) * g.edges_w[i];
                }
            }
            return cut / 2;
        }

        void recurse(const CSRGraph &g,
                     PartitionManager &pm,
                     partition_t k,
                     partition_t total_k,
                     weight_t global_weight,
                     f64 imbalance,
                     const RecursiveBisectionConfiguration &config) {
            if (g.n == 0) return;

            if (k == 1) {
                HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "recurse_base_case");
                pm.reset_weights();
                for (vertex_t v = 0; v < g.n; ++v) {
                    pm.set(v, g.v_weights[v], 0);
                }
                return;
            }

            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "find_components");
            partition_t k_left = k / 2;
            partition_t k_right = k - k_left;

            weight_t avg_weight = std::ceil((1.0 + imbalance) * ((f64) global_weight / (f64) total_k));
            weight_t lmax_left = avg_weight * (weight_t) k_left;
            weight_t lmax_right = avg_weight * (weight_t) k_right;

            // --- Phase 1: Greedy Placement ---
            std::vector<u8> side(g.n, 2); // 2 means unassigned
            weight_t w0_base = 0;
            weight_t w1_base = 0;
            std::vector<std::vector<vertex_t> > components_to_split;

            std::vector<std::vector<vertex_t> > components;
            std::vector<weight_t> component_weights;
            std::vector<s8> visited(g.n, 0);

            for (vertex_t i = 0; i < g.n; ++i) {
                if (!visited[i]) {
                    components.emplace_back();
                    auto &comp = components.back();
                    weight_t w = 0;
                    std::queue<vertex_t> q;
                    q.push(i);
                    visited[i] = 1;
                    while (!q.empty()) {
                        vertex_t u = q.front();
                        q.pop();
                        comp.push_back(u);
                        w += g.v_weights[u];
                        for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                            vertex_t v = g.edges_v[j];
                            if (!visited[v]) {
                                visited[v] = 1;
                                q.push(v);
                            }
                        }
                    }
                    component_weights.push_back(w);
                }
            }

            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "assign_components");
            std::vector<size_t> sorted_indices(components.size());
            std::iota(sorted_indices.begin(), sorted_indices.end(), 0);
            std::sort(sorted_indices.begin(), sorted_indices.end(), [&](size_t i, size_t j) {
                return component_weights[i] > component_weights[j];
            });

            for (size_t idx: sorted_indices) {
                weight_t w = component_weights[idx];
                bool can_fit_left = (w0_base + w <= lmax_left);
                bool can_fit_right = (w1_base + w <= lmax_right);

                if (can_fit_left || can_fit_right) {
                    u8 target_side = 2;
                    if (can_fit_left && can_fit_right) {
                        target_side = ((f64) (w0_base + w) / lmax_left <= (f64) (w1_base + w) / lmax_right) ? 0 : 1;
                    } else if (can_fit_left) {
                        target_side = 0;
                    } else {
                        target_side = 1;
                    }

                    for (vertex_t v: components[idx]) side[v] = target_side;
                    if (target_side == 0) w0_base += w;
                    else w1_base += w;
                } else {
                    components_to_split.push_back(components[idx]);
                }
            }

            // --- Phase 2: Iterative Splitting ---
            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "split_components");
            for (const auto &component: components_to_split) {
                CSRGraph sub_g;
                TranslationTable<vertex_t> sub_to_g;
                std::vector<u8> mask(g.n, 0);
                for (vertex_t v: component) mask[v] = 1;
                SubgraphExtractor::extract(g, mask, 1, sub_g, sub_to_g);

                weight_t sub_lmax_0 = lmax_left - w0_base;
                weight_t sub_lmax_1 = lmax_right - w1_base;

                std::vector<u8> sub_side;
                weight_t sub_cut = 0; // This is ignored for now, can be used for logging
                bool sub_balanced = false;

                perform_single_bisection(sub_g, sub_lmax_0, sub_lmax_1, config, sub_side, sub_cut, sub_balanced);

                weight_t w0_sub = 0;
                weight_t w1_sub = 0;
                for (vertex_t v = 0; v < sub_g.n; ++v) {
                    side[sub_to_g.get_o(v)] = sub_side[v];
                    if (sub_side[v] == 0) w0_sub += sub_g.v_weights[v];
                    else w1_sub += sub_g.v_weights[v];
                }
                w0_base += w0_sub;
                w1_base += w1_sub;
            }

            if (config.use_full_refine) {
                perform_full_refinement(g, side, lmax_left, lmax_right, config);
            }

            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "extract_subgraphs");
            // --- Final recursion step ---
            CSRGraph left_g, right_g;
            TranslationTable<vertex_t> left_to_g, right_to_g;
            SubgraphExtractor::extract(g, side, 0, left_g, left_to_g);
            SubgraphExtractor::extract(g, side, 1, right_g, right_to_g);

            PartitionManager left_pm;
            left_pm.initialize(left_g.n, k_left, left_g.g_weight);
            recurse(left_g, left_pm, k_left, total_k, global_weight, imbalance, config);

            PartitionManager right_pm;
            right_pm.initialize(right_g.n, k_right, right_g.g_weight);
            recurse(right_g, right_pm, k_right, total_k, global_weight, imbalance, config);

            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "merge_results");
            // Merge back into pm
            pm.reset_weights();
            for (vertex_t v = 0; v < left_g.n; ++v) {
                vertex_t u = left_to_g.get_o(v);
                pm.set(u, left_g.v_weights[v], left_pm[v]);
            }
            for (vertex_t v = 0; v < right_g.n; ++v) {
                vertex_t u = right_to_g.get_o(v);
                pm.set(u, right_g.v_weights[v], right_pm[v] + k_left);
            }

            if (config.use_full_refine) {
                AlignedArray<weight_t> kway_lmax;
                kway_lmax.initialize(k);
                for (partition_t i = 0; i < k; ++i) {
                    kway_lmax[i] = avg_weight;
                }
                refine_pm(g, pm, kway_lmax, config);
            }
        }

    public:
        void partition(const CSRGraph &g,
                       PartitionManager &p_manager,
                       partition_t k,
                       u64 seed,
                       f64 imbalance,
                       RecursiveBisectionConfiguration &config) {
            rng = RandomEngine(seed);
            original_k = k;

            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "allocate");
            dist.assign(g.n, s32(-1));
            pq.initialize(g.n);

            if (p_manager.k != k || p_manager.n != g.n) {
                p_manager.initialize(g.n, k, g.g_weight);
            } else {
                p_manager.reset_weights();
            }

            recurse(g, p_manager, k, k, g.g_weight, imbalance, config);
        }
    };
} // namespace HeiProMap

#endif // HEIPROMAP_RECURSIVE_BISECTION_H
