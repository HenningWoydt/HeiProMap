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
#include "../datastructures/subgraph_extractor.h"
#include "../definitions.h"
#include "../utility/random_engine.h"
#include "../utility/indexed_max_heap.h"
#include "../utility/translation_table.h"

#include "../refinement/label_propagation_refinement.h"
#include "../refinement/quotient_graph_refinement.h"
#include "../refinement/flow_based_refinement.h"

namespace HeiProMap {
    enum struct BisectionMethod {
        BFS,
        GGG,
        HYBRID
    };

    struct RecursiveBisectionConfiguration {
        u64 kappa = 1;
        bool use_full_refine = false;
        BisectionMethod method = BisectionMethod::BFS;
        LabelPropagationConfiguration lp_config = LabelPropagationConfiguration("LP");
        QuotientGraphRefinementConfiguration qg_config = QuotientGraphRefinementConfiguration("QG");
        FlowBasedRefinementConfiguration flow_config = FlowBasedRefinementConfiguration("Flow");

        RecursiveBisectionConfiguration() {
            lp_config.enabled = true;
            lp_config.max_iteration = 5;

            qg_config.enabled = true;
            qg_config.max_iteration = 1;
            qg_config.alpha = 5.0;
            qg_config.min_n_steps = 3;
            qg_config.use_preemptive_exit = true;

            // enable flow based refinement
            flow_config.enabled = false;
            flow_config.max_global_iteration = 1;
            flow_config.max_local_iteration = 1;
            flow_config.alpha = 1.0;
            flow_config.alpha_upper_bound = 64.0;
            flow_config.alpha_modifier = 2.0;
            flow_config.use_closed_vertex_set = true;
            flow_config.closed_vertex_sets_repeats = 500;
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

        RandomEngine *rng = nullptr;

        // Refinement instances
        LabelPropagationRefinement lp_refine;
        QuotientGraphRefinement qg_refine;
        FlowBasedRefinement flow_refine;

        void refine_pm(const CSRGraph &g,
                       PartitionManager &pm,
                       const AlignedArray<weight_t> &lmax_constraints,
                       const RecursiveBisectionConfiguration &config) {
            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "refine_pm");
            if (pm.k <= 1) return;

            bv_manager_t bv;
            bv.initialize(g.n, pm.k);
            q_graph_t qg;
            qg.initialize(pm.k);
            block_conn_t bc;
            bc.initialize(g.n, g.m, pm.k);

            for (vertex_t u = 0; u < g.n; ++u) {
                bc.begin_vertex(g, u);
                partition_t u_id = pm[u];
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    vertex_t v = g.edges_v[i];
                    partition_t v_id = pm[v];
                    bc.add_connection(u, v_id, g.edges_w[i]);
                    if (u_id != v_id) {
                        bv.add(u, u_id);
                        if (u < v) qg.add_edge(u_id, v_id, g.edges_w[i]);
                    }
                }
            }

            d_oracle_t do_oracle;
            std::vector<partition_t> k_vec(1, pm.k);
            std::vector<weight_t> dist_vec(1, 1);
            do_oracle.initialize(k_vec, dist_vec);

            // Run refinements
            if (config.lp_config.enabled) {
                lp_refine.initialize(g.n, g.m, pm.k, 1, 0, config.lp_config);
                lp_refine.refine(const_cast<CSRGraph &>(g), do_oracle, bv, pm, qg, bc, lmax_constraints, g.uniform_v_weights, g.uniform_e_weights);
            }
            if (config.qg_config.enabled) {
                qg_refine.initialize(g.n, g.m, pm.k, 1, 0, config.qg_config);
                qg_refine.refine(const_cast<CSRGraph &>(g), do_oracle, bv, pm, qg, bc, lmax_constraints, g.uniform_v_weights, g.uniform_e_weights);
            }
            if (config.flow_config.enabled) {
                flow_refine.initialize(g.n, g.m, pm.k, 1, 0, config.flow_config);
                flow_refine.refine(const_cast<CSRGraph &>(g), do_oracle, bv, pm, qg, bc, lmax_constraints, g.uniform_v_weights, g.uniform_e_weights);
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

        void assign_unassigned_components(const CSRGraph &g, std::vector<u8> &side, weight_t &w0, weight_t &w1, weight_t lmax_left, weight_t lmax_right) {
            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "assign_unassigned");
            std::vector<std::vector<vertex_t>> components;
            
            for (vertex_t u = 0; u < g.n; ++u) {
                if (side[u] == 2 && dist[u] == -1) {
                    components.emplace_back();
                    auto &comp = components.back();
                    
                    std::queue<vertex_t> q;
                    q.push(u);
                    dist[u] = -2;
                    
                    while (!q.empty()) {
                        vertex_t curr = q.front();
                        q.pop();
                        comp.push_back(curr);
                        
                        for (size_t i = g.neighborhoods[curr]; i < g.neighborhoods[curr + 1]; ++i) {
                            vertex_t v = g.edges_v[i];
                            if (side[v] == 2 && dist[v] == -1) {
                                dist[v] = -2;
                                q.push(v);
                            }
                        }
                    }
                }
            }
            
            for (const auto &comp : components) {
                for (vertex_t v : comp) dist[v] = -1;
                
                weight_t comp_weight = 0;
                for (vertex_t v : comp) comp_weight += g.v_weights[v];
                
                u8 s = 0;
                if (w0 + comp_weight <= lmax_left && w1 + comp_weight <= lmax_right) {
                    s = ((f64) w0 / (f64) lmax_left <= (f64) w1 / (f64) lmax_right) ? u8(0) : u8(1);
                } else if (w0 + comp_weight <= lmax_left) {
                    s = 0;
                } else if (w1 + comp_weight <= lmax_right) {
                    s = 1;
                } else {
                    s = ((f64) (w0 + comp_weight) / (f64) lmax_left <= (f64) (w1 + comp_weight) / (f64) lmax_right) ? u8(0) : u8(1);
                }
                
                for (vertex_t v : comp) side[v] = s;
                if (s == 0) w0 += comp_weight;
                else w1 += comp_weight;
            }
        }

        void bisect_bfs(const CSRGraph &g,
                        weight_t lmax_left,
                        weight_t lmax_right,
                        std::vector<u8> &side) {
            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "bisect_bfs");
            if (g.n == 0) return;

            // ---- Furthest-point seed selection ----
            vertex_t temp_s = rng->get_u64() % g.n;
            vertex_t seed_0 = bfs_furthest(g, temp_s);
            vertex_t seed_1 = bfs_furthest(g, seed_0);

            if (seed_1 == seed_0 && g.n > 1) {
                seed_1 = (seed_0 == 0) ? 1 : 0;
            }

            // ---- Interleaved BFS growing ----
            side.assign(g.n, 2);
            side[seed_0] = 0;
            side[seed_1] = 1;
            weight_t w0 = g.v_weights[seed_0];
            weight_t w1 = g.v_weights[seed_1];

            if (seed_0 == seed_1) {
                side[seed_0] = 0;
                w0 = g.v_weights[seed_0];
                w1 = 0;
            }

            std::queue<vertex_t> q0, q1;
            q0.push(seed_0);
            if (seed_1 != seed_0) { q1.push(seed_1); }

            while (!q0.empty() || !q1.empty()) {
                const bool can_grow_0 = !q0.empty() && w0 < lmax_left;
                const bool can_grow_1 = !q1.empty() && w1 < lmax_right;

                if (!can_grow_0 && !can_grow_1) { break; }

                f64 fill_0 = (f64) w0 / (f64) lmax_left;
                f64 fill_1 = (f64) w1 / (f64) lmax_right;
                const bool grow_left = (can_grow_0 && can_grow_1)
                                           ? (fill_0 <= fill_1)
                                           : can_grow_0;

                std::queue<vertex_t> &q = grow_left ? q0 : q1;
                const u8 s = grow_left ? u8(0) : u8(1);
                weight_t &w = grow_left ? w0 : w1;
                const weight_t target_w = grow_left ? lmax_left : lmax_right;

                vertex_t u = q.front();
                q.pop();

                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    vertex_t v = g.edges_v[i];
                    if (side[v] == 2) {
                        if (w + g.v_weights[v] > target_w) continue;
                        side[v] = s;
                        w += g.v_weights[v];
                        q.push(v);
                    }
                }
            }

            assign_unassigned_components(g, side, w0, w1, lmax_left, lmax_right);
        }

        void bisect_ggg(const CSRGraph &g,
                        weight_t lmax_left,
                        weight_t lmax_right,
                        std::vector<u8> &side) {
            if (g.n == 0) return;
            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "bisect_ggg");

            pq.clear();
            side.assign(g.n, 2);

            vertex_t seed = rng->get_u64() % g.n;
            side[seed] = 0;
            weight_t w0 = g.v_weights[seed];

            for (size_t i = g.neighborhoods[seed]; i < g.neighborhoods[seed + 1]; ++i) {
                vertex_t v = g.edges_v[i];
                if (side[v] == 2) {
                    pq.push_increment(v, g.edges_w[i]);
                }
            }

            while (!pq.empty() && w0 < lmax_left) {
                vertex_t u = pq.top_key();
                pq.pop();

                if (w0 + g.v_weights[u] > lmax_left) continue;

                side[u] = 0;
                w0 += g.v_weights[u];

                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    vertex_t v = g.edges_v[i];
                    if (side[v] == 2) {
                        pq.push_increment(v, g.edges_w[i]);
                    }
                }
            }

            weight_t w1 = 0;
            assign_unassigned_components(g, side, w0, w1, lmax_left, lmax_right);
        }

        weight_t compute_bisection_cut(const CSRGraph &g, const std::vector<u8> &side) {
            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "compute_bisection_cut");
            weight_t cut = 0;
            for (vertex_t u = 0; u < g.n; ++u) {
                u8 u_side = side[u];
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    vertex_t v = g.edges_v[i];
                    if (side[v] != u_side) {
                        cut += g.edges_w[i];
                    }
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

            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "recurse_overhead");
            partition_t k_left = k / 2;
            partition_t k_right = k - k_left;

            f64 max_avg_weight = (1.0 + imbalance) * ((f64) global_weight / (f64) total_k);
            weight_t lmax_left = (weight_t) std::ceil(max_avg_weight * (f64) k_left);
            weight_t lmax_right = (weight_t) std::ceil(max_avg_weight * (f64) k_right);

            weight_t best_bisect_cut = std::numeric_limits<weight_t>::max();
            bool best_is_balanced = false;
            std::vector<u8> best_side;

            auto perform_trials = [&](BisectionMethod m) {
                for (u64 trial = 0; trial < config.kappa; ++trial) {
                    std::vector<u8> current_side;
                    if (m == BisectionMethod::BFS) {
                        bisect_bfs(g, lmax_left, lmax_right, current_side);
                    } else {
                        bisect_ggg(g, lmax_left, lmax_right, current_side);
                    }

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
                        best_side = std::move(current_side);
                    }
                }
            };

            if (config.method == BisectionMethod::HYBRID) {
                perform_trials(BisectionMethod::BFS);
                perform_trials(BisectionMethod::GGG);
            } else {
                perform_trials(config.method);
            }

            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "recurse_overhead");
            CSRGraph left_g, right_g;
            TranslationTable<vertex_t> left_to_g, right_to_g;
            SubgraphExtractor::extract(g, best_side, 0, left_g, left_to_g);
            SubgraphExtractor::extract(g, best_side, 1, right_g, right_to_g);

            PartitionManager left_pm;
            left_pm.initialize(left_g.n, k_left, left_g.g_weight);
            recurse(left_g, left_pm, k_left, total_k, global_weight, imbalance, config);

            PartitionManager right_pm;
            right_pm.initialize(right_g.n, k_right, right_g.g_weight);
            recurse(right_g, right_pm, k_right, total_k, global_weight, imbalance, config);

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
                weight_t lmax_single = (weight_t) std::ceil(max_avg_weight);
                for (partition_t i = 0; i < k; ++i) {
                    kway_lmax[i] = lmax_single;
                }
                refine_pm(g, pm, kway_lmax, config);
            }
        }

    public:
        void partition(const CSRGraph &g,
                       PartitionManager &out_pm,
                       partition_t k,
                       u64 seed,
                       f64 imbalance = 0.03,
                       const RecursiveBisectionConfiguration &config = RecursiveBisectionConfiguration()) {
            RandomEngine rand_engine(seed);
            rng = &rand_engine;

            HEIPROMAP_PROFILE_SCOPE("partition", "RecursiveBisectionPartitioner", "allocate");
            dist.assign(g.n, s32(-1));
            pq.initialize(g.n);
            
            if (out_pm.k != k || out_pm.n != g.n) {
                out_pm.initialize(g.n, k, g.g_weight);
            } else {
                out_pm.reset_weights();
            }

            recurse(g, out_pm, k, k, g.g_weight, imbalance, config);
        }
    };
} // namespace HeiProMap

#endif // HEIPROMAP_RECURSIVE_BISECTION_H

