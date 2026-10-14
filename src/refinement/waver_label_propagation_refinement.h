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

#ifndef HEIPROMAP_WAVER_LABEL_PROPAGATION_REFINEMENT_H
#define HEIPROMAP_WAVER_LABEL_PROPAGATION_REFINEMENT_H

#include <cmath>
#include <limits>
#include <vector>

#include <omp.h>

#include "../definitions.h"
#include "../datastructures/block_conn.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../datastructures/distance_1_matching.h"
#include "../datastructures/active_block_scheduling.h"
#include "../datastructures/block_adjacency.h"
#include "../utility/aligned_array.h"
#include "../utility/profiler.h"
#include "../utility/random_engine.h"

namespace HeiProMap {

    class WaverLabelPropagationConfiguration {
    public:
        explicit WaverLabelPropagationConfiguration(std::string t_name) : name(std::move(t_name)) {}

        std::string name;
        bool enabled = false;
        u64 max_iteration = 25;

        bool force_parallel_alg = false;
        bool use_active_scheduling = true;
        bool use_static_matchings = true;

        u64 max_matching_rounds = 512 * 512;
        u64 min_matching_threshold = 8;

        // four-term cost model parameters
        f64 alpha_end  = 0.0;   // alpha_send + alpha_recv
        f64 beta_end   = 0.0;   // beta_send  + beta_recv
        f64 alpha_link = 0.0;
        f64 beta_link  = 0.0;
    };

    template<bool LARGE_K>
    class WaverLabelPropagationRefinement {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
        u64 m_threads = 1;

        std::vector<RandomEngine> rnd_engines;
        WaverLabelPropagationConfiguration config = WaverLabelPropagationConfiguration("waver_label_propagation");

        ActiveBlockScheduling active_block_scheduling;
        Distance1Matching<LARGE_K> d1_matcher;
        BlockAdjacency block_adjacency;

        static inline __attribute__((always_inline)) void apply_move(graph_t &g,
                                                                     bv_manager_t &bv_manager,
                                                                     p_manager_t &p_manager,
                                                                     const vertex_t u,
                                                                     const weight_t u_weight,
                                                                     const partition_t from_id,
                                                                     const partition_t to_id) {
            bv_manager.move(g, p_manager, u, from_id, to_id);
            p_manager.move_serial(u, u_weight, from_id, to_id);
        }

        // Four-term cost delta for moving u from old_id to new_id.
        // Returns -Δcost; positive means the move reduces total cost.
        //
        // Δcost = alpha_end  * Δ|E_π|
        //       + beta_end   * ΔΣS
        //       + alpha_link * ΔΣh
        //       + beta_link  * ΔΣh·S
        //
        // Δ|E_π| approximation: alpha_end is paid once per communicating (src,dst)
        // PE pair (due to message coalescing), NOT once per edge. We approximate
        // Δ|E_π| using boolean indicators over u's neighbors:
        //   - any neighbor in old_id → u's A-A edges become B-A cuts → B-A pair appears: -alpha_end
        //   - any neighbor in new_id → u's A-B cuts become B-B internal → A-B pair may disappear: +alpha_end
        // This charges alpha_end at most once per direction across the A-B boundary,
        // replacing the previous per-edge overcounting.
        //
        // ΔΣS, ΔΣh, ΔΣh·S: per-edge (correct — each edge contributes independently).
        template<bool t_uniform_e_weights, typename DistanceOracleT, typename QGraphT>
        inline __attribute__((always_inline)) f64 compute_cost_delta(const graph_t &g,
                                                                     const vertex_t u,
                                                                     const partition_t old_id,
                                                                     const partition_t new_id,
                                                                     const p_manager_t &p_manager,
                                                                     DistanceOracleT &d_oracle,
                                                                     const QGraphT &q_graph,
                                                                     const block_conn_t &block_conn) const {
            const f64 dist_old_new = static_cast<f64>(d_oracle.get(old_id, new_id));
            f64 delta = 0.0;

            // alpha_end: exact Δ|E_π| via block_conn (u's per-block weights) and q_graph (block-pair totals).
            // u moves from old_id (A) to new_id (B). Positive delta = cost reduction.
            for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) {
                const partition_t C       = block_conn.get_id(i);
                const weight_t    u_to_C  = block_conn.get_w(i);

                if (C == old_id) {
                    // A-A internal → B-A cut; new pair if B-A didn't previously exist
                    if (q_graph.get_weight(new_id, old_id) == 0) {
                        delta -= config.alpha_end;
                    }
                } else if (C == new_id) {
                    // A-B cut → B-B internal; pair disappears iff u was the sole A-B connection
                    if (q_graph.get_weight(old_id, new_id) == u_to_C) {
                        delta += config.alpha_end;
                    }
                } else {
                    // A-C cut → B-C cut
                    if (q_graph.get_weight(old_id, C) == u_to_C) {
                        delta += config.alpha_end;  // A-C pair disappears
                    }
                    if (q_graph.get_weight(new_id, C) == 0) {
                        delta -= config.alpha_end;  // B-C pair appears
                    }
                }
            }

            // beta_end, alpha_link, beta_link: per-edge terms
            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t    v    = g.edges_v[i];
                const partition_t v_id = p_manager[v];
                const f64         w    = t_uniform_e_weights ? 1.0 : static_cast<f64>(g.edge_w(i));

                if (v_id == old_id) {
                    delta -= config.beta_end * w
                           + (config.alpha_link + config.beta_link * w) * dist_old_new;
                } else if (v_id == new_id) {
                    delta += config.beta_end * w
                           + (config.alpha_link + config.beta_link * w) * dist_old_new;
                } else {
                    const f64 d_old = static_cast<f64>(d_oracle.get(v_id, old_id));
                    const f64 d_new = static_cast<f64>(d_oracle.get(v_id, new_id));
                    delta += (config.alpha_link + config.beta_link * w) * (d_old - d_new);
                }
            }
            return delta;
        }

    public:
        WaverLabelPropagationRefinement() = default;
        ~WaverLabelPropagationRefinement() = default;

        void set_pool(MemoryPool *pool) {
            active_block_scheduling.set_pool(pool);
            d1_matcher.set_pool(pool);
            block_adjacency.set_pool(pool);
        }

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const u64 t_threads,
                        const u64 t_seed,
                        const WaverLabelPropagationConfiguration &t_config,
                        MemoryPool *pool = nullptr) {
            set_pool(pool);
            m_n = t_n;
            m_m = t_m;
            m_k = t_k;
            m_threads = t_threads;
            config = t_config;

            rnd_engines.resize(m_threads);
            for (u64 t = 0; t < m_threads; ++t) {
                rnd_engines[t] = RandomEngine(t_seed + t);
            }

            active_block_scheduling.initialize(m_k, pool);
            d1_matcher.initialize(m_k, pool);
            block_adjacency.initialize(m_k, m_threads, pool);
        }

        template<typename DistanceOracleT, typename QGraphT>
        void refine(graph_t &g,
                    DistanceOracleT &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    QGraphT &q_graph,
                    block_conn_t &block_conn) {
            if (g.uniform_v_weights && g.uniform_e_weights) {
                refine_impl<true, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn);
            } else if (g.uniform_v_weights) {
                refine_impl<true, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn);
            } else if (g.uniform_e_weights) {
                refine_impl<false, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn);
            } else {
                refine_impl<false, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn);
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights, typename DistanceOracleT, typename QGraphT>
        void refine_impl(graph_t &g,
                         DistanceOracleT &d_oracle,
                         bv_manager_t &bv_manager,
                         p_manager_t &p_manager,
                         QGraphT &q_graph,
                         block_conn_t &block_conn) {
            if (config.force_parallel_alg || m_threads > 1) {
                refine_impl_parallel<t_uniform_v_weights, t_uniform_e_weights>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn);
            } else {
                refine_impl_serial<t_uniform_v_weights, t_uniform_e_weights>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn);
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights, typename DistanceOracleT, typename QGraphT>
        void refine_impl_parallel(graph_t &g,
                                  DistanceOracleT &d_oracle,
                                  bv_manager_t &bv_manager,
                                  p_manager_t &p_manager,
                                  QGraphT &q_graph,
                                  block_conn_t &block_conn) {
            active_block_scheduling.reset(m_k);

            std::vector<std::vector<std::pair<partition_t, partition_t>>> static_matchings;
            std::vector<std::pair<partition_t, partition_t>> matching;
            std::vector<std::vector<vertex_t>> thread_local_boundary(m_threads);

            for (u64 iteration = 0; iteration < config.max_iteration; ++iteration) {
                HEIPROMAP_PROFILE_SCOPE("refinement", "WaverLPRefinement", "recompute_adjacency");
                block_adjacency.compute(g, p_manager, bv_manager, m_threads);

                d1_matcher.reset_used_edges();
                bool positive_move_occurred = false;

                if (config.use_static_matchings) {
                    HEIPROMAP_PROFILE_SCOPE("refinement", "WaverLPRefinement", "matching");
                    d1_matcher.compute_static_matchings(block_adjacency, p_manager, active_block_scheduling.active_this_round, static_matchings, config.max_matching_rounds, config.min_matching_threshold);

                    if (static_matchings.empty()) { break; }

                    for (size_t i = 0; i < static_matchings.size(); ++i) {
                        matching.assign(static_matchings[i].begin(), static_matchings[i].end());
                        if (matching.empty()) { continue; }

                        HEIPROMAP_PROFILE_SCOPE("refinement", "WaverLPRefinement", "process");

                        #pragma omp parallel for num_threads(m_threads) schedule(dynamic)
                        for (size_t j = 0; j < matching.size(); ++j) {
                            const partition_t A = matching[j].first;
                            const partition_t B = matching[j].second;

                            const u64 tid = omp_get_thread_num();
                            RandomEngine &rng = rnd_engines[tid];

                            auto &local_boundary = thread_local_boundary[tid];
                            local_boundary.clear();
                            local_boundary.insert(local_boundary.end(), bv_manager.boundary(A).begin(), bv_manager.boundary(A).end());
                            local_boundary.insert(local_boundary.end(), bv_manager.boundary(B).begin(), bv_manager.boundary(B).end());

                            bool moved_in_pair = false;

                            for (const vertex_t u : local_boundary) {
                                const partition_t u_id = p_manager[u];
                                if (u_id != A && u_id != B) { continue; }

                                const partition_t target_id = (u_id == A) ? B : A;
                                const weight_t u_weight = t_uniform_v_weights ? 1 : g.v_weight(u);

                                if (p_manager.get_bweight(target_id) + u_weight > p_manager.lmax[target_id]) { continue; }

                                const f64 delta = compute_cost_delta<t_uniform_e_weights>(g, u, u_id, target_id, p_manager, d_oracle, q_graph, block_conn);

                                if (delta > 0.0 || (delta == 0.0 && rng.get_f32() < 0.5f)) {
                                    apply_move(g, bv_manager, p_manager, u, u_weight, u_id, target_id);
                                    if (delta > 0.0) { moved_in_pair = true; }
                                }
                            }

                            if (moved_in_pair) {
                                if (config.use_active_scheduling) {
                                    active_block_scheduling.activate(A, B);
                                }
                                #pragma omp atomic write
                                positive_move_occurred = true;
                            }
                        }
                    }
                } else {
                    std::vector<std::pair<partition_t, partition_t>> dyn_matching;

                    HEIPROMAP_PROFILE_SCOPE("refinement", "WaverLPRefinement", "matching");
                    d1_matcher.find_matching(q_graph, p_manager, active_block_scheduling.active_this_round, dyn_matching);

                    if (dyn_matching.empty()) { break; }

                    bool found_matching = true;
                    while (found_matching) {
                        HEIPROMAP_PROFILE_SCOPE("refinement", "WaverLPRefinement", "process");

                        #pragma omp parallel for num_threads(m_threads) schedule(dynamic)
                        for (size_t i = 0; i < dyn_matching.size(); ++i) {
                            const partition_t A = dyn_matching[i].first;
                            const partition_t B = dyn_matching[i].second;

                            const u64 tid = omp_get_thread_num();
                            RandomEngine &rng = rnd_engines[tid];

                            auto &local_boundary = thread_local_boundary[tid];
                            local_boundary.clear();
                            local_boundary.insert(local_boundary.end(), bv_manager.boundary(A).begin(), bv_manager.boundary(A).end());
                            local_boundary.insert(local_boundary.end(), bv_manager.boundary(B).begin(), bv_manager.boundary(B).end());

                            bool moved_in_pair = false;

                            for (const vertex_t u : local_boundary) {
                                const partition_t u_id = p_manager[u];
                                if (u_id != A && u_id != B) { continue; }

                                const partition_t target_id = (u_id == A) ? B : A;
                                const weight_t u_weight = t_uniform_v_weights ? 1 : g.v_weight(u);

                                if (p_manager.get_bweight(target_id) + u_weight > p_manager.lmax[target_id]) { continue; }

                                const f64 delta = compute_cost_delta<t_uniform_e_weights>(g, u, u_id, target_id, p_manager, d_oracle, q_graph, block_conn);

                                if (delta > 0.0 || (delta == 0.0 && rng.get_f32() < 0.5f)) {
                                    apply_move(g, bv_manager, p_manager, u, u_weight, u_id, target_id);
                                    if (delta > 0.0) { moved_in_pair = true; }
                                }
                            }

                            if (moved_in_pair) {
                                if (config.use_active_scheduling) {
                                    active_block_scheduling.activate(A, B);
                                }
                                #pragma omp atomic write
                                positive_move_occurred = true;
                            }
                        }

                        if constexpr (LARGE_K) {
                            if (dyn_matching.size() < config.min_matching_threshold) { break; }
                        }

                        HEIPROMAP_PROFILE_SCOPE("refinement", "WaverLPRefinement", "matching");
                        found_matching = d1_matcher.find_matching(block_adjacency, p_manager, active_block_scheduling.active_this_round, dyn_matching);
                    }
                }

                if (!positive_move_occurred) { break; }

                if (config.use_active_scheduling) {
                    active_block_scheduling.next_round();
                }
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights, typename DistanceOracleT, typename QGraphT>
        void refine_impl_serial(graph_t &g,
                                DistanceOracleT &d_oracle,
                                bv_manager_t &bv_manager,
                                p_manager_t &p_manager,
                                QGraphT &q_graph,
                                block_conn_t &block_conn) {
            RandomEngine &random_engine = rnd_engines[0];

            bool positive_move_occurred = true;
            for (u64 iteration = 0; iteration < config.max_iteration && positive_move_occurred; ++iteration) {
                positive_move_occurred = false;

                HEIPROMAP_PROFILE_SCOPE("refinement", "WaverLPRefinement", "process_vertices");
                for (vertex_t u = 0; u < g.n; ++u) {
                    if (!bv_manager.is_boundary(u)) { continue; }

                    const weight_t u_weight = t_uniform_v_weights ? 1 : g.v_weight(u);
                    const partition_t u_id  = p_manager[u];

                    partition_t best_id    = NO_ID;
                    f64         best_delta = -std::numeric_limits<f64>::max();
                    f32         counter    = 0.0f;

                    partition_t neighbor_blocks[256];
                    size_t n_neighbor_blocks = 0;

                    for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                        const partition_t v_id = p_manager[g.edges_v[i]];
                        if (v_id == u_id) { continue; }

                        bool found = false;
                        for (size_t j = 0; j < n_neighbor_blocks; ++j) {
                            if (neighbor_blocks[j] == v_id) { found = true; break; }
                        }
                        if (!found && n_neighbor_blocks < 256) {
                            neighbor_blocks[n_neighbor_blocks++] = v_id;
                        }
                    }

                    for (size_t bi = 0; bi < n_neighbor_blocks; ++bi) {
                        const partition_t id = neighbor_blocks[bi];
                        if (p_manager.get_bweight(id) + u_weight > p_manager.lmax[id]) { continue; }

                        const f64 delta = compute_cost_delta<t_uniform_e_weights>(g, u, u_id, id, p_manager, d_oracle, q_graph, block_conn);

                        if (delta > best_delta) {
                            best_id    = id;
                            best_delta = delta;
                            counter    = 1.0f;
                        } else if (delta == best_delta) {
                            counter += 1.0f;
                            if (random_engine.get_f32() < 1.0f / counter) {
                                best_id = id;
                            }
                        }
                    }

                    if (best_delta > 0.0 || (best_delta == 0.0 && random_engine.get_f32() < 0.5f)) {
                        q_graph.move(g, p_manager, u, u_id, best_id);
                        block_conn.move(g, u, u_id, best_id);
                        apply_move(g, bv_manager, p_manager, u, u_weight, u_id, best_id);
                        positive_move_occurred |= (best_delta > 0.0);
                    }
                }
            }
        }
    };

}

#endif //HEIPROMAP_WAVER_LABEL_PROPAGATION_REFINEMENT_H
