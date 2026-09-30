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

#ifndef HEIPROMAP_LABEL_PROPAGATION_REFINEMENT_H
#define HEIPROMAP_LABEL_PROPAGATION_REFINEMENT_H

#include <cmath>
#include <limits>
#include <vector>

#include <omp.h>

#include "../definitions.h"
#include "../datastructures/block_conn.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/csr_graph.h"
#include "../distance_oracles/distance_oracle.h"
#include "../datastructures/partition_manager.h"
#include "../datastructures/quotient_graph.h"
#include "../datastructures/distance_3_matching.h"
#include "../datastructures/active_block_scheduling.h"
#include "../utility/aligned_array.h"
#include "../utility/profiler.h"
#include "../utility/qap.h"
#include "../utility/random_engine.h"

namespace HeiProMap {
    class LabelPropagationConfiguration {
    public:
        explicit LabelPropagationConfiguration(std::string t_name) : name(std::move(t_name)) {
        }

        std::string name;
        bool enabled = false;
        u64 max_iteration = 25; // how many iterations to run the algorithm at most

        bool force_parallel_alg = false;
        bool use_edge_cut = true;
        bool use_active_scheduling = true;
        u64 min_matching_threshold = 8;
        bool use_static_matchings = true;
    };

    template<bool LARGE_K>
    class LabelPropagationRefinement {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
        u64 m_threads = 1;

        std::vector<RandomEngine> rnd_engines;
        LabelPropagationConfiguration config = LabelPropagationConfiguration("default_label_propagation");

        ActiveBlockScheduling active_block_scheduling;
        Distance3Matching<LARGE_K> d3_matcher;

        template<typename QGraphT>
        static inline __attribute__((always_inline)) void apply_move(graph_t &g,
                                                                     bv_manager_t &bv_manager,
                                                                     p_manager_t &p_manager,
                                                                     QGraphT &q_graph,
                                                                     block_conn_t &block_conn,
                                                                     const vertex_t u,
                                                                     const weight_t u_weight,
                                                                     const partition_t from_id,
                                                                     const partition_t to_id) {
            bv_manager.move(g, p_manager, u, from_id, to_id);
            q_graph.move(g, p_manager, u, from_id, to_id);
            block_conn.move(g, u, from_id, to_id);
            p_manager.move_serial(u, u_weight, from_id, to_id);
        }

    public:
        LabelPropagationRefinement() = default;

        ~LabelPropagationRefinement() = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const u64 t_threads,
                        const u64 t_seed,
                        const LabelPropagationConfiguration &t_config) {
            m_n = t_n;
            m_m = t_m;
            m_k = t_k;
            m_threads = t_threads;

            config = t_config;

            rnd_engines.resize(m_threads);
            for (u64 t = 0; t < m_threads; ++t) {
                rnd_engines[t] = RandomEngine(t_seed + t);
            }

            active_block_scheduling.initialize(m_k);
            d3_matcher.initialize(m_k);
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

            std::vector<std::vector<std::pair<partition_t, partition_t> > > static_matchings;
            std::vector<std::pair<partition_t, partition_t> > matching;

            for (u64 iteration = 0; iteration < config.max_iteration; ++iteration) {
                HEIPROMAP_PROFILE_SCOPE("refinement", "LabelPropagationRefinement", "reset_used_edges");
                d3_matcher.reset_used_edges();

                bool positive_move_occurred = false;

                if (config.use_static_matchings) {
                    HEIPROMAP_PROFILE_SCOPE("refinement", "LabelPropagationRefinement", "matching");
                    d3_matcher.compute_static_matchings(q_graph, p_manager, active_block_scheduling.active_this_round, static_matchings, 512, config.min_matching_threshold);

                    if (static_matchings.empty()) {
                        break;
                    }

                    for (size_t i = 0; i < static_matchings.size(); ++i) {
                        HEIPROMAP_PROFILE_SCOPE("refinement", "LabelPropagationRefinement", "filter_matching");

                        const auto &raw_batch = static_matchings[i];

                        if (i > 0) {
                            d3_matcher.filter_valid_matching(q_graph, raw_batch, matching);
                        }

                        if (matching.empty()) {
                            continue;
                        }

                        HEIPROMAP_PROFILE_SCOPE("refinement", "LabelPropagationRefinement", "process");
                        // Collect changed block IDs from parallel loop
                        std::vector<partition_t> changed_list;

                        #pragma omp parallel for num_threads(m_threads) schedule(dynamic)
                        for (size_t j = 0; j < matching.size(); ++j) {
                            const partition_t A = matching[j].first;
                            const partition_t B = matching[j].second;

                            const u64 tid = omp_get_thread_num();
                            RandomEngine &rng = rnd_engines[tid];

                            // Copy boundary vertices of A and B to local vector to avoid modification issues
                            const size_t size_A = bv_manager.size(A);
                            const size_t size_B = bv_manager.size(B);
                            std::vector<vertex_t> local_boundary;
                            local_boundary.reserve(size_A + size_B);
                            for (size_t idx = 0; idx < size_A; ++idx) {
                                local_boundary.push_back(bv_manager.get(A, idx));
                            }
                            for (size_t idx = 0; idx < size_B; ++idx) {
                                local_boundary.push_back(bv_manager.get(B, idx));
                            }

                            const bool last_level_pair = config.use_edge_cut && d_oracle.last_level_pair(A, B);
                            bool moved_in_pair = false;

                            // Refine vertices sequentially within matched pair (A, B)
                            for (const vertex_t u: local_boundary) {
                                const partition_t u_id = p_manager[u];
                                if (u_id != A && u_id != B) { continue; }

                                const partition_t target_id = (u_id == A) ? B : A;
                                const weight_t u_weight = t_uniform_v_weights ? 1 : g.v_weights[u];

                                if (p_manager.get_bweight(target_id) + u_weight > p_manager.lmax[target_id]) { continue; }

                                weight_t qap_delta;
                                if (last_level_pair) {
                                    qap_delta = get_u_edge_cut_delta_t<t_uniform_e_weights>(g, u, u_id, target_id, p_manager, block_conn);
                                } else {
                                    qap_delta = get_u_qap_delta_t<t_uniform_e_weights>(g, u, u_id, target_id, p_manager, d_oracle, block_conn);
                                }

                                if (qap_delta > 0 || (qap_delta == 0 && rng.get_f32() < 0.5f)) {
                                    apply_move(g, bv_manager, p_manager, q_graph, block_conn, u, u_weight, u_id, target_id);
                                    if (qap_delta > 0) {
                                        moved_in_pair = true;
                                    }
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
                    HEIPROMAP_PROFILE_SCOPE("refinement", "LabelPropagationRefinement", "matching");
                    std::vector<std::pair<partition_t, partition_t> > dyn_matching;
                    bool found_matching = d3_matcher.find_matching(q_graph, p_manager, active_block_scheduling.active_this_round, dyn_matching);

                    if (!found_matching) {
                        break;
                    }

                    while (found_matching) {
                        HEIPROMAP_PROFILE_SCOPE("refinement", "LabelPropagationRefinement", "process");

                        #pragma omp parallel for num_threads(m_threads) schedule(dynamic)
                        for (size_t i = 0; i < dyn_matching.size(); ++i) {
                            const partition_t A = dyn_matching[i].first;
                            const partition_t B = dyn_matching[i].second;

                            const u64 tid = omp_get_thread_num();
                            RandomEngine &rng = rnd_engines[tid];

                            // Copy boundary vertices of A and B to local vector to avoid modification issues
                            const size_t size_A = bv_manager.size(A);
                            const size_t size_B = bv_manager.size(B);
                            std::vector<vertex_t> local_boundary;
                            local_boundary.reserve(size_A + size_B);
                            for (size_t idx = 0; idx < size_A; ++idx) {
                                local_boundary.push_back(bv_manager.get(A, idx));
                            }
                            for (size_t idx = 0; idx < size_B; ++idx) {
                                local_boundary.push_back(bv_manager.get(B, idx));
                            }

                            const bool last_level_pair = config.use_edge_cut && d_oracle.last_level_pair(A, B);
                            bool moved_in_pair = false;

                            // Refine vertices sequentially within matched pair (A, B)
                            for (const vertex_t u: local_boundary) {
                                const partition_t u_id = p_manager[u];
                                if (u_id != A && u_id != B) { continue; }

                                const partition_t target_id = (u_id == A) ? B : A;
                                const weight_t u_weight = t_uniform_v_weights ? 1 : g.v_weights[u];

                                if (p_manager.get_bweight(target_id) + u_weight > p_manager.lmax[target_id]) { continue; }

                                weight_t qap_delta;
                                if (last_level_pair) {
                                    qap_delta = get_u_edge_cut_delta_t<t_uniform_e_weights>(g, u, u_id, target_id, p_manager, block_conn);
                                } else {
                                    qap_delta = get_u_qap_delta_t<t_uniform_e_weights>(g, u, u_id, target_id, p_manager, d_oracle, block_conn);
                                }

                                if (qap_delta > 0 || (qap_delta == 0 && rng.get_f32() < 0.5f)) {
                                    apply_move(g, bv_manager, p_manager, q_graph, block_conn, u, u_weight, u_id, target_id);
                                    if (qap_delta > 0) {
                                        moved_in_pair = true;
                                    }
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
                            if (dyn_matching.size() < config.min_matching_threshold) {
                                break;
                            }
                        }

                        HEIPROMAP_PROFILE_SCOPE("refinement", "LabelPropagationRefinement", "matching");
                        found_matching = d3_matcher.find_matching(q_graph, p_manager, active_block_scheduling.active_this_round, dyn_matching);
                    }
                }

                if (!positive_move_occurred) {
                    break;
                }

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

                HEIPROMAP_PROFILE_SCOPE("refinement", "LabelPropagationRefinement", "process_vertices");
                for (vertex_t u = 0; u < g.n; ++u) {
                    if (!bv_manager.is_boundary(u)) { continue; }

                    const weight_t u_weight = t_uniform_v_weights ? 1 : g.v_weights[u];
                    const partition_t u_id = p_manager[u];

                    partition_t best_id = NO_ID;
                    weight_t best_qap_delta = -std::numeric_limits<weight_t>::max();
                    f32 counter = 0.0f;

                    for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) {
                        const partition_t id = block_conn.get_id(i);
                        const weight_t v_id_weight = p_manager.get_bweight(id);

                        if (id == u_id || v_id_weight + u_weight > p_manager.lmax[id]) { continue; }

                        const weight_t qap_delta = get_u_qap_delta_t<t_uniform_e_weights>(g, u, u_id, id, p_manager, d_oracle, block_conn);
                        if (qap_delta > best_qap_delta) {
                            best_id = id;
                            best_qap_delta = qap_delta;
                            counter = 1.0f;
                        } else if (qap_delta == best_qap_delta) {
                            counter += 1.0f;
                            if (random_engine.get_f32() < 1.0f / counter) {
                                best_id = id;
                            }
                        }
                    }

                    if (best_qap_delta > 0 || (best_qap_delta == 0 && random_engine.get_f32() < 0.5f)) {
                        apply_move(g, bv_manager, p_manager, q_graph, block_conn, u, u_weight, u_id, best_id);
                        positive_move_occurred |= (best_qap_delta > 0);
                    }
                }
            }
        }
    };
}

#endif //HEIPROMAP_LABEL_PROPAGATION_REFINEMENT_H
