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
#include "../datastructures/distance_oracle.h"
#include "../datastructures/partition_manager.h"
#include "../datastructures/quotient_graph.h"
#include "../utility/aligned_array.h"
#include "../utility/profiler.h"
#include "../utility/qap.h"
#include "../utility/random_engine.h"

namespace HeiProMap {
    class LabelPropagationConfiguration final {
    public:
        explicit LabelPropagationConfiguration(const std::string &t_name) {
            name = t_name;
        }

        std::string name;
        bool enabled = false;
        u64 max_iteration = 25; // how many iterations to run the algorithm at most

        bool use_parallel_alg = false;
    };

    class LabelPropagationRefinement final {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
        u64 m_threads = 1;

        AlignedArray<vertex_t> curr_boundary;
        size_t curr_boundary_size = 0;

        AlignedArray<partition_t> blocks;
        AlignedArray<weight_t> blocks_qap_delta;
        size_t blocks_size = 0;

        std::vector<RandomEngine> rnd_engines;
        const LabelPropagationConfiguration *config = nullptr;

        AlignedArray<u8> active_this_round;
        AlignedArray<u8> used_this_round;

    public:
        LabelPropagationRefinement() = default;

        ~LabelPropagationRefinement() = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const u64 t_threads,
                        const u64 seed,
                        const LabelPropagationConfiguration &i_config) {
            m_n = t_n;
            m_m = t_m;
            m_k = t_k;
            m_threads = t_threads;

            config = &i_config;

            curr_boundary.initialize(m_n);
            curr_boundary_size = 0;

            blocks.initialize(m_k);
            blocks_qap_delta.initialize(m_k);
            blocks_size = 0;

            rnd_engines.resize(m_threads);
            for (u64 t = 0; t < m_threads; ++t) {
                rnd_engines[t] = RandomEngine(seed + t);
            }

            active_this_round.initialize(m_k);
            used_this_round.initialize(m_k * m_k);
        }

        void refine(graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph,
                    block_conn_t &block_conn,
                    const AlignedArray<weight_t> &lmax_constraints) {
            if (g.uniform_v_weights && g.uniform_e_weights) refine_impl<true, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
            else if (g.uniform_v_weights) refine_impl<true, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
            else if (g.uniform_e_weights) refine_impl<false, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
            else refine_impl<false, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void refine_impl(graph_t &g,
                         d_oracle_t &d_oracle,
                         bv_manager_t &bv_manager,
                         p_manager_t &p_manager,
                         q_graph_t &q_graph,
                         block_conn_t &block_conn,
                         const AlignedArray<weight_t> &lmax_constraints) {
            if (config->use_parallel_alg) {
                refine_impl_parallel<t_uniform_v_weights, t_uniform_e_weights>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
            } else {
                refine_impl_serial<t_uniform_v_weights, t_uniform_e_weights>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void refine_impl_parallel(graph_t &g,
                                  d_oracle_t &d_oracle,
                                  bv_manager_t &bv_manager,
                                  p_manager_t &p_manager,
                                  q_graph_t &q_graph,
                                  block_conn_t &block_conn,
                                  const AlignedArray<weight_t> &lmax_constraints) {
            bool positive_move_occurred = true;
            for (u64 iteration = 0; iteration < config->max_iteration && positive_move_occurred; ++iteration) {
                positive_move_occurred = false;

                HEIPROMAP_PROFILE_SCOPE("refinement", "LabelPropagationRefinement", "pick_vertices_parallel");
                active_this_round.initialize(m_k, 1);
                std::fill_n(used_this_round.get_ptr(), m_k * m_k, 0);

                std::vector<std::pair<partition_t, partition_t>> matching;
                bool found_matching = q_graph.find_distance_3_matching(active_this_round, used_this_round, matching);

                while (found_matching) {
                    #pragma omp parallel for num_threads(m_threads) schedule(dynamic)
                    for (size_t i = 0; i < matching.size(); ++i) {
                        partition_t A = matching[i].first;
                        partition_t B = matching[i].second;

                        u64 tid = omp_get_thread_num();
                        RandomEngine &rng = rnd_engines[tid];

                        // Copy the boundary vertices of A and B to a local vector to avoid modification issues
                        std::vector<vertex_t> local_boundary;
                        local_boundary.reserve(bv_manager.size(A) + bv_manager.size(B));
                        for (size_t idx = 0; idx < bv_manager.size(A); ++idx) {
                            local_boundary.push_back(bv_manager.get(A, idx));
                        }
                        for (size_t idx = 0; idx < bv_manager.size(B); ++idx) {
                            local_boundary.push_back(bv_manager.get(B, idx));
                        }

                        // Refine vertices sequentially within matched pair (A, B)
                        for (vertex_t u : local_boundary) {
                            partition_t u_id = p_manager[u];
                            if (u_id != A && u_id != B) { continue; }

                            partition_t target_id = (u_id == A) ? B : A;
                            weight_t u_weight = t_uniform_v_weights ? 1 : g.v_weights[u];

                            if (p_manager.get_bweight(target_id) + u_weight > lmax_constraints[target_id]) { continue; }

                            weight_t qap_delta = get_u_qap_delta_t<t_uniform_e_weights>(g, u, u_id, target_id, p_manager, d_oracle, block_conn);

                            if (qap_delta > 0 || (qap_delta == 0 && rng.get_f32() < 0.5)) {
                                bv_manager.move(g, p_manager, u, u_id, target_id);
                                q_graph.move(g, p_manager, u, u_id, target_id);
                                block_conn.move(g, u, u_id, target_id);
                                p_manager.move_serial(u, u_weight, u_id, target_id);
                                #pragma omp atomic write
                                positive_move_occurred = true;
                            }
                        }
                    }

                    found_matching = q_graph.find_distance_3_matching(active_this_round, used_this_round, matching);
                }
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void refine_impl_serial(graph_t &g,
                                d_oracle_t &d_oracle,
                                bv_manager_t &bv_manager,
                                p_manager_t &p_manager,
                                q_graph_t &q_graph,
                                block_conn_t &block_conn,
                                const AlignedArray<weight_t> &lmax_constraints) {
            RandomEngine &random_engine = rnd_engines[0];

            bool positive_move_occurred = true;
            for (u64 iteration = 0; iteration < config->max_iteration && positive_move_occurred; ++iteration) {
                positive_move_occurred = false;

                HEIPROMAP_PROFILE_SCOPE("refinement", "LabelPropagationRefinement", "process_vertices");
                // for (size_t j = 0; j < curr_boundary_size; ++j) {
                for (vertex_t u = 0; u < g.n; ++u) {
                    // vertex_t u = curr_boundary[j];

                    if (!bv_manager.is_boundary(u)) { continue; }

                    weight_t u_weight = t_uniform_v_weights ? 1 : g.v_weights[u];
                    partition_t u_id = p_manager[u];

                    partition_t best_id = NO_ID;
                    weight_t best_qap_delta = -std::numeric_limits<weight_t>::max();
                    f32 counter = 0;

                    for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) {
                        partition_t id = block_conn.get_id(i);
                        weight_t v_id_weight = p_manager.get_bweight(id);

                        if (id == u_id) { continue; }
                        if (v_id_weight + u_weight > lmax_constraints[id]) { continue; }

                        weight_t qap_delta = get_u_qap_delta_t<t_uniform_e_weights>(g, u, u_id, id, p_manager, d_oracle, block_conn);
                        if (qap_delta > best_qap_delta) {
                            best_id = id;
                            best_qap_delta = qap_delta;
                            counter = 1.0;
                        } else if (qap_delta == best_qap_delta) {
                            counter += 1.0;
                            if (random_engine.get_f32() < 1.0f / counter) { best_id = id; }
                        }
                    }

                    if (best_qap_delta > 0 || (best_qap_delta == 0 && random_engine.get_f32() < 0.5)) {
                        bv_manager.move(g, p_manager, u, u_id, best_id);
                        q_graph.move(g, p_manager, u, u_id, best_id);
                        block_conn.move(g, u, u_id, best_id);
                        p_manager.move_serial(u, u_weight, u_id, best_id);
                        positive_move_occurred |= best_qap_delta > 0;
                    }
                }
            }
        }
    };
}

#endif //HEIPROMAP_LABEL_PROPAGATION_REFINEMENT_H
