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

#include <omp.h>

#include "ISerialRefiner.h"
#include "../definitions.h"
#include "../utility/random_engine.h"

namespace HeiProMap {
    class LabelPropagationConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit LabelPropagationConfiguration(const std::string &t_name) : ISerialRefinerConfiguration(t_name) {
        }

        u64 max_iteration = 25; // how many iterations to run the algorithm at most
    };

    class LabelPropagationRefinement final : public ISerialRefiner {
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

        // per-vertex proposed move for parallel path
        struct ProposedMove {
            partition_t best_id;
            weight_t best_qap_delta;
            weight_t u_weight;
        };
        std::vector<ProposedMove> proposed_moves;

    public:
        weight_t min_improvement = 0;

    public:
        LabelPropagationRefinement() = default;

        ~LabelPropagationRefinement() override = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const u64 t_threads,
                        const u64 seed,
                        const ISerialRefinerConfiguration &i_config) override {
            m_n = t_n;
            m_m = t_m;
            m_k = t_k;
            m_threads = t_threads;

            config = dynamic_cast<const LabelPropagationConfiguration *>(&i_config);

            curr_boundary.initialize(m_n);
            curr_boundary_size = 0;

            blocks.initialize(m_k);
            blocks_qap_delta.initialize(m_k);
            blocks_size = 0;

            rnd_engines.resize(m_threads);
            for (u64 t = 0; t < m_threads; ++t) {
                rnd_engines[t] = RandomEngine(seed + t);
            }

            proposed_moves.resize(m_n);
        }

        void refine(graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph,
                    block_conn_t &block_conn,
                    f64 imbalance,
                    bool uniform_v_weights,
                    bool uniform_e_weights) override {
            if (uniform_v_weights && uniform_e_weights)      refine_impl<true, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, imbalance);
            else if (uniform_v_weights)                      refine_impl<true, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, imbalance);
            else if (uniform_e_weights)                      refine_impl<false, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, imbalance);
            else                                             refine_impl<false, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, imbalance);
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void refine_impl(graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph,
                    block_conn_t &block_conn,
                    f64 imbalance) {
            ScopedTimer _t("refinement", "LabelPropagationRefinement", "refine");
            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            RandomEngine &random_engine = rnd_engines[0];

            bool positive_move_occurred = true;
            for (u64 iteration = 0; iteration < config->max_iteration && positive_move_occurred; ++iteration) {
                positive_move_occurred = false;

                {
                    curr_boundary_size = 0;
                    for (partition_t id = 0; id < bv_manager.get_k(); ++id) {
                        for (size_t i = 0; i < bv_manager.size(id); ++i) { const vertex_t u = bv_manager.get(id, i);
                            {
                                curr_boundary[curr_boundary_size++] = u;
                            }
                        }
                    }
                    fast_shuffle_unchecked(curr_boundary.get_ptr(), curr_boundary.get_ptr() + curr_boundary_size, random_engine.generator);
                }

                if (m_threads > 1) {
                    // Phase 1: compute best moves in parallel (read-only on shared state)
                    #pragma omp parallel for num_threads(m_threads) schedule(dynamic, 64)
                    for (size_t j = 0; j < curr_boundary_size; ++j) {
                        vertex_t u = curr_boundary[j];
                        proposed_moves[j].best_id = NO_ID;

                        if (!bv_manager.is_boundary(u)) { continue; }

                        weight_t u_weight = g.v_weights[u];
                        partition_t u_id = p_manager[u];

                        partition_t best_id = NO_ID;
                        weight_t best_qap_delta = min_improvement;
                        f32 counter = 0;

                        u64 tid = omp_get_thread_num();
                        RandomEngine &rng = rnd_engines[tid];

                        for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) { const partition_t id = block_conn.get_id(i);
                            {
                                if (id == u_id) { continue; }
                                weight_t v_id_weight = p_manager.get_bweight(id);
                                if (v_id_weight + u_weight <= lmax) {
                                    weight_t qap_delta = get_u_qap_delta_t<t_uniform_e_weights>(g, u, u_id, id, p_manager, d_oracle, block_conn);
                                    if (qap_delta > best_qap_delta) {
                                        best_id = id;
                                        best_qap_delta = qap_delta;
                                        counter = 1.0;
                                    } else if (qap_delta == best_qap_delta) {
                                        counter += 1.0;
                                        if (rng.get_f32() < 1.0f / counter) { best_id = id; }
                                    }
                                }
                            }
                        }

                        proposed_moves[j] = {best_id, best_qap_delta, u_weight};
                    }

                    // Phase 2: apply moves sequentially, re-checking balance
                    for (size_t j = 0; j < curr_boundary_size; ++j) {
                        auto [best_id, best_qap_delta, u_weight] = proposed_moves[j];
                        if (best_id == NO_ID) { continue; }

                        vertex_t u = curr_boundary[j];
                        partition_t u_id = p_manager[u];

                        if (u_id == best_id) { continue; }
                        if (p_manager.get_bweight(best_id) + u_weight > lmax) { continue; }

                        if (best_qap_delta >= min_improvement || random_engine.get_f32() < 0.5) {
                            bv_manager.move(g, p_manager, u, u_id, best_id);
                            q_graph.move(g, p_manager, u, u_id, best_id);
                            block_conn.move(g, u, u_id, best_id);
                            p_manager.move(u, u_weight, u_id, best_id);
                            positive_move_occurred |= best_qap_delta > 0;
                        }
                    }
                } else {
                    // Serial path
                    for (size_t j = 0; j < curr_boundary_size; ++j) {
                        vertex_t u = curr_boundary[j];

                        if (!bv_manager.is_boundary(u)) { continue; }

                        weight_t u_weight = g.v_weights[u];
                        partition_t u_id = p_manager[u];

                        partition_t best_id = NO_ID;
                        weight_t best_qap_delta = min_improvement;
                        f32 counter = 0;

                        for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) { const partition_t id = block_conn.get_id(i);
                            {
                                if (id == u_id) { continue; }
                                weight_t v_id_weight = p_manager.get_bweight(id);
                                if (v_id_weight + u_weight <= lmax) {
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
                            }
                        }

                        if (best_id != NO_ID && (best_qap_delta >= min_improvement || random_engine.get_f32() < 0.5)) {
                            bv_manager.move(g, p_manager, u, u_id, best_id);
                            q_graph.move(g, p_manager, u, u_id, best_id);
                            block_conn.move(g, u, u_id, best_id);
                            p_manager.move(u, u_weight, u_id, best_id);
                            positive_move_occurred |= best_qap_delta > 0;
                        }
                    }
                }
            }
        }
    };
}

#endif //HEIPROMAP_LABEL_PROPAGATION_REFINEMENT_H
