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

        AlignedArray<vertex_t> curr_boundary;
        size_t curr_boundary_size = 0;

        AlignedArray<partition_t> blocks;
        AlignedArray<s64> blocks_qap_delta;
        size_t blocks_size = 0;

        RandomEngine random_engine = RandomEngine(0);
        const LabelPropagationConfiguration *config = nullptr;

    public:
        weight_t min_improvement = 0;

    public:
        LabelPropagationRefinement() = default;

        ~LabelPropagationRefinement() override = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const ISerialRefinerConfiguration &i_config) override {
            m_n = t_n;
            m_m = t_m;
            m_k = t_k;

            config = dynamic_cast<const LabelPropagationConfiguration *>(&i_config);

            curr_boundary.initialize(m_n);
            curr_boundary_size = 0;

            blocks.initialize(m_k);
            blocks_qap_delta.initialize(m_k);
            blocks_size = 0;
        }

        void refine(graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph,
                    block_conn_t &block_conn,
                    f64 imbalance) override {
            ScopedTimer _t("refinement", "LabelPropagationRefinement", "refine");
            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            bool positive_move_occurred = true;
            for (u64 iteration = 0; iteration < config->max_iteration && positive_move_occurred; ++iteration) {
                positive_move_occurred = false;

                {
                    // ScopedTimer _t("refinement", "LabelPropagationRefinement", "get_boundary");
                    curr_boundary_size = 0;
                    for (partition_t id = 0; id < bv_manager.get_k(); ++id) {
                        forall_bv_id_iu(bv_manager, id, i, u)
                            {
                                curr_boundary[curr_boundary_size++] = u;
                            }
                        endfor
                    }
                    fast_shuffle_unchecked(curr_boundary.get_ptr(), curr_boundary.get_ptr() + curr_boundary_size, random_engine.generator);
                    // std::shuffle(curr_boundary.get_ptr(), curr_boundary.get_ptr() + curr_boundary_size, random_engine.generator);
                }

                for (size_t j = 0; j < curr_boundary_size; ++j) {
                    vertex_t u = curr_boundary[j];

                    if (!bv_manager.is_boundary(u)) { continue; }

                    weight_t u_weight = g.v_weights[u];
                    partition_t u_id = p_manager[u];

                    // make the move that reduces qap the most
                    partition_t best_id = NO_ID;
                    weight_t best_qap_delta = min_improvement;
                    f32 counter = 0;

                    {
                        // ScopedTimer _t("refinement", "LabelPropagationRefinement", "find_move");

                        forall_bc_ui_id(block_conn, u, i, id)
                            {
                                if (id == u_id) { continue; }

                                weight_t v_id_weight = p_manager.get_bweight(id);

                                if (v_id_weight + u_weight <= lmax) {
                                    s64 qap_delta = get_u_qap_delta(g, u, u_id, id, p_manager, d_oracle, block_conn);

                                    if (qap_delta > best_qap_delta) {
                                        best_id = id;
                                        best_qap_delta = qap_delta;
                                        counter = 1.0;
                                    } else if (qap_delta == best_qap_delta) {
                                        counter += 1.0;
                                        // choose with probability 1/counter as it ensures uniform distribution
                                        if (random_engine.get_f32() < 1.0f / counter) { best_id = id; }
                                    }
                                }
                            }
                        endfor
                    }

                    if (best_id != NO_ID && (best_qap_delta >= min_improvement || random_engine.get_f32() < 0.5)) {
                        // choose if positive, if 0-gain choose 50% of the time
                        // ScopedTimer _t("refinement", "LabelPropagationRefinement", "make_move");

                        bv_manager.move(g, p_manager, u, u_id, best_id);
                        q_graph.move(g, p_manager, u, u_id, best_id);
                        block_conn.move(g, u, u_id, best_id);
                        p_manager.move(u, u_weight, u_id, best_id);
                        positive_move_occurred |= best_qap_delta > 0;
                    }
                }
            }
        }
    };
}

#endif //HEIPROMAP_LABEL_PROPAGATION_REFINEMENT_H
