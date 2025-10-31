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

#include "../definitions.h"
#include "../utility/JSON_utils.h"
#include "../utility/random_engine.h"

namespace HeiProMap {
    class LabelPropagationConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit LabelPropagationConfiguration(const std::string &t_name) : ISerialRefinerConfiguration(t_name) {}

        u64 max_iteration = 25; // how many iterations to run the algorithm at most
    };

    class LabelPropagationRefinement final : public ISerialRefiner {
        vertex_t                 m_n         = 0;
        vertex_t                 m_m         = 0;
        partition_t              m_k         = 0;
        f64                      m_imbalance = 0.0;
        weight_t                 m_lmax      = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t>    m_distance;

        AlignedArray<u32> vertex_used;
        u32               vertex_marker = 0;

        AlignedArray<u32> block_used;
        u32               block_marker = 0;

        AlignedArray<vertex_t> curr_boundary;
        size_t                 curr_boundary_size = 0;

        AlignedArray<partition_t> blocks;
        AlignedArray<s64>         blocks_qap_delta;
        size_t                    blocks_size = 0;

        RandomEngine                        *random_engine = nullptr;
        const LabelPropagationConfiguration *config        = nullptr;

    public:
        LabelPropagationRefinement() = default;

        ~LabelPropagationRefinement() override = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const f64 t_imbalance,
                        const weight_t t_lmax,
                        const std::vector<partition_t> &t_hierarchy,
                        const std::vector<weight_t> &t_distance,
                        RandomEngine &t_random_engine,
                        const ISerialRefinerConfiguration &i_config) override {
            m_n         = t_n;
            m_m         = t_m;
            m_k         = t_k;
            m_imbalance = t_imbalance;
            m_lmax      = t_lmax;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;

            random_engine = &t_random_engine;
            config        = dynamic_cast<const LabelPropagationConfiguration *>(&i_config);

            vertex_used.initialize(m_n, 0);
            block_used.initialize(m_m, 0);

            curr_boundary.initialize(m_n);
            curr_boundary_size = 0;

            blocks.initialize(m_k);
            blocks_qap_delta.initialize(m_k);
            blocks_size = 0;
        }

        void refine([[maybe_unused]] const u64 level,
                    [[maybe_unused]] const u64 max_level,
                    graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph) override {
            bool     positive_move_occurred = true;
            for (u64 iteration              = 0; iteration < config->max_iteration && positive_move_occurred; ++iteration) {
                positive_move_occurred = false;

                curr_boundary_size = 0;
                for (partition_t id = 0; id < bv_manager.get_k(); ++id) {
                    forall_bv_id_iu(bv_manager, id, i, u)
                    {
                        curr_boundary[curr_boundary_size++] = u;
                    }
                    endfor
                }
                std::shuffle(curr_boundary.get_ptr(), curr_boundary.get_ptr() + curr_boundary_size, random_engine->generator);

                vertex_marker += 1;
                for (size_t j = 0; j < curr_boundary_size; ++j) {
                    vertex_t u = curr_boundary[j];

                    if (vertex_used[u] == vertex_marker) { continue; }
                    if (!bv_manager.is_boundary(u)) { continue; }

                    weight_t    u_weight = g.weight(u);
                    partition_t u_id     = p_manager[u];

                    // make the move that reduces qap the most
                    partition_t best_id        = u_id;
                    weight_t    best_id_weight = 0;
                    s64         best_qap_delta = -1;
                    f32         counter        = 0;

                    block_marker += 1;
                    block_used[u_id] = block_marker;
                    forall_guiv(g, u, i, v)
                        {
                            partition_t v_id        = p_manager[v];
                            weight_t    v_id_weight = p_manager.get_bweight(v_id);

                            if (block_used[v_id] != block_marker) {
                                if (v_id_weight + u_weight <= m_lmax) {
                                    s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);

                                    if (qap_delta > best_qap_delta || (qap_delta == best_qap_delta && v_id_weight < best_id_weight)) {
                                        best_id        = v_id;
                                        best_id_weight = v_id_weight;
                                        best_qap_delta = qap_delta;
                                        counter        = 1.0;
                                    } else if (qap_delta == best_qap_delta && qap_delta != -1) {
                                        counter += 1.0;
                                        // choose with probability 1/counter as it ensures uniform distribution
                                        if (random_engine->get_f32() < 1.0f / counter) {
                                            best_id = v_id;
                                        }
                                    }
                                }
                                block_used[v_id] = block_marker;
                            }
                        }
                    endfor

                    if (best_id != u_id) {
                        // choose if positive, if 0-gain choose 50% of the time
                        if (best_qap_delta > 0 || random_engine->get_f32() < 0.5) {
                            bv_manager.move(g, p_manager, u, u_id, best_id);
                            q_graph.move(g, p_manager, u, u_id, best_id);
                            p_manager.move(u, u_weight, u_id, best_id);
                            positive_move_occurred |= best_qap_delta > 0;
                        }
                    }
                    vertex_used[u] = vertex_marker;
                }
            }
        }

        void refine_layer([[maybe_unused]] const u64 level,
                          [[maybe_unused]] const u64 max_level,
                          [[maybe_unused]] graph_t &g,
                          [[maybe_unused]] d_oracle_t &d_oracle,
                          [[maybe_unused]] bv_manager_t &bv_manager,
                          [[maybe_unused]] p_manager_t &p_manager,
                          [[maybe_unused]] q_graph_t &q_graph,
                          [[maybe_unused]] size_t layer) override {}
    };
}

#endif //HEIPROMAP_LABEL_PROPAGATION_REFINEMENT_H
