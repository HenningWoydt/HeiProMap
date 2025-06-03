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

#include "../../commons/definitions.h"
#include "../../commons/JSON_utils.h"
#include "../../commons/random_engine.h"
#include "../../commons/statistic_collector.h"

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

        AlignedArray <u32> vertex_used;
        u32                vertex_marker = 0;

        AlignedArray <u32> block_used;
        u32                block_marker = 0;

        AlignedArray <vertex_t> curr_boundary;
        size_t                  curr_boundary_size = 0;

        AlignedArray <partition_t> blocks;
        AlignedArray <s64>         blocks_qap_delta;
        size_t                     blocks_size = 0;

        RandomEngine                        *random_engine    = nullptr;
        const LabelPropagationConfiguration *config           = nullptr;
        StatisticCollector                  *m_stat_collector = nullptr;

        METRICS(f64 global_time              = 0;)
        METRICS(f64 global_time_get_boundary = 0.0;)
        METRICS(f64 global_time_iterate      = 0.0;)

        METRICS(s64 global_qap_delta     = 0;)
        METRICS(u64 global_n_pos_moves   = 0;)
        METRICS(u64 global_n_0gain_moves = 0;)

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
                        const ISerialRefinerConfiguration &i_config,
                        StatisticCollector &t_stat_collect) override {
            m_n         = t_n;
            m_m         = t_m;
            m_k         = t_k;
            m_imbalance = t_imbalance;
            m_lmax      = t_lmax;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;

            random_engine    = &t_random_engine;
            config           = dynamic_cast<const LabelPropagationConfiguration *>(&i_config);
            m_stat_collector = &t_stat_collect;

            vertex_used.initialize(m_n, 0);
            block_used.initialize(m_m, 0);

            curr_boundary.initialize(m_n);
            curr_boundary_size = 0;

            blocks.initialize(m_k);
            blocks_qap_delta.initialize(m_k);
            blocks_size = 0;
        }

        void refine(const u64 level,
                    const u64 max_level,
                    graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph) override {
            METRICS(std::vector<f64> iteration_time);
            METRICS(std::vector<f64> iteration_time_get_boundary);
            METRICS(std::vector<f64> iteration_time_iterate);

            METRICS(f64 level_time = 0.0);
            METRICS(f64 level_time_get_boundary = 0.0);
            METRICS(f64 level_time_iterate = 0.0);

            METRICS(std::vector<s64> iteration_qap_delta);
            METRICS(std::vector<u64> iteration_n_pos_moves);
            METRICS(std::vector<u64> iteration_n_0gain_moves);

            METRICS(s64 level_qap_delta = 0);
            METRICS(u64 level_n_pos_moves = 0);
            METRICS(u64 level_n_0gain_moves = 0);

            bool     positive_move_occurred = true;
            for (u64 iteration              = 0; iteration < config->max_iteration && positive_move_occurred; ++iteration) {
                positive_move_occurred = false;

                METRICS(s64 temp_qap_delta = 0);
                METRICS(u64 temp_n_pos_moves = 0);
                METRICS(u64 temp_n_0gain_moves = 0);
                METRICS(auto sp_iteration = std::chrono::high_resolution_clock::now());
                METRICS(auto sp_get_boundary = std::chrono::high_resolution_clock::now());

                curr_boundary_size = 0;
                forall_bv_iu(bv_manager, i, u)
                    {
                        curr_boundary[curr_boundary_size++] = u;
                    }
                endfor
                std::shuffle(curr_boundary.get_ptr(), curr_boundary.get_ptr() + curr_boundary_size, random_engine->gen);

                METRICS(auto ep_get_boundary = std::chrono::high_resolution_clock::now());
                METRICS(auto sp_iterate = std::chrono::high_resolution_clock::now());

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

                            METRICS(temp_qap_delta += best_qap_delta * (best_qap_delta > 0));
                            METRICS(temp_n_pos_moves += (best_qap_delta > 0));
                            METRICS(temp_n_0gain_moves += (best_qap_delta == 0));
                        }
                    }
                    vertex_used[u] = vertex_marker;
                }

#if COLLECT_METRICS
                auto ep_iterate   = std::chrono::high_resolution_clock::now();
                auto ep_iteration = std::chrono::high_resolution_clock::now();

                f64 t_iteration    = get_seconds(sp_iteration, ep_iteration);
                f64 t_get_boundary = get_seconds(sp_get_boundary, ep_get_boundary);
                f64 t_iterate      = get_seconds(sp_iterate, ep_iterate);

                iteration_time.push_back(t_iteration);
                iteration_time_get_boundary.push_back(t_get_boundary);
                iteration_time_iterate.push_back(t_iterate);

                level_time += t_iteration;
                level_time_get_boundary += t_get_boundary;
                level_time_iterate += t_iterate;

                global_time += t_iteration;
                global_time_get_boundary += t_get_boundary;
                global_time_iterate += t_iterate;

                temp_qap_delta *= 2;

                iteration_qap_delta.push_back(temp_qap_delta);
                iteration_n_pos_moves.push_back(temp_n_pos_moves);
                iteration_n_0gain_moves.push_back(temp_n_0gain_moves);

                level_qap_delta += temp_qap_delta;
                level_n_pos_moves += temp_n_pos_moves;
                level_n_0gain_moves += temp_n_0gain_moves;

                global_qap_delta += temp_qap_delta;
                global_n_pos_moves += temp_n_pos_moves;
                global_n_0gain_moves += temp_n_0gain_moves;
#endif
            }
        }

        JSONString get_stats() override { return {}; };
    };
}

#endif //HEIPROMAP_LABEL_PROPAGATION_REFINEMENT_H
