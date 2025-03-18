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

#ifndef HEIPROMAP_TWO_VERTEX_LABEL_PROPAGATION_REFINEMENT_H
#define HEIPROMAP_TWO_VERTEX_LABEL_PROPAGATION_REFINEMENT_H

#include <iostream>

#include "../../commons/definitions.h"
#include "../../commons/utils.h"
#include "../interfaces/ISerialRefiner.h"

namespace HeiProMap {
    class TwoVertexLabelPropagationConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit TwoVertexLabelPropagationConfiguration(const std::string &t_name) : ISerialRefinerConfiguration(t_name) {}
        u64 last_n_levels = 2;
        u64 max_iteration = 25; // how many iterations to run the algorithm at most
    };

    class TwoVertexLabelPropagationRefinement final : public ISerialRefiner {
        vertex_t                 m_n    = 0;
        vertex_t                 m_m    = 0;
        partition_t              m_k    = 0;
        weight_t                 m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t>    m_distance;
        u64                      m_seed = 0;

        u32 *vertex_used = nullptr;
        u32 vertex_marker = 0;

        u32 *block_used = nullptr;
        u32 block_marker = 0;

        vertex_t *curr_boundary = nullptr;
        size_t curr_boundary_size = 0;

        partition_t *u_move_ids = nullptr;
        size_t u_move_ids_size = 0;

        partition_t *v_move_ids = nullptr;
        size_t v_move_ids_size = 0;

        RandomEngine                                 *random_engine    = nullptr;
        const TwoVertexLabelPropagationConfiguration *config           = nullptr;
        StatisticCollector                           *m_stat_collector = nullptr;

        METRICS(f64 global_time              = 0;)
        METRICS(f64 global_time_get_boundary = 0.0;)
        METRICS(f64 global_time_iterate      = 0.0;)

        METRICS(s64 global_qap_delta     = 0;)
        METRICS(u64 global_n_pos_moves   = 0;)
        METRICS(u64 global_n_0gain_moves = 0;)

    public:
        TwoVertexLabelPropagationRefinement() = default;

        ~TwoVertexLabelPropagationRefinement() override {
            free(vertex_used);
            free(block_used);
            free(curr_boundary);
            free(u_move_ids);
            free(v_move_ids);
        }

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const weight_t t_lmax,
                        const std::vector<partition_t> &t_hierarchy,
                        const std::vector<weight_t> &t_distance,
                        RandomEngine &t_random_engine,
                        const ISerialRefinerConfiguration &i_config,
                        StatisticCollector &t_stat_collect) override {
            m_n         = t_n;
            m_m         = t_m;
            m_k         = t_k;
            m_lmax      = t_lmax;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;

            random_engine    = &t_random_engine;
            config           = dynamic_cast<const TwoVertexLabelPropagationConfiguration *>(&i_config);
            m_stat_collector = &t_stat_collect;

            vertex_t    m_n_64 = round_up_64(m_n);
            partition_t m_k_64 = round_up_64(m_k);

            vertex_used = (u32 *) aligned_alloc(64, m_n_64 * sizeof(u32));
            std::fill_n(vertex_used, m_n_64, vertex_marker);

            block_used = (u32 *) aligned_alloc(64, m_k_64 * sizeof(u32));
            std::fill_n(block_used, m_k_64, block_marker);

            curr_boundary      = (vertex_t *) aligned_alloc(64, m_n_64 * sizeof(vertex_t));
            curr_boundary_size = 0;

            u_move_ids      = (vertex_t *) aligned_alloc(64, m_k_64 * sizeof(vertex_t));
            u_move_ids_size = 0;

            v_move_ids      = (vertex_t *) aligned_alloc(64, m_k_64 * sizeof(vertex_t));
            v_move_ids_size = 0;
        }

        void refine(const u64 level,
                    const u64 max_level,
                    const graph_t &g,
                    const d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph) override {
            if (level + config->last_n_levels < max_level) {
                return;
            }

            METRICS(std::vector<f64> iteration_time);
            METRICS(std::vector<f64> iteration_time_get_boundary);
            METRICS(std::vector<f64> iteration_time_iterate);

            METRICS(f64 level_time              = 0.0);
            METRICS(f64 level_time_get_boundary = 0.0);
            METRICS(f64 level_time_iterate      = 0.0);

            METRICS(std::vector<s64> iteration_qap_delta);
            METRICS(std::vector<u64> iteration_n_pos_moves);
            METRICS(std::vector<u64> iteration_n_0gain_moves);

            METRICS(s64 level_qap_delta     = 0);
            METRICS(u64 level_n_pos_moves   = 0);
            METRICS(u64 level_n_0gain_moves = 0);

            bool     move_occurred = true;
            for (u64 iteration     = 0; iteration < config->max_iteration && move_occurred; ++iteration) {
                move_occurred = false;

                METRICS(s64  temp_qap_delta     = 0);
                METRICS(u64  temp_n_pos_moves   = 0);
                METRICS(u64  temp_n_0gain_moves = 0);
                METRICS(auto sp_iteration       = std::chrono::high_resolution_clock::now());
                METRICS(auto sp_get_boundary    = std::chrono::high_resolution_clock::now());

                curr_boundary_size = 0;
                forall_bv_iu(bv_manager, i, u)
                    {
                        curr_boundary[curr_boundary_size++] = u;
                    }
                endfor
                std::shuffle(curr_boundary, curr_boundary + curr_boundary_size, random_engine->gen);

                METRICS(auto ep_get_boundary = std::chrono::high_resolution_clock::now());
                METRICS(auto sp_iterate      = std::chrono::high_resolution_clock::now());

                vertex_marker += 1;
                for (size_t i = 0; i < curr_boundary_size; ++i) {
                    vertex_t u = curr_boundary[i];
                    if (vertex_used[u] == vertex_marker) { continue; } // vertex was used
                    if (!bv_manager.is_boundary(u)) { continue; } // vertex is not boundary

                    weight_t    u_weight = g.weight(u);
                    partition_t u_id     = p_manager[u];

                    // get all connected partitions to u
                    block_marker += 1;
                    block_used[u_id] = block_marker;
                    u_move_ids_size = 0;
                    forall_guiv(g, u, i, neighbor)
                        {
                            partition_t neighbor_id = p_manager[neighbor];
                            if (block_used[neighbor_id] != block_marker) {
                                u_move_ids[u_move_ids_size++] = neighbor_id;
                                block_used[neighbor_id]       = block_marker;
                            }
                        }
                    endfor

                    partition_t best_u_move_id = 0;
                    vertex_t    best_v         = 0;
                    partition_t best_v_id      = 0;
                    weight_t    best_v_weight  = 0;
                    partition_t best_v_move_id = 0;
                    s64         best_qap_delta = -1;

                    // check all neighbors v
                    forall_guiv(g, u, i, v)
                        {
                            if (vertex_used[v] == vertex_marker) { continue; } // vertex was used
                            if (!bv_manager.is_boundary(v)) { continue; } // vertex is not boundary

                            weight_t    v_weight = g.weight(v);
                            partition_t v_id     = p_manager[v];

                            // get all connected partitions to v
                            block_marker += 1;
                            block_used[v_id] = block_marker;
                            v_move_ids_size = 0;
                            forall_guiv(g, v, i, neighbor)
                                {
                                    partition_t neighbor_id = p_manager[neighbor];

                                    if (block_used[neighbor_id] != block_marker) {
                                        v_move_ids[v_move_ids_size++] = neighbor_id;
                                        block_used[neighbor_id]       = block_marker;
                                    }
                                }
                            endfor

                            // check if moving u to u_ids and v to v_ids simultaneously would improve the score
                            for (size_t j = 0; j < u_move_ids_size; ++j) {
                                for (size_t l = 0; l < v_move_ids_size; ++l) {
                                    partition_t u_move_id = u_move_ids[j];
                                    partition_t v_move_id = v_move_ids[l];

                                    weight_t u_move_id_weight = p_manager.get_bweight(u_move_id);
                                    weight_t v_move_id_weight = p_manager.get_bweight(v_move_id);

                                    if (u_move_id == v_id && u_move_id_weight + u_weight - v_weight > m_lmax) { continue; }
                                    if (u_move_id != v_id && u_move_id_weight + u_weight > m_lmax) { continue; }
                                    if (v_move_id == u_id && v_move_id_weight + v_weight - u_weight > m_lmax) { continue; }
                                    if (v_move_id != u_id && v_move_id_weight + v_weight > m_lmax) { continue; }
                                    if (u_move_id == v_move_id && u_move_id_weight + u_weight + v_weight > m_lmax) { continue; }

                                    // no overloading is happening, now compute the qap_delta
                                    s64 qap_delta = get_qap_delta(g, u, u_id, u_move_id, v, v_id, v_move_id, p_manager, d_oracle);

                                    if (qap_delta > best_qap_delta) {
                                        best_u_move_id = u_move_id;
                                        best_v         = v;
                                        best_v_id      = v_id;
                                        best_v_weight  = v_weight;
                                        best_v_move_id = v_move_id;
                                        best_qap_delta = qap_delta;
                                    }
                                }
                            }
                        }
                    endfor

                    if (best_qap_delta > 0) {
                        bv_manager.move(g, p_manager, u, u_id, best_u_move_id);
                        q_graph.move(g, p_manager, u, u_id, best_u_move_id);
                        p_manager.move(u, u_weight, u_id, best_u_move_id);

                        bv_manager.move(g, p_manager, best_v, best_v_id, best_v_move_id);
                        q_graph.move(g, p_manager, best_v, best_v_id, best_v_move_id);
                        p_manager.move(best_v, best_v_weight, best_v_id, best_v_move_id);

                        vertex_used[best_v] = vertex_marker;
                        move_occurred = true;

                        METRICS(temp_qap_delta += best_qap_delta * (best_qap_delta > 0));
                        METRICS(temp_n_pos_moves += (best_qap_delta > 0));
                        METRICS(temp_n_0gain_moves += (best_qap_delta == 0));
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

#endif //HEIPROMAP_TWO_VERTEX_LABEL_PROPAGATION_REFINEMENT_H
