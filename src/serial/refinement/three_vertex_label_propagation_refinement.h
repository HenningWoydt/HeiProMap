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

#ifndef HEIPROMAP_THREE_VERTEX_LABEL_PROPAGATION_REFINEMENT_H
#define HEIPROMAP_THREE_VERTEX_LABEL_PROPAGATION_REFINEMENT_H

#include <iostream>

#include "../../commons/definitions.h"
#include "../../commons/utils.h"
#include "../interfaces/ISerialRefiner.h"

namespace HeiProMap {
    class ThreeVertexLabelPropagationConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit ThreeVertexLabelPropagationConfiguration(const std::string& t_name) : ISerialRefinerConfiguration(t_name) {}
        u64 last_n_levels = 2;
        u64 max_iteration = 25; // how many iterations to run the algorithm at most
    };

    class ThreeVertexLabelPropagationRefinement final : public ISerialRefiner {
        vertex_t m_n    = 0;
        vertex_t m_m    = 0;
        partition_t m_k = 0;
        f64 m_imbalance = 0.0;
        weight_t m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;
        u64 m_seed = 0;

        u32* vertex_used  = nullptr;
        u32 vertex_marker = 0;

        u32* block_used  = nullptr;
        u32 block_marker = 0;

        vertex_t* curr_boundary   = nullptr;
        size_t curr_boundary_size = 0;

        partition_t* u_move_ids = nullptr;
        size_t u_move_ids_size  = 0;

        partition_t* v_move_ids = nullptr;
        size_t v_move_ids_size  = 0;

        RandomEngine* random_engine                            = nullptr;
        const ThreeVertexLabelPropagationConfiguration* config = nullptr;
        StatisticCollector* m_stat_collector                   = nullptr;

        METRICS(f64 global_time = 0;)
        METRICS(f64 global_time_get_boundary = 0.0;)
        METRICS(f64 global_time_iterate = 0.0;)

        METRICS(s64 global_qap_delta = 0;)
        METRICS(u64 global_n_pos_moves = 0;)
        METRICS(u64 global_n_0gain_moves = 0;)

    public:
        ThreeVertexLabelPropagationRefinement() = default;

        ~ThreeVertexLabelPropagationRefinement() override {
            free(vertex_used);
            free(block_used);
            free(curr_boundary);
            free(u_move_ids);
            free(v_move_ids);
        }

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const f64 t_imbalance,
                        const weight_t t_lmax,
                        const std::vector<partition_t>& t_hierarchy,
                        const std::vector<weight_t>& t_distance,
                        RandomEngine& t_random_engine,
                        const ISerialRefinerConfiguration& i_config,
                        StatisticCollector& t_stat_collect) override {
            m_n         = t_n;
            m_m         = t_m;
            m_k         = t_k;
            m_lmax      = t_lmax;
            m_imbalance = t_imbalance;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;

            random_engine    = &t_random_engine;
            config           = dynamic_cast<const ThreeVertexLabelPropagationConfiguration*>(&i_config);
            m_stat_collector = &t_stat_collect;

            vertex_t m_n_64    = round_up_64(m_n);
            partition_t m_k_64 = round_up_64(m_k);

            vertex_used = (u32*)aligned_alloc(64, m_n_64 * sizeof(u32));
            std::fill_n(vertex_used, m_n_64, vertex_marker);

            block_used = (u32*)aligned_alloc(64, m_k_64 * sizeof(u32));
            std::fill_n(block_used, m_k_64, block_marker);

            curr_boundary      = (vertex_t*)aligned_alloc(64, m_n_64 * sizeof(vertex_t));
            curr_boundary_size = 0;

            u_move_ids      = (vertex_t*)aligned_alloc(64, m_k_64 * sizeof(vertex_t));
            u_move_ids_size = 0;

            v_move_ids      = (vertex_t*)aligned_alloc(64, m_k_64 * sizeof(vertex_t));
            v_move_ids_size = 0;
        }

        void refine(const u64 level,
                    const u64 max_level,
                    const graph_t& g,
                    d_oracle_t& d_oracle,
                    bv_manager_t& bv_manager,
                    p_manager_t& p_manager,
                    q_graph_t& q_graph) override {
            if (level + config->last_n_levels < max_level) {
                return;
            }

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

            bool move_occurred = true;
            for (u64 iteration = 0; iteration < config->max_iteration && move_occurred; ++iteration) {
                move_occurred = false;

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
                std::shuffle(curr_boundary, curr_boundary + curr_boundary_size, random_engine->gen);

                METRICS(auto ep_get_boundary = std::chrono::high_resolution_clock::now());
                METRICS(auto sp_iterate = std::chrono::high_resolution_clock::now());

                vertex_marker += 1;
                for (size_t i = 0; i < curr_boundary_size; ++i) {
                    vertex_t u = curr_boundary[i];
                    if (vertex_used[u] == vertex_marker) { continue; }
                    if (!bv_manager.is_boundary(u)) { continue; }

                    // collect all vertices u, uu, uuu
                    std::vector<vertex_t> vertices;
                    vertices.push_back(u);

                    std::vector<u8> inserted(g.get_n(), 0);
                    inserted[u] = 1;

                    forall_guiv(g, u, j, uu)
                        {
                            if (uu >= u) { continue; }
                            if (inserted[uu] == 0) {
                                vertices.push_back(uu);
                                inserted[uu] = 1;
                            }
                            forall_guiv(g, uu, k, uuu)
                                {
                                    if (uuu >= uu || uuu >= u) { continue; }
                                    if (inserted[uuu] == 0) {
                                        vertices.push_back(uuu);
                                        inserted[uuu] = 1;
                                    }
                                }
                            endfor
                        }
                    endfor

                    // collect for all vertices the connected partitions
                    std::vector<u8> id_inserted(m_k);
                    std::vector<std::vector<partition_t>> partitions(vertices.size());

                    for (size_t j = 0; j < vertices.size(); ++j) {
                        vertex_t v       = vertices[j];
                        partition_t v_id = p_manager[v];
                        std::fill(id_inserted.begin(), id_inserted.end(), 0);

                        forall_guiv(g, v, k, vv)
                            {
                                partition_t vv_id = p_manager[vv];
                                if (v_id == vv_id || id_inserted[vv_id] == 1) { continue; }

                                partitions[j].push_back(vv_id);
                                id_inserted[vv_id] = 1;
                            }
                        endfor
                    }

                    vertex_t best_v, best_vv, best_vvv;
                    weight_t best_v_weight, best_vv_weight, best_vvv_weight;
                    partition_t best_v_id, best_vv_id, best_vvv_id;
                    partition_t best_new_v_id, best_new_vv_id, best_new_vvv_id;
                    s64 best_qap_delta = -1;

                    // search for the best triple combination that does not overload
                    for (size_t j = 0; j < vertices.size(); ++j) {
                        for (size_t k = j + 1; k < vertices.size(); ++k) {
                            for (size_t l = k + 1; l < vertices.size(); ++l) {
                                // get the three vertices
                                vertex_t v   = vertices[j];
                                vertex_t vv  = vertices[k];
                                vertex_t vvv = vertices[l];

                                weight_t v_weight   = g.weight(v);
                                weight_t vv_weight  = g.weight(vv);
                                weight_t vvv_weight = g.weight(vvv);

                                partition_t v_id   = p_manager[v];
                                partition_t vv_id  = p_manager[vv];
                                partition_t vvv_id = p_manager[vvv];

                                for (partition_t new_v_id : partitions[j]) {
                                    for (partition_t new_vv_id : partitions[k]) {
                                        for (partition_t new_vvv_id : partitions[l]) {
                                            // make temporary moves
                                            p_manager.move(v, v_weight, v_id, new_v_id);
                                            p_manager.move(vv, vv_weight, vv_id, new_vv_id);
                                            p_manager.move(vvv, vvv_weight, vvv_id, new_vvv_id);

                                            bool overloaded = p_manager.get_bweight(new_v_id) > m_lmax || p_manager.get_bweight(new_vv_id) > m_lmax || p_manager.get_bweight(new_vvv_id) > m_lmax;

                                            // revert the moves
                                            p_manager.move(v, v_weight, new_v_id, v_id);
                                            p_manager.move(vv, vv_weight, new_vv_id, vv_id);
                                            p_manager.move(vvv, vvv_weight, new_vvv_id, vvv_id);

                                            if (overloaded) { continue; }

                                            // get the qap delta
                                            s64 qap_delta = get_qap_delta(g, v, vv, vvv, v_id, vv_id, vvv_id, new_v_id, new_vv_id, new_vvv_id, p_manager, d_oracle);

                                            if (qap_delta > best_qap_delta) {
                                                best_qap_delta = qap_delta;

                                                best_v = v;
                                                best_vv = vv;
                                                best_vvv = vvv;

                                                best_v_weight = v_weight;
                                                best_vv_weight = vv_weight;
                                                best_vvv_weight = vvv_weight;

                                                best_v_id = v_id;
                                                best_vv_id = vv_id;
                                                best_vvv_id = vvv_id;

                                                best_new_v_id = new_v_id;
                                                best_new_vv_id = new_vv_id;
                                                best_new_vvv_id = new_vvv_id;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    if (best_qap_delta > 0) {
                        bv_manager.move(g, p_manager, best_v, best_v_id, best_new_v_id);
                        q_graph.move(g, p_manager, best_v, best_v_id, best_new_v_id);
                        p_manager.move(best_v, best_v_weight, best_v_id, best_new_v_id);

                        bv_manager.move(g, p_manager, best_vv, best_vv_id, best_new_vv_id);
                        q_graph.move(g, p_manager, best_vv, best_vv_id, best_new_vv_id);
                        p_manager.move(best_vv, best_vv_weight, best_vv_id, best_new_vv_id);

                        bv_manager.move(g, p_manager, best_vvv, best_vvv_id, best_new_vvv_id);
                        q_graph.move(g, p_manager, best_vvv, best_vvv_id, best_new_vvv_id);
                        p_manager.move(best_vvv, best_vvv_weight, best_vvv_id, best_new_vvv_id);

                        vertex_used[best_v] = vertex_marker;
                        vertex_used[best_vv] = vertex_marker;
                        vertex_used[best_vvv] = vertex_marker;
                        move_occurred = true;

                        METRICS(temp_qap_delta += best_qap_delta * (best_qap_delta > 0));
                        METRICS(temp_n_pos_moves += (best_qap_delta > 0));
                        METRICS(temp_n_0gain_moves += (best_qap_delta == 0));
                    }
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

#endif //HEIPROMAP_THREE_VERTEX_LABEL_PROPAGATION_REFINEMENT_H
