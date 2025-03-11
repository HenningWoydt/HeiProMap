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

#ifndef HEIPROMAP_MULTI_TRY_FM_REFINEMENT_H
#define HEIPROMAP_MULTI_TRY_FM_REFINEMENT_H

#include <algorithm>

#include "k_way_fm_refinement_Faraj20.h"
#include "../datastructures/functions.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/qap.h"
#include "../../commons/utils.h"

namespace HeiProMap {
    class MultiTryFmRefinementConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit MultiTryFmRefinementConfiguration(const std::string &t_name) : ISerialRefinerConfiguration(t_name) {}

        u64 max_iteration = 1;
        f64 alpha         = 10.0;
        f64 beta          = 1.0;
    };

    /*
    class MultiTryFMRefinement final : public ISerialRefiner {
    private:
        vertex_t                 m_n    = 0;
        vertex_t                 m_m    = 0;
        partition_t              m_k    = 0;
        weight_t                 m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t>    m_distance;

        u32 *vertex_used = nullptr;
        u32 vertex_mark  = 0;

        u32 *block_used = nullptr;
        u32 block_mark  = 0;

        vertex_t *curr_boundary     = nullptr;
        size_t   curr_boundary_size = 0;

        KWayFMMove *moves     = nullptr;
        size_t     moves_size = 0;

        // priority queues
        IndexedMaxHeap<KWayFMMove> boundary_vertices;

        RandomEngine                            *random_engine    = nullptr;
        const MultiTryFmRefinementConfiguration *config           = nullptr;
        StatisticCollector                      *m_stat_collector = nullptr;

        METRICS(std::vector<std::vector<f64>> iteration_time;)
        METRICS(std::vector<std::vector<f64>> iteration_time_get_boundary;)
        METRICS(std::vector<std::vector<f64>> iteration_time_initialize;)
        METRICS(std::vector<std::vector<f64>> iteration_time_queue;)
        METRICS(std::vector<std::vector<f64>> iteration_time_moves;)
        METRICS(std::vector<std::vector<s64>> iteration_qap_delta;)

    public:
        MultiTryFMRefinement() = default;

        ~MultiTryFMRefinement() override {
            free(vertex_used);
            free(block_used);
            free(curr_boundary);
            free(moves);
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
            config           = dynamic_cast<const MultiTryFmRefinementConfiguration *>(&i_config);
            m_stat_collector = &t_stat_collect;

            vertex_t    m_n_64 = round_up_64(m_n);
            partition_t m_k_64 = round_up_64(m_k);

            vertex_used = (u32 *) aligned_alloc(64, sizeof(u32) * m_n_64);
            block_used  = (u32 *) aligned_alloc(64, sizeof(u32) * m_k_64);

            curr_boundary      = (vertex_t *) aligned_alloc(64, sizeof(vertex_t) * m_n_64);
            curr_boundary_size = 0;

            moves      = (KWayFMMove *) aligned_alloc(64, sizeof(KWayFMMove) * m_n_64);
            moves_size = 0;

            boundary_vertices.initialize(m_n);
        }

        void refine(const u64 level,
                    const graph_t &g,
                    const av_manager_t &av_manager,
                    const d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph) override {
            METRICS(std::vector<f64> time;)
            METRICS(std::vector<f64> time_get_boundary;)
            METRICS(std::vector<f64> time_initialize;)
            METRICS(std::vector<f64> time_queue;)
            METRICS(std::vector<f64> time_moves;)
            METRICS(std::vector<s64> qap_delta;)

            f64 alpha = config->alpha;
            f64 beta  = std::log(av_manager.size());

            bool     positive_move_occurred = true;
            for (u64 iteration              = 0; iteration < config->max_iteration && positive_move_occurred; ++iteration) {
                METRICS(time.push_back(0.0);)
                METRICS(time_get_boundary.push_back(0.0);)
                METRICS(time_initialize.push_back(0.0);)
                METRICS(time_queue.push_back(0.0);)
                METRICS(time_moves.push_back(0.0);)
                METRICS(qap_delta.push_back(0);)

                METRICS_TIME(sp)

                positive_move_occurred = false;

                METRICS_TIME(sp_get_boundary)
                curr_boundary_size = 0;
                forall_bv_iu(bv_manager, i, u)
                    {
                        curr_boundary[curr_boundary_size++] = u;
                    }
                endfor
                std::shuffle(curr_boundary, curr_boundary + curr_boundary_size, random_engine->gen);

                METRICS_TIME(ep_get_boundary)

                METRICS(time_get_boundary.back() += get_seconds(sp_get_boundary, ep_get_boundary);)

                vertex_mark += 1;
                for (size_t ii = 0; ii < curr_boundary_size; ++ii) {
                    vertex_t u = curr_boundary[ii];
                    if (vertex_used[u] == vertex_mark) { continue; }

                    METRICS_TIME(sp_initialize)

                    boundary_vertices.clear();

                    // insert u into the priority queue
                    partition_t u_id     = p_manager[u];
                    weight_t    u_weight = g.get_weight(u);

                    // find all connected partitions to u
                    block_mark += 1;
                    partition_t best_id;
                    s64         best_qap_delta = -std::numeric_limits<s64>::max();
                    forall_guiv(g, u, i, v)
                        {
                            partition_t v_id = p_manager[v];
                            if (v_id != u_id && block_used[v_id] != block_mark && p_manager.get_bweight(v_id) + u_weight <= m_lmax) {
                                s64 u_qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                                block_used[v_id] = block_mark;

                                if (u_qap_delta > best_qap_delta) {
                                    best_qap_delta = u_qap_delta;
                                    best_id        = v_id;
                                }
                            }
                        }
                    endfor
                    if (best_qap_delta != -std::numeric_limits<s64>::max()) {
                        boundary_vertices.push(u, KWayFMMove(u, u_id, best_id, best_qap_delta));
                    }

                    // insert all neighbors of u that are boundary into the queue
                    forall_guiv(g, u, i, neighbor)
                        {
                            if (vertex_used[neighbor] == vertex_mark) { continue; }
                            if (!bv_manager.is_boundary(neighbor)) { continue; }

                            partition_t neighbor_id     = p_manager[neighbor];
                            weight_t    neighbor_weight = g.get_weight(neighbor);

                            block_mark += 1;
                            best_qap_delta              = -std::numeric_limits<s64>::max();
                            forall_guiv(g, neighbor, j, v)
                                {
                                    partition_t v_id = p_manager[v];
                                    if (v_id != neighbor_id && block_used[v_id] != block_mark && p_manager.get_bweight(v_id) + neighbor_weight <= m_lmax) {
                                        s64 u_qap_delta = get_u_qap_delta(g, neighbor, neighbor_id, v_id, p_manager, d_oracle);
                                        block_used[v_id] = block_mark;

                                        if (u_qap_delta > best_qap_delta) {
                                            best_qap_delta = u_qap_delta;
                                            best_id        = v_id;
                                        }
                                    }
                                }
                            endfor
                            if (best_qap_delta != -std::numeric_limits<s64>::max()) {
                                boundary_vertices.push(neighbor, KWayFMMove(neighbor, neighbor_id, best_id, best_qap_delta));
                            }
                        }
                    endfor

                    METRICS_TIME(ep_initialize)
                    METRICS(time_initialize.back() += get_seconds(sp_initialize, ep_initialize);)

                    // process the queue
                    moves_size = 0;
                    size_t best_idx      = 0;
                    s64    curr_qap_gain = 0;
                    s64    max_qap_gain  = 0;

                    f64         steps_since_last_improvement = 0.0;
                    f64         qap_gain_mean                = 0.0;
                    f64         qap_gain_var                 = 0.0;

                    METRICS_TIME(sp_queue)

                    while (!boundary_vertices.empty()) {
                        const KWayFMMove move = boundary_vertices.top();
                        boundary_vertices.pop();

                        vertex_t    vertex        = move.u;
                        partition_t vertex_id     = p_manager[vertex];
                        weight_t    vertex_weight = g.get_weight(vertex);
                        partition_t move_id       = move.to_move_id;

                        if (!is_connected_to(g, p_manager, vertex, move_id)) { continue; }
                        if (p_manager.get_bweight(move_id) + vertex_weight > m_lmax) { continue; }

                        moves[moves_size++] = move;
                        curr_qap_gain += move.qap_delta;
                        if (curr_qap_gain >= max_qap_gain) {
                            best_idx     = moves_size;
                            max_qap_gain = curr_qap_gain;

                            steps_since_last_improvement = 0.0;
                            qap_gain_mean                = 0.0;
                            qap_gain_var                 = 0.0;
                        }

                        // make move in structures
                        p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                        vertex_used[vertex] = vertex_mark;

                        steps_since_last_improvement += 1.0;
                        f64 new_qap_gain_mean = qap_gain_mean + ((f64) move.qap_delta - qap_gain_mean) / steps_since_last_improvement;
                        f64 new_qap_gain_var  = (qap_gain_var + ((f64) move.qap_delta - qap_gain_mean) * ((f64) move.qap_delta - new_qap_gain_mean)) / steps_since_last_improvement;

                        qap_gain_mean = new_qap_gain_mean;
                        qap_gain_var  = new_qap_gain_var;

                        if (steps_since_last_improvement > 2.0 && steps_since_last_improvement * qap_gain_mean * qap_gain_mean > alpha * qap_gain_var + beta) { break; }

                        // we have to push or update the neighbors that were not moved already
                        forall_guiv(g, vertex, i, neighbor)
                            {
                                if (vertex_used[neighbor] == vertex_mark) { continue; }
                                if (!is_boundary(g, p_manager, neighbor)) { continue; }

                                partition_t neighbor_id     = p_manager[neighbor];
                                weight_t    neighbor_weight = g.get_weight(neighbor);

                                block_mark += 1;
                                best_qap_delta              = -std::numeric_limits<s64>::max();
                                forall_guiv(g, neighbor, j, v)
                                    {
                                        partition_t v_id = p_manager[v];
                                        if (v_id != neighbor_id && block_used[v_id] != block_mark && p_manager.get_bweight(v_id) + neighbor_weight <= m_lmax) {
                                            s64 v_qap_delta = get_u_qap_delta(g, neighbor, neighbor_id, v_id, p_manager, d_oracle);
                                            block_used[v_id] = block_mark;

                                            if (v_qap_delta > best_qap_delta) {
                                                best_qap_delta = v_qap_delta;
                                                best_id        = v_id;
                                            }
                                        }
                                    }
                                endfor
                                if (best_qap_delta != -std::numeric_limits<s64>::max()) {
                                    boundary_vertices.push_update(neighbor, KWayFMMove(neighbor, neighbor_id, best_id, best_qap_delta));
                                }
                            }
                        endfor
                    }

                    METRICS_TIME(ep_queue)
                    METRICS(time_queue.back() += get_seconds(sp_queue, ep_queue);)

                    METRICS_TIME(sp_moves)
                    // revert all moves in partitioning manager
                    for (size_t i                            = 0; i < moves_size; i++) {
                        vertex_t    vertex        = moves[moves_size - 1 - i].u;
                        weight_t    vertex_weight = g.get_weight(vertex);
                        partition_t vertex_id     = moves[moves_size - 1 - i].to_move_id;
                        partition_t move_id       = moves[moves_size - 1 - i].u_id;

                        p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                    }

                    // make all moves to best index
                    for (size_t i = 0; i < best_idx; ++i) {
                        vertex_t    vertex        = moves[i].u;
                        weight_t    vertex_weight = g.get_weight(vertex);
                        partition_t vertex_id     = moves[i].u_id;
                        partition_t move_id       = moves[i].to_move_id;

                        bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                        q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                        p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                    }

                    positive_move_occurred |= max_qap_gain > 0;

                    METRICS_TIME(ep_moves)
                    METRICS(time_moves.back() += get_seconds(sp_moves, ep_moves);)
                    METRICS(qap_delta.back() += max_qap_gain;)
                }

                METRICS_TIME(ep)
                METRICS(time.back() += get_seconds(sp, ep);)
            }
            METRICS(iteration_time.push_back(time);)
            METRICS(iteration_time_get_boundary.push_back(time_get_boundary);)
            METRICS(iteration_time_initialize.push_back(time_initialize);)
            METRICS(iteration_time_queue.push_back(time_queue);)
            METRICS(iteration_time_moves.push_back(time_moves);)
            METRICS(iteration_qap_delta.push_back(qap_delta);)
        }

        JSONString get_stats() override {
            std::string stats = "{ \n";
#if COLLECT_METRICS
            std::vector<f64> level_time(iteration_time.size(), 0.0);
            std::vector<f64> level_time_get_boundary(iteration_time.size(), 0.0);
            std::vector<f64> level_time_initialize(iteration_time.size(), 0.0);
            std::vector<f64> level_time_queue(iteration_time.size(), 0.0);
            std::vector<f64> level_time_moves(iteration_time.size(), 0.0);
            std::vector<s64> level_qap_delta(iteration_time.size(), 0);

            for (size_t i = 0; i < iteration_time.size(); ++i) {
                level_time[i]              = sum<f64>(iteration_time[i]);
                level_time_get_boundary[i] = sum<f64>(iteration_time_get_boundary[i]);
                level_time_initialize[i]   = sum<f64>(iteration_time_initialize[i]);
                level_time_queue[i]        = sum<f64>(iteration_time_queue[i]);
                level_time_moves[i]        = sum<f64>(iteration_time_moves[i]);
                level_qap_delta[i]         = sum<s64>(iteration_qap_delta[i]);
            }

            f64 global_time              = sum<f64>(level_time);
            f64 global_time_get_boundary = sum<f64>(level_time_get_boundary);
            f64 global_time_initialize   = sum<f64>(level_time_initialize);
            f64 global_time_queue        = sum<f64>(level_time_queue);
            f64 global_time_moves        = sum<f64>(level_time_moves);
            s64 global_qap_delta         = sum<s64>(level_qap_delta);

            stats += to_JSON_MACRO(global_time);
            stats += to_JSON_MACRO(global_time_get_boundary);
            stats += to_JSON_MACRO(global_time_initialize);
            stats += to_JSON_MACRO(global_time_queue);
            stats += to_JSON_MACRO(global_time_moves);
            stats += to_JSON_MACRO(global_qap_delta);
            stats += to_JSON_MACRO(level_time);
            stats += to_JSON_MACRO(level_time_get_boundary);
            stats += to_JSON_MACRO(level_time_initialize);
            stats += to_JSON_MACRO(level_time_queue);
            stats += to_JSON_MACRO(level_time_moves);
            stats += to_JSON_MACRO(level_qap_delta);
            stats += to_JSON_MACRO(iteration_time);
            stats += to_JSON_MACRO(iteration_time_get_boundary);
            stats += to_JSON_MACRO(iteration_time_initialize);
            stats += to_JSON_MACRO(iteration_time_queue);
            stats += to_JSON_MACRO(iteration_time_moves);
            stats += to_JSON_MACRO(iteration_qap_delta);
#endif
            stats.pop_back();
            stats.pop_back();
            stats += "\n}";

            JSONString json_stats;
            json_stats.s = stats;
            return json_stats;
        }
    };
     */
}

#endif //HEIPROMAP_MULTI_TRY_FM_REFINEMENT_H
