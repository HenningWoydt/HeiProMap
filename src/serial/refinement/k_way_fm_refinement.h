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

#ifndef HEIPROMAP_K_WAY_FM_REFINEMENT_H
#define HEIPROMAP_K_WAY_FM_REFINEMENT_H

#include <queue>

#include "../../commons/utils.h"
#include "../datastructures/functions.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/qap.h"

namespace HeiProMap {
    class KWayFMRefinementConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit KWayFMRefinementConfiguration(const std::string& t_name) : ISerialRefinerConfiguration(t_name) {}

        u64 max_iteration = 1; // how many iterations to run the algorithm at most
        f64 alpha         = 10.0;
        f64 beta          = 1.0;
    };

    class KWayFMRefinement final : public ISerialRefiner {
    private:
        vertex_t m_n    = 0;
        vertex_t m_m    = 0;
        partition_t m_k = 0;
        f64 m_imbalance = 0.0;
        weight_t m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;

        u32* vertex_used  = nullptr;
        u32 vertex_marker = 0;

        u32* block_used  = nullptr;
        u32 block_marker = 0;

        std::priority_queue<KWayFMMove> heap;

        Move* moves       = nullptr;
        size_t moves_size = 0;

        RandomEngine* random_engine                 = nullptr;
        const KWayFMRefinementConfiguration* config = nullptr;
        StatisticCollector* m_stat_collector        = nullptr;

    public:
        KWayFMRefinement() = default;

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
            m_imbalance = t_imbalance;
            m_lmax      = t_lmax;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;

            random_engine    = &t_random_engine;
            config           = dynamic_cast<const KWayFMRefinementConfiguration*>(&i_config);
            m_stat_collector = &t_stat_collect;

            vertex_t t_n_64 = round_up_64(t_n);
            vertex_t t_k_64 = round_up_64(t_k);

            vertex_used = (u32*)aligned_alloc(64, t_n_64 * sizeof(u32));
            std::fill_n(vertex_used, t_n_64, 0);
            vertex_marker = 0;

            block_used = (u32*)aligned_alloc(64, t_k_64 * sizeof(u32));
            std::fill_n(block_used, t_k_64, 0);
            block_marker = 0;

            moves      = (Move*)aligned_alloc(64, t_n_64 * sizeof(Move));
            moves_size = 0;
        }

        void refine(const u64 level,
                    const u64 max_level,
                    const graph_t& g,
                    const d_oracle_t& d_oracle,
                    bv_manager_t& bv_manager,
                    p_manager_t& p_manager,
                    q_graph_t& q_graph) override {
            f64 alpha = config->alpha;
            f64 beta  = std::log(g.get_n());

            for (u64 iteration = 0; iteration < config->max_iteration; ++iteration) {
                heap = std::priority_queue<KWayFMMove>();

                // insert all boundary vertices
                forall_bv_iu(bv_manager, j, u)
                    {
                        partition_t u_id  = p_manager[u];
                        weight_t u_weight = g.weight(u);

                        // find all connected partitions to u
                        block_marker += 1;
                        forall_guiv(g, u, i, v)
                            {
                                partition_t v_id = p_manager[v];
                                if (v_id == u_id) { continue; }
                                if (block_used[v_id] == block_marker) { continue; }
                                if (p_manager.get_bweight(v_id) + u_weight > m_lmax) { continue; }

                                s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                                heap.emplace(u, u_id, v_id, qap_delta);

                                block_used[v_id] = block_marker;
                            }
                        endfor
                    }
                endfor

                moves_size        = 0;
                size_t best_idx   = 0;
                s64 curr_qap_gain = 0;
                s64 max_qap_gain  = 0;

                f64 steps_since_last_improvement = 0.0;
                f64 qap_gain_mean                = 0.0;
                f64 qap_gain_var                 = 0.0;

                vertex_marker += 1;
                while (!heap.empty()) {
                    const KWayFMMove move = heap.top();
                    heap.pop();

                    vertex_t vertex        = move.u;
                    partition_t vertex_id  = p_manager[vertex];
                    weight_t vertex_weight = g.weight(vertex);
                    partition_t move_id    = move.to_move_id;

                    if (vertex_used[vertex] == vertex_marker) { continue; }
                    if (p_manager.get_bweight(move_id) + vertex_weight > m_lmax) { continue; }

                    s64 temp_qap_delta = get_u_qap_delta(g, vertex, vertex_id, move_id, p_manager, d_oracle);
                    if (temp_qap_delta != move.qap_delta) { continue; }

                    moves[moves_size++] = Move(vertex, vertex_id, move_id);
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
                    vertex_used[vertex] = vertex_marker;

                    steps_since_last_improvement += 1.0;
                    f64 new_qap_gain_mean = qap_gain_mean + ((f64)move.qap_delta - qap_gain_mean) / steps_since_last_improvement;
                    f64 new_qap_gain_var  = (qap_gain_var + ((f64)move.qap_delta - qap_gain_mean) * ((f64)move.qap_delta - new_qap_gain_mean)) / steps_since_last_improvement;

                    qap_gain_mean = new_qap_gain_mean;
                    qap_gain_var  = new_qap_gain_var;

                    if (steps_since_last_improvement > 2.0 && steps_since_last_improvement * qap_gain_mean * qap_gain_mean > alpha * qap_gain_var + beta) { break; }

                    // we have to push or update the neighbors that were not moved already
                    forall_guiv(g, vertex, i, neighbor)
                        {
                            if (vertex_used[neighbor] == vertex_marker) { continue; }
                            if (!is_boundary(g, p_manager, neighbor)) { continue; }

                            partition_t neighbor_id  = p_manager[neighbor];
                            weight_t neighbor_weight = g.weight(neighbor);

                            block_marker += 1;
                            forall_guiv(g, neighbor, j, v)
                                {
                                    partition_t v_id = p_manager[v];
                                    if (v_id == neighbor_id) { continue; }
                                    if (block_used[v_id] == block_marker) { continue; }
                                    if (p_manager.get_bweight(v_id) + neighbor_weight > m_lmax) { continue; }

                                    s64 qap_delta = get_u_qap_delta(g, neighbor, neighbor_id, v_id, p_manager, d_oracle);
                                    heap.emplace(neighbor, neighbor_id, v_id, qap_delta);

                                    block_used[v_id] = block_marker;
                                }
                            endfor
                        }
                    endfor
                }

                // revert all moves in partitioning manager
                for (size_t i = 0; i < moves_size; i++) {
                    vertex_t vertex        = moves[moves_size - 1 - i].u;
                    weight_t vertex_weight = g.weight(vertex);
                    partition_t vertex_id  = moves[moves_size - 1 - i].to_move_id;
                    partition_t move_id    = moves[moves_size - 1 - i].u_id;

                    p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                }

                // make all moves to best index
                for (size_t i = 0; i < best_idx; ++i) {
                    vertex_t vertex        = moves[i].u;
                    weight_t vertex_weight = g.weight(vertex);
                    partition_t vertex_id  = moves[i].u_id;
                    partition_t move_id    = moves[i].to_move_id;

                    bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                    q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                    p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                }
            }
        }

        JSONString get_stats() override { return {}; };
    };
}

#endif //HEIPROMAP_K_WAY_FM_REFINEMENT_H
