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

#ifndef HEIPROMAP_QUOTIENT_GRAPH_REFINEMENT_H
#define HEIPROMAP_QUOTIENT_GRAPH_REFINEMENT_H

#include "quotient_graph_refinement_Faraj20.h"
#include "../../commons/utils.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/qap.h"

namespace HeiProMap {
    struct PairWeight {
        partition_t id1;
        partition_t id2;
        weight_t distance;

        bool operator<(const PairWeight& p) const {
            return distance < p.distance;
        }

        bool operator>(const PairWeight& p) const {
            return distance > p.distance;
        }
    };

    class QuotientGraphRefinementConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit QuotientGraphRefinementConfiguration(const std::string& t_name) : ISerialRefinerConfiguration(t_name) {}
        u64 max_iteration = 1;
        f64 alpha = 100.0;
        f64 beta  = 1.0;
    };

    class QuotientGraphRefinement final : public ISerialRefiner {
        vertex_t m_n    = 0;
        vertex_t m_m    = 0;
        partition_t m_k = 0;
        weight_t m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;

        // priority queues
        IndexedMaxHeap<s64> boundary_vertices_u;
        IndexedMaxHeap<s64> boundary_vertices_v;

        // store change
        vertex_t* moves   = nullptr;
        size_t moves_size = 0;
        s64 curr_qap_gain = 0;
        s64 max_qap_gain  = 0;
        size_t best_idx   = 0;

        // active block scheduling
        u8* active_this_round = nullptr;
        u8* active_next_round = nullptr;
        PairWeight* pairs     = nullptr;
        size_t pairs_size     = 0;

        // store which vertices have been moved
        u32* vertex_used = nullptr;
        u32 vertex_mark  = 0;

        RandomEngine* random_engine                        = nullptr;
        const QuotientGraphRefinementConfiguration* config = nullptr;
        StatisticCollector* m_stat_collector               = nullptr;

        METRICS(std::vector<std::vector<f64>> iteration_time;)
        METRICS(std::vector<std::vector<f64>> iteration_time_get_pairs;)
        METRICS(std::vector<std::vector<f64>> iteration_time_initialize;)
        METRICS(std::vector<std::vector<f64>> iteration_time_queue;)
        METRICS(std::vector<std::vector<f64>> iteration_time_moves;)
        METRICS(std::vector<std::vector<s64>> iteration_qap_delta;)

    public:
        QuotientGraphRefinement() = default;

        ~QuotientGraphRefinement() override {
            free(vertex_used);
            free(moves);
            free(active_this_round);
            free(active_next_round);
            free(pairs);
        }

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const weight_t t_lmax,
                        const std::vector<partition_t>& t_hierarchy,
                        const std::vector<weight_t>& t_distance,
                        RandomEngine& t_random_engine,
                        const ISerialRefinerConfiguration& i_config,
                        StatisticCollector& t_stat_collect) override {
            vertex_t t_n_64        = round_up_64(t_n);
            partition_t t_k_64     = round_up_64(t_k);
            partition_t t_k_t_k_64 = round_up_64(t_k * t_k);

            m_n         = t_n;
            m_m         = t_m;
            m_k         = t_k;
            m_lmax      = t_lmax;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;

            random_engine    = &t_random_engine;
            config           = dynamic_cast<const QuotientGraphRefinementConfiguration*>(&i_config);
            m_stat_collector = &t_stat_collect;

            // priority queues
            boundary_vertices_u.initialize(m_n);
            boundary_vertices_v.initialize(m_n);

            vertex_mark = 0;
            vertex_used = (u32*)aligned_alloc(64, t_n_64 * sizeof(u32));
            std::fill_n(vertex_used, t_n_64, vertex_mark);

            moves      = (vertex_t*)aligned_alloc(64, t_n_64 * sizeof(vertex_t));
            moves_size = 0;

            // active block scheduling
            active_this_round = (u8*)aligned_alloc(64, t_k_64 * sizeof(u8));
            active_next_round = (u8*)aligned_alloc(64, t_k_64 * sizeof(u8));
            pairs             = (PairWeight*)aligned_alloc(64, t_k_t_k_64 * sizeof(PairWeight));
            pairs_size        = 0;
        }

        void refine(const u64 level,
                    const u64 max_level,
                    const graph_t& g,
                    const d_oracle_t& d_oracle,
                    bv_manager_t& bv_manager,
                    p_manager_t& p_manager,
                    q_graph_t& q_graph) override {
            METRICS(iteration_time.emplace_back();)
            METRICS(iteration_time_get_pairs.emplace_back();)
            METRICS(iteration_time_initialize.emplace_back();)
            METRICS(iteration_time_queue.emplace_back();)
            METRICS(iteration_time_moves.emplace_back();)
            METRICS(iteration_qap_delta.emplace_back();)

            f64 alpha = config->alpha;
            f64 beta  = std::log(g.get_n());

            std::fill_n(active_this_round, m_k, 1);
            std::fill_n(active_next_round, m_k, 0);

            u64 iteration = 0;
            s64 last_qap_gain = std::numeric_limits<s64>::max();
            bool one_pair_active = true;
            while (one_pair_active && iteration < config->max_iteration && last_qap_gain > 100) {
                METRICS(iteration_time.back().push_back(0.0);)
                METRICS(iteration_time_get_pairs.back().push_back(0.0);)
                METRICS(iteration_time_initialize.back().push_back(0.0);)
                METRICS(iteration_time_queue.back().push_back(0.0);)
                METRICS(iteration_time_moves.back().push_back(0.0);)
                METRICS(iteration_qap_delta.back().push_back(0);)

                METRICS_TIME(sp)

                iteration += 1;
                last_qap_gain = 0;
                one_pair_active = false;

                METRICS_TIME(sp_get_pairs)

                // determine all pairs in the quotient graph
                pairs_size = 0;
                for (partition_t u_id = 0; u_id < m_k; ++u_id) {
                    for (partition_t v_id = u_id + 1; v_id < m_k; ++v_id) {
                        if (q_graph.has_edge(u_id, v_id) && (active_this_round[u_id] || active_this_round[v_id])) {
                            pairs[pairs_size++] = {u_id, v_id, d_oracle.get(u_id, v_id)};
                        }
                    }
                }
                std::sort(pairs, pairs + pairs_size, std::greater<>());

                METRICS_TIME(ep_get_pairs)
                METRICS(iteration_time_get_pairs.back().back() += get_seconds(sp_get_pairs, ep_get_pairs);)

                // for each pair do fm refinement
                for (size_t j = 0; j < pairs_size; ++j) {
                    auto [u_id, v_id, distance] = pairs[j];

                    METRICS_TIME(sp_initialize)

                    size_t max_n_swaps = 0;

                    // add all boundary vertices with gain
                    boundary_vertices_u.clear();
                    boundary_vertices_v.clear();
                    vertex_mark += 1;
                    forall_bv_id_iu(bv_manager, u_id, k, u)
                        {
                            forall_guiv(g, u, i, v)
                                {
                                    if (p_manager[v] == v_id) {
                                        // u is connected to block v_id
                                        if (vertex_used[u] != vertex_mark) {
                                            s64 qap_delta_u = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                                            boundary_vertices_u.push(u, qap_delta_u);
                                            vertex_used[u] = vertex_mark;
                                            max_n_swaps += 1;
                                        }

                                        if (vertex_used[v] != vertex_mark) {
                                            s64 qap_delta_v = get_u_qap_delta(g, v, v_id, u_id, p_manager, d_oracle);
                                            boundary_vertices_v.push(v, qap_delta_v);
                                            vertex_used[v] = vertex_mark;
                                            max_n_swaps += 1;
                                        }
                                    }
                                }
                            endfor
                        }
                    endfor

                    METRICS_TIME(ep_initialize)
                    METRICS(iteration_time_initialize.back().back() += get_seconds(sp_initialize, ep_initialize);)

                    // start executing moves based on the TopGain method
                    vertex_mark += 1;
                    moves_size    = 0;
                    best_idx      = 0;
                    curr_qap_gain = 0;
                    max_qap_gain  = 0;

                    f64 steps_since_last_improvement = 0.0;
                    f64 qap_gain_mean                = 0.0;
                    f64 qap_gain_var                 = 0.0;

                    METRICS_TIME(sp_queue)

                    while ((!boundary_vertices_u.empty() || !boundary_vertices_v.empty()) && moves_size < max_n_swaps) {
                        // determine from which block to choose
                        bool choose_u = true;
                        // 1. if one block is empty, then choose the other one
                        if (boundary_vertices_u.empty() || boundary_vertices_v.empty()) {
                            choose_u = boundary_vertices_v.empty();
                        } else {
                            // 2. choose the block with greater gain and randomly if even
                            if (boundary_vertices_v.top() > boundary_vertices_u.top()) {
                                choose_u = false;
                            } else if (boundary_vertices_v.top() == boundary_vertices_u.top()) {
                                choose_u = random_engine->get_f32() < 0.5;
                            }

                            // 3. if one block is overloaded, choose the larger one, if both same sizes, then randomly
                            weight_t u_id_weight = p_manager.get_bweight(u_id);
                            weight_t v_id_weight = p_manager.get_bweight(v_id);

                            if (u_id_weight > m_lmax && u_id_weight > v_id_weight) { choose_u = true; }
                            if (v_id_weight > m_lmax && v_id_weight > u_id_weight) { choose_u = false; }
                            if (u_id_weight > m_lmax && v_id_weight > m_lmax && u_id_weight == v_id_weight) { choose_u = random_engine->get_f32() < 0.5; }
                        }

                        // choose the priority queue
                        IndexedMaxHeap<s64>& boundary_vertices = choose_u ? boundary_vertices_u : boundary_vertices_v;

                        vertex_t vertex           = boundary_vertices.top_key();
                        s64 qap_delta             = boundary_vertices.top();
                        weight_t vertex_weight    = g.weight(vertex);
                        partition_t vertex_id     = choose_u ? u_id : v_id;
                        partition_t move_id       = choose_u ? v_id : u_id;
                        weight_t partition_weight = p_manager.get_bweight(move_id);
                        boundary_vertices.pop();

                        // move the vertex
                        moves[moves_size++] = vertex;
                        curr_qap_gain += qap_delta;
                        if (curr_qap_gain >= max_qap_gain && partition_weight + vertex_weight <= m_lmax) {
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
                        f64 new_qap_gain_mean = qap_gain_mean + ((f64)qap_delta - qap_gain_mean) / steps_since_last_improvement;
                        f64 new_qap_gain_var  = (qap_gain_var + ((f64)qap_delta - qap_gain_mean) * ((f64)qap_delta - new_qap_gain_mean)) / steps_since_last_improvement;

                        qap_gain_mean = new_qap_gain_mean;
                        qap_gain_var  = new_qap_gain_var;

                        if (steps_since_last_improvement > 2.0 && steps_since_last_improvement * qap_gain_mean * qap_gain_mean > alpha * qap_gain_var + beta) { break; }

                        // we have to push or update the neighbors that were not moved already
                        forall_guiv(g, vertex, i, neighbor)
                            {
                                if (vertex_used[neighbor] == vertex_mark) { continue; }

                                partition_t neighbor_id = p_manager[neighbor];

                                if (neighbor_id != u_id && neighbor_id != v_id) { continue; }

                                partition_t new_id = neighbor_id == vertex_id ? move_id : vertex_id;

                                bool is_connected_to_new_id;
                                s64 new_qap_delta = get_u_qap_delta_and_is_connected_to(g, neighbor, neighbor_id, new_id, is_connected_to_new_id, p_manager, d_oracle);

                                if (!is_connected_to_new_id) { continue; }

                                if (neighbor_id == u_id) {
                                    boundary_vertices_u.push_update(neighbor, new_qap_delta);
                                } else {
                                    boundary_vertices_v.push_update(neighbor, new_qap_delta);
                                }
                            }
                        endfor

                        // remove vertex from u if it is not boundary
                        while (!boundary_vertices_u.empty() && !is_connected_to(g, p_manager, boundary_vertices_u.top_key(), v_id)) { boundary_vertices_u.pop(); }

                        // remove vertex from v if it is not boundary
                        while (!boundary_vertices_v.empty() && !is_connected_to(g, p_manager, boundary_vertices_v.top_key(), u_id)) { boundary_vertices_v.pop(); }
                    }

                    METRICS_TIME(ep_queue)
                    METRICS(iteration_time_queue.back().back() += get_seconds(sp_queue, ep_queue);)

                    METRICS_TIME(sp_moves)
                    // revert all moves in partitioning manager
                    for (size_t i = 0; i < moves_size; i++) {
                        vertex_t vertex        = moves[moves_size - 1 - i];
                        weight_t vertex_weight = g.weight(vertex);
                        partition_t vertex_id  = p_manager[vertex];
                        partition_t move_id    = u_id == vertex_id ? v_id : u_id;

                        p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                    }

                    // make all moves to best index
                    for (size_t i = 0; i < best_idx; ++i) {
                        vertex_t vertex        = moves[i];
                        weight_t vertex_weight = g.weight(vertex);
                        partition_t vertex_id  = p_manager[vertex];
                        partition_t move_id    = u_id == vertex_id ? v_id : u_id;

                        bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                        q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                        p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                    }
                    METRICS_TIME(ep_moves)
                    METRICS(iteration_time_moves.back().back() += get_seconds(sp_moves, ep_moves);)

                    if (max_qap_gain > 0) {
                        active_next_round[u_id] = 1;
                        active_next_round[v_id] = 1;
                        one_pair_active         = true;
                    }
                    last_qap_gain += max_qap_gain;

                    METRICS(iteration_qap_delta.back().back() += max_qap_gain;)
                }
                std::swap(active_this_round, active_next_round);
                std::fill_n(active_next_round, m_k, 0);

                METRICS_TIME(ep)
                METRICS(iteration_time.back().back() += get_seconds(sp, ep);)
            }
        }

        JSONString get_stats() override {
            std::string stats = "{ \n";
#if COLLECT_METRICS
            std::vector<f64> level_time(iteration_time.size(), 0.0);
            std::vector<f64> level_time_get_pairs(iteration_time.size(), 0.0);
            std::vector<f64> level_time_initialize(iteration_time.size(), 0.0);
            std::vector<f64> level_time_queue(iteration_time.size(), 0.0);
            std::vector<f64> level_time_moves(iteration_time.size(), 0.0);
            std::vector<s64> level_qap_delta(iteration_time.size(), 0);

            for (size_t i = 0; i < iteration_time.size(); ++i) {
                level_time[i]            = sum<f64>(iteration_time[i]);
                level_time_get_pairs[i]  = sum<f64>(iteration_time_get_pairs[i]);
                level_time_initialize[i] = sum<f64>(iteration_time_initialize[i]);
                level_time_queue[i]      = sum<f64>(iteration_time_queue[i]);
                level_time_moves[i]      = sum<f64>(iteration_time_moves[i]);
                level_qap_delta[i]       = sum<s64>(iteration_qap_delta[i]);
            }

            f64 global_time            = sum<f64>(level_time);
            f64 global_time_get_pairs  = sum<f64>(level_time_get_pairs);
            f64 global_time_initialize = sum<f64>(level_time_initialize);
            f64 global_time_queue      = sum<f64>(level_time_queue);
            f64 global_time_moves      = sum<f64>(level_time_moves);
            s64 global_qap_delta       = sum<s64>(level_qap_delta);

            stats += to_JSON_MACRO(global_time);
            stats += to_JSON_MACRO(global_time_get_pairs);
            stats += to_JSON_MACRO(global_time_initialize);
            stats += to_JSON_MACRO(global_time_queue);
            stats += to_JSON_MACRO(global_time_moves);
            stats += to_JSON_MACRO(global_qap_delta);
            stats += to_JSON_MACRO(level_time);
            stats += to_JSON_MACRO(level_time_get_pairs);
            stats += to_JSON_MACRO(level_time_initialize);
            stats += to_JSON_MACRO(level_time_queue);
            stats += to_JSON_MACRO(level_time_moves);
            stats += to_JSON_MACRO(level_qap_delta);
            stats += to_JSON_MACRO(iteration_time);
            stats += to_JSON_MACRO(iteration_time_get_pairs);
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
}

#endif //HEIPROMAP_QUOTIENT_GRAPH_REFINEMENT_H
