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

#ifndef HEIPROMAP_DEEP_QUOTIENT_GRAPH_REFINEMENT_H
#define HEIPROMAP_DEEP_QUOTIENT_GRAPH_REFINEMENT_H

#include "ISerialDeepRefiner.h"
#include "../../../commons/utils.h"
#include "../../datastructures/indexed_max_heap.h"
#include "../../utility/functions.h"
#include "../../utility/qap.h"

namespace HeiProMap {
    class DeepQuotientGraphRefinementConfiguration final : public ISerialDeepRefinerConfiguration {
    public:
        explicit DeepQuotientGraphRefinementConfiguration(const std::string &t_name) : ISerialDeepRefinerConfiguration(t_name) {}

        u64 max_iteration = 5;
        f64 alpha = 100.0;
        f64 beta_factor = 1.0;

        f64 alpha_edge_cut = 100.0;
        f64 beta_factor_edge_cut = 1.0;
    };

    class DeepQuotientGraphRefinement final : public ISerialDeepRefiner {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
        u64 m_threads = 1;

        struct thread_info {
            // priority queues
            IndexedMaxHeap<s64> boundary_vertices_u;
            IndexedMaxHeap<s64> boundary_vertices_v;

            // store change
            AlignedArray<vertex_t> moves;
            size_t moves_size = 0;
            s64 curr_qap_gain = 0;
            s64 max_qap_gain = 0;
            s64 curr_edge_cut_gain = 0;
            s64 max_edge_cut_gain = 0;
            size_t best_idx = 0;

            // store which vertices have been moved
            AlignedArray<u32> vertex_used;
            u32 vertex_mark = 0;

            RandomEngine random_engine;
        };

        std::vector<thread_info> thread_infos;
        std::mutex mutex;
        f64 time = 0;

        // active block scheduling
        AlignedArray<u8> active_this_round;
        AlignedArray<u8> active_next_round;

        const DeepQuotientGraphRefinementConfiguration *config = nullptr;

    public:
        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const f64 t_imbalance,
                        const u64 t_threads,
                        const std::vector<partition_t> &t_hierarchy,
                        const std::vector<weight_t> &t_distance,
                        RandomEngine &t_random_engine,
                        const ISerialDeepRefinerConfiguration &i_config) override {
            m_n = t_n;
            m_m = t_m;
            m_k = t_k;
            m_threads = t_threads;

            config = dynamic_cast<const DeepQuotientGraphRefinementConfiguration *>(&i_config);

            thread_infos.resize(m_threads);
            for (size_t i = 0; i < m_threads; ++i) {
                thread_infos[i].boundary_vertices_u.initialize(m_n);
                thread_infos[i].boundary_vertices_v.initialize(m_n);

                thread_infos[i].vertex_mark = 0;
                thread_infos[i].vertex_used.initialize(m_n, 0);

                thread_infos[i].moves.initialize(m_n);
                thread_infos[i].moves_size = 0;

                thread_infos[i].random_engine = RandomEngine(t_random_engine.get_u32());
            }

            // active block scheduling
            active_this_round.initialize(m_k);
            active_next_round.initialize(m_k);
        }

        void refine(const u64 level,
                    const u64 max_level,
                    const graph_t &g,
                    deep_d_oracle_t &d_oracle,
                    deep_bv_manager_t &bv_manager,
                    deep_p_manager_t &p_manager,
                    deep_q_graph_t &q_graph) override {
            active_this_round.initialize(m_k, 1);
            active_next_round.initialize(m_k, 0);

            u64 iteration = 0;
            while (iteration < config->max_iteration) {
                iteration += 1;

                // determine all pairs in the quotient graph
                std::vector<std::vector<std::pair<partition_t, partition_t>>> matchings = q_graph.get_distance_3_matchings(active_this_round);
                for (auto &matching: matchings) {
#pragma omp parallel for num_threads(m_threads) schedule(dynamic)
                    for (auto [u_id, v_id]: matching) {
                        u64 thread_id = omp_get_thread_num();
                        if (d_oracle.last_level_pair(u_id, v_id)) {
                            refine_blocks_edge_cut(g, bv_manager, p_manager, q_graph, u_id, v_id, thread_id);
                        } else {
                            refine_blocks(g, d_oracle, bv_manager, p_manager, q_graph, u_id, v_id, thread_id);
                        }
                    }
                }

                std::swap(active_this_round, active_next_round);
                active_next_round.initialize(m_k, 0);
            }
        }

        void refine_blocks(const graph_t &g,
                           deep_d_oracle_t &d_oracle,
                           deep_bv_manager_t &bv_manager,
                           deep_p_manager_t &p_manager,
                           deep_q_graph_t &q_graph,
                           partition_t u_id,
                           partition_t v_id,
                           u64 thread_id) {
            f64 alpha = config->alpha;
            f64 beta = std::log(g.get_n()) * config->beta_factor;

            // get data for thread
            IndexedMaxHeap<s64> &boundary_vertices_u = thread_infos[thread_id].boundary_vertices_u;
            IndexedMaxHeap<s64> &boundary_vertices_v = thread_infos[thread_id].boundary_vertices_v;

            // store change
            AlignedArray<vertex_t> &moves = thread_infos[thread_id].moves;
            size_t &moves_size = thread_infos[thread_id].moves_size;
            s64 &curr_qap_gain = thread_infos[thread_id].curr_qap_gain;
            s64 &max_qap_gain = thread_infos[thread_id].max_qap_gain;
            size_t &best_idx = thread_infos[thread_id].best_idx;

            // store which vertices have been moved
            AlignedArray<u32> &vertex_used = thread_infos[thread_id].vertex_used;
            u32 &vertex_mark = thread_infos[thread_id].vertex_mark;

            RandomEngine &random_engine = thread_infos[thread_id].random_engine;

            // add all boundary vertices with gain
            boundary_vertices_u.clear();
            boundary_vertices_v.clear();
            vertex_mark += 1;

            size_t max_n_swaps = 0;
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

            // start executing moves based on the TopGain method
            vertex_mark += 1;
            moves_size = 0;
            best_idx = 0;
            curr_qap_gain = 0;
            max_qap_gain = 0;

            f64 steps_since_last_improvement = 0.0;
            f64 qap_gain_mean = 0.0;
            f64 qap_gain_var = 0.0;

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
                        choose_u = random_engine.get_f32() < 0.5;
                    }

                    // 3. if one block is overloaded, choose the larger one, if both same sizes, then randomly
                    weight_t u_id_weight = p_manager.get_bweight(u_id);
                    weight_t u_id_lmax = p_manager.get_lmax(u_id);
                    weight_t v_id_weight = p_manager.get_bweight(v_id);
                    weight_t v_id_lmax = p_manager.get_lmax(v_id);

                    if (u_id_weight > u_id_lmax && u_id_weight > v_id_weight) { choose_u = true; }
                    if (v_id_weight > v_id_lmax && v_id_weight > u_id_weight) { choose_u = false; }
                    if (u_id_weight > u_id_lmax && v_id_weight > v_id_lmax && u_id_weight == v_id_weight) { choose_u = random_engine.get_f32() < 0.5; }
                }

                // choose the priority queue
                IndexedMaxHeap<s64> &boundary_vertices = choose_u ? boundary_vertices_u : boundary_vertices_v;

                vertex_t vertex = boundary_vertices.top_key();
                s64 qap_delta = boundary_vertices.top();
                weight_t vertex_weight = g.weight(vertex);
                partition_t vertex_id = choose_u ? u_id : v_id;
                partition_t move_id = choose_u ? v_id : u_id;
                weight_t move_lmax = choose_u ? p_manager.get_lmax(v_id) : p_manager.get_lmax(u_id);
                weight_t stay_lmax = choose_u ? p_manager.get_lmax(u_id) : p_manager.get_lmax(v_id);
                boundary_vertices.pop();

                // move the vertex
                moves[moves_size++] = vertex;
                curr_qap_gain += qap_delta;
                if (curr_qap_gain >= max_qap_gain && p_manager.get_bweight(move_id) + vertex_weight <= move_lmax && p_manager.get_bweight(vertex_id) - vertex_weight <= stay_lmax) {
                    best_idx = moves_size;
                    max_qap_gain = curr_qap_gain;

                    steps_since_last_improvement = 0.0;
                    qap_gain_mean = 0.0;
                    qap_gain_var = 0.0;
                }

                // make move in structures
                p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                vertex_used[vertex] = vertex_mark;

                steps_since_last_improvement += 1.0;
                f64 new_qap_gain_mean = qap_gain_mean + ((f64) qap_delta - qap_gain_mean) / steps_since_last_improvement;
                f64 new_qap_gain_var = (qap_gain_var + ((f64) qap_delta - qap_gain_mean) * ((f64) qap_delta - new_qap_gain_mean)) / steps_since_last_improvement;

                qap_gain_mean = new_qap_gain_mean;
                qap_gain_var = new_qap_gain_var;

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

            // revert all moves in partitioning manager
            for (size_t i = 0; i < moves_size; i++) {
                vertex_t vertex = moves[moves_size - 1 - i];
                weight_t vertex_weight = g.weight(vertex);
                partition_t vertex_id = p_manager[vertex];
                partition_t move_id = u_id == vertex_id ? v_id : u_id;

                p_manager.move(vertex, vertex_weight, vertex_id, move_id);
            }

            // make all moves to best index
            for (size_t i = 0; i < best_idx; ++i) {
                vertex_t vertex = moves[i];
                weight_t vertex_weight = g.weight(vertex);
                partition_t vertex_id = p_manager[vertex];
                partition_t move_id = u_id == vertex_id ? v_id : u_id;

                bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                p_manager.move(vertex, vertex_weight, vertex_id, move_id);
            }

            if (max_qap_gain > 0) {
                active_next_round[u_id] = 1;
                active_next_round[v_id] = 1;
            }
        }

        void refine_blocks_edge_cut(const graph_t &g,
                                    deep_bv_manager_t &bv_manager,
                                    deep_p_manager_t &p_manager,
                                    deep_q_graph_t &q_graph,
                                    partition_t u_id,
                                    partition_t v_id,
                                    u64 thread_id) {
            f64 alpha = config->alpha_edge_cut;
            f64 beta = std::log(g.get_n()) * config->beta_factor_edge_cut;

            // get data for thread
            IndexedMaxHeap<s64> &boundary_vertices_u = thread_infos[thread_id].boundary_vertices_u;
            IndexedMaxHeap<s64> &boundary_vertices_v = thread_infos[thread_id].boundary_vertices_v;

            // store change
            AlignedArray<vertex_t> &moves = thread_infos[thread_id].moves;
            size_t &moves_size = thread_infos[thread_id].moves_size;
            s64 &curr_edge_cut_gain = thread_infos[thread_id].curr_edge_cut_gain;
            s64 &max_edge_cut_gain = thread_infos[thread_id].max_edge_cut_gain;
            size_t &best_idx = thread_infos[thread_id].best_idx;

            // store which vertices have been moved
            AlignedArray<u32> &vertex_used = thread_infos[thread_id].vertex_used;
            u32 &vertex_mark = thread_infos[thread_id].vertex_mark;

            RandomEngine &random_engine = thread_infos[thread_id].random_engine;

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
                                    s64 edge_cut_delta_u = get_u_edge_cut_delta(g, u, u_id, v_id, p_manager);
                                    boundary_vertices_u.push(u, edge_cut_delta_u);
                                    vertex_used[u] = vertex_mark;
                                    max_n_swaps += 1;
                                }

                                if (vertex_used[v] != vertex_mark) {
                                    s64 edge_cut_delta_v = get_u_edge_cut_delta(g, v, v_id, u_id, p_manager);
                                    boundary_vertices_v.push(v, edge_cut_delta_v);
                                    vertex_used[v] = vertex_mark;
                                    max_n_swaps += 1;
                                }
                            }
                        }
                    endfor
                }
            endfor

            // start executing moves based on the TopGain method
            vertex_mark += 1;
            moves_size = 0;
            best_idx = 0;
            curr_edge_cut_gain = 0;
            max_edge_cut_gain = 0;

            f64 steps_since_last_improvement = 0.0;
            f64 edge_cut_gain_mean = 0.0;
            f64 edge_cut_gain_var = 0.0;

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
                        choose_u = random_engine.get_f32() < 0.5;
                    }

                    // 3. if one block is overloaded, choose the larger one, if both same sizes, then randomly
                    weight_t u_id_weight = p_manager.get_bweight(u_id);
                    weight_t u_id_lmax = p_manager.get_lmax(u_id);
                    weight_t v_id_weight = p_manager.get_bweight(v_id);
                    weight_t v_id_lmax = p_manager.get_lmax(v_id);

                    if (u_id_weight > u_id_lmax && u_id_weight > v_id_weight) { choose_u = true; }
                    if (v_id_weight > v_id_lmax && v_id_weight > u_id_weight) { choose_u = false; }
                    if (u_id_weight > u_id_lmax && v_id_weight > v_id_lmax && u_id_weight == v_id_weight) { choose_u = random_engine.get_f32() < 0.5; }
                }

                // choose the priority queue
                IndexedMaxHeap<s64> &boundary_vertices = choose_u ? boundary_vertices_u : boundary_vertices_v;

                vertex_t vertex = boundary_vertices.top_key();
                s64 edge_cut_delta = boundary_vertices.top();
                weight_t vertex_weight = g.weight(vertex);
                partition_t vertex_id = choose_u ? u_id : v_id;
                partition_t move_id = choose_u ? v_id : u_id;
                weight_t move_lmax = choose_u ? p_manager.get_lmax(v_id) : p_manager.get_lmax(u_id);
                weight_t stay_lmax = choose_u ? p_manager.get_lmax(u_id) : p_manager.get_lmax(v_id);
                boundary_vertices.pop();

                // move the vertex
                moves[moves_size++] = vertex;
                curr_edge_cut_gain += edge_cut_delta;
                if (curr_edge_cut_gain >= max_edge_cut_gain && p_manager.get_bweight(move_id) + vertex_weight <= move_lmax && p_manager.get_bweight(vertex_id) - vertex_weight <= stay_lmax) {
                    best_idx = moves_size;
                    max_edge_cut_gain = curr_edge_cut_gain;

                    steps_since_last_improvement = 0.0;
                    edge_cut_gain_mean = 0.0;
                    edge_cut_gain_var = 0.0;
                }

                // make move in structures
                p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                vertex_used[vertex] = vertex_mark;

                steps_since_last_improvement += 1.0;
                f64 new_edge_cut_gain_mean = edge_cut_gain_mean + ((f64) edge_cut_delta - edge_cut_gain_mean) / steps_since_last_improvement;
                f64 new_edge_cut_gain_var = (edge_cut_gain_var + ((f64) edge_cut_delta - edge_cut_gain_mean) * ((f64) edge_cut_delta - new_edge_cut_gain_mean)) / steps_since_last_improvement;

                edge_cut_gain_mean = new_edge_cut_gain_mean;
                edge_cut_gain_var = new_edge_cut_gain_var;

                if (steps_since_last_improvement > 2.0 && steps_since_last_improvement * edge_cut_gain_mean * edge_cut_gain_mean > alpha * edge_cut_gain_var + beta) { break; }

                // we have to push or update the neighbors that were not moved already
                forall_guiv(g, vertex, i, neighbor)
                    {
                        if (vertex_used[neighbor] == vertex_mark) { continue; }

                        partition_t neighbor_id = p_manager[neighbor];

                        if (neighbor_id != u_id && neighbor_id != v_id) { continue; }

                        partition_t new_id = neighbor_id == vertex_id ? move_id : vertex_id;

                        bool is_connected_to_new_id;
                        s64 new_edge_cut_delta = get_u_edge_cut_delta_and_is_connected_to(g, neighbor, neighbor_id, new_id, is_connected_to_new_id, p_manager);

                        if (!is_connected_to_new_id) { continue; }

                        if (neighbor_id == u_id) {
                            boundary_vertices_u.push_update(neighbor, new_edge_cut_delta);
                        } else {
                            boundary_vertices_v.push_update(neighbor, new_edge_cut_delta);
                        }
                    }
                endfor

                // remove vertex from u if it is not boundary
                while (!boundary_vertices_u.empty() && !is_connected_to(g, p_manager, boundary_vertices_u.top_key(), v_id)) { boundary_vertices_u.pop(); }

                // remove vertex from v if it is not boundary
                while (!boundary_vertices_v.empty() && !is_connected_to(g, p_manager, boundary_vertices_v.top_key(), u_id)) { boundary_vertices_v.pop(); }
            }

            // revert all moves in partitioning manager
            for (size_t i = 0; i < moves_size; i++) {
                vertex_t vertex = moves[moves_size - 1 - i];
                weight_t vertex_weight = g.weight(vertex);
                partition_t vertex_id = p_manager[vertex];
                partition_t move_id = u_id == vertex_id ? v_id : u_id;

                p_manager.move(vertex, vertex_weight, vertex_id, move_id);
            }

            // make all moves to best index
            // mutex.lock();
            for (size_t i = 0; i < best_idx; ++i) {
                vertex_t vertex = moves[i];
                weight_t vertex_weight = g.weight(vertex);
                partition_t vertex_id = p_manager[vertex];
                partition_t move_id = u_id == vertex_id ? v_id : u_id;

                bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                p_manager.move(vertex, vertex_weight, vertex_id, move_id);
            }
            // mutex.unlock();

            if (max_edge_cut_gain > 0) {
                active_next_round[u_id] = 1;
                active_next_round[v_id] = 1;
            }
        }
    };
}

#endif //HEIPROMAP_DEEP_QUOTIENT_GRAPH_REFINEMENT_H
