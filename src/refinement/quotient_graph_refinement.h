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

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include <omp.h>

#include "../definitions.h"
#include "../datastructures/block_conn.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/distance_oracle.h"
#include "../datastructures/partition_manager.h"
#include "../datastructures/quotient_graph.h"
#include "../utility/aligned_array.h"
#include "../utility/profiler.h"
#include "../utility/random_engine.h"
#include "../utility/utils.h"
#include "../utility/qap.h"
#include "../utility/indexed_max_heap.h"
#include "../utility/functions.h"

namespace HeiProMap {
    class QuotientGraphRefinementConfiguration final {
    public:
        explicit QuotientGraphRefinementConfiguration(const std::string &t_name) {
            name = t_name;
        }

        std::string name;
        bool enabled = false;
        u64 max_iteration = 1;
        u64 min_n_steps = 4;
        f64 alpha = 100.0;
        f64 beta = 1.0;
        bool use_active_scheduling = true;
        bool use_preemptive_exit = true;
    };

    class QuotientGraphRefinement final {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
        u64 m_threads = 1;

        // active block scheduling
        AlignedArray<u8> active_this_round;
        AlignedArray<u8> active_next_round;

        std::vector<std::vector<vertex_t>> moves_vec;

        // store which vertices have been moved
        AlignedArray<u32> vertex_used;
        u32 global_vertex_mark = 0;

        std::vector<IndexedMaxHeap<weight_t> > boundary_vertices_u_vec;
        std::vector<IndexedMaxHeap<weight_t> > boundary_vertices_v_vec;

        std::vector<RandomEngine> rnd_engines;
        AlignedArray<u8> used_this_round;

        const QuotientGraphRefinementConfiguration *config = nullptr;

    public:
        QuotientGraphRefinement() = default;

        ~QuotientGraphRefinement() = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const u64 t_threads,
                        const u64 seed,
                        const QuotientGraphRefinementConfiguration &i_config) {
            HEIPROMAP_PROFILE_SCOPE("misc", "QuotientGraphRefinement", "initialize");

            m_n = t_n;
            m_m = t_m;
            m_k = t_k;
            m_threads = t_threads;

            config = &i_config;

            global_vertex_mark = 0;
            vertex_used.initialize(m_n, 0);

            // active block scheduling
            active_this_round.initialize(m_k);
            active_next_round.initialize(m_k);

            for (size_t i = 0; i < m_threads; ++i) {
                moves_vec.emplace_back();
                boundary_vertices_u_vec.emplace_back();
                boundary_vertices_v_vec.emplace_back();
            }

            rnd_engines.resize(m_threads);
            for (u64 t = 0; t < m_threads; ++t) {
                rnd_engines[t] = RandomEngine(seed + t);
            }

            used_this_round.initialize(m_k * m_k);

            #pragma omp parallel for num_threads(m_threads)
            for (size_t i = 0; i < m_threads; ++i) {
                boundary_vertices_u_vec[i].initialize(m_n);
                boundary_vertices_v_vec[i].initialize(m_n);
            }
        }

        void refine(graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph,
                    block_conn_t &block_conn,
                    const AlignedArray<weight_t> &lmax_constraints) {
            if (g.uniform_v_weights && g.uniform_e_weights) refine_impl<true, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
            else if (g.uniform_v_weights) refine_impl<true, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
            else if (g.uniform_e_weights) refine_impl<false, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
            else refine_impl<false, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void refine_impl(graph_t &g,
                         d_oracle_t &d_oracle,
                         bv_manager_t &bv_manager,
                         p_manager_t &p_manager,
                         q_graph_t &q_graph,
                         block_conn_t &block_conn,
                         const AlignedArray<weight_t> &lmax_constraints) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "QuotientGraphRefinement", "init_block_scheduling");
            active_this_round.initialize(m_k, 1);
            active_next_round.initialize(m_k, 0);

            std::vector<std::pair<partition_t, partition_t> > matching;

            for (u64 iteration = 0; iteration < config->max_iteration; ++iteration) {
                HEIPROMAP_PROFILE_SCOPE("refinement", "QuotientGraphRefinement", "reset_used_edges");
                std::fill_n(used_this_round.get_ptr(), m_k * m_k, 0);

                HEIPROMAP_PROFILE_SCOPE("refinement", "QuotientGraphRefinement", "matching");
                bool found_matching = false;
                if (m_threads > 1) {
                    found_matching = q_graph.find_distance_3_matching(active_this_round, used_this_round, matching);
                } else {
                    found_matching = q_graph.find_all_pairs(active_this_round, used_this_round, matching);
                }

                while (found_matching) {
                    // Pre-assign unique marks for each thread
                    u32 base_mark = global_vertex_mark + 1;
                    global_vertex_mark += static_cast<u32>(matching.size());

                    #pragma omp parallel for num_threads(m_threads) schedule(dynamic)
                    for (size_t i = 0; i < matching.size(); ++i) {
                        partition_t u_id = matching[i].first;
                        partition_t v_id = matching[i].second;

                        u64 tid = omp_get_thread_num();
                        u32 mark = base_mark + static_cast<u32>(i);

                        if (d_oracle.last_level_pair(u_id, v_id)) {
                            refine_blocks_edge_cut<t_uniform_v_weights, t_uniform_e_weights>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, u_id, v_id, moves_vec[tid], lmax_constraints, boundary_vertices_u_vec[tid], boundary_vertices_v_vec[tid], mark, rnd_engines[tid]);
                        } else {
                            refine_blocks<t_uniform_v_weights, t_uniform_e_weights>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, u_id, v_id, moves_vec[tid], lmax_constraints, boundary_vertices_u_vec[tid], boundary_vertices_v_vec[tid], mark, rnd_engines[tid]);
                        }
                    }

                    HEIPROMAP_PROFILE_SCOPE("refinement", "QuotientGraphRefinement", "matching");
                    if (m_threads > 1) {
                        found_matching = q_graph.find_distance_3_matching(active_this_round, used_this_round, matching);
                    } else {
                        found_matching = q_graph.find_all_pairs(active_this_round, used_this_round, matching);
                    }
                }

                HEIPROMAP_PROFILE_SCOPE("refinement", "QuotientGraphRefinement", "swap");
                std::swap(active_this_round, active_next_round);
                active_next_round.initialize(m_k, 0);
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void refine_blocks(const graph_t &g,
                           d_oracle_t &d_oracle,
                           bv_manager_t &bv_manager,
                           p_manager_t &p_manager,
                           q_graph_t &q_graph,
                           block_conn_t &block_conn,
                           partition_t u_id,
                           partition_t v_id,
                           std::vector<vertex_t> &moves,
                           const AlignedArray<weight_t> &lmax_constraints,
                           IndexedMaxHeap<weight_t> &boundary_vertices_u,
                           IndexedMaxHeap<weight_t> &boundary_vertices_v,
                           vertex_t vertex_mark,
                           RandomEngine &random_engine) {
            f64 alpha = config->alpha * (f64) d_oracle.get(u_id, v_id);
            f64 beta = std::log(g.n) * (f64) d_oracle.get(u_id, v_id);

            HEIPROMAP_PROFILE_SCOPE("refinement", "QuotientGraphRefinement", "initial_qap");
            // add all boundary vertices with gain
            u64 n_init_moves = 0;
            boundary_vertices_u.clear();
            boundary_vertices_v.clear();
            for (size_t j = 0; j < bv_manager.size(u_id); ++j) {
                const vertex_t u = bv_manager.get(u_id, j);
                for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) {
                    if (block_conn.get_id(i) == v_id) {
                        weight_t qap_delta_u = get_u_qap_delta_t<t_uniform_e_weights>(g, u, u_id, v_id, p_manager, d_oracle, block_conn);
                        boundary_vertices_u.push(u, qap_delta_u);
                        n_init_moves += 1;
                        break;
                    }
                }
            }

            for (size_t j = 0; j < bv_manager.size(v_id); ++j) {
                const vertex_t v = bv_manager.get(v_id, j);
                for (size_t i = block_conn.start(v); i < block_conn.end(v); ++i) {
                    if (block_conn.get_id(i) == u_id) {
                        weight_t qap_delta_v = get_u_qap_delta_t<t_uniform_e_weights>(g, v, v_id, u_id, p_manager, d_oracle, block_conn);
                        boundary_vertices_v.push(v, qap_delta_v);
                        n_init_moves += 1;
                        break;
                    }
                }
            }

            // store change
            weight_t curr_qap_gain = 0;
            weight_t max_qap_gain = 0;
            // weight_t curr_edge_cut_gain = 0;
            // weight_t max_edge_cut_gain = 0;
            size_t best_idx = 0;

            u64 steps_since_last_improvement = 0;
            f64 qap_gain_mean = 0.0;
            f64 qap_gain_var = 0.0;

            weight_t u_id_weight_initial = p_manager.get_bweight(u_id);
            weight_t v_id_weight_initial = p_manager.get_bweight(v_id);
            bool best_is_balanced = (u_id_weight_initial <= lmax_constraints[u_id] && v_id_weight_initial <= lmax_constraints[v_id]);

            moves.clear();

            HEIPROMAP_PROFILE_SCOPE("refinement", "QuotientGraphRefinement", "process_queue");
            while ((!boundary_vertices_u.empty() || !boundary_vertices_v.empty()) && moves.size() < n_init_moves) {
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
                    weight_t v_id_weight = p_manager.get_bweight(v_id);

                    if (u_id_weight > lmax_constraints[u_id] && u_id_weight > v_id_weight) { choose_u = true; }
                    if (v_id_weight > lmax_constraints[v_id] && v_id_weight > u_id_weight) { choose_u = false; }
                    if (u_id_weight > lmax_constraints[u_id] && v_id_weight > lmax_constraints[v_id] && u_id_weight == v_id_weight) { choose_u = random_engine.get_f32() < 0.5; }
                }

                // choose the priority queue
                IndexedMaxHeap<weight_t> &boundary_vertices = choose_u ? boundary_vertices_u : boundary_vertices_v;

                vertex_t vertex = boundary_vertices.top_key();
                weight_t qap_delta = boundary_vertices.top();
                weight_t vertex_weight = t_uniform_v_weights ? 1 : g.v_weights[vertex];
                partition_t vertex_id = choose_u ? u_id : v_id;
                partition_t move_id = choose_u ? v_id : u_id;
                boundary_vertices.pop();

                // move the vertex
                moves.push_back(vertex);
                curr_qap_gain += qap_delta;

                weight_t current_move_id_weight = p_manager.get_bweight(move_id) + vertex_weight;
                weight_t current_vertex_id_weight = p_manager.get_bweight(vertex_id) - vertex_weight;
                bool current_is_balanced = (current_move_id_weight <= lmax_constraints[move_id] && current_vertex_id_weight <= lmax_constraints[vertex_id]);

                bool update_best = false;
                if (current_is_balanced) {
                    if (!best_is_balanced || curr_qap_gain >= max_qap_gain) {
                        update_best = true;
                    }
                } else {
                    if (!best_is_balanced) {
                        weight_t move_id_weight_initial = move_id == u_id ? u_id_weight_initial : v_id_weight_initial;
                        if (curr_qap_gain >= max_qap_gain && current_move_id_weight <= std::max(move_id_weight_initial, lmax_constraints[move_id])) {
                            update_best = true;
                        }
                    }
                }

                if (update_best) {
                    best_idx = moves.size();
                    max_qap_gain = curr_qap_gain;
                    best_is_balanced = current_is_balanced;

                    steps_since_last_improvement = 0;
                    qap_gain_mean = 0.0;
                    qap_gain_var = 0.0;
                }

                // make move in structures
                p_manager.move_serial(vertex, vertex_weight, vertex_id, move_id);
                vertex_used[vertex] = vertex_mark;

                steps_since_last_improvement += 1;
                f64 new_qap_gain_mean = qap_gain_mean + ((f64) qap_delta - qap_gain_mean) / (f64) steps_since_last_improvement;
                f64 new_qap_gain_var = (qap_gain_var + ((f64) qap_delta - qap_gain_mean) * ((f64) qap_delta - new_qap_gain_mean)) / (f64) steps_since_last_improvement;

                qap_gain_mean = new_qap_gain_mean;
                qap_gain_var = new_qap_gain_var;

                if (config->use_preemptive_exit) {
                    if (steps_since_last_improvement > config->min_n_steps && (f64) steps_since_last_improvement * qap_gain_mean * qap_gain_mean > alpha * qap_gain_var + beta) {
                        break;
                    }
                }

                // we have to push or update the neighbors that were not moved already
                for (size_t i = g.neighborhoods[vertex]; i < g.neighborhoods[vertex + 1]; ++i) {
                    const vertex_t neighbor = g.edges_v[i];

                    if (vertex_used[neighbor] == vertex_mark) { continue; }

                    partition_t neighbor_id = p_manager[neighbor];

                    if (neighbor_id != u_id && neighbor_id != v_id) { continue; }

                    partition_t new_id = neighbor_id == vertex_id ? move_id : vertex_id;

                    bool is_connected_to_new_id;
                    weight_t new_qap_delta = get_u_qap_delta_and_is_connected_to_t<t_uniform_e_weights>(g, neighbor, neighbor_id, new_id, is_connected_to_new_id, p_manager, d_oracle);

                    if (!is_connected_to_new_id) { continue; }

                    if (neighbor_id == u_id) {
                        boundary_vertices_u.push_update(neighbor, new_qap_delta);
                    } else {
                        boundary_vertices_v.push_update(neighbor, new_qap_delta);
                    }
                }

                // remove vertex from u if it is not boundary
                while (!boundary_vertices_u.empty() && !is_connected_to(g, p_manager, boundary_vertices_u.top_key(), v_id)) { boundary_vertices_u.pop(); }

                // remove vertex from v if it is not boundary
                while (!boundary_vertices_v.empty() && !is_connected_to(g, p_manager, boundary_vertices_v.top_key(), u_id)) { boundary_vertices_v.pop(); }
            }


            HEIPROMAP_PROFILE_SCOPE("refinement", "QuotientGraphRefinement", "make_moves");
            // revert all moves in partitioning manager
            for (size_t i = 0; i < moves.size(); i++) {
                vertex_t vertex = moves[moves.size() - 1 - i];
                weight_t vertex_weight = t_uniform_v_weights ? 1 : g.v_weights[vertex];
                partition_t vertex_id = p_manager[vertex];
                partition_t move_id = u_id == vertex_id ? v_id : u_id;
                vertex_used[vertex] = global_vertex_mark - 1;

                p_manager.move_serial(vertex, vertex_weight, vertex_id, move_id);
            }

            // make all moves to best index
            for (size_t i = 0; i < best_idx; ++i) {
                vertex_t vertex = moves[i];
                weight_t vertex_weight = t_uniform_v_weights ? 1 : g.v_weights[vertex];
                partition_t vertex_id = p_manager[vertex];
                partition_t move_id = u_id == vertex_id ? v_id : u_id;
                vertex_used[vertex] = vertex_mark;

                bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                block_conn.move(g, vertex, vertex_id, move_id);
                p_manager.move_serial(vertex, vertex_weight, vertex_id, move_id);
            }

            if (max_qap_gain > 0) {
                active_next_round[u_id] = 1;
                active_next_round[v_id] = 1;
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void refine_blocks_edge_cut(const graph_t &g,
                                    d_oracle_t &d_oracle,
                                    bv_manager_t &bv_manager,
                                    p_manager_t &p_manager,
                                    q_graph_t &q_graph,
                                    block_conn_t &block_conn,
                                    partition_t u_id,
                                    partition_t v_id,
                                    std::vector<vertex_t> &moves,
                                    const AlignedArray<weight_t> &lmax_constraints,
                                    IndexedMaxHeap<weight_t> &boundary_vertices_u,
                                    IndexedMaxHeap<weight_t> &boundary_vertices_v,
                                    vertex_t vertex_mark,
                                    RandomEngine &random_engine) {
            f64 alpha = config->alpha * (f64) d_oracle.get(u_id, v_id);
            f64 beta = std::log(g.n) * (f64) d_oracle.get(u_id, v_id);

            HEIPROMAP_PROFILE_SCOPE("refinement", "QuotientGraphRefinement", "initial_edge_cut");
            // add all boundary vertices with gain
            u64 n_init_moves = 0;
            boundary_vertices_u.clear();
            boundary_vertices_v.clear();
            for (size_t j = 0; j < bv_manager.size(u_id); ++j) {
                const vertex_t u = bv_manager.get(u_id, j);
                for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) {
                    if (block_conn.get_id(i) == v_id) {
                        weight_t qap_delta_u = get_u_edge_cut_delta_t<t_uniform_e_weights>(g, u, u_id, v_id, p_manager, block_conn);
                        boundary_vertices_u.push(u, qap_delta_u);
                        n_init_moves += 1;
                        break;
                    }
                }
            }

            for (size_t j = 0; j < bv_manager.size(v_id); ++j) {
                const vertex_t v = bv_manager.get(v_id, j);
                for (size_t i = block_conn.start(v); i < block_conn.end(v); ++i) {
                    if (block_conn.get_id(i) == u_id) {
                        weight_t qap_delta_v = get_u_edge_cut_delta_t<t_uniform_e_weights>(g, v, v_id, u_id, p_manager, block_conn);
                        boundary_vertices_v.push(v, qap_delta_v);
                        n_init_moves += 1;
                        break;
                    }
                }
            }

            // store change
            weight_t curr_qap_gain = 0;
            weight_t max_qap_gain = 0;
            size_t best_idx = 0;

            u64 steps_since_last_improvement = 0;
            f64 qap_gain_mean = 0.0;
            f64 qap_gain_var = 0.0;

            weight_t u_id_weight_initial = p_manager.get_bweight(u_id);
            weight_t v_id_weight_initial = p_manager.get_bweight(v_id);
            bool best_is_balanced = (u_id_weight_initial <= lmax_constraints[u_id] && v_id_weight_initial <= lmax_constraints[v_id]);

            moves.clear();

            HEIPROMAP_PROFILE_SCOPE("refinement", "QuotientGraphRefinement", "process_queue_edge_cut");
            while ((!boundary_vertices_u.empty() || !boundary_vertices_v.empty()) && moves.size() < n_init_moves) {
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
                    weight_t v_id_weight = p_manager.get_bweight(v_id);

                    if (u_id_weight > lmax_constraints[u_id] && u_id_weight > v_id_weight) { choose_u = true; }
                    if (v_id_weight > lmax_constraints[v_id] && v_id_weight > u_id_weight) { choose_u = false; }
                    if (u_id_weight > lmax_constraints[u_id] && v_id_weight > lmax_constraints[v_id] && u_id_weight == v_id_weight) { choose_u = random_engine.get_f32() < 0.5; }
                }

                // choose the priority queue
                IndexedMaxHeap<weight_t> &boundary_vertices = choose_u ? boundary_vertices_u : boundary_vertices_v;

                vertex_t vertex = boundary_vertices.top_key();
                weight_t qap_delta = boundary_vertices.top();
                weight_t vertex_weight = t_uniform_v_weights ? 1 : g.v_weights[vertex];
                partition_t vertex_id = choose_u ? u_id : v_id;
                partition_t move_id = choose_u ? v_id : u_id;
                boundary_vertices.pop();

                // move the vertex
                moves.push_back(vertex);
                curr_qap_gain += qap_delta;

                weight_t current_move_id_weight = p_manager.get_bweight(move_id) + vertex_weight;
                weight_t current_vertex_id_weight = p_manager.get_bweight(vertex_id) - vertex_weight;
                bool current_is_balanced = (current_move_id_weight <= lmax_constraints[move_id] && current_vertex_id_weight <= lmax_constraints[vertex_id]);

                bool update_best = false;
                if (current_is_balanced) {
                    if (!best_is_balanced || curr_qap_gain >= max_qap_gain) {
                        update_best = true;
                    }
                } else {
                    if (!best_is_balanced) {
                        weight_t move_id_weight_initial = move_id == u_id ? u_id_weight_initial : v_id_weight_initial;
                        if (curr_qap_gain >= max_qap_gain && current_move_id_weight <= std::max(move_id_weight_initial, lmax_constraints[move_id])) {
                            update_best = true;
                        }
                    }
                }

                if (update_best) {
                    best_idx = moves.size();
                    max_qap_gain = curr_qap_gain;
                    best_is_balanced = current_is_balanced;

                    steps_since_last_improvement = 0;
                    qap_gain_mean = 0.0;
                    qap_gain_var = 0.0;
                }

                // make move in structures
                p_manager.move_serial(vertex, vertex_weight, vertex_id, move_id);
                vertex_used[vertex] = vertex_mark;

                steps_since_last_improvement += 1;
                f64 new_qap_gain_mean = qap_gain_mean + ((f64) qap_delta - qap_gain_mean) / (f64) steps_since_last_improvement;
                f64 new_qap_gain_var = (qap_gain_var + ((f64) qap_delta - qap_gain_mean) * ((f64) qap_delta - new_qap_gain_mean)) / (f64) steps_since_last_improvement;

                qap_gain_mean = new_qap_gain_mean;
                qap_gain_var = new_qap_gain_var;

                if (config->use_preemptive_exit) {
                    if (steps_since_last_improvement > config->min_n_steps && (f64) steps_since_last_improvement * qap_gain_mean * qap_gain_mean > alpha * qap_gain_var + beta) {
                        break;
                    }
                }

                // we have to push or update the neighbors that were not moved already
                for (size_t i = g.neighborhoods[vertex]; i < g.neighborhoods[vertex + 1]; ++i) {
                    const vertex_t neighbor = g.edges_v[i];

                    if (vertex_used[neighbor] == vertex_mark) { continue; }

                    partition_t neighbor_id = p_manager[neighbor];

                    if (neighbor_id != u_id && neighbor_id != v_id) { continue; }

                    partition_t new_id = neighbor_id == vertex_id ? move_id : vertex_id;

                    bool is_connected_to_new_id;
                    weight_t new_qap_delta = get_u_edge_cut_delta_and_is_connected_to_t<t_uniform_e_weights>(g, neighbor, neighbor_id, new_id, is_connected_to_new_id, p_manager);

                    if (!is_connected_to_new_id) { continue; }

                    if (neighbor_id == u_id) {
                        boundary_vertices_u.push_update(neighbor, new_qap_delta);
                    } else {
                        boundary_vertices_v.push_update(neighbor, new_qap_delta);
                    }
                }

                // remove vertex from u if it is not boundary
                while (!boundary_vertices_u.empty() && !is_connected_to(g, p_manager, boundary_vertices_u.top_key(), v_id)) { boundary_vertices_u.pop(); }

                // remove vertex from v if it is not boundary
                while (!boundary_vertices_v.empty() && !is_connected_to(g, p_manager, boundary_vertices_v.top_key(), u_id)) { boundary_vertices_v.pop(); }
            }

            HEIPROMAP_PROFILE_SCOPE("refinement", "QuotientGraphRefinement", "make_moves_edge_cut");
            // revert all moves in partitioning manager
            for (size_t i = 0; i < moves.size(); i++) {
                vertex_t vertex = moves[moves.size() - 1 - i];
                weight_t vertex_weight = t_uniform_v_weights ? 1 : g.v_weights[vertex];
                partition_t vertex_id = p_manager[vertex];
                partition_t move_id = u_id == vertex_id ? v_id : u_id;
                vertex_used[vertex] = global_vertex_mark - 1;

                p_manager.move_serial(vertex, vertex_weight, vertex_id, move_id);
            }

            // make all moves to best index
            for (size_t i = 0; i < best_idx; ++i) {
                vertex_t vertex = moves[i];
                weight_t vertex_weight = t_uniform_v_weights ? 1 : g.v_weights[vertex];
                partition_t vertex_id = p_manager[vertex];
                partition_t move_id = u_id == vertex_id ? v_id : u_id;
                vertex_used[vertex] = vertex_mark;

                bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                block_conn.move(g, vertex, vertex_id, move_id);
                p_manager.move_serial(vertex, vertex_weight, vertex_id, move_id);
            }

            if (max_qap_gain > 0) {
                active_next_round[u_id] = 1;
                active_next_round[v_id] = 1;
            }
        }
    };
}

#endif //HEIPROMAP_QUOTIENT_GRAPH_REFINEMENT_H
