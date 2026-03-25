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

#include "../utility/utils.h"
#include "ISerialRefiner.h"
#include "../utility/qap.h"
#include "../datastructures/indexed_max_heap.h"
#include "../utility/functions.h"

namespace HeiProMap {
    struct PairWeight {
        partition_t id1;
        partition_t id2;
        weight_t distance;

        bool operator<(const PairWeight &p) const {
            return distance < p.distance;
        }

        bool operator>(const PairWeight &p) const {
            return distance > p.distance;
        }
    };

    class QuotientGraphRefinementConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit QuotientGraphRefinementConfiguration(const std::string &t_name) : ISerialRefinerConfiguration(t_name) {
        }

        u64 max_iteration = 1;
        u64 min_n_steps = 4;
        f64 alpha = 100.0;
        f64 beta = 1.0;
        bool use_active_scheduling = true;
        bool use_preemptive_exit = true;
    };

    class QuotientGraphRefinement final : public ISerialRefiner {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
        u64 m_threads = 1;

        // active block scheduling
        AlignedArray<u8> active_this_round;
        AlignedArray<u8> active_next_round;
        AlignedArray<PairWeight> pairs;
        size_t pairs_size = 0;

        // store which vertices have been moved
        AlignedArray<u32> vertex_used;
        u32 global_vertex_mark = 0;

        std::vector<IndexedMaxHeap<s64> > boundary_vertices_u_vec;
        std::vector<IndexedMaxHeap<s64> > boundary_vertices_v_vec;

        const QuotientGraphRefinementConfiguration *config = nullptr;

    public:
        QuotientGraphRefinement() = default;

        ~QuotientGraphRefinement() override = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const u64 t_threads,
                        const ISerialRefinerConfiguration &i_config) override {
            ScopedTimer _t("io", "QuotientGraphRefinement", "initialize");

            m_n = t_n;
            m_m = t_m;
            m_k = t_k;
            m_threads = t_threads;

            config = dynamic_cast<const QuotientGraphRefinementConfiguration *>(&i_config);

            global_vertex_mark = 0;
            vertex_used.initialize(m_n, 0);

            // active block scheduling
            active_this_round.initialize(m_k);
            active_next_round.initialize(m_k);
            pairs.initialize(m_k * m_k);
            pairs_size = 0;

            for (size_t i = 0; i < m_threads; ++i) {
                boundary_vertices_u_vec.push_back(IndexedMaxHeap<s64>());
                boundary_vertices_v_vec.push_back(IndexedMaxHeap<s64>());
            }

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
                    f64 imbalance) override {
            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            RandomEngine random_engine = RandomEngine(0);

            active_this_round.initialize(m_k, 1);
            active_next_round.initialize(m_k, 0);

            for (u64 iteration = 0; iteration < config->max_iteration; ++iteration) {
                // determine all pairs in the quotient graph
                {
                    ScopedTimer _t("refinement", "QuotientGraphRefinement", "get_pairs");

                    pairs_size = 0;

                    for (partition_t u_id = 0; u_id < m_k; ++u_id) {
                        for (partition_t v_id = u_id + 1; v_id < m_k; ++v_id) {
                            if (!q_graph.has_edge(u_id, v_id)) continue;

                            if (config->use_active_scheduling && !(active_this_round[u_id] || active_this_round[v_id])) {
                                continue;
                            }

                            pairs[pairs_size++] = {u_id, v_id, d_oracle.get(u_id, v_id)};
                        }
                    }
                    if (pairs_size == 0) { return; }
                }
                //
                {
                    ScopedTimer _t("refinement", "QuotientGraphRefinement", "shuffle");

                    fast_shuffle_unchecked(pairs.get_ptr(), pairs.get_ptr() + pairs_size, random_engine.generator);
                }

                for (size_t j = 0; j < pairs_size; ++j) {
                    auto [u_id, v_id, distance] = pairs[j];
                    // if (d_oracle.last_level_pair(u_id, v_id)) {
                    //     refine_blocks_edge_cut(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, u_id, v_id);
                    // } else {

                    // priority queues
                    IndexedMaxHeap<s64> &boundary_vertices_u = boundary_vertices_u_vec[0];
                    IndexedMaxHeap<s64> &boundary_vertices_v = boundary_vertices_v_vec[0];

                    global_vertex_mark += 1;
                    refine_blocks(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, u_id, v_id, lmax, boundary_vertices_u, boundary_vertices_v, global_vertex_mark, random_engine);
                    // }
                }

                //
                {
                    ScopedTimer _t("refinement", "QuotientGraphRefinement", "swap");

                    std::swap(active_this_round, active_next_round);
                    active_next_round.initialize(m_k, 0);
                }
            }

            /*
            AlignedArray<u8> used_this_round;
            //
            {
                ScopedTimer _t("refinement", "QuotientGraphRefinement", "allocate");

                used_this_round.initialize(m_k * m_k);
            }

            std::vector<std::pair<partition_t, partition_t> > matching;

            for (u64 iteration = 0; iteration < config->max_iteration; ++iteration) {
                //
                {
                    ScopedTimer _t("refinement", "QuotientGraphRefinement", "reset_used_edges");

                    std::fill_n(used_this_round.get_ptr(), m_k * m_k, 0);
                }

                bool found_matching = false;
                //
                {
                    ScopedTimer _t("refinement", "QuotientGraphRefinement", "matching");
                    found_matching = q_graph.find_distance_3_matching(active_this_round, used_this_round, matching);
                }

                while (found_matching) {
                    #pragma omp parallel for num_threads(m_threads) schedule(dynamic)
                    for (size_t i = 0; i < matching.size(); ++i) {
                        partition_t u_id = matching[i].first;
                        partition_t v_id = matching[i].second;

                        u64 thread_id = omp_get_thread_num();

                        // priority queues
                        IndexedMaxHeap<s64> &boundary_vertices_u = boundary_vertices_u_vec[thread_id];
                        IndexedMaxHeap<s64> &boundary_vertices_v = boundary_vertices_v_vec[thread_id];

                        refine_blocks(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, u_id, v_id, lmax, boundary_vertices_u, boundary_vertices_v, global_vertex_mark);
                    }

                    //
                    {
                        ScopedTimer _t("refinement", "QuotientGraphRefinement", "matching");
                        found_matching = q_graph.find_distance_3_matching(active_this_round, used_this_round, matching);
                    }
                }

                // swap active
                {
                    ScopedTimer _t("refinement", "QuotientGraphRefinement", "swap_active");
                    std::swap(active_this_round, active_next_round);
                    active_next_round.initialize(m_k, 0);
                }
            }
            */
        }

        void refine_blocks(const graph_t &g,
                           d_oracle_t &d_oracle,
                           bv_manager_t &bv_manager,
                           p_manager_t &p_manager,
                           q_graph_t &q_graph,
                           block_conn_t &block_conn,
                           partition_t u_id,
                           partition_t v_id,
                           weight_t lmax,
                           IndexedMaxHeap<s64> &boundary_vertices_u,
                           IndexedMaxHeap<s64> &boundary_vertices_v,
                           vertex_t vertex_mark,
                           RandomEngine &random_engine) {
            f64 alpha = config->alpha;
            f64 beta = std::log(g.n);

            size_t max_n_swaps = 0;
            //
            {
                ScopedTimer _t("refinement", "QuotientGraphRefinement", "initial_qap");

                // add all boundary vertices with gain
                boundary_vertices_u.clear();
                boundary_vertices_v.clear();
                forall_bv_id_iu(bv_manager, u_id, j, u)
                    {
                        forall_bc_ui_id(block_conn, u, i, id)
                            {
                                if (id == v_id) {
                                    weight_t qap_delta_u = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle, block_conn);
                                    boundary_vertices_u.push(u, qap_delta_u);
                                    break;
                                }
                            }
                        endfor
                    }
                endfor

                forall_bv_id_iu(bv_manager, v_id, j, v)
                    {
                        forall_bc_ui_id(block_conn, v, i, id)
                            {
                                if (id == u_id) {
                                    s64 qap_delta_v = get_u_qap_delta(g, v, v_id, u_id, p_manager, d_oracle, block_conn);
                                    boundary_vertices_v.push(v, qap_delta_v);
                                    break;
                                }
                            }
                        endfor
                    }
                endfor
            }

            max_n_swaps = boundary_vertices_u.size() + boundary_vertices_v.size();

            // store change
            s64 curr_qap_gain = 0;
            s64 max_qap_gain = 0;
            // s64 curr_edge_cut_gain = 0;
            // s64 max_edge_cut_gain = 0;
            size_t best_idx = 0;

            u64 steps_since_last_improvement = 0;
            f64 qap_gain_mean = 0.0;
            f64 qap_gain_var = 0.0;

            std::vector<vertex_t> moves;
            //
            {
                ScopedTimer _t("refinement", "QuotientGraphRefinement", "process_queue");

                while ((!boundary_vertices_u.empty() || !boundary_vertices_v.empty()) && moves.size() < max_n_swaps) {
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

                        if (u_id_weight > lmax && u_id_weight > v_id_weight) { choose_u = true; }
                        if (v_id_weight > lmax && v_id_weight > u_id_weight) { choose_u = false; }
                        if (u_id_weight > lmax && v_id_weight > lmax && u_id_weight == v_id_weight) { choose_u = random_engine.get_f32() < 0.5; }
                    }

                    // choose the priority queue
                    IndexedMaxHeap<s64> &boundary_vertices = choose_u ? boundary_vertices_u : boundary_vertices_v;

                    vertex_t vertex = boundary_vertices.top_key();
                    s64 qap_delta = boundary_vertices.top();
                    weight_t vertex_weight = g.v_weights[vertex];
                    partition_t vertex_id = choose_u ? u_id : v_id;
                    partition_t move_id = choose_u ? v_id : u_id;
                    boundary_vertices.pop();

                    // move the vertex
                    moves.push_back(vertex);
                    curr_qap_gain += qap_delta;
                    if (curr_qap_gain >= max_qap_gain && p_manager.get_bweight(move_id) + vertex_weight <= lmax && p_manager.get_bweight(vertex_id) - vertex_weight <= lmax) {
                        best_idx = moves.size();
                        max_qap_gain = curr_qap_gain;

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
            }
            //
            {
                ScopedTimer _t("refinement", "QuotientGraphRefinement", "make_moves");

                // revert all moves in partitioning manager
                for (size_t i = 0; i < moves.size(); i++) {
                    vertex_t vertex = moves[moves.size() - 1 - i];
                    weight_t vertex_weight = g.v_weights[vertex];
                    partition_t vertex_id = p_manager[vertex];
                    partition_t move_id = u_id == vertex_id ? v_id : u_id;
                    vertex_used[vertex] = vertex_mark - 1;

                    p_manager.move_serial(vertex, vertex_weight, vertex_id, move_id);
                }

                // make all moves to best index
                for (size_t i = 0; i < best_idx; ++i) {
                    vertex_t vertex = moves[i];
                    weight_t vertex_weight = g.v_weights[vertex];
                    partition_t vertex_id = p_manager[vertex];
                    partition_t move_id = u_id == vertex_id ? v_id : u_id;
                    vertex_used[vertex] = vertex_mark;

                    bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                    q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                    block_conn.move(g, vertex, vertex_id, move_id);
                    p_manager.move_serial(vertex, vertex_weight, vertex_id, move_id);
                }
            }

            if (max_qap_gain > 0) {
                active_next_round[u_id] = 1;
                active_next_round[v_id] = 1;
            }
        }

        /*
        void refine_blocks_edge_cut([[maybe_unused]] const u64 level,
                                    [[maybe_unused]] const u64 max_level,
                                    const graph_t &g,
                                    [[maybe_unused]] const d_oracle_t &d_oracle,
                                    bv_manager_t &bv_manager,
                                    p_manager_t &p_manager,
                                    q_graph_t &q_graph,
                                    block_conn_t &block_conn,
                                    partition_t u_id,
                                    partition_t v_id,
                                    weight_t lmax) {
            f64 alpha = config->alpha;
            f64 beta = std::log(g.n);

            size_t max_n_swaps = 0;

            // priority queues
            IndexedMaxHeap<s64> boundary_vertices_u;
            IndexedMaxHeap<s64> boundary_vertices_v;
            boundary_vertices_u.initialize(m_n);
            boundary_vertices_v.initialize(m_n);

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

            u64 steps_since_last_improvement = 0;
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
                    weight_t v_id_weight = p_manager.get_bweight(v_id);

                    if (u_id_weight > lmax && u_id_weight > v_id_weight) { choose_u = true; }
                    if (v_id_weight > lmax && v_id_weight > u_id_weight) { choose_u = false; }
                    if (u_id_weight > lmax && v_id_weight > lmax && u_id_weight == v_id_weight) { choose_u = random_engine.get_f32() < 0.5; }
                }

                // choose the priority queue
                IndexedMaxHeap<s64> &boundary_vertices = choose_u ? boundary_vertices_u : boundary_vertices_v;

                vertex_t vertex = boundary_vertices.top_key();
                s64 edge_cut_delta = boundary_vertices.top();
                weight_t vertex_weight = g.v_weights[vertex];
                partition_t vertex_id = choose_u ? u_id : v_id;
                partition_t move_id = choose_u ? v_id : u_id;
                boundary_vertices.pop();

                // move the vertex
                moves[moves_size++] = vertex;
                curr_edge_cut_gain += edge_cut_delta;
                if (curr_edge_cut_gain >= max_edge_cut_gain && p_manager.get_bweight(move_id) + vertex_weight <= lmax && p_manager.get_bweight(vertex_id) - vertex_weight <= lmax) {
                    best_idx = moves_size;
                    max_edge_cut_gain = curr_edge_cut_gain;

                    steps_since_last_improvement = 0;
                    edge_cut_gain_mean = 0.0;
                    edge_cut_gain_var = 0.0;
                }

                // make move in structures
                p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                vertex_used[vertex] = vertex_mark;

                steps_since_last_improvement += 1;
                f64 new_edge_cut_gain_mean = edge_cut_gain_mean + ((f64) edge_cut_delta - edge_cut_gain_mean) / (f64) steps_since_last_improvement;
                f64 new_edge_cut_gain_var = (edge_cut_gain_var + ((f64) edge_cut_delta - edge_cut_gain_mean) * ((f64) edge_cut_delta - new_edge_cut_gain_mean)) / (f64) steps_since_last_improvement;

                edge_cut_gain_mean = new_edge_cut_gain_mean;
                edge_cut_gain_var = new_edge_cut_gain_var;

                if (steps_since_last_improvement > config->min_n_steps && (f64) steps_since_last_improvement * edge_cut_gain_mean * edge_cut_gain_mean > alpha * edge_cut_gain_var + beta) { break; }

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
                weight_t vertex_weight = g.v_weights[vertex];
                partition_t vertex_id = p_manager[vertex];
                partition_t move_id = u_id == vertex_id ? v_id : u_id;
                vertex_used[vertex] = vertex_mark - 1;

                p_manager.move(vertex, vertex_weight, vertex_id, move_id);
            }

            // make all moves to best index
            for (size_t i = 0; i < best_idx; ++i) {
                vertex_t vertex = moves[i];
                weight_t vertex_weight = g.v_weights[vertex];
                partition_t vertex_id = p_manager[vertex];
                partition_t move_id = u_id == vertex_id ? v_id : u_id;
                vertex_used[vertex] = vertex_mark;

                bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                block_conn.move(g, vertex, vertex_id, move_id);
                p_manager.move(vertex, vertex_weight, vertex_id, move_id);
            }

            if (max_edge_cut_gain > 0) {
                active_next_round[u_id] = 1;
                active_next_round[v_id] = 1;
            }
        }
        */
    };
}

#endif //HEIPROMAP_QUOTIENT_GRAPH_REFINEMENT_H
