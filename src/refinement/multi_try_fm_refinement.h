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

#include "../utility/indexed_update_heap.h"
#include "../utility/utils.h"
#include "../utility/functions.h"
#include "ISerialRefiner.h"
#include "../utility/qap.h"
#include "../utility/random_engine.h"

namespace HeiProMap {
    class MultiTryFmRefinementConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit MultiTryFmRefinementConfiguration(const std::string &t_name) : ISerialRefinerConfiguration(t_name) {
        }

        u64 max_iteration = 1;
        f64 alpha = 10.0;
        f64 beta = 1.0;
        u64 min_n_steps = 2;
    };

    class MultiTryFMRefinement final : public ISerialRefiner {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;

        AlignedArray<u32> vertex_used;
        u32 vertex_mark = 0;

        AlignedArray<u32> block_used;
        u32 block_mark = 0;

        AlignedArray<vertex_t> curr_boundary;
        size_t curr_boundary_size = 0;

        AlignedArray<Move> moves;
        size_t moves_size = 0;

        // priority queues
        IndexedUpdateHeap heap;

        RandomEngine random_engine;
        const MultiTryFmRefinementConfiguration *config = nullptr;

    public:
        MultiTryFMRefinement() = default;

        ~MultiTryFMRefinement() override = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const u64 t_threads,
                        const u64 seed,
                        const ISerialRefinerConfiguration &i_config) override {
            m_n = t_n;
            m_m = t_m;
            m_k = t_k;

            config = dynamic_cast<const MultiTryFmRefinementConfiguration *>(&i_config);

            vertex_used.initialize(m_n, 0);
            vertex_mark = 0;

            block_used.initialize(m_k, 0);
            block_mark = 0;

            curr_boundary.initialize(m_n);
            curr_boundary_size = 0;

            moves.initialize(m_n);
            moves_size = 0;

            heap.initialize(m_n);

            random_engine = RandomEngine(seed);
        }

        void refine(graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph,
                    block_conn_t &block_conn,
                    f64 imbalance,
                    bool uniform_v_weights,
                    bool uniform_e_weights) override {
            if (uniform_v_weights && uniform_e_weights)      refine_impl<true, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, imbalance);
            else if (uniform_v_weights)                      refine_impl<true, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, imbalance);
            else if (uniform_e_weights)                      refine_impl<false, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, imbalance);
            else                                             refine_impl<false, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, imbalance);
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void refine_impl(graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph,
                    block_conn_t &block_conn,
                    f64 imbalance) {
            f64 alpha = config->alpha;
            f64 beta = std::log(g.n);
            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            bool positive_move_occurred = true;
            for (u64 iteration = 0; iteration < config->max_iteration && positive_move_occurred; ++iteration) {
                positive_move_occurred = false;
                // collect boundary
                {
                    ScopedTimer _t("refinement", "MultiTryFMRefinement", "collect_boundary");
                    curr_boundary_size = 0;
                    for (partition_t id = 0; id < m_k; ++id) {
                        for (size_t i = 0; i < bv_manager.size(id); ++i) { const vertex_t u = bv_manager.get(id, i);
                            {
                                curr_boundary[curr_boundary_size++] = u;
                            }
                        }
                    }
                    std::shuffle(curr_boundary.get_ptr(), curr_boundary.get_ptr() + curr_boundary_size, random_engine.generator);
                }

                vertex_mark += 1;
                for (size_t ii = 0; ii < curr_boundary_size; ++ii) {
                    vertex_t u = curr_boundary[ii];
                    if (vertex_used[u] == vertex_mark) { continue; }
                    if (!bv_manager.is_boundary(u)) { continue; }

                    heap.clear();

                    // insert u into the priority queue
                    partition_t u_id = p_manager[u];
                    weight_t u_weight = g.v_weights[u];

                    // initial qap for all blocks
                    {
                        ScopedTimer _t("refinement", "MultiTryFMRefinement", "initial_block_qap");

                        // find all connected partitions to u
                        partition_t best_v_id = NO_ID;
                        weight_t best_qap_delta = -std::numeric_limits<weight_t>::max();

                        for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) { const partition_t id = block_conn.get_id(i);
                            {
                                if (id == u_id) { continue; }
                                if (p_manager.get_bweight(id) + u_weight > lmax) { continue; }

                                weight_t qap_delta = get_u_qap_delta(g, u, u_id, id, p_manager, d_oracle, block_conn);

                                if (qap_delta > best_qap_delta) {
                                    best_qap_delta = qap_delta;
                                    best_v_id = id;
                                }
                            }
                        }

                        if (best_v_id != NO_ID) { heap.push(u, best_v_id, best_qap_delta); }
                    }

                    // insert all neighbors of u that are boundary into the queue
                    {
                        ScopedTimer _t("refinement", "MultiTryFMRefinement", "initial_boundary_qap");

                        for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) { const vertex_t neighbor = g.edges_v[i];
                            {
                                if (vertex_used[neighbor] == vertex_mark) { continue; }
                                if (!bv_manager.is_boundary(neighbor)) { continue; }

                                partition_t neighbor_id = p_manager[neighbor];
                                weight_t neighbor_weight = g.v_weights[neighbor];

                                // find all connected partitions to neighbor
                                partition_t best_v_id = NO_ID;
                                weight_t best_qap_delta = -std::numeric_limits<weight_t>::max();

                                for (size_t j = block_conn.start(neighbor); j < block_conn.end(neighbor); ++j) { const partition_t id = block_conn.get_id(j);
                                    {
                                        if (id == neighbor_id) { continue; }
                                        if (p_manager.get_bweight(id) + neighbor_weight > lmax) { continue; }

                                        weight_t u_qap_delta = get_u_qap_delta(g, neighbor, neighbor_id, id, p_manager, d_oracle, block_conn);

                                        if (u_qap_delta > best_qap_delta) {
                                            best_qap_delta = u_qap_delta;
                                            best_v_id = id;
                                        }
                                    }
                                }

                                if (best_v_id != NO_ID) { heap.push(neighbor, best_v_id, best_qap_delta); }
                            }
                        }
                    }

                    if (heap.empty()) { continue; }

                    moves_size = 0;
                    size_t best_idx = 0;
                    weight_t max_qap_gain = 0;
                    // process the queue
                    {
                        ScopedTimer _t("refinement", "MultiTryFMRefinement", "process_queue");

                        weight_t curr_qap_gain = 0;
                        u64 steps_since_last_improvement = 0;
                        f64 qap_gain_mean = 0.0;
                        f64 qap_gain_var = 0.0;

                        while (!heap.empty()) {
                            vertex_t vertex = heap.top_u();
                            partition_t vertex_id = p_manager[vertex];
                            weight_t vertex_weight = g.v_weights[vertex];
                            partition_t move_id = heap.top_id();
                            weight_t move_qap_delta = heap.top_qap_delta();
                            heap.pop();

                            if (p_manager.get_bweight(move_id) + vertex_weight > lmax) { continue; }

                            moves[moves_size++] = Move(vertex, vertex_id, move_id);
                            curr_qap_gain += move_qap_delta;
                            if (curr_qap_gain >= max_qap_gain) {
                                best_idx = moves_size;
                                max_qap_gain = curr_qap_gain;

                                steps_since_last_improvement = 0;
                                qap_gain_mean = 0.0;
                                qap_gain_var = 0.0;
                            }

                            // make move in structures
                            p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                            vertex_used[vertex] = vertex_mark;

                            steps_since_last_improvement += 1;
                            f64 new_qap_gain_mean = qap_gain_mean + ((f64) move_qap_delta - qap_gain_mean) / (f64) steps_since_last_improvement;
                            f64 new_qap_gain_var = (qap_gain_var + ((f64) move_qap_delta - qap_gain_mean) * ((f64) move_qap_delta - new_qap_gain_mean)) / (f64) steps_since_last_improvement;

                            qap_gain_mean = new_qap_gain_mean;
                            qap_gain_var = new_qap_gain_var;

                            if (steps_since_last_improvement > config->min_n_steps && (f64) steps_since_last_improvement * qap_gain_mean * qap_gain_mean > alpha * qap_gain_var + beta) { break; }

                            // we have to push or update the neighbors that were not moved already
                            for (size_t i = g.neighborhoods[vertex]; i < g.neighborhoods[vertex + 1]; ++i) { const vertex_t neighbor = g.edges_v[i];
                                {
                                    if (vertex_used[neighbor] == vertex_mark) { continue; }
                                    if (!is_boundary(g, p_manager, neighbor)) {
                                        if (heap.entry_exists(neighbor)) { heap.remove(neighbor); }
                                        continue;
                                    }

                                    partition_t neighbor_id = p_manager[neighbor];
                                    weight_t neighbor_weight = g.v_weights[neighbor];

                                    partition_t best_v_id = NO_ID;
                                    weight_t best_qap_delta = -std::numeric_limits<weight_t>::max();

                                    block_mark += 1;
                                    for (size_t j = g.neighborhoods[neighbor]; j < g.neighborhoods[neighbor + 1]; ++j) { const vertex_t v = g.edges_v[j];
                                        {
                                            partition_t v_id = p_manager[v];
                                            if (v_id == neighbor_id) { continue; }
                                            if (block_used[v_id] == block_mark) { continue; }
                                            if (p_manager.get_bweight(v_id) + neighbor_weight > lmax) { continue; }

                                            weight_t v_qap_delta = get_u_qap_delta(g, neighbor, neighbor_id, v_id, p_manager, d_oracle);
                                            block_used[v_id] = block_mark;

                                            if (v_qap_delta > best_qap_delta) {
                                                best_qap_delta = v_qap_delta;
                                                best_v_id = v_id;
                                            }
                                        }
                                    }

                                    if (best_v_id != NO_ID) { heap.push_update(neighbor, best_v_id, best_qap_delta); }
                                }
                            }
                        }
                    }

                    // revert all moves in partitioning manager
                    {
                        ScopedTimer _t("refinement", "MultiTryFMRefinement", "revert_moves");

                        for (size_t i = 0; i < moves_size; i++) {
                            vertex_t vertex = moves[moves_size - 1 - i].u;
                            weight_t vertex_weight = g.v_weights[vertex];
                            partition_t vertex_id = moves[moves_size - 1 - i].to_move_id;
                            partition_t move_id = moves[moves_size - 1 - i].u_id;
                            vertex_used[vertex] = vertex_mark - 1;

                            p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                        }
                    }

                    // make all moves to best index
                    {
                        ScopedTimer _t("refinement", "MultiTryFMRefinement", "make_moves");

                        for (size_t i = 0; i < best_idx; ++i) {
                            vertex_t vertex = moves[i].u;
                            weight_t vertex_weight = g.v_weights[vertex];
                            partition_t vertex_id = moves[i].u_id;
                            partition_t move_id = moves[i].to_move_id;
                            vertex_used[vertex] = vertex_mark;

                            bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                            q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                            block_conn.move(g, vertex, vertex_id, move_id);
                            p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                        }
                    }

                    positive_move_occurred |= max_qap_gain > 0;
                }
            }
        }

        void refine_new(graph_t &g,
                        d_oracle_t &d_oracle,
                        bv_manager_t &bv_manager,
                        p_manager_t &p_manager,
                        q_graph_t &q_graph,
                        block_conn_t &block_conn,
                        f64 imbalance) {
            f64 alpha = config->alpha;
            f64 beta = std::log(g.n);
            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            bool positive_move_occurred = true;
            for (u64 iteration = 0; iteration < config->max_iteration && positive_move_occurred; ++iteration) {
                positive_move_occurred = false;
                // collect boundary
                {
                    ScopedTimer _t("refinement", "MultiTryFMRefinement", "collect_boundary");
                    curr_boundary_size = 0;
                    for (partition_t id = 0; id < m_k; ++id) {
                        for (size_t i = 0; i < bv_manager.size(id); ++i) { const vertex_t u = bv_manager.get(id, i);
                            {
                                curr_boundary[curr_boundary_size++] = u;
                            }
                        }
                    }
                    std::shuffle(curr_boundary.get_ptr(), curr_boundary.get_ptr() + curr_boundary_size, random_engine.generator);
                }

                vertex_mark += 1;
                for (size_t ii = 0; ii < curr_boundary_size; ++ii) {
                    vertex_t u = curr_boundary[ii];
                    if (vertex_used[u] == vertex_mark) { continue; }
                    if (!bv_manager.is_boundary(u)) { continue; }

                    heap.clear();

                    // insert u into the priority queue
                    partition_t u_id = p_manager[u];
                    weight_t u_weight = g.v_weights[u];

                    // find all connected partitions to u
                    partition_t best_v_id = 0;
                    weight_t best_qap_delta = -std::numeric_limits<weight_t>::max();

                    // initial qap for all blocks
                    {
                        ScopedTimer _t("refinement", "MultiTryFMRefinement", "initial_block_qap");

                        for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) { const partition_t id = block_conn.get_id(i);
                            {
                                if (id == u_id) { continue; }
                                if (p_manager.get_bweight(id) + u_weight > lmax) { continue; }

                                weight_t qap_delta = get_u_qap_delta(g, u, u_id, id, p_manager, d_oracle, block_conn);

                                if (qap_delta > best_qap_delta) {
                                    best_qap_delta = qap_delta;
                                    best_v_id = id;
                                }
                            }
                        }
                        if (best_qap_delta != -std::numeric_limits<weight_t>::max()) {
                            heap.push(u, best_v_id, best_qap_delta);
                        }
                    }

                    // insert all neighbors of u that are boundary into the queue
                    {
                        ScopedTimer _t("refinement", "MultiTryFMRefinement", "initial_boundary_qap");

                        for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) { const vertex_t neighbor = g.edges_v[i];
                            {
                                if (vertex_used[neighbor] == vertex_mark) { continue; }
                                if (!bv_manager.is_boundary(neighbor)) { continue; }

                                partition_t neighbor_id = p_manager[neighbor];
                                weight_t neighbor_weight = g.v_weights[neighbor];

                                best_qap_delta = -std::numeric_limits<weight_t>::max();
                                block_mark += 1;

                                for (size_t j = block_conn.start(neighbor); j < block_conn.end(neighbor); ++j) { const partition_t id = block_conn.get_id(j);
                                    {
                                        if (id == neighbor_id) { continue; }
                                        if (p_manager.get_bweight(id) + neighbor_weight > lmax) { continue; }

                                        weight_t u_qap_delta = get_u_qap_delta(g, neighbor, neighbor_id, id, p_manager, d_oracle, block_conn);

                                        if (u_qap_delta > best_qap_delta) {
                                            best_qap_delta = u_qap_delta;
                                            best_v_id = id;
                                        }
                                    }
                                }

                                if (best_qap_delta != -std::numeric_limits<weight_t>::max()) {
                                    heap.push(neighbor, best_v_id, best_qap_delta);
                                }
                            }
                        }
                    }

                    if (heap.empty()) { continue; }

                    moves_size = 0;
                    size_t best_idx = 0;
                    weight_t max_qap_gain = 0;
                    // process the queue
                    {
                        ScopedTimer _t("refinement", "MultiTryFMRefinement", "process_queue");

                        weight_t curr_qap_gain = 0;
                        u64 steps_since_last_improvement = 0;
                        f64 qap_gain_mean = 0.0;
                        f64 qap_gain_var = 0.0;

                        while (!heap.empty()) {
                            vertex_t vertex = heap.top_u();
                            partition_t vertex_id = p_manager[vertex];
                            weight_t vertex_weight = g.v_weights[vertex];
                            partition_t move_id = heap.top_id();
                            weight_t move_qap_delta = heap.top_qap_delta();
                            heap.pop();

                            if (p_manager.get_bweight(move_id) + vertex_weight > lmax) { continue; }

                            moves[moves_size++] = Move(vertex, vertex_id, move_id);
                            curr_qap_gain += move_qap_delta;
                            if (curr_qap_gain > max_qap_gain) {
                                best_idx = moves_size;
                                max_qap_gain = curr_qap_gain;

                                steps_since_last_improvement = 0;
                                qap_gain_mean = 0.0;
                                qap_gain_var = 0.0;
                            }

                            // make move in structures
                            bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                            q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                            block_conn.move(g, vertex, vertex_id, move_id);
                            p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                            vertex_used[vertex] = vertex_mark;

                            steps_since_last_improvement += 1;
                            f64 new_qap_gain_mean = qap_gain_mean + ((f64) move_qap_delta - qap_gain_mean) / (f64) steps_since_last_improvement;
                            f64 new_qap_gain_var = (qap_gain_var + ((f64) move_qap_delta - qap_gain_mean) * ((f64) move_qap_delta - new_qap_gain_mean)) / (f64) steps_since_last_improvement;

                            qap_gain_mean = new_qap_gain_mean;
                            qap_gain_var = new_qap_gain_var;

                            if (steps_since_last_improvement > config->min_n_steps && (f64) steps_since_last_improvement * qap_gain_mean * qap_gain_mean > alpha * qap_gain_var + beta) { break; }

                            // we have to push or update the neighbors that were not moved already
                            for (size_t i = g.neighborhoods[vertex]; i < g.neighborhoods[vertex + 1]; ++i) { const vertex_t neighbor = g.edges_v[i];
                                {
                                    if (vertex_used[neighbor] == vertex_mark) { continue; }
                                    if (!bv_manager.is_boundary(neighbor)) {
                                        if (heap.entry_exists(neighbor)) { heap.remove(neighbor); }
                                        continue;
                                    }

                                    partition_t neighbor_id = p_manager[neighbor];
                                    weight_t neighbor_weight = g.v_weights[neighbor];

                                    best_qap_delta = -std::numeric_limits<weight_t>::max();

                                    for (size_t j = block_conn.start(neighbor); j < block_conn.end(neighbor); ++j) { const partition_t id = block_conn.get_id(j);
                                        {
                                            if (id == neighbor_id) { continue; }
                                            if (p_manager.get_bweight(id) + neighbor_weight > lmax) { continue; }

                                            weight_t v_qap_delta = get_u_qap_delta(g, neighbor, neighbor_id, id, p_manager, d_oracle, block_conn);

                                            if (v_qap_delta > best_qap_delta) {
                                                best_qap_delta = v_qap_delta;
                                                best_v_id = id;
                                            }
                                        }
                                    }


                                    if (best_qap_delta != -std::numeric_limits<weight_t>::max()) { heap.push_update(neighbor, best_v_id, best_qap_delta); }
                                }
                            }
                        }
                    }

                    // revert all moves in partitioning manager
                    {
                        ScopedTimer _t("refinement", "MultiTryFMRefinement", "revert_moves");

                        for (size_t i = moves_size; i > best_idx; --i) {
                            vertex_t vertex = moves[i - 1].u;
                            weight_t vertex_weight = g.v_weights[vertex];
                            partition_t vertex_id = moves[i - 1].to_move_id;
                            partition_t move_id = moves[i - 1].u_id;

                            bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                            q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                            block_conn.move(g, vertex, vertex_id, move_id);
                            p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                        }
                    }

                    positive_move_occurred |= max_qap_gain > 0;
                }
            }
        }
    };
}

#endif //HEIPROMAP_MULTI_TRY_FM_REFINEMENT_H
