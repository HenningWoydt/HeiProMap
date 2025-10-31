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

namespace HeiProMap {
    class MultiTryFmRefinementConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit MultiTryFmRefinementConfiguration(const std::string &t_name) : ISerialRefinerConfiguration(t_name) {
        }

        u64 max_iteration = 1;
        f64 alpha = 10.0;
        f64 beta = 1.0;
    };

    class MultiTryFMRefinement final : public ISerialRefiner {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
        f64 m_imbalance = 0.0;
        weight_t m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;

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

        RandomEngine *random_engine = nullptr;
        const MultiTryFmRefinementConfiguration *config = nullptr;

    public:
        MultiTryFMRefinement() = default;

        ~MultiTryFMRefinement() override = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const f64 t_imbalance,
                        const weight_t t_lmax,
                        const std::vector<partition_t> &t_hierarchy,
                        const std::vector<weight_t> &t_distance,
                        RandomEngine &t_random_engine,
                        const ISerialRefinerConfiguration &i_config) override {
            m_n = t_n;
            m_m = t_m;
            m_k = t_k;
            m_imbalance = t_imbalance;
            m_lmax = t_lmax;
            m_hierarchy = t_hierarchy;
            m_distance = t_distance;

            random_engine = &t_random_engine;
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
        }

        void refine([[maybe_unused]] const u64 level,
                    [[maybe_unused]] const u64 max_level,
                    graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph) override {
            ScopedTimer _t("refinement", "MultiTryFMRefinement", "refine");

            f64 alpha = config->alpha;
            f64 beta = std::log(g.get_n());

            bool positive_move_occurred = true;
            for (u64 iteration = 0; iteration < config->max_iteration && positive_move_occurred; ++iteration) {
                positive_move_occurred = false;

                curr_boundary_size = 0;
                for (partition_t id = 0; id < bv_manager.get_k(); ++id) {
                    forall_bv_id_iu(bv_manager, id, i, u) {
                            curr_boundary[curr_boundary_size++] = u;
                        }
                    endfor
                }
                std::shuffle(curr_boundary.get_ptr(), curr_boundary.get_ptr() + curr_boundary_size, random_engine->generator);


                vertex_mark += 1;
                for (size_t ii = 0; ii < curr_boundary_size; ++ii) {
                    vertex_t u = curr_boundary[ii];
                    if (vertex_used[u] == vertex_mark) { continue; }
                    if (!bv_manager.is_boundary(u)) { continue; }

                    heap.clear();

                    s64 best_initial_qap = -std::numeric_limits<s64>::max();

                    // insert u into the priority queue
                    partition_t u_id = p_manager[u];
                    weight_t u_weight = g.weight(u);

                    // find all connected partitions to u
                    partition_t best_v_id = 0;
                    s64 best_qap_delta = -std::numeric_limits<s64>::max();

                    block_mark += 1;
                    forall_guiv(g, u, i, v) {
                            partition_t v_id = p_manager[v];
                            if (v_id == u_id) { continue; }
                            if (block_used[v_id] == block_mark) { continue; }
                            if (p_manager.get_bweight(v_id) + u_weight > m_lmax) { continue; }

                            s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                            block_used[v_id] = block_mark;

                            if (qap_delta > best_qap_delta) {
                                best_qap_delta = qap_delta;
                                best_v_id = v_id;
                            }
                        }
                    endfor
                    best_initial_qap = std::max(best_initial_qap, best_qap_delta);
                    if (best_qap_delta != -std::numeric_limits<s64>::max()) {
                        heap.push(u, best_v_id, best_qap_delta);
                    }

                    // insert all neighbors of u that are boundary into the queue
                    forall_guiv(g, u, i, neighbor) {
                            if (vertex_used[neighbor] == vertex_mark) { continue; }
                            if (!bv_manager.is_boundary(neighbor)) { continue; }

                            partition_t neighbor_id = p_manager[neighbor];
                            weight_t neighbor_weight = g.weight(neighbor);

                            best_qap_delta = -std::numeric_limits<s64>::max();
                            block_mark += 1;

                            forall_guiv(g, neighbor, j, v) {
                                    partition_t v_id = p_manager[v];
                                    if (v_id == neighbor_id) { continue; }
                                    if (block_used[v_id] == block_mark) { continue; }
                                    if (p_manager.get_bweight(v_id) + neighbor_weight > m_lmax) { continue; }

                                    s64 u_qap_delta = get_u_qap_delta(g, neighbor, neighbor_id, v_id, p_manager, d_oracle);
                                    block_used[v_id] = block_mark;

                                    if (u_qap_delta > best_qap_delta) {
                                        best_qap_delta = u_qap_delta;
                                        best_v_id = v_id;
                                    }
                                }
                            endfor

                            best_initial_qap = std::max(best_initial_qap, best_qap_delta);
                            if (best_qap_delta != -std::numeric_limits<s64>::max()) {
                                heap.push(neighbor, best_v_id, best_qap_delta);
                            }
                        }
                    endfor

                    if (heap.empty()) { continue; }

                    // process the queue
                    moves_size = 0;
                    size_t best_idx = 0;
                    s64 curr_qap_gain = 0;
                    s64 max_qap_gain = 0;

                    f64 steps_since_last_improvement = 0.0;
                    f64 qap_gain_mean = 0.0;
                    f64 qap_gain_var = 0.0;

                    while (!heap.empty()) {
                        vertex_t vertex = heap.top_u();
                        partition_t vertex_id = p_manager[vertex];
                        weight_t vertex_weight = g.weight(vertex);
                        partition_t move_id = heap.top_id();
                        s64 move_qap_delta = heap.top_qap_delta();
                        heap.pop();

                        if (p_manager.get_bweight(move_id) + vertex_weight > m_lmax) { continue; }

                        moves[moves_size++] = Move(vertex, vertex_id, move_id);
                        curr_qap_gain += move_qap_delta;
                        if (curr_qap_gain > max_qap_gain) {
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
                        f64 new_qap_gain_mean = qap_gain_mean + ((f64) move_qap_delta - qap_gain_mean) / steps_since_last_improvement;
                        f64 new_qap_gain_var = (qap_gain_var + ((f64) move_qap_delta - qap_gain_mean) * ((f64) move_qap_delta - new_qap_gain_mean)) / steps_since_last_improvement;

                        qap_gain_mean = new_qap_gain_mean;
                        qap_gain_var = new_qap_gain_var;

                        if (steps_since_last_improvement > 2.0 && steps_since_last_improvement * qap_gain_mean * qap_gain_mean > alpha * qap_gain_var + beta) { break; }

                        // we have to push or update the neighbors that were not moved already
                        forall_guiv(g, vertex, i, neighbor) {
                                if (vertex_used[neighbor] == vertex_mark) { continue; }
                                if (!is_boundary(g, p_manager, neighbor)) { continue; }

                                partition_t neighbor_id = p_manager[neighbor];
                                weight_t neighbor_weight = g.weight(neighbor);

                                best_qap_delta = -std::numeric_limits<s64>::max();

                                block_mark += 1;
                                forall_guiv(g, neighbor, j, v) {
                                        partition_t v_id = p_manager[v];
                                        if (v_id == neighbor_id) { continue; }
                                        if (block_used[v_id] == block_mark) { continue; }
                                        if (p_manager.get_bweight(v_id) + neighbor_weight > m_lmax) { continue; }

                                        s64 v_qap_delta = get_u_qap_delta(g, neighbor, neighbor_id, v_id, p_manager, d_oracle);
                                        block_used[v_id] = block_mark;

                                        if (v_qap_delta > best_qap_delta) {
                                            best_qap_delta = v_qap_delta;
                                            best_v_id = v_id;
                                        }
                                    }
                                endfor

                                if (best_qap_delta != -std::numeric_limits<s64>::max()) {
                                    heap.push_update(neighbor, best_v_id, best_qap_delta);
                                }
                            }
                        endfor
                    }

                    // revert all moves in partitioning manager
                    for (size_t i = 0; i < moves_size; i++) {
                        vertex_t vertex = moves[moves_size - 1 - i].u;
                        weight_t vertex_weight = g.weight(vertex);
                        partition_t vertex_id = moves[moves_size - 1 - i].to_move_id;
                        partition_t move_id = moves[moves_size - 1 - i].u_id;

                        p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                    }

                    // make all moves to best index
                    for (size_t i = 0; i < best_idx; ++i) {
                        vertex_t vertex = moves[i].u;
                        weight_t vertex_weight = g.weight(vertex);
                        partition_t vertex_id = moves[i].u_id;
                        partition_t move_id = moves[i].to_move_id;

                        bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                        q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                        p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                    }

                    positive_move_occurred |= max_qap_gain > 0;
                }
            }
        }

        void refine_layer([[maybe_unused]] const u64 level,
                          [[maybe_unused]] const u64 max_level,
                          [[maybe_unused]] graph_t &g,
                          [[maybe_unused]] d_oracle_t &d_oracle,
                          [[maybe_unused]] bv_manager_t &bv_manager,
                          [[maybe_unused]] p_manager_t &p_manager,
                          [[maybe_unused]] q_graph_t &q_graph,
                          [[maybe_unused]] size_t layer) override {
        }
    };
}

#endif //HEIPROMAP_MULTI_TRY_FM_REFINEMENT_H
