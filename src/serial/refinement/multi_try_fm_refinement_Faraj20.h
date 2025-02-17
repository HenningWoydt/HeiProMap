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

#ifndef HEIPROMAP_MULTI_TRY_FM_REFINEMENT_FARAJ20_H
#define HEIPROMAP_MULTI_TRY_FM_REFINEMENT_FARAJ20_H

#include <algorithm>
#include <random>

#include "../datastructures/distance_oracle.h"
#include "../datastructures/functions.h"
#include "../interfaces/ISerialActiveVertexManager.h"
#include "../interfaces/ISerialBoundaryVertexManager.h"
#include "../interfaces/ISerialQuotientGraph.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/qap.h"
#include "../utility/utils.h"

namespace HeiProMap {
    struct MultiTryFmRefinementFaraj20Configuration {
        u64 max_iteration = 2;
        f64 alpha         = 1000000.0;
        f64 beta          = 1.0;
    };

    /**
     * Executes Multi-Try FM Refinement as described in
     * > Marcelo Fonseca Faraj, Alexander van der Grinten, Henning Meyerhenke, Jesper Larsson Träff, and Christian Schulz.
     * > High-quality Hierarchical Process Mapping.
     * > In 18th International Symposium on Experimental Algorithms, SEA 2020, June 16-18, 2020, Catania, Italy, volume 160 of LIPIcs, pages 4:1–4:15.
     * > Schloss Dagstuhl - Leibniz-Zentrum für Informatik, 2020.
     */
    class MultiTryFMRefinementFaraj20 final : public ISerialRefiner {
    private:
        std::vector<partition_t> hierarchy;
        std::vector<weight_t> distance;
        partition_t k = 0;
        weight_t lmax = 0;

        std::vector<s32> used;
        s32 mark = -1;

        std::random_device rd;
        std::mt19937 gen;
        std::uniform_real_distribution<float> dis;

        KWayFMPriorityQueue queue;

    public:
        MultiTryFMRefinementFaraj20() : gen(rd()), dis(0.0f, 1.0f) {}

        void initialize(const vertex_t n,
                        std::vector<partition_t>& t_hierarchy,
                        std::vector<weight_t>& t_distance,
                        const weight_t t_lmax) override {
            hierarchy = t_hierarchy;
            distance  = t_distance;
            k         = prod<partition_t>(hierarchy);
            lmax      = t_lmax;

            used.resize(n, -1);

            queue = KWayFMPriorityQueue(n);
        }

        template <typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle, typename TSerialQuotientGraph>
        void refine([[maybe_unused]] MultiTryFmRefinementFaraj20Configuration& config,
                    [[maybe_unused]] TSerialGraph& g,
                    [[maybe_unused]] TSerialActiveVertexManager& av_manager,
                    [[maybe_unused]] TSerialBoundaryVertexManager& bv_manager,
                    [[maybe_unused]] TSerialPartitionManager& p_manager,
                    [[maybe_unused]] TSerialDistanceOracle& d_oracle,
                    [[maybe_unused]] TSerialQuotientGraph& q_graph) {
            static_assert(std::is_base_of_v<ISerialGraph, TSerialGraph>, "TSerialGraph must inherit from ISerialGraph");
            static_assert(std::is_base_of_v<ISerialActiveVertexManager, TSerialActiveVertexManager>, "TSerialActiveVertexManager must inherit from ISerialActiveVertexManager");
            static_assert(std::is_base_of_v<ISerialBoundaryVertexManager, TSerialBoundaryVertexManager>, "TSerialBoundaryVertexManager must inherit from ISerialBoundaryVertexManager");
            static_assert(std::is_base_of_v<ISerialPartitionManager, TSerialPartitionManager>, "TSerialPartitionManager must inherit from ISerialPartitionManager");
            static_assert(std::is_base_of_v<ISerialDistanceOracle, TSerialDistanceOracle>, "TSerialDistanceOracle must inherit from ISerialDistanceOracle");
            static_assert(std::is_base_of_v<ISerialQuotientGraph, TSerialQuotientGraph>, "TSerialQuotientGraph must inherit from ISerialQuotientGraph");
            u32 counter = 0;
            std::vector<u32> found_ids_mark(k, counter);

            std::vector<KWayFMMove> moves;

            config.beta = std::log(g.get_n());

            std::vector<vertex_t> curr_boundary;

            for (u64 iteration = 0; iteration < config.max_iteration; ++iteration) {
                auto sp = std::chrono::high_resolution_clock::now();

                u64 iteration_qap_gain = 0;
                u64 iteration_n_moves = 0;
                u64 iteration_n_queue_push = 0;

                mark += 1;

                for (vertex_t u : bv_manager) { curr_boundary.push_back(u); } // get the list
                std::shuffle(curr_boundary.begin(), curr_boundary.end(), gen);

                while (!curr_boundary.empty()) {
                    vertex_t u = curr_boundary.back();
                    curr_boundary.pop_back();

                    if (used[u] == mark) {
                        continue;
                    }

                    queue.clear();

                    // insert u into the priority queue
                    partition_t u_id  = p_manager[u];
                    weight_t u_weight = g.get_weight(u);

                    // find all connected partitions to u
                    counter += 1;
                    bool one_id_is_valid = false;
                    for (const auto& [v, w] : g[u]) {
                        // for (size_t idx = 0; idx < g.size(u); ++idx) {
                        // const vertex_t v = g.neighbor(u, idx);
                        // const weight_t w = g.get_weight(u, idx);
                        partition_t v_id = p_manager[v];
                        if (v_id != u_id && found_ids_mark[v_id] != counter && p_manager.get_bweight(v_id) + u_weight <= lmax) {
                            found_ids_mark[v_id] = counter;
                            s64 qap_delta        = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);

                            queue.push(u, u_id, v_id, qap_delta);
                            one_id_is_valid = true;
                        }
                    }
                    if (one_id_is_valid) {
                        queue.sort(u);
                    }

                    // insert all neighbor of u that are boundary into the queue
                    for (const auto& [neighbor, w] : g[u]) {
                        // for (size_t idx = 0; idx < g.size(vertex); ++idx) {
                        // const vertex_t neighbor = g.neighbor(vertex, idx);
                        // const weight_t w = g.get_weight(vertex, idx);
                        if (used[neighbor] == mark || !is_boundary(g, p_manager, neighbor)) {
                            continue;
                        }
                        partition_t neighbor_id = p_manager[neighbor];

                        counter += 1;
                        one_id_is_valid = false;
                        for (const auto& [v, w] : g[neighbor]) {
                            partition_t v_id = p_manager[v];
                            if (v_id != neighbor_id && found_ids_mark[v_id] != counter) {
                                found_ids_mark[v_id] = counter;
                                s64 qap_delta        = get_u_qap_delta(g, neighbor, neighbor_id, v_id, p_manager, d_oracle);

                                queue.push(neighbor, neighbor_id, v_id, qap_delta);
                                one_id_is_valid = true;
                            }
                        }
                        if (one_id_is_valid) {
                            queue.sort(neighbor);
                        }
                    }

                    // process the queue
                    moves.clear();
                    size_t best_idx   = 0;
                    s64 curr_qap_gain = 0;
                    s64 max_qap_gain  = 0;

                    u64 steps_since_last_improvement = 0;
                    f64 qap_gain_mean                = 0.0;
                    f64 qap_gain_var                 = 1.0;

                    while (!queue.empty()) {
                        const KWayFMMove move = queue.top();
                        queue.pop();

                        vertex_t vertex        = move.u;
                        weight_t vertex_weight = g.get_weight(vertex);
                        partition_t vertex_id  = p_manager[vertex];
                        partition_t move_id    = move.to_move_id;

                        if (vertex_id != move.u_id || used[vertex] == mark || p_manager.get_bweight(move_id) + vertex_weight > lmax) {
                            // if vertex_id and old vertex_id don't match
                            // vertex was already used,
                            // moving would overload
                            continue;
                        }

                        used[vertex] = mark;

                        // make move in structures
                        p_manager.move(vertex, vertex_weight, vertex_id, move_id);

                        moves.push_back(move);
                        curr_qap_gain += move.qap_delta;
                        if (curr_qap_gain > max_qap_gain && !p_manager.is_overloaded()) {
                            best_idx     = moves.size();
                            max_qap_gain = curr_qap_gain;

                            steps_since_last_improvement = 0;
                            qap_gain_mean                = 0.0;
                            qap_gain_var                 = 1.0;
                        }

                        steps_since_last_improvement += 1;
                        f64 new_qap_gain_mean = qap_gain_mean + ((f64)move.qap_delta - qap_gain_mean) / (f64)steps_since_last_improvement;
                        f64 new_qap_gain_var  = qap_gain_var + ((f64)move.qap_delta - qap_gain_mean) * ((f64)move.qap_delta - new_qap_gain_mean);

                        qap_gain_mean = new_qap_gain_mean;
                        qap_gain_var  = new_qap_gain_var;

                        if (steps_since_last_improvement > 3 && (f64)steps_since_last_improvement * qap_gain_mean * qap_gain_mean > config.alpha * qap_gain_var + config.beta) {
                            std::cout << "Stop on random walk: " << steps_since_last_improvement << " " << qap_gain_mean << " " << qap_gain_var << std::endl;
                            break;
                        }

                        // we have to push or update the neighbors that were not moved already
                        for (const auto& [neighbor, w] : g[vertex]) {
                            // for (size_t idx = 0; idx < g.size(vertex); ++idx) {
                            // const vertex_t neighbor = g.neighbor(vertex, idx);
                            // const weight_t w = g.get_weight(vertex, idx);
                            if (used[neighbor] == mark || !is_boundary(g, p_manager, vertex)) {
                                continue;
                            }
                            partition_t neighbor_id = p_manager[neighbor];

                            counter += 1;
                            one_id_is_valid = false;
                            for (const auto& [v, w] : g[neighbor]) {
                                partition_t v_id = p_manager[v];
                                if (v_id != neighbor_id && found_ids_mark[v_id] != counter) {
                                    found_ids_mark[v_id] = counter;
                                    s64 qap_delta        = get_u_qap_delta(g, neighbor, neighbor_id, v_id, p_manager, d_oracle);

                                    queue.push(neighbor, neighbor_id, v_id, qap_delta);
                                    one_id_is_valid = true;
                                }
                            }
                            if (one_id_is_valid) {
                                queue.sort(neighbor);
                            }
                        }
                    }

                    // revert all moves in partitioning manager
                    for (size_t i = 0; i < moves.size(); i++) {
                        vertex_t vertex        = moves[moves.size() - 1 - i].u;
                        weight_t vertex_weight = g.get_weight(vertex);
                        partition_t vertex_id  = moves[moves.size() - 1 - i].to_move_id;
                        partition_t move_id    = moves[moves.size() - 1 - i].u_id;

                        p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                    }

                    // make all moves to best index
                    for (size_t i = 0; i < best_idx; ++i) {
                        vertex_t vertex        = moves[i].u;
                        weight_t vertex_weight = g.get_weight(vertex);
                        partition_t vertex_id  = moves[i].u_id;
                        partition_t move_id    = moves[i].to_move_id;

                        bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                        q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                        p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                    }

                    iteration_qap_gain += max_qap_gain;
                    iteration_n_moves += moves.size();
                    iteration_n_queue_push += queue.push_operations;
                }

                auto ep = std::chrono::high_resolution_clock::now();
                std::cout << "iteration: " << iteration << " gain: " << iteration_qap_gain << " push ops: " << iteration_n_queue_push << " moves: " << iteration_n_moves << " time: " << get_seconds(sp, ep) << std::endl;
            }
        }
    };
}

#endif //HEIPROMAP_MULTI_TRY_FM_REFINEMENT_FARAJ20_H
