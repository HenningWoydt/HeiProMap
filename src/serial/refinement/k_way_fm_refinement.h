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
#include <random>

#include "k_way_fm_refinement_Faraj20.h"
#include "../datastructures/distance_oracle.h"
#include "../datastructures/functions.h"
#include "../interfaces/ISerialActiveVertexManager.h"
#include "../interfaces/ISerialBoundaryVertexManager.h"
#include "../interfaces/ISerialQuotientGraph.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/qap.h"
#include "../utility/utils.h"

namespace HeiProMap {
    struct KWayFMRefinementConfiguration {
        u64 max_iteration = 1; // how many iterations to run the algorithm at most
        f64 alpha         = 1000.0;
        f64 beta          = 1.0;
    };

    class KWayFMRefinement final : public ISerialRefiner {
    private:
        vertex_t                 m_n    = 0;
        vertex_t                 m_m    = 0;
        partition_t              m_k    = 0;
        weight_t                 m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t>    m_distance;
        u64                      m_seed = 0;

        std::vector<u32> vertex_used;
        u32              vertex_mark = 0;

        std::vector<u32> block_used;
        u32              block_marker = 0;

        // indexed max heaps
        // KWayFMPriorityQueue queue;

        std::priority_queue<KWayFMMove> prio_queue;

        std::mt19937                          gen;
        std::uniform_real_distribution<float> dis;

    public:
        KWayFMRefinement() = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const weight_t t_lmax,
                        const std::vector<partition_t> &t_hierarchy,
                        const std::vector<weight_t> &t_distance,
                        const u64 t_seed) override {
            m_n         = t_n;
            m_m         = t_m;
            m_k         = t_k;
            m_lmax      = t_lmax;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;
            m_seed      = t_seed;

            vertex_used.resize(t_n, 0);
            block_used.resize(t_n, 0);

            // queue = KWayFMPriorityQueue(t_n);

            gen.seed(m_seed);
            dis = std::uniform_real_distribution<float>(0.0f, 1.0f);
        }

        template<typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle, typename TSerialQuotientGraph>
        void refine([[maybe_unused]] KWayFMRefinementConfiguration &config,
                    [[maybe_unused]] TSerialGraph &g,
                    [[maybe_unused]] TSerialActiveVertexManager &av_manager,
                    [[maybe_unused]] TSerialBoundaryVertexManager &bv_manager,
                    [[maybe_unused]] TSerialPartitionManager &p_manager,
                    [[maybe_unused]] TSerialDistanceOracle &d_oracle,
                    [[maybe_unused]] TSerialQuotientGraph &q_graph) {
            static_assert(std::is_base_of_v<ISerialGraph, TSerialGraph>, "TSerialGraph must inherit from ISerialGraph");
            static_assert(std::is_base_of_v<ISerialActiveVertexManager, TSerialActiveVertexManager>, "TSerialActiveVertexManager must inherit from ISerialActiveVertexManager");
            static_assert(std::is_base_of_v<ISerialBoundaryVertexManager, TSerialBoundaryVertexManager>, "TSerialBoundaryVertexManager must inherit from ISerialBoundaryVertexManager");
            static_assert(std::is_base_of_v<ISerialPartitionManager, TSerialPartitionManager>, "TSerialPartitionManager must inherit from ISerialPartitionManager");
            static_assert(std::is_base_of_v<ISerialDistanceOracle, TSerialDistanceOracle>, "TSerialDistanceOracle must inherit from ISerialDistanceOracle");
            static_assert(std::is_base_of_v<ISerialQuotientGraph, TSerialQuotientGraph>, "TSerialQuotientGraph must inherit from ISerialQuotientGraph");

            std::vector<KWayFMMove> moves;

            config.beta = std::log(g.get_n());

            // std::cout << "alpha = " << config.alpha << std::endl;
            // std::cout << "beta  = " << config.beta << std::endl;

            for (u64 iteration = 0; iteration < config.max_iteration; ++iteration) {
                auto sp    = std::chrono::high_resolution_clock::now();

                // queue.clear();
                vertex_mark += 1;
                prio_queue = std::priority_queue<KWayFMMove>();

                // boost::heap::priority_queue<KWayFMMove> boost_prio_queue;

                // insert all boundary vertices
                for (vertex_t u: bv_manager) {
                    partition_t u_id     = p_manager[u];
                    weight_t    u_weight = g.get_weight(u);

                    // find all connected partitions to u
                    block_marker += 1;
                    bool one_id_is_valid = false;
                    for (const auto [v, w]: g[u]) {
                        partition_t v_id = p_manager[v];
                        if (v_id == u_id) { continue; }
                        if (block_used[v_id] == block_marker) { continue; }
                        if (p_manager.get_bweight(v_id) + u_weight > m_lmax) { continue; }

                        s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                        // queue.push(u, u_id, v_id, qap_delta);
                        prio_queue.emplace(u, u_id, v_id, qap_delta);
                        // boost_prio_queue.emplace(u, u_id, v_id, qap_delta);
                        one_id_is_valid = true;

                        block_used[v_id] = block_marker;
                    }
                    if (one_id_is_valid) {
                        // queue.sort(u);
                    }
                }

                moves.clear();
                size_t best_idx      = 0;
                s64    curr_qap_gain = 0;
                s64    max_qap_gain  = 0;

                u64 steps_since_last_improvement = 0;
                f64 qap_gain_mean                = 0.0;
                f64 qap_gain_var                 = 1.0;

                while (!prio_queue.empty()) {
                    const KWayFMMove move = prio_queue.top();
                    prio_queue.pop();

                    vertex_t vertex = move.u;
                    if (vertex_used[vertex] == vertex_mark) { continue; }

                    partition_t move_id = move.to_move_id;
                    if (!is_connected_to(g, p_manager, vertex, move_id)) { continue; }

                    partition_t vertex_id     = p_manager[vertex];
                    weight_t    vertex_weight = g.get_weight(vertex);
                    if (p_manager.get_bweight(move_id) + vertex_weight > m_lmax) { continue; }

                    s64 curr_qap_delta = get_u_qap_delta(g, vertex, vertex_id, move_id, p_manager, d_oracle);
                    if (curr_qap_delta != move.qap_delta) { continue; }

                    moves.push_back(move);
                    curr_qap_gain += move.qap_delta;
                    if (curr_qap_gain >= max_qap_gain) {
                        best_idx     = moves.size();
                        max_qap_gain = curr_qap_gain;

                        steps_since_last_improvement = 0;
                        qap_gain_mean                = 0.0;
                        qap_gain_var                 = 1.0;
                    }

                    // make move in structures
                    p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                    vertex_used[vertex] = vertex_mark;

                    steps_since_last_improvement += 1;
                    f64 new_qap_gain_mean = qap_gain_mean + ((f64) move.qap_delta - qap_gain_mean) / (f64) steps_since_last_improvement;
                    f64 new_qap_gain_var  = (qap_gain_var + ((f64) move.qap_delta - qap_gain_mean) * ((f64) move.qap_delta - new_qap_gain_mean)) / (f64) steps_since_last_improvement;

                    qap_gain_mean = new_qap_gain_mean;
                    qap_gain_var  = new_qap_gain_var;

                    if (steps_since_last_improvement > 3 && (f64) steps_since_last_improvement * qap_gain_mean * qap_gain_mean > config.alpha * qap_gain_var + config.beta) {
                        // std::cout << "Stop on random walk: " << steps_since_last_improvement << " " << qap_gain_mean << " " << qap_gain_var << std::endl;
                        break;
                    }


                    // we have to push or update the neighbors that were not moved already
                    for (const auto [neighbor, _]: g[vertex]) {
                        if (vertex_used[neighbor] == vertex_mark) { continue; }
                        if (!is_boundary(g, p_manager, neighbor)) { continue; }

                        partition_t neighbor_id     = p_manager[neighbor];
                        weight_t    neighbor_weight = g.get_weight(neighbor);

                        block_marker += 1;
                        bool one_id_is_valid = false;
                        for (const auto [v, w]: g[neighbor]) {
                            partition_t v_id = p_manager[v];
                            if (v_id == neighbor_id) { continue; }
                            if (block_used[v_id] == block_marker) { continue; }
                            if (p_manager.get_bweight(v_id) + neighbor_weight > m_lmax) { continue; }

                            s64 qap_delta = get_u_qap_delta(g, neighbor, neighbor_id, v_id, p_manager, d_oracle);
                            // queue.push(neighbor, neighbor_id, v_id, qap_delta);
                            prio_queue.emplace(neighbor, neighbor_id, v_id, qap_delta);
                            // boost_prio_queue.emplace(neighbor, neighbor_id, v_id, qap_delta);
                            one_id_is_valid = true;
                            block_used[v_id] = block_marker;
                        }
                        if (one_id_is_valid) {
                            // queue.sort(neighbor);
                        }
                    }
                }

                // revert all moves in partitioning manager
                for (size_t i = 0; i < moves.size(); i++) {
                    vertex_t    vertex        = moves[moves.size() - 1 - i].u;
                    weight_t    vertex_weight = g.get_weight(vertex);
                    partition_t vertex_id     = moves[moves.size() - 1 - i].to_move_id;
                    partition_t move_id       = moves[moves.size() - 1 - i].u_id;

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

                auto ep = std::chrono::high_resolution_clock::now();
                // std::cout << "iteration: " << iteration << " best_idx: " << best_idx << " new_gain: " << max_qap_gain << " max moves: " << moves.size() << " time: " << get_seconds(sp, ep) << std::endl;
            }
        }
    };
}

#endif //HEIPROMAP_K_WAY_FM_REFINEMENT_H
