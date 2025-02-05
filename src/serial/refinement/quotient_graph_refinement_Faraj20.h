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

#ifndef HEIPROMAP_QUOTIENT_GRAPH_REFINEMENT_FARAJ20_H
#define HEIPROMAP_QUOTIENT_GRAPH_REFINEMENT_FARAJ20_H

#include <algorithm>
#include <queue>
#include <random>

#include "../datastructures/distance_oracle.h"
#include "../datastructures/indexed_max_heap.h"
#include "../interfaces/ISerialQuotientGraph.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/qap.h"
#include "../utility/utils.h"

namespace HeiProMap {
    /**
     * Executes quotient graph refinement as described in
     * > Marcelo Fonseca Faraj, Alexander van der Grinten, Henning Meyerhenke, Jesper Larsson Träff, and Christian Schulz.
     * > High-quality Hierarchical Process Mapping.
     * > In 18th International Symposium on Experimental Algorithms, SEA 2020, June 16-18, 2020, Catania,Italy, volume 160 of LIPIcs, pages 4:1–4:15.
     * > Schloss Dagstuhl - Leibniz-Zentrum für Informatik, 2020.
     */
    class QuotientGraphRefinementFaraj20 final : public ISerialRefiner {
    private:
        std::vector<partition_t> hierarchy;
        std::vector<weight_t>    distance;
        partition_t              k    = 0;
        weight_t                 lmax = 0;

        // indexed max heaps
        IndexedMaxHeap<s64>   boundary_vertices_u;
        IndexedMaxHeap<s64>   boundary_vertices_v;
        std::vector<vertex_t> moves;
        std::vector<s64>      curr_qap_gain;
        s64                   max_qap_gain = 0;

        // active block scheduling
        std::vector<bool> active_this_round;
        std::vector<bool> active_next_round;

        std::vector<s32> used;
        s32              mark = -1;

        std::random_device                    rd;
        std::mt19937                          gen;
        std::uniform_real_distribution<float> dis;

    public:
        QuotientGraphRefinementFaraj20() : gen(rd()), dis(0.0f, 1.0f) {}

        void initialize(const vertex_t n,
                        std::vector<partition_t> &t_hierarchy,
                        std::vector<weight_t> &t_distance,
                        const weight_t t_lmax) override {
            hierarchy = t_hierarchy;
            distance  = t_distance;
            k         = prod<partition_t>(hierarchy);
            lmax      = t_lmax;

            used.resize(n, -1);

            // indexed max heaps
            boundary_vertices_u = IndexedMaxHeap<s64>(n);
            boundary_vertices_v = IndexedMaxHeap<s64>(n);

            // active block scheduling
            active_this_round.resize(k);
            active_next_round.resize(k);
        }

        /*
        template <typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle, typename TSerialQuotientGraph>
        void refine([[maybe_unused]] TSerialGraph& g,
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

            bool global_moved         = true;
            u64 global_max_iterations = 2;
            for (u64 global_i = 0; global_i < global_max_iterations && global_moved; ++global_i) {
                global_moved = false;

                for (partition_t u_id = 0; u_id < k; ++u_id) {
                    for (partition_t v_id = u_id + 1; v_id < k; ++v_id) {
                        if (!q_graph.has_edge(u_id, v_id)) {
                            // no boundary between u_id and v_id
                            continue;
                        }

                        bool local_moved         = true;
                        u64 local_max_iterations = 3;
                        for (u64 local_i = 0; local_i < local_max_iterations && local_moved; ++local_i) {
                            local_moved = false;

                            // add all boundary vertices with qap
                            boundary_vertices.clear();
                            for (vertex_t u : bv_manager[u_id]) {
                                for (const auto& [v, w] : g[u]) {
                                    if (p_manager[v] == v_id) {
                                        // u is connected to block v_id
                                        if (!boundary_vertices.entry_exists(u)) {
                                            s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                                            boundary_vertices.push(u, qap_delta);
                                        }

                                        // v is connected to block u_id
                                        if (!boundary_vertices.entry_exists(v)) {
                                            s64 qap_delta = get_u_qap_delta(g, v, v_id, u_id, p_manager, d_oracle);
                                            boundary_vertices.push(v, qap_delta);
                                        }
                                    }
                                }
                            }

                            // make moves
                            mark += 1;
                            moves.clear();
                            curr_qap_gain.clear();
                            curr_qap_gain.push_back(0);
                            max_qap_gain = 0;
                            while (!boundary_vertices.empty()) {
                                vertex_t vertex        = boundary_vertices.top_key();
                                partition_t vertex_id  = p_manager[vertex];
                                partition_t move_id    = u_id == vertex_id ? v_id : u_id;
                                weight_t vertex_weight = g.get_weight(vertex);
                                s64 qap_delta          = boundary_vertices.top();
                                bool overloads         = p_manager.get_bweight(move_id) + vertex_weight > lmax;

                                boundary_vertices.pop();
                                used[vertex] = mark;

                                if (overloads || !bv_manager.is_boundary(vertex)) {
                                    // if the move overloads the block, then do not move
                                    // if the vertex is not boundary anymore, then do not move
                                    continue;
                                }

                                // move the vertex
                                moves.push_back(vertex);
                                curr_qap_gain.push_back(curr_qap_gain.back() + qap_delta);
                                max_qap_gain = std::max(max_qap_gain, curr_qap_gain.back());

                                // make move in structures
                                bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                                // q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                                p_manager.move(vertex, vertex_weight, vertex_id, move_id);

                                // we have to push or update the neighbors that were not moved already
                                for (const auto& [neighbor, w] : g[vertex]) {
                                    if (used[neighbor] == mark ||
                                        (p_manager[neighbor] != vertex_id && p_manager[neighbor] != move_id) ||
                                        !bv_manager.is_boundary(neighbor)) {
                                        continue;
                                    }

                                    partition_t old_id = p_manager[neighbor];
                                    partition_t new_id = old_id == vertex_id ? move_id : vertex_id;

                                    s64 new_qap_delta = get_u_qap_delta(g, neighbor, old_id, new_id, p_manager, d_oracle);
                                    boundary_vertices.push_update(neighbor, new_qap_delta);
                                }
                            }

                            // revert to state with best qap delta
                            while (!curr_qap_gain.empty() && curr_qap_gain.back() != max_qap_gain) {
                                vertex_t vertex        = moves.back();
                                weight_t vertex_weight = g.get_weight(vertex);
                                partition_t vertex_id  = p_manager[vertex];
                                partition_t move_id    = u_id == vertex_id ? v_id : u_id;

                                bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                                // q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                                p_manager.move(vertex, vertex_weight, vertex_id, move_id);

                                moves.pop_back();
                                curr_qap_gain.pop_back();
                            }

                            // keep quotient graph up to date
                            for (vertex_t vertex : moves) {
                                partition_t vertex_id = p_manager[vertex];
                                partition_t move_id   = u_id == vertex_id ? v_id : u_id;
                                q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                            }

                            if (!moves.empty()) {
                                local_moved = true;
                            }
                        }
                        global_moved |= local_moved;
                    }
                }
            }
        }
         */

        template<typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle, typename TSerialQuotientGraph>
        void refine([[maybe_unused]] TSerialGraph &g,
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

            std::fill(active_this_round.begin(), active_this_round.end(), true);
            std::fill(active_next_round.begin(), active_next_round.end(), false);
            u64 max_iterations = 1000;
            u64 iteration      = 0;

            while (std::any_of(active_this_round.begin(), active_this_round.end(), [](bool value) { return value; }) && iteration < max_iterations) {
                iteration++;

                struct Pair {
                    partition_t u_id;
                    partition_t v_id;

                    Pair(partition_t t_u_id, partition_t t_v_id) {
                        u_id = t_u_id;
                        v_id = t_v_id;
                    }
                };
                std::vector<Pair> pairs;

                for (partition_t u_id = 0; u_id < k; ++u_id) {
                    for (partition_t v_id = u_id + 1; v_id < k; ++v_id) {
                        if (!q_graph.has_edge(u_id, v_id) || (!active_this_round[u_id] && !active_this_round[v_id])) {
                            // no boundary between u_id and v_id
                            active_next_round[u_id] = false;
                            active_next_round[v_id] = false;
                        } else {
                            pairs.push_back(Pair(u_id, v_id));
                        }
                    }
                }

                auto rng = std::default_random_engine{};
                std::shuffle(std::begin(pairs), std::end(pairs), rng);
                for (auto p: pairs) {
                    vertex_t u_id = p.u_id;
                    vertex_t v_id = p.v_id;

                    // add all boundary vertices with gain
                    boundary_vertices_u.clear();
                    boundary_vertices_v.clear();
                    for (vertex_t u: bv_manager[u_id]) {
                        for (const auto &[v, w]: g[u]) {
                            if (p_manager[v] == v_id) {
                                // u is connected to block v_id
                                if (!boundary_vertices_u.entry_exists(u)) {
                                    s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                                    boundary_vertices_u.push(u, qap_delta);
                                }

                                // v is connected to block u_id
                                if (!boundary_vertices_v.entry_exists(v)) {
                                    s64 qap_delta = get_u_qap_delta(g, v, v_id, u_id, p_manager, d_oracle);
                                    boundary_vertices_v.push(v, qap_delta);
                                }
                            }
                        }
                    }

                    // start executing moves based on the TopGain method
                    mark += 1;
                    moves.clear();
                    curr_qap_gain.clear();
                    curr_qap_gain.push_back(0);
                    max_qap_gain = 0;
                    while (!boundary_vertices_u.empty() || !boundary_vertices_v.empty()) {

                        // determine from which block to choose
                        bool choose_u = true;
                        // 1. if one block is empty, then choose the other one
                        if (boundary_vertices_u.empty() || boundary_vertices_v.empty()) {
                            if (boundary_vertices_u.empty()) {
                                choose_u = false;
                            }
                        } else {
                            // 2. choose the block with greater gain and randomly if even
                            if (boundary_vertices_v.top() > boundary_vertices_u.top()) {
                                choose_u = false;
                            } else if (boundary_vertices_v.top() == boundary_vertices_u.top()) {
                                choose_u = dis(gen) < 0.5;
                            }

                            // 3. if one block is overloaded, choose the larger one, if both same sizes, then randomly
                            if (p_manager.get_bweight(u_id) > lmax && p_manager.get_bweight(u_id) > p_manager.get_bweight(v_id)) { choose_u = true; }
                            if (p_manager.get_bweight(v_id) > lmax && p_manager.get_bweight(v_id) > p_manager.get_bweight(u_id)) { choose_u = false; }
                            if (p_manager.get_bweight(v_id) > lmax && p_manager.get_bweight(u_id) && p_manager.get_bweight(v_id) == p_manager.get_bweight(u_id)) { choose_u = dis(gen) < 0.5; }
                        }


                        IndexedMaxHeap<s64> &boundary_vertices = choose_u ? boundary_vertices_u : boundary_vertices_v;
                        vertex_t            vertex             = boundary_vertices.top_key();
                        partition_t         vertex_id          = p_manager[vertex];
                        partition_t         move_id            = u_id == vertex_id ? v_id : u_id;
                        weight_t            vertex_weight      = g.get_weight(vertex);
                        s64                 qap_delta          = boundary_vertices.top();
                        bool                overloads          = p_manager.get_bweight(move_id) + vertex_weight > lmax;

                        boundary_vertices.pop();

                        if (overloads || !bv_manager.is_boundary(vertex) || used[vertex] == mark) {
                            // if the move overloads the block, then do not move
                            // if the vertex is not boundary anymore, then do not move
                            // if the vertex was already used, then do not move
                            continue;
                        }

                        // move the vertex
                        moves.push_back(vertex);
                        curr_qap_gain.push_back(curr_qap_gain.back() + qap_delta);
                        max_qap_gain = std::max(max_qap_gain, curr_qap_gain.back());
                        used[vertex] = mark;

                        // make move in structures
                        bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                        // q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                        p_manager.move(vertex, vertex_weight, vertex_id, move_id);

                        // we have to push or update the neighbors that were not moved already
                        for (const auto &[neighbor, w]: g[vertex]) {
                            if (used[neighbor] == mark ||
                                (p_manager[neighbor] != vertex_id && p_manager[neighbor] != move_id) ||
                                !bv_manager.is_boundary(neighbor)) {
                                continue;
                            }

                            partition_t old_id = p_manager[neighbor];
                            partition_t new_id = old_id == vertex_id ? move_id : vertex_id;

                            s64 new_qap_delta = get_u_qap_delta(g, neighbor, old_id, new_id, p_manager, d_oracle);
                            if (old_id == u_id) {
                                boundary_vertices_u.push_update(neighbor, new_qap_delta);
                            } else {
                                boundary_vertices_v.push_update(neighbor, new_qap_delta);
                            }
                        }
                    }

                    // revert to state with best qap delta
                    while (!curr_qap_gain.empty() && curr_qap_gain.back() != max_qap_gain) {
                        vertex_t    vertex        = moves.back();
                        weight_t    vertex_weight = g.get_weight(vertex);
                        partition_t vertex_id     = p_manager[vertex];
                        partition_t move_id       = u_id == vertex_id ? v_id : u_id;

                        bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                        // q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                        p_manager.move(vertex, vertex_weight, vertex_id, move_id);

                        moves.pop_back();
                        curr_qap_gain.pop_back();
                    }

                    // keep quotient graph up to date
                    for (vertex_t vertex: moves) {
                        partition_t vertex_id = p_manager[vertex];
                        partition_t move_id   = u_id == vertex_id ? v_id : u_id;
                        q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                    }

                    if (max_qap_gain > 0) {
                        active_next_round[u_id] = true;
                        active_next_round[v_id] = true;
                    }
                }
                active_this_round.swap(active_next_round);
                std::fill(active_next_round.begin(), active_next_round.end(), false);
            }
        }

    };
}

#endif //HEIPROMAP_QUOTIENT_GRAPH_REFINEMENT_FARAJ20_H
