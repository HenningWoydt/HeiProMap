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

#include "../datastructures/indexed_max_heap.h"
#include "../interfaces/ISerialQuotientGraph.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/qap.h"
#include "../utility/utils.h"

namespace HeiProMap {
    class QuotientGraphRefinement final : public ISerialRefiner {
    private:
        std::vector<partition_t> hierarchy;
        std::vector<weight_t> distance;
        partition_t k = 0;
        weight_t lmax = 0;

        // indexed max heaps
        IndexedMaxHeap<s64> boundary_vertices;
        std::vector<vertex_t> moves;
        std::vector<s64> curr_qap_gain;
        s64 max_qap_gain = 0;

        std::vector<s32> used;
        s32 mark = -1;

    public:
        QuotientGraphRefinement() = default;

        void initialize(const vertex_t n,
                        std::vector<partition_t>& t_hierarchy,
                        std::vector<weight_t>& t_distance,
                        const weight_t t_lmax) override {
            hierarchy = t_hierarchy;
            distance  = t_distance;
            k         = prod<partition_t>(hierarchy);
            lmax      = t_lmax;

            used.resize(n, -1);

            // indexed max heaps
            boundary_vertices = IndexedMaxHeap<s64>(n);
        }

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
    };
}

#endif //HEIPROMAP_QUOTIENT_GRAPH_REFINEMENT_H
