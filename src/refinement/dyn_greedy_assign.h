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

#ifndef HEIPROMAP_DYN_GREEDY_ASSIGN_H
#define HEIPROMAP_DYN_GREEDY_ASSIGN_H

#include <vector>
#include <limits>
#include <algorithm>

#include "../definitions.h"
#include "../datastructures/dyn_graph.h"
#include "../datastructures/partition_manager.h"
#include "../distance_oracles/matrix_distance_oracle.h"

namespace HeiProMap {

    class DynGreedyAssign {
    public:
        /**
         * Greedily assigns the vertices in new_vertices to blocks that minimize communication
         * cost to their already placed neighbors while respecting max_block_weight.
         * Vertices are placed in BFS order starting from vertices that have a placed neighbor, so a
         * new vertex never sees an unplaced neighbor as a member of some block. Components without
         * any placed vertex start in the lightest block.
         * New vertices are not counted in p_manager before they are placed here.
         * Every placed vertex pays new_vertex_unit_cost per unit of weight.
         */
        static weight_t assign(const DynGraph &g,
                               PartitionManager &p_manager,
                               const DistanceOracle &oracle,
                               const std::vector<std::pair<vertex_t, weight_t> > &new_vertices,
                               weight_t max_block_weight = 0,
                               weight_t new_vertex_unit_cost = 0) {
            weight_t placement_cost = 0;

            std::vector<char> unplaced(g.n, 0);
            std::vector<vertex_t> todo;
            for (const auto &entry: new_vertices) {
                vertex_t v = entry.first;
                if (!g.vertex_exists(v) || unplaced[v]) continue;
                unplaced[v] = 1;
                todo.push_back(v);
            }
            if (todo.empty()) return 0;

            auto has_placed_neighbor = [&](vertex_t v) {
                for (const auto &edge: g.neighbors[v]) {
                    if (!unplaced[edge.u]) return true;
                }
                return false;
            };

            auto place = [&](vertex_t v) {
                weight_t v_weight = g.v_weights[v];

                std::vector<partition_t> candidates;
                for (const auto &edge: g.neighbors[v]) {
                    if (!unplaced[edge.u]) candidates.push_back(p_manager[edge.u]);
                }
                std::sort(candidates.begin(), candidates.end());
                candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
                if (candidates.empty()) {
                    for (partition_t b = 0; b < p_manager.k; ++b) candidates.push_back(b);
                }

                partition_t best_block = candidates[0];
                weight_t min_comm_cost = std::numeric_limits<weight_t>::max();
                partition_t min_weight_block = candidates[0];
                weight_t min_target_weight = std::numeric_limits<weight_t>::max();

                for (partition_t target_block: candidates) {
                    weight_t target_weight = p_manager.bweights[target_block];
                    if (target_weight < min_target_weight) {
                        min_target_weight = target_weight;
                        min_weight_block = target_block;
                    }
                    if (target_weight + v_weight > p_manager.get_lmax(target_block)) continue;

                    weight_t comm_cost = 0;
                    for (const auto &edge: g.neighbors[v]) {
                        if (unplaced[edge.u]) continue;
                        comm_cost += edge.w * oracle.get(p_manager[edge.u], target_block);
                    }

                    // Ties go to the lighter block.
                    if (comm_cost < min_comm_cost ||
                        (comm_cost == min_comm_cost && target_weight < p_manager.bweights[best_block])) {
                        min_comm_cost = comm_cost;
                        best_block = target_block;
                    }
                }

                if (min_comm_cost == std::numeric_limits<weight_t>::max()) {
                    // All candidates exceeded capacity, fall back to the lightest block
                    best_block = min_weight_block;
                }

                p_manager.set(v, v_weight, best_block);
                unplaced[v] = 0;
                placement_cost += v_weight * new_vertex_unit_cost;
            };

            std::vector<vertex_t> queue;
            std::vector<char> queued(g.n, 0);
            size_t next_seed = 0;
            size_t remaining = todo.size();

            auto enqueue = [&](vertex_t v) {
                if (unplaced[v] && !queued[v]) {
                    queued[v] = 1;
                    queue.push_back(v);
                }
            };

            for (vertex_t v: todo) {
                if (has_placed_neighbor(v)) enqueue(v);
            }

            size_t head = 0;
            while (remaining > 0) {
                if (head == queue.size()) {
                    // No vertex adjacent to a placed one is left, start a new component.
                    while (next_seed < todo.size() && !unplaced[todo[next_seed]]) next_seed++;
                    enqueue(todo[next_seed]);
                }
                vertex_t v = queue[head++];
                place(v);
                remaining--;
                for (const auto &edge: g.neighbors[v]) enqueue(edge.u);
            }

            return placement_cost;
        }
    };

}

#endif // HEIPROMAP_DYN_GREEDY_ASSIGN_H
