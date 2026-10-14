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
         * Greedily assigns newly added/unassigned vertices in g.new_vertices
         * to blocks that minimize communication cost while respecting max_block_weight.
         */
        static weight_t assign(const DynGraph &g,
                               PartitionManager &p_manager,
                               const DistanceOracle &oracle,
                               weight_t max_block_weight = 0) {
            weight_t migration_cost = 0;
            for (vertex_t v: g.new_vertices) {
                if (!g.vertex_exists(v)) continue;

                partition_t current_block = p_manager[v];
                weight_t v_weight = g.v_weights[v];

                std::vector<partition_t> candidates;
                candidates.push_back(0);
                for (const auto &edge: g.neighbors[v]) {
                    candidates.push_back(p_manager[edge.u]);
                }
                std::sort(candidates.begin(), candidates.end());
                candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());

                partition_t best_block = current_block;
                weight_t min_comm_cost = std::numeric_limits<weight_t>::max();
                partition_t min_weight_block = current_block;
                weight_t min_target_weight = std::numeric_limits<weight_t>::max();

                for (partition_t target_block: candidates) {
                    weight_t current_target_weight = p_manager.bweights[target_block];
                    if (current_target_weight < min_target_weight) {
                        min_target_weight = current_target_weight;
                        min_weight_block = target_block;
                    }
                    if (target_block != current_block && current_target_weight + v_weight > p_manager.get_lmax(target_block)) {
                        continue;
                    }

                    weight_t current_comm_cost = 0;
                    for (const auto &edge: g.neighbors[v]) {
                        partition_t nb_block = p_manager[edge.u];
                        current_comm_cost += edge.w * oracle.get(nb_block, target_block);
                    }

                    if (current_comm_cost < min_comm_cost) {
                        min_comm_cost = current_comm_cost;
                        best_block = target_block;
                    }
                }

                if (min_comm_cost == std::numeric_limits<weight_t>::max()) {
                    // All candidates exceeded capacity, fall back to candidate with minimum load
                    best_block = min_weight_block;
                }

                if (best_block != current_block) {
                    migration_cost += v_weight * oracle.get(current_block, best_block);
                    p_manager.move_serial(v, v_weight, current_block, best_block);
                }
            }
            return migration_cost;
        }
    };

}

#endif // HEIPROMAP_DYN_GREEDY_ASSIGN_H
