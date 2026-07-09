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

#ifndef HEIPROMAP_GREEDY_PARTITIONER_H
#define HEIPROMAP_GREEDY_PARTITIONER_H

#include <algorithm>
#include <limits>
#include <numeric>
#include <vector>

#include "../definitions.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../datastructures/distance_oracle.h"
#include "../utility/profiler.h"
#include "../utility/random_engine.h"

namespace HeiProMap {
    class UniformDistanceOracle {
    public:
        UniformDistanceOracle(const partition_t k) : m_k(k) {
        }

        partition_t get_k() const { return m_k; }

        weight_t get(const partition_t u_id, const partition_t v_id) const {
            return (u_id == v_id) ? 0 : 1;
        }

    private:
        partition_t m_k;
    };

    /**
     * Greedy algorithm that assigns each vertex of the graph one after the other
     * always to the best fitting block.
     * Best fitting is defined as the block that minimizes the QAP cost
     * (communication cost with already assigned neighbors) while respecting
     * the maximum block weight constraint.
     */
    template<typename DistanceOracleT>
    inline void greedy_partition(const graph_t &g,
                                 const DistanceOracleT &d_oracle,
                                 const f64 imbalance,
                                 const u64 seed,
                                 PartitionManager &p_manager) {
        HEIPROMAP_PROFILE_SCOPE("partitioning", "greedy_partitioner", "greedy_partition");

        const vertex_t n = g.n;
        const partition_t k = d_oracle.get_k();
        const weight_t total_weight = g.g_weight;
        const weight_t max_block_weight = (weight_t) ((1.0 + imbalance) * (total_weight / k));

        p_manager.reset_weights();
        for (vertex_t u = 0; u < n; ++u) {
            p_manager.partition[u] = NO_ID;
        }

        std::vector<vertex_t> vertices(n);
        std::iota(vertices.begin(), vertices.end(), 0);

        if (seed != 0) {
            RandomEngine re(seed);
            std::shuffle(vertices.begin(), vertices.end(), re.generator);
        }

        for (const vertex_t u: vertices) {
            const weight_t u_weight = g.v_weights[u];

            partition_t best_block = NO_ID;
            weight_t min_cost = std::numeric_limits<weight_t>::max();

            for (partition_t block = 0; block < k; ++block) {
                if (p_manager.bweights[block] + u_weight > max_block_weight) {
                    continue;
                }

                weight_t current_cost = 0;
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    const vertex_t v = g.edges_v[i];
                    const partition_t v_block = p_manager.partition[v];
                    if (v_block != NO_ID) {
                        current_cost += g.edges_w[i] * d_oracle.get(block, v_block);
                    }
                }

                if (current_cost < min_cost) {
                    min_cost = current_cost;
                    best_block = block;
                } else if (current_cost == min_cost) {
                    // Tie-breaking: prefer blocks with less weight
                    if (best_block == NO_ID || p_manager.bweights[block] < p_manager.bweights[best_block]) {
                        best_block = block;
                    }
                }
            }

            // If no block fits (should not happen if max_block_weight is reasonable and u_weight <= max_block_weight),
            // pick the one with the least weight.
            if (best_block == NO_ID) {
                best_block = 0;
                weight_t min_weight = p_manager.bweights[0];
                for (partition_t block = 1; block < k; ++block) {
                    if (p_manager.bweights[block] < min_weight) {
                        min_weight = p_manager.bweights[block];
                        best_block = block;
                    }
                }
            }

            p_manager.set(u, u_weight, best_block);
        }
    }
}

#endif // HEIPROMAP_GREEDY_PARTITIONER_H
