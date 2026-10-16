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
#include "../distance_oracles/matrix_distance_oracle.h"
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
                                 PartitionManager &p_manager,
                                 AlignedArray<vertex_t> &vertices) {
        const vertex_t n = g.n;
        const partition_t k = d_oracle.get_k();
        const weight_t total_weight = g.g_weight;
        const weight_t max_block_weight = (weight_t) ((1.0 + imbalance) * (total_weight / k));

        p_manager.reset_weights();
        for (vertex_t u = 0; u < n; ++u) {
            p_manager.partition[u] = NO_ID;
        }

        for (vertex_t u = 0; u < n; ++u) {
            vertices[u] = u;
        }

        if (seed != 0) {
            RandomEngine re(seed);
            std::shuffle(vertices.get_ptr(), vertices.get_ptr() + n, re.generator);
        }

        for (vertex_t i = 0; i < n; ++i) {
            const vertex_t u = vertices[i];
            const weight_t u_weight = g.v_weight(u);

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
                        current_cost += g.edge_w(i) * d_oracle.get(block, v_block);
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

    template<typename DistanceOracleT>
    inline void compute_best_move(const graph_t &g,
                                  const DistanceOracleT &d_oracle,
                                  const PartitionManager &p_manager,
                                  const weight_t max_block_weight,
                                  const vertex_t u,
                                  weight_t &gain,
                                  partition_t &target) {
        const weight_t u_weight = g.v_weight(u);
        const partition_t old_block = p_manager[u];

        weight_t old_cost = 0;
        for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
            old_cost += g.edge_w(i) * d_oracle.get(old_block, p_manager[g.edges_v[i]]);
        }

        gain = 0;
        target = old_block;

        for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
            const partition_t cand = p_manager[g.edges_v[i]];
            if (cand == old_block) continue;
            if (p_manager.bweights[cand] + u_weight > max_block_weight) continue;

            weight_t cand_cost = 0;
            for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                cand_cost += g.edge_w(j) * d_oracle.get(cand, p_manager[g.edges_v[j]]);
            }

            weight_t cand_gain = old_cost - cand_cost;
            if (cand_gain > gain || (cand_gain == gain && cand_gain > 0 && p_manager.bweights[cand] < p_manager.bweights[target])) {
                gain = cand_gain;
                target = cand;
            }
        }
    }

    template<typename DistanceOracleT>
    inline void greedy_refine(const graph_t &g,
                              const DistanceOracleT &d_oracle,
                              const f64 imbalance,
                              PartitionManager &p_manager) {
        const vertex_t n = g.n;
        const weight_t max_block_weight = (weight_t) ((1.0 + imbalance) * (g.g_weight / d_oracle.get_k()));

        struct PQEntry {
            weight_t gain;
            vertex_t vertex;
            bool operator<(const PQEntry &o) const { return gain < o.gain; }
        };

        std::priority_queue<PQEntry> pq;
        std::vector<weight_t> cached_gain(n, 0);
        std::vector<partition_t> cached_target(n);
        std::vector<bool> in_queue(n, false);

        for (vertex_t u = 0; u < n; ++u) {
            weight_t gain; partition_t target;
            compute_best_move(g, d_oracle, p_manager, max_block_weight, u, gain, target);
            cached_gain[u] = gain;
            cached_target[u] = target;
            if (gain > 0) {
                pq.push({gain, u});
                in_queue[u] = true;
            }
        }

        while (!pq.empty()) {
            auto [gain, u] = pq.top();
            pq.pop();
            in_queue[u] = false;

            if (gain != cached_gain[u]) continue;

            weight_t fresh_gain; partition_t fresh_target;
            compute_best_move(g, d_oracle, p_manager, max_block_weight, u, fresh_gain, fresh_target);
            cached_gain[u] = fresh_gain;
            cached_target[u] = fresh_target;

            if (fresh_gain != gain) {
                if (fresh_gain > 0) {
                    pq.push({fresh_gain, u});
                    in_queue[u] = true;
                }
                continue;
            }

            if (fresh_gain <= 0) continue;

            const partition_t old_block = p_manager[u];
            const weight_t u_weight = g.v_weight(u);
            if (p_manager.bweights[fresh_target] + u_weight > max_block_weight) continue;
            p_manager.move(u, u_weight, old_block, fresh_target);

            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t v = g.edges_v[i];
                weight_t v_gain; partition_t v_target;
                compute_best_move(g, d_oracle, p_manager, max_block_weight, v, v_gain, v_target);
                cached_gain[v] = v_gain;
                cached_target[v] = v_target;
                if (v_gain > 0 && !in_queue[v]) {
                    pq.push({v_gain, v});
                    in_queue[v] = true;
                }
            }
        }
    }
}

#endif //HEIPROMAP_GREEDY_PARTITIONER_H
