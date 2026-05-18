/*******************************************************************************
 * MIT License
 *
 * This file is part of Dyn-HeiProMap.
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

#ifndef DYN_HEIPROMAP_LABEL_PROPAGATION_REFINEMENT_H
#define DYN_HEIPROMAP_LABEL_PROPAGATION_REFINEMENT_H

#include <vector>
#include <random>
#include <algorithm>
#include <map>

#include "../datastructures/dyn_graph.h"
#include "../datastructures/quotient_graph.h"
#include "../utility/distance_oracle.h"
#include "../utility/configuration.h"
#include "../../utility/profiler.h"

namespace HeiProMap {
    class LabelPropagationRefinement {
    public:
        LabelPropagationRefinement() = default;

        weight_t refine(const DynGraph &g,
                        std::vector<partition_t> &partition,
                        const DistanceOracle &oracle,
                        const Configuration &config,
                        const std::vector<vertex_t> &initial_active_vertices,
                        const std::vector<vertex_t> &new_vertices,
                        u32 num_iterations,
                        QuotientGraph &q) {
            ScopedTimer _t("refinement", "label_propagation", "refine");

            if (partition.empty() && initial_active_vertices.empty() && new_vertices.empty()) return 0;

            u64 num_blocks = 1;
            for (auto h: config.hierarchy) num_blocks *= h;

            if (partition.size() < g.n) partition.resize(g.n, 0);

            std::vector<weight_t> block_weights(num_blocks, 0);
            for (vertex_t v = 0; v < g.n; ++v) {
                if (g.vertex_exists(v)) {
                    block_weights[partition[v]] += g.v_weights[v];
                }
            }

            weight_t total_weight = g.g_weight;
            weight_t max_block_weight = (weight_t) ((1.0 + config.imbalance) * ((f64) total_weight / num_blocks));

            weight_t total_migration_cost = 0;

            // 1. Label Propagation starting with initial_active_vertices
            std::vector<vertex_t> active_vertices = initial_active_vertices;
            // Ensure new vertices are also in the active set for refinement
            for (vertex_t v: new_vertices) {
                active_vertices.push_back(v);
            }
            std::sort(active_vertices.begin(), active_vertices.end());
            active_vertices.erase(std::unique(active_vertices.begin(), active_vertices.end()), active_vertices.end());

            std::vector<u8> is_active(g.n, 0);
            for (vertex_t v: active_vertices) {
                if (v < g.n) is_active[v] = 1;
            }

            std::mt19937 rng(config.seed);

            for (u32 iter = 0; iter < num_iterations; ++iter) {
                if (active_vertices.empty()) break;

                std::shuffle(active_vertices.begin(), active_vertices.end(), rng);

                std::vector<vertex_t> next_active;
                std::vector<u8> next_is_active(g.n, 0);

                u32 moves = 0;
                for (vertex_t v: active_vertices) {
                    is_active[v] = 0; // No longer active in current round
                    if (!g.vertex_exists(v)) continue;

                    partition_t current_block = partition[v];
                    weight_t v_weight = g.v_weights[v];

                    // Collect candidate blocks from neighbors
                    std::vector<partition_t> candidates;
                    candidates.push_back(current_block);
                    for (const auto &edge: g.neighbors[v]) {
                        candidates.push_back(edge.u < partition.size() ? partition[edge.u] : 0);
                    }
                    std::sort(candidates.begin(), candidates.end());
                    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());

                    partition_t best_block = current_block;
                    weight_t min_comm_cost = std::numeric_limits<weight_t>::max();

                    for (partition_t target_block: candidates) {
                        if (target_block != current_block && block_weights[target_block] + v_weight > max_block_weight) {
                            continue;
                        }

                        weight_t current_comm_cost = 0;
                        for (const auto &edge: g.neighbors[v]) {
                            partition_t nb_block = edge.u < partition.size() ? partition[edge.u] : 0;
                            current_comm_cost += edge.w * oracle.query(nb_block, target_block);
                        }

                        if (current_comm_cost < min_comm_cost) {
                            min_comm_cost = current_comm_cost;
                            best_block = target_block;
                        } else if (current_comm_cost == min_comm_cost && target_block == current_block) {
                            best_block = target_block;
                        }
                    }

                    if (best_block != current_block) {
                        total_migration_cost += v_weight * oracle.query(current_block, best_block);
                        q.move_vertex(v, current_block, best_block, g, partition);
                        block_weights[current_block] -= v_weight;
                        block_weights[best_block] += v_weight;
                        partition[v] = best_block;
                        moves++;

                        // Activate self and neighbors for next round
                        if (!next_is_active[v]) {
                            next_is_active[v] = 1;
                            next_active.push_back(v);
                        }
                        for (const auto &edge: g.neighbors[v]) {
                            if (!next_is_active[edge.u]) {
                                next_is_active[edge.u] = 1;
                                next_active.push_back(edge.u);
                            }
                        }
                    }
                }

                active_vertices = std::move(next_active);
                is_active = std::move(next_is_active);

                if (moves == 0) break;
            }
            return total_migration_cost;
        }
    };
}

#endif //DYN_HEIPROMAP_LABEL_PROPAGATION_REFINEMENT_H
