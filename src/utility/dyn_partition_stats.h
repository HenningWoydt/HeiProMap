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

#ifndef HEIPROMAP_DYN_PARTITION_STATS_H
#define HEIPROMAP_DYN_PARTITION_STATS_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

#include "../definitions.h"
#include "../datastructures/dyn_graph.h"
#include "../datastructures/partition_manager.h"
#include "../distance_oracles/matrix_distance_oracle.h"
#include "../utility/qap.h"

namespace HeiProMap {
    // new_partition can be any container indexable by vertex (vector or AlignedArray).
    template<typename NewPartition>
    inline weight_t calculate_migration_cost(const DynGraph &g,
                                             const DistanceOracle &oracle,
                                             const AlignedArray<partition_t> &old_partition,
                                             const NewPartition &new_partition) {
        // New vertices have no previous location, so only old vertices pay migration.
        std::vector<char> is_new = g.new_vertex_mask();
        weight_t total_migration_cost = 0;
        for (vertex_t v = 0; v < g.n; ++v) {
            if (is_new[v]) continue;
            if (old_partition[v] != new_partition[v]) {
                total_migration_cost += g.v_weights[v] * oracle.get(old_partition[v], new_partition[v]);
            }
        }
        return total_migration_cost;
    }

    inline weight_t calculate_edge_cut(const DynGraph &g, const PartitionManager &p_manager) {
        weight_t edge_cut = 0;
        for (vertex_t u = 0; u < g.n; ++u) {
            for (const auto &edge: g.neighbors[u]) {
                if (p_manager[u] != p_manager[edge.u]) {
                    edge_cut += edge.w;
                }
            }
        }
        return edge_cut / 2;
    }

    struct ComponentStats {
        u64 total = 0;
        u64 max_per_block = 0;
        // Component counts by size: 1, 2-9, 10-99, 100 and more vertices.
        u64 by_size[4] = {0, 0, 0, 0};
    };

    // Counts the connected components of each block's induced subgraph via BFS.
    inline ComponentStats calculate_connected_components(const DynGraph &g,
                                                         const PartitionManager &p_manager) {
        ComponentStats stats;
        std::vector<char> visited(g.n, 0);
        std::vector<u64> components_per_block(p_manager.k, 0);
        std::vector<vertex_t> queue;

        for (vertex_t start = 0; start < g.n; ++start) {
            if (!g.vertex_exists(start) || visited[start]) continue;

            partition_t block = p_manager[start];
            components_per_block[block]++;
            stats.total++;

            visited[start] = 1;
            queue.clear();
            queue.push_back(start);
            for (size_t head = 0; head < queue.size(); ++head) {
                vertex_t u = queue[head];
                for (const auto &edge: g.neighbors[u]) {
                    vertex_t v = edge.u;
                    if (visited[v] || p_manager[v] != block) continue;
                    visited[v] = 1;
                    queue.push_back(v);
                }
            }

            size_t size = queue.size();
            if (size == 1) stats.by_size[0]++;
            else if (size < 10) stats.by_size[1]++;
            else if (size < 100) stats.by_size[2]++;
            else stats.by_size[3]++;
        }

        for (u64 c: components_per_block) {
            stats.max_per_block = std::max(stats.max_per_block, c);
        }
        return stats;
    }

    inline void print_dyn_stats(const DynGraph &g,
                      const PartitionManager &p_manager,
                      const DistanceOracle &oracle,
                      weight_t total_migration_cost_from_start,
                      f64 imbalance = 0.03,
                      weight_t total_comm_cost_from_start = 0) {
        std::cout << "--------------------------------------------------" << std::endl;
        std::cout << "Graph Statistics:" << std::endl;
        u64 active_n = 0;
        for (vertex_t v = 0; v < g.n; ++v) {
            if (g.vertex_exists(v)) active_n++;
        }
        std::cout << "  Total Nodes (n): " << g.n << std::endl;
        std::cout << "  Curr. Nodes (n): " << active_n << std::endl;
        std::cout << "  Curr. Edges (m): " << g.m / 2 << std::endl;
        std::cout << "  Curr. Weight:    " << g.g_weight << std::endl;

        if (p_manager.n > 0) {
            weight_t edge_cut = calculate_edge_cut(g, p_manager);
            weight_t max_block_weight = p_manager.max_weight();
            weight_t lmax = (p_manager.k > 0) ? (weight_t) std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k)) : 0;
            double balance = (p_manager.k > 0 && g.g_weight > 0) ? (double) max_block_weight / ((double) g.g_weight / (double) p_manager.k) : 0.0;

            std::cout << "Partition Statistics:" << std::endl;
            std::cout << "  Curr. Edge Cut:       " << edge_cut << std::endl;
            std::cout << "  Curr. Comm Cost:      " << get_qap(g, p_manager.partition.get_vector(), oracle) << std::endl;
            std::cout << "  Curr. Max Block Weight: " << max_block_weight << std::endl;
            std::cout << "  Curr. LMax (Allowed):   " << lmax << std::endl;
            std::cout << "  Curr. Balance (max/avg): " << balance << std::endl;

            ComponentStats components = calculate_connected_components(g, p_manager);
            std::cout << "  Curr. Components (total): " << components.total << std::endl;
            std::cout << "  Curr. Components (max/block): " << components.max_per_block << std::endl;
            std::cout << "  Curr. Components by size (1 / 2-9 / 10-99 / 100+): " << components.by_size[0] << " / " << components.by_size[1] << " / " << components.by_size[2] << " / " << components.by_size[3] << std::endl;

            std::cout << "  Total Migration Cost: " << total_migration_cost_from_start << std::endl;
            std::cout << "  Total Comm Cost: " << total_comm_cost_from_start << std::endl;
        }

        std::cout << "--------------------------------------------------" << std::endl;
    }
}

#endif // HEIPROMAP_DYN_PARTITION_STATS_H
