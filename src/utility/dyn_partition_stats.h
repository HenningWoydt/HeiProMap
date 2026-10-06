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

    struct KTStat {
        u32 calls = 0;
        double total_ms = 0;

        void add(double ms) {
            calls++;
            total_ms += ms;
        }

        double avg() const {
            return calls > 0 ? total_ms / calls : 0.0;
        }
    };

    class DynPartitionStats {
    public:
        static weight_t calculate_migration_cost(const DynGraph &g,
                                                const DistanceOracle &oracle,
                                                const std::vector<partition_t> &old_partition,
                                                const std::vector<partition_t> &new_partition) {
            weight_t total_migration_cost = 0;
            vertex_t num_vertices = std::min((vertex_t) old_partition.size(), (vertex_t) new_partition.size());
            for (vertex_t v = 0; v < num_vertices; ++v) {
                if (old_partition[v] != new_partition[v]) {
                    total_migration_cost += g.v_weights[v] * oracle.get(old_partition[v], new_partition[v]);
                }
            }
            return total_migration_cost;
        }

        static u64 count_moved_vertices(const std::vector<partition_t> &old_partition,
                                        const std::vector<partition_t> &new_partition) {
            u64 moved_count = 0;
            vertex_t num_vertices = std::min((vertex_t) old_partition.size(), (vertex_t) new_partition.size());
            for (vertex_t v = 0; v < num_vertices; ++v) {
                if (old_partition[v] != new_partition[v]) {
                    moved_count++;
                }
            }
            return moved_count;
        }

        static weight_t calculate_edge_cut(const DynGraph &g, const PartitionManager &p_manager) {
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

        static void print(const DynGraph &g,
                          const PartitionManager &p_manager,
                          const DistanceOracle &oracle,
                          const std::vector<partition_t> &previous_partition,
                          const std::vector<partition_t> &initial_partition,
                          weight_t total_migration_cost_from_start,
                          const std::map<std::string, KTStat> &command_stats,
                          f64 imbalance = 0.03) {
            std::cout << "--------------------------------------------------" << std::endl;
            std::cout << "Graph Statistics:" << std::endl;
            std::cout << "  Nodes (n):      " << g.n << std::endl;
            std::cout << "  Edges (m):      " << g.m << std::endl;
            std::cout << "  Total Weight:   " << g.g_weight << std::endl;

            if (p_manager.n > 0) {
                weight_t edge_cut = calculate_edge_cut(g, p_manager);
                weight_t max_block_weight = p_manager.max_weight();
                weight_t lmax = (p_manager.k > 0)
                    ? (weight_t) std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k))
                    : 0;
                double balance = (p_manager.k > 0 && g.g_weight > 0)
                    ? (double) max_block_weight / ((double) g.g_weight / (double) p_manager.k)
                    : 0.0;

                std::cout << "Partition Statistics:" << std::endl;
                std::cout << "  Edge Cut:       " << edge_cut << std::endl;
                std::cout << "  Comm Cost:      " << get_qap(g, p_manager.partition.get_vector(), oracle) << std::endl;
                std::cout << "  Max Block Weight: " << max_block_weight << std::endl;
                std::cout << "  LMax (Allowed):   " << lmax << std::endl;
                std::cout << "  Balance (max/avg): " << balance << std::endl;

                if (!previous_partition.empty()) {
                    std::cout << "  Migration Cost (prev): " << calculate_migration_cost(g, oracle, previous_partition, p_manager.partition.get_vector()) << std::endl;
                    std::cout << "  Moved Vertices (prev): " << count_moved_vertices(previous_partition, p_manager.partition.get_vector()) << std::endl;
                }
                if (!initial_partition.empty()) {
                    std::cout << "  Migration Cost (init): " << calculate_migration_cost(g, oracle, initial_partition, p_manager.partition.get_vector()) << std::endl;
                    std::cout << "  Moved Vertices (init): " << count_moved_vertices(initial_partition, p_manager.partition.get_vector()) << std::endl;
                }
                std::cout << "  Cumulative Migration Cost: " << total_migration_cost_from_start << std::endl;
            }

            std::cout << "--------------------------------------------------" << std::endl;
            std::cout << "Command Timings (User):" << std::endl;
            std::cout << "  " << std::left << std::setw(20) << "Command" << std::right << std::setw(8) << "Calls" << std::setw(15) << "Total (ms)" << std::setw(15) << "Avg (ms)" << std::endl;
            for (const auto &[name, stat]: command_stats) {
                std::cout << "  " << std::left << std::setw(20) << name << std::right << std::setw(8) << stat.calls << std::setw(15) << std::fixed << std::setprecision(3) << stat.total_ms << std::setw(15) << stat.avg() << std::endl;
            }
            std::cout << "--------------------------------------------------" << std::endl;
        }
    };

}

#endif // HEIPROMAP_DYN_PARTITION_STATS_H
