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

#ifndef HEIPROMAP_DYN_SOLVER_H
#define HEIPROMAP_DYN_SOLVER_H

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>

#include "dyn_graph.h"
#include "csr_graph.h"
#include "quotient_graph.h"
#include "partition_manager.h"
#include "boundary_vertex_manger.h"
#include "block_conn.h"
#include "../partitioning/heipromap_partition.h"
#include "../refinement/label_propagation_refinement.h"
#include "../DynHeiProMap_configuration.h"
#include "../utility/profiler.h"
#include "distance_oracle.h"
#include "../utility/hungarian.h"
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
            return calls > 0 ? total_ms / calls : 0;
        }
    };

    class DynHeiProMapSolver {
        DynGraph g;
        
        PartitionManager p_manager;
        BoundaryVertexManager bv_manager;
        BlockConn b_conn;

        std::vector<partition_t> previous_partition;
        std::vector<partition_t> initial_partition;

        QuotientGraph q;
        QuotientGraph initial_q;

        DynConfiguration config;
        DistanceOracle oracle;
        std::map<std::string, KTStat> command_stats;

        weight_t total_migration_cost_from_start = 0;

        void exec_file(const std::string &file_path) {
            HEIPROMAP_PROFILE_SCOPE("solver", "exec_file", "exec_file");
            std::ifstream file(file_path);
            if (!file.is_open()) {
                std::cerr << "Cannot open file: " << file_path << std::endl;
                return;
            }

            vertex_t n_before = g.n;
            vertex_t m_before = g.m;
            size_t dirty_before = g.dirty_list.size();

            std::string line;
            u32 lines_processed = 0;
            while (std::getline(file, line)) {
                if (line.empty() || line[0] == '#') continue;
                process_command(line);
                lines_processed++;
            }

            std::cout << "Finished reading " << file_path << ":" << std::endl;
            std::cout << "  Processed: " << lines_processed << " commands" << std::endl;
            std::cout << "  Vertices:  " << g.n << " (+" << (g.n - n_before) << ")" << std::endl;
            std::cout << "  Edges:     " << g.m / 2 << " (+" << (g.m - m_before) / 2 << ")" << std::endl;
            std::cout << "  Dirty:     " << g.dirty_list.size() << " (+" << (g.dirty_list.size() - dirty_before) << ")"
                    << std::endl;
        }

    public:
        void process_command(const std::string &line) {
            std::istringstream iss(line);
            std::string cmd;
            iss >> cmd;
            if (cmd.empty()) return;

            auto start = std::chrono::steady_clock::now();

            if (cmd == "+v") {
                vertex_t v;
                weight_t w = 1;
                iss >> v >> w;
                g.add_vertex(v, w);
            } else if (cmd == "+e") {
                vertex_t u, v;
                weight_t w = 1;
                iss >> u >> v >> w;
                g.add_edge(u, v, w);
            } else if (cmd == "-v") {
                vertex_t v;
                weight_t w;
                iss >> v;
                if (iss >> w) {
                    g.remove_vertex(v, w);
                } else {
                    g.remove_vertex(v);
                }
            } else if (cmd == "-e") {
                vertex_t u, v;
                weight_t w;
                iss >> u >> v;
                if (iss >> w) {
                    g.remove_edge(u, v, w);
                } else {
                    g.remove_edge(u, v);
                }
            } else if (cmd == "vw") {
                vertex_t v;
                weight_t w;
                iss >> v >> w;
                g.set_vertex_weight(v, w);
            } else if (cmd == "ew") {
                vertex_t u, v;
                weight_t w;
                iss >> u >> v >> w;
                g.set_edge_weight(u, v, w);
            } else if (cmd == "exec") {
                std::string file_path;
                iss >> file_path;
                exec_file(file_path);
            } else if (cmd == "set-hierarchy") {
                std::string h_str, d_str;
                iss >> h_str >> d_str;

                config.hierarchy.clear();
                std::istringstream h_iss(h_str);
                std::string token;
                while (std::getline(h_iss, token, ':')) {
                    config.hierarchy.push_back(std::stoull(token));
                }

                config.distance.clear();
                std::istringstream d_iss(d_str);
                while (std::getline(d_iss, token, ':')) {
                    config.distance.push_back(std::stoull(token));
                }

                oracle.initialize(config.hierarchy, config.distance);
                
                u64 num_blocks = 1;
                for (auto h : config.hierarchy) num_blocks *= h;
                
                // Reset managers if hierarchy changes
                p_manager.initialize(g.n, (partition_t)num_blocks, g.g_weight);
                bv_manager.initialize(g.n, (partition_t)num_blocks);
                b_conn.initialize(g.n, g.m, (partition_t)num_blocks);
                q.initialize((partition_t)num_blocks);
                initial_q.initialize((partition_t)num_blocks);
                
                std::cout << "Hierarchy set to " << h_str << " with distances " << d_str << std::endl;
            } else if (cmd == "partition") {
                std::string config_str;
                iss >> config_str;
                run_partition(config_str);
            } else if (cmd == "save-partition") {
                std::string file_path;
                iss >> file_path;

                std::ofstream file(file_path);
                if (!file.is_open()) {
                    std::cerr << "Cannot open file: " << file_path << std::endl;
                    return;
                }

                for (vertex_t v = 0; v < p_manager.n; ++v) {
                    file << p_manager[v] << "\n";
                }
                std::cout << "Partition saved to " << file_path << "\n";
            } else if (cmd == "save-graph") {
                std::string file_path;
                iss >> file_path;
                g.save_to_metis(file_path);
                std::cout << "Graph saved to " << file_path << "\n";
            } else if (cmd == "refine-fast") {
                u32 num_iterations = 5;
                iss >> num_iterations;
                run_refine_fast(num_iterations);
            } else if (cmd == "autorefine") {
                double threshold_percent = 10.0;
                std::string config_str = "fast";
                iss >> threshold_percent >> config_str;

                double current_dirty_percent = (g.n > 0) ? (100.0 * g.dirty_list_partition.size() / g.n) : 0.0;
                if (current_dirty_percent >= threshold_percent) {
                    std::cout << "Autorefine: " << current_dirty_percent << "% dirty (>= " << threshold_percent <<
                            "%), partitioning (" << config_str << ")..." << std::endl;
                    run_partition(config_str);
                } else {
                    u32 iters = 5;
                    if (config_str == "eco") iters = 10;
                    else if (config_str == "strong") iters = 20;

                    std::cout << "Autorefine: " << current_dirty_percent << "% dirty (< " << threshold_percent <<
                            "%), refining (" << iters << " iterations)..." << std::endl;
                    run_refine_fast(iters);
                }
            } else if (cmd == "no-refinement") {
                // Dummy command, does nothing.
            } else if (cmd == "stats") {
                std::cout << "--------------------------------------------------" << std::endl;
                std::cout << "Graph Statistics:" << std::endl;
                std::cout << "  Nodes (n):      " << g.n << std::endl;
                std::cout << "  Edges (m):      " << g.m << std::endl;
                std::cout << "  Total Weight:   " << g.g_weight << std::endl;
                if (p_manager.n > 0) {
                    weight_t edge_cut = 0;
                    for (vertex_t u = 0; u < g.n; ++u) {
                        for (const auto &edge: g.neighbors[u]) {
                            if (p_manager[u] != p_manager[edge.u]) {
                                edge_cut += edge.w;
                            }
                        }
                    }
                    edge_cut /= 2;
                    weight_t max_block_weight = p_manager.max_weight();
                    double balance = (double) max_block_weight / ((double)g.g_weight / (double)p_manager.k);

                    std::cout << "Partition Statistics:" << std::endl;
                    std::cout << "  Edge Cut:       " << edge_cut << std::endl;
                    std::cout << "  Comm Cost:      " << get_qap(g, p_manager.partition.get_vector(), oracle) << std::endl;
                    std::cout << "  Max Block (LMax): " << max_block_weight << std::endl;
                    std::cout << "  Balance (max/avg): " << balance << std::endl;

                    if (!previous_partition.empty() && previous_partition.size() == p_manager.n) {
                        std::cout << "  Migration Cost (prev): " << calculate_migration_cost(
                            previous_partition, p_manager.partition.get_vector()) << std::endl;
                        std::cout << "  Moved Vertices (prev): " << count_moved_vertices(previous_partition, p_manager.partition.get_vector())
                                << std::endl;
                    }
                    if (!initial_partition.empty() && initial_partition.size() <= p_manager.n) {
                        std::cout << "  Migration Cost (init): " << calculate_migration_cost(
                            initial_partition, p_manager.partition.get_vector()) << std::endl;
                        std::cout << "  Moved Vertices (init): " << count_moved_vertices(initial_partition, p_manager.partition.get_vector())
                                << std::endl;
                    }
                    std::cout << "  Cumulative Migration Cost: " << total_migration_cost_from_start << std::endl;
                }
                std::cout << "--------------------------------------------------" << std::endl;
                std::cout << "Command Timings (User):" << std::endl;
                std::cout << "  " << std::left << std::setw(20) << "Command"
                        << std::right << std::setw(8) << "Calls"
                        << std::setw(15) << "Total (ms)"
                        << std::setw(15) << "Avg (ms)" << std::endl;
                for (const auto &[name, stat]: command_stats) {
                    std::cout << "  " << std::left << std::setw(20) << name
                            << std::right << std::setw(8) << stat.calls
                            << std::setw(15) << std::fixed << std::setprecision(3) << stat.total_ms
                            << std::setw(15) << stat.avg() << std::endl;
                }
                std::cout << "--------------------------------------------------" << std::endl;
            }

            auto end = std::chrono::steady_clock::now();
            double ms = std::chrono::duration<double, std::milli>(end - start).count();

            // Only measure significant commands for the user
            command_stats[cmd].add(ms);
        }

        explicit DynHeiProMapSolver(DynConfiguration &t_config) : config(t_config) {
            oracle.initialize(t_config.hierarchy, t_config.distance);
            
            u64 num_blocks = 1;
            for (auto h : config.hierarchy) num_blocks *= h;
            
            p_manager.initialize(0, (partition_t)num_blocks, 0);
            bv_manager.initialize(0, (partition_t)num_blocks);
            b_conn.initialize(0, 0, (partition_t)num_blocks);
            q.initialize((partition_t)num_blocks);
            initial_q.initialize((partition_t)num_blocks);

            g.dirty_callback = [this](vertex_t v) {
                // Placeholder for incremental manager updates if needed
            };
        }

        weight_t calculate_migration_cost(const std::vector<partition_t> &old_partition,
                                          const std::vector<partition_t> &new_partition) {
            weight_t total_migration_cost = 0;
            vertex_t num_vertices = std::min((vertex_t)old_partition.size(), (vertex_t)new_partition.size());
            for (vertex_t v = 0; v < num_vertices; ++v) {
                if (old_partition[v] != new_partition[v]) {
                    total_migration_cost += g.v_weights[v] * oracle.get(old_partition[v], new_partition[v]);
                }
            }
            return total_migration_cost;
        }

        u64 count_moved_vertices(const std::vector<partition_t> &old_partition,
                                 const std::vector<partition_t> &new_partition) {
            u64 moved_count = 0;
            vertex_t num_vertices = std::min((vertex_t)old_partition.size(), (vertex_t)new_partition.size());
            for (vertex_t v = 0; v < num_vertices; ++v) {
                if (old_partition[v] != new_partition[v]) {
                    moved_count++;
                }
            }
            return moved_count;
        }

        void align_partitions_hierarchically(std::vector<partition_t> &new_partition) {
            if (previous_partition.empty()) return;

            bool all_zeros = true;
            for (auto p: previous_partition) {
                if (p != 0) {
                    all_zeros = false;
                    break;
                }
            }
            if (all_zeros) return;

            int num_levels = (int) config.hierarchy.size();
            std::vector<partition_t> D(num_levels + 1, 1);
            for (int i = 0; i < num_levels; ++i) D[i + 1] = D[i] * config.hierarchy[i];

            std::vector<std::map<partition_t, std::vector<int> > > inv_m(num_levels);
            run_alignment_recursive(num_levels - 1, 0, 0, new_partition, D, inv_m);

            for (vertex_t v = 0; v < g.n; ++v) {
                partition_t old_new_id = new_partition[v];
                partition_t aligned_id = 0;
                for (int l = num_levels - 1; l >= 0; --l) {
                    partition_t prefix_new = old_new_id / D[l + 1];
                    int comp_new = (int) ((old_new_id / D[l]) % config.hierarchy[l]);
                    int comp_old = inv_m[l][prefix_new][comp_new];
                    aligned_id = aligned_id * config.hierarchy[l] + comp_old;
                }
                new_partition[v] = aligned_id;
            }
        }

    private:
        void rebuild_managers() {
            u64 num_blocks = p_manager.k;
            p_manager.initialize(g.n, (partition_t)num_blocks, g.g_weight);
            bv_manager.initialize(g.n, (partition_t)num_blocks);
            b_conn.initialize(g.n, g.m, (partition_t)num_blocks);
            q.initialize((partition_t)num_blocks);

            // This is a full rebuild, so we just iterate over the current graph
            for (vertex_t u = 0; u < g.n; ++u) {
                if (!g.vertex_exists(u)) continue;
                // Note: partition was updated elsewhere, we need to populate managers
                // But in partition() or refine(), managers are usually built from scratch anyway
            }
        }

        void rebuild_q(QuotientGraph &target_q) {
            u64 num_blocks = p_manager.k;
            target_q.initialize((partition_t)num_blocks);
            for (vertex_t u = 0; u < g.n; ++u) {
                if (!g.vertex_exists(u)) continue;
                partition_t u_id = p_manager[u];
                for (const auto &edge: g.neighbors[u]) {
                    if (u < edge.u) {
                        partition_t v_id = p_manager[edge.u];
                        if (u_id != v_id) {
                            target_q.add_edge(u_id, v_id, edge.w);
                        }
                    }
                }
            }
        }

        void run_alignment_recursive(int level, partition_t p_old, partition_t p_new,
                                     std::vector<partition_t> &new_partition,
                                     const std::vector<partition_t> &D,
                                     std::vector<std::map<partition_t, std::vector<int> > > &inv_m) {
            if (level < 0) return;

            std::vector<int> m = match_blocks_general(level, p_old, p_new, new_partition, D);

            std::vector<int> inv(m.size());
            for (size_t i = 0; i < m.size(); ++i) inv[m[i]] = (int) i;
            inv_m[level][p_new] = inv;

            int count = (int) config.hierarchy[level];
            for (int i = 0; i < count; ++i) {
                run_alignment_recursive(level - 1, p_old * count + i, p_new * count + m[i], new_partition, D, inv_m);
            }
        }

        std::vector<int> match_blocks_general(int level, partition_t p_old, partition_t p_new,
                                              const std::vector<partition_t> &new_partition,
                                              const std::vector<partition_t> &D) {
            int count = (int) config.hierarchy[level];
            std::vector<std::vector<weight_t> > cost(count, std::vector<weight_t>(count, 0));
            int num_levels = (int) config.hierarchy.size();

            for (vertex_t v = 0; v < g.n; ++v) {
                partition_t id_old = previous_partition[v];
                partition_t id_new = new_partition[v];

                if (level < num_levels - 1) {
                    if (id_old / D[level + 1] != p_old || id_new / D[level + 1] != p_new) continue;
                }

                cost[(id_old / D[level]) % count][(id_new / D[level]) % count] += g.v_weights[v];
            }

            weight_t max_w = 0;
            for (int i = 0; i < count; ++i)
                for (int j = 0; j < count; ++j) max_w = std::max(max_w, cost[i][j]);

            std::vector<std::vector<weight_t> > neg_cost(count, std::vector<weight_t>(count));
            for (int i = 0; i < count; ++i)
                for (int j = 0; j < count; ++j) neg_cost[i][j] = max_w - cost[i][j];

            return solve_hungarian(neg_cost);
        }

        void run_partition(const std::string &config_str) {
            HEIPROMAP_PROFILE_SCOPE("solver", "run_partition", "run_partition");

            std::cout << "Partitioning with config: " << config_str << " (HeiProMap multisection)" << std::endl;
            
            std::vector<partition_t> new_partition;
            heipromap_partition(g, config.hierarchy, config.distance, config.imbalance, config.seed, config.n_threads, config_str, new_partition);

            u64 num_blocks = p_manager.k;

            if (initial_partition.empty()) {
                // First partition
                p_manager.initialize(g.n, (partition_t)num_blocks, g.g_weight);
                for (vertex_t u = 0; u < g.n; ++u) {
                    p_manager.set(u, g.v_weights[u], new_partition[u]);
                }
                
                initial_partition = new_partition;
                previous_partition = new_partition;

                rebuild_q(q);
                rebuild_q(initial_q);

                std::cout << "Comm Cost: " << get_qap(g, p_manager.partition.get_vector(), oracle) << "\n";
                std::cout << "Migration Cost: Initialized baselines.\n";
            } else {
                // Subsequent partition
                std::vector<partition_t> reference = p_manager.partition.get_vector();
                if (reference.size() < g.n) reference.resize(g.n, 0);

                weight_t migration_before = calculate_migration_cost(reference, new_partition);
                align_partitions_hierarchically(new_partition);
                weight_t migration_after = calculate_migration_cost(reference, new_partition);

                previous_partition = reference; 
                
                p_manager.initialize(g.n, (partition_t)num_blocks, g.g_weight);
                for (vertex_t u = 0; u < g.n; ++u) {
                    p_manager.set(u, g.v_weights[u], new_partition[u]);
                }
                
                total_migration_cost_from_start += migration_after;

                rebuild_q(q);

                std::cout << "Comm Cost: " << get_qap(g, p_manager.partition.get_vector(), oracle) << "\n";
                std::cout << "Migration Cost: Before=" << migration_before << " After=" << migration_after <<
                        " Improvement=" << (migration_before - migration_after) << "\n";
            }

            g.clear_dirty_status();
            g.clear_dirty_status_partition();
            g.clear_new_vertices();
            std::cout << "PARTITION_DONE\n" << std::flush;
        }

        void run_refine_fast(u32 num_iterations) {
            HEIPROMAP_PROFILE_SCOPE("solver", "refine-fast", "refine-fast");

            if (p_manager.n < g.n) {
                // Should not happen if incremental updates are implemented, 
                // but for now we re-initialize if needed.
                p_manager.initialize(g.n, p_manager.k, g.g_weight);
                // Need to recover partition data... this highlights why we need incremental updates.
            }

            size_t num_dirty = g.dirty_list.size();
            size_t num_new = g.new_vertices.size();
            
            weight_t comm_cost_before = get_qap(g, p_manager.partition.get_vector(), oracle);
            double balance_before = (double)p_manager.max_weight() / ((double)g.g_weight / (double)p_manager.k);

            u64 num_blocks = p_manager.k;
            weight_t allowed_max_block_weight = (weight_t) ((1.0 + config.imbalance) * ((f64) g.g_weight / (f64) num_blocks));

            // 1. Greedy assignment for new vertices
            if (!g.new_vertices.empty()) {
                total_migration_cost_from_start += greedy_assign(allowed_max_block_weight);
            }

            std::vector<partition_t> partition_before = p_manager.partition.get_vector();

            // 2. Incremental Refinement
            {
                // Create temporary CSRGraph (Still needed because refinement algorithms expect CSRGraph)
                ::HeiProMap::CSRGraph csr_g(g.n, g.m, g.g_weight);
                for (vertex_t u = 0; u < g.n; ++u) {
                    csr_g.v_weights[u] = g.v_weights[u];
                    csr_g.neighborhoods[u + 1] = csr_g.neighborhoods[u] + g.neighbors[u].size();
                    for (size_t i = 0; i < g.neighbors[u].size(); ++i) {
                        csr_g.edges_v[csr_g.neighborhoods[u] + i] = g.neighbors[u][i].u;
                        csr_g.edges_w[csr_g.neighborhoods[u] + i] = g.neighbors[u][i].w;
                    }
                }

                // Initialize managers (fully for now, later incrementally)
                bv_manager.initialize(csr_g.n, (partition_t)num_blocks);
                b_conn.initialize(csr_g.n, csr_g.m, (partition_t)num_blocks);
                b_conn.reset_build();
                q.initialize((partition_t)num_blocks);

                // Populate structures
                for (vertex_t u = 0; u < csr_g.n; ++u) {
                    b_conn.begin_vertex(csr_g, u);
                    partition_t u_id = p_manager[u];
                    for (size_t i = csr_g.neighborhoods[u]; i < csr_g.neighborhoods[u + 1]; ++i) {
                        vertex_t v = csr_g.edges_v[i];
                        weight_t w = csr_g.edges_w[i];
                        partition_t v_id = p_manager[v];
                        b_conn.add_connection(u, v_id, w);
                        if (u_id != v_id) {
                            bv_manager.add(u, u_id);
                            if (u < v) q.add_edge(u_id, v_id, w);
                        }
                    }
                }

                // Initialize and run refinement
                ::HeiProMap::LabelPropagationConfiguration lp_config("Label Propagation");
                lp_config.enabled = true;
                lp_config.max_iteration = num_iterations;

                ::HeiProMap::LabelPropagationRefinement lp_refine;
                lp_refine.initialize(csr_g.n, csr_g.m, (partition_t)num_blocks, (u32)config.n_threads, (u32)config.seed, lp_config);

                AlignedArray<weight_t> lmax_constraints;
                lmax_constraints.initialize((partition_t)num_blocks);
                weight_t lmax = std::ceil((1.0 + config.imbalance) * ((f64) csr_g.g_weight / (f64) num_blocks));
                for (partition_t i = 0; i < (partition_t)num_blocks; ++i) {
                    lmax_constraints[i] = lmax;
                }

                lp_refine.refine(csr_g, oracle, bv_manager, p_manager, q, b_conn, lmax_constraints);
            }

            total_migration_cost_from_start += calculate_migration_cost(partition_before, p_manager.partition.get_vector());

            g.clear_dirty_status();
            g.clear_new_vertices();

            weight_t comm_cost_after = get_qap(g, p_manager.partition.get_vector(), oracle);
            double balance_after = (double)p_manager.max_weight() / ((double)g.g_weight / (double)p_manager.k);

            std::cout << "Incremental Refinement done (" << num_dirty << " dirty, " << num_new << " new vertices):" <<
                    std::endl;
            std::cout << "  Comm Cost: Before=" << comm_cost_before << " After=" << comm_cost_after << " (Diff=" << (
                s64) comm_cost_before - (s64) comm_cost_after << ")" << std::endl;
            std::cout << "  Balance:   Before=" << balance_before << " After=" << balance_after << std::endl;
        }

        weight_t greedy_assign(weight_t max_block_weight) {
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

                partition_t best_block = 0;
                weight_t min_comm_cost = std::numeric_limits<weight_t>::max();

                for (partition_t target_block: candidates) {
                    weight_t current_target_weight = p_manager.bweights[target_block];
                    if (target_block != current_block && current_target_weight + v_weight > max_block_weight) {
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

                if (best_block != current_block) {
                    migration_cost += v_weight * oracle.get(current_block, best_block);
                    p_manager.move_serial(v, v_weight, current_block, best_block);
                }
            }
            return migration_cost;
        }
    };

    inline void interactive_mode(DynHeiProMapSolver &solver) {
        std::string line;
        while (std::getline(std::cin, line)) {
            if (line == "quit") break;
            solver.process_command(line);
        }
    }

    inline void execute_commands(DynHeiProMapSolver &solver, const std::vector<std::string> &commands) {
        for (const auto &line: commands) {
            solver.process_command(line);
        }
    }
}

#endif // HEIPROMAP_DYN_SOLVER_H
