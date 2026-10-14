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

#include "../datastructures/dyn_graph.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/quotient_graph.h"
#include "../datastructures/partition_manager.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/block_conn.h"
#include "../init_partitioning/heipromap_partition.h"
#include "../refinement/label_propagation_refinement.h"
#include "../refinement/dyn_greedy_assign.h"
#include "../rebalance/rebalancer.h"
#include "../configuration/DynHeiProMap_configuration.h"
#include "../utility/profiler.h"
#include "../distance_oracles/matrix_distance_oracle.h"
#include "../utility/dyn_partition_stats.h"
#include "partition_aligner.h"

namespace HeiProMap {
    class DynHeiProMapSolver {
        DynGraph g;

        PartitionManager p_manager;
        BoundaryVertexManager bv_manager;
        BlockConn b_conn;
        Rebalancer rebalancer;

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
            std::cout << "  Dirty:     " << g.dirty_list.size() << " (+" << (g.dirty_list.size() - dirty_before) << ")" << std::endl;
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
                    config.distance.push_back((weight_t) std::stoull(token));
                }

                oracle.initialize(config.hierarchy, config.distance);

                u64 num_blocks = 1;
                for (auto h: config.hierarchy) num_blocks *= h;

                // Reset managers if hierarchy changes
                p_manager.initialize(g.n, num_blocks, g.g_weight);
                p_manager.activate_all_blocks();
                bv_manager.initialize(g.n, num_blocks);
                b_conn.initialize(g.n, g.m, num_blocks);
                q.initialize(num_blocks);
                initial_q.initialize(num_blocks);

                std::cout << "Hierarchy set to " << h_str << " with distances " << d_str << std::endl;
            } else if (cmd == "partition") {
                std::string config_str;
                iss >> config_str;
                run_partition(config_str);
                print_stats();
            } else if (cmd == "save-partition" || cmd == "save_partition") {
                std::string file_path;
                iss >> file_path;
                std::ofstream file(file_path);
                if (file.is_open()) {
                    for (vertex_t u = 0; u < g.n; ++u) {
                        file << p_manager[u] << "\n";
                    }
                    std::cout << "Partition saved to " << file_path << std::endl;
                } else {
                    std::cerr << "Cannot open file: " << file_path << std::endl;
                }
            } else if (cmd == "load-partition" || cmd == "load_partition") {
                std::string file_path;
                iss >> file_path;
                std::ifstream file(file_path);
                if (file.is_open()) {
                    std::vector<partition_t> loaded_part;
                    partition_t p;
                    while (file >> p) {
                        loaded_part.push_back(p);
                    }
                    if (loaded_part.size() == g.n) {
                        for (vertex_t u = 0; u < g.n; ++u) {
                            p_manager.partition[u] = loaded_part[u];
                        }
                        rebuild_managers();
                        std::cout << "Partition loaded from " << file_path << std::endl;
                    } else {
                        std::cerr << "Loaded partition size mismatch: " << loaded_part.size() << " vs " << g.n <<
                                std::endl;
                    }
                } else {
                    std::cerr << "Cannot open file: " << file_path << std::endl;
                }
            } else if (cmd == "refine-fast") {
                u32 iters = 2;
                iss >> iters;
                run_refine_fast(iters);
                print_stats();
            } else if (cmd == "autorefine") {
                f64 threshold_percent = 10.0;
                std::string config_str = "fast";
                iss >> threshold_percent >> config_str;

                f64 current_dirty_percent = (g.n > 0) ? (100.0 * (f64) g.dirty_list_partition.size() / (f64) g.n) : 0.0;
                if (current_dirty_percent >= threshold_percent) {
                    std::cout << "Autorefine: " << current_dirty_percent << "% dirty (>= " << threshold_percent << "%), partitioning (" << config_str << ")..." << std::endl;
                    run_partition(config_str);
                } else {
                    u32 iters = 5;
                    if (config_str == "eco") iters = 10;
                    else if (config_str == "strong") iters = 20;

                    std::cout << "Autorefine: " << current_dirty_percent << "% dirty (< " << threshold_percent << "%), refining (" << iters << " iterations)..." << std::endl;
                    run_refine_fast(iters);
                }
                print_stats();
            } else if (cmd == "no-refinement") {
                // Dummy command, does nothing.
            } else if (cmd == "stats") {
                print_stats();
            }

            auto end = std::chrono::steady_clock::now();
            f64 ms = std::chrono::duration<f64, std::milli>(end - start).count();

            // Only measure significant commands for the user
            command_stats[cmd].add(ms);
        }

        void print_stats() const {
            DynPartitionStats::print(g, p_manager, oracle, previous_partition,
                                     initial_partition, total_migration_cost_from_start,
                                     command_stats, config.imbalance);
        }

        explicit DynHeiProMapSolver(const DynConfiguration &t_config) : config(t_config) {
            oracle.initialize(t_config.hierarchy, t_config.distance);

            u64 num_blocks = 1;
            for (auto h: config.hierarchy) num_blocks *= h;

            p_manager.initialize(0, num_blocks, 0);
            bv_manager.initialize(0, num_blocks);
            b_conn.initialize(0, 0, num_blocks);
            q.initialize(num_blocks);
            initial_q.initialize(num_blocks);

            g.dirty_callback = [this](vertex_t v) {
                if (v >= p_manager.n) {
                    p_manager.resize(v + 1, 0);
                }
            };
        }

        weight_t calculate_migration_cost(const std::vector<partition_t> &old_partition,
                                          const std::vector<partition_t> &new_partition) const {
            return DynPartitionStats::calculate_migration_cost(g, oracle, old_partition, new_partition);
        }

        static u64 count_moved_vertices(const std::vector<partition_t> &old_partition,
                                        const std::vector<partition_t> &new_partition) {
            return DynPartitionStats::count_moved_vertices(old_partition, new_partition);
        }

        void align_partitions_hierarchically(const std::vector<partition_t> &ref_partition,
                                             std::vector<partition_t> &new_partition) const {
            PartitionAligner::align(g, ref_partition, new_partition, config.hierarchy);
        }

    private:
        void rebuild_managers() {
            u64 num_blocks = p_manager.k;
            p_manager.initialize(g.n, num_blocks, g.g_weight);
            p_manager.activate_all_blocks();
            bv_manager.initialize(g.n, num_blocks);
            b_conn.initialize(g.n, g.m, num_blocks);
            q.initialize(num_blocks);

            // This is a full rebuild, so we just iterate over the current graph
            for (vertex_t u = 0; u < g.n; ++u) {
                if (!g.vertex_exists(u)) continue;
                // Note: partition was updated elsewhere, we need to populate managers
                // But in partition() or refine(), managers are usually built from scratch anyway
            }
        }

        void rebuild_q(QuotientGraph &target_q) const {
            u64 num_blocks = p_manager.k;
            target_q.initialize(num_blocks);
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

    public:
        void run_partition(const std::string &config_str) {
            HEIPROMAP_PROFILE_SCOPE("solver", "run_partition", "run_partition");

            std::cout << "Partitioning with config: " << config_str << " (HeiProMap multisection)" << std::endl;

            std::vector<partition_t> new_partition;
            heipromap_partition(g, config.hierarchy, config.distance, config.imbalance, config.seed, config.n_threads, config_str, new_partition);

            u64 num_blocks = p_manager.k;

            if (initial_partition.empty()) {
                // First partition
                p_manager.initialize(g.n, num_blocks, g.g_weight);
                p_manager.activate_all_blocks();
                p_manager.reset_weights();
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
                align_partitions_hierarchically(reference, new_partition);
                weight_t migration_after = calculate_migration_cost(reference, new_partition);

                previous_partition = reference;

                p_manager.initialize(g.n, num_blocks, g.g_weight);
                p_manager.activate_all_blocks();
                p_manager.reset_weights();
                for (vertex_t u = 0; u < g.n; ++u) {
                    p_manager.set(u, g.v_weights[u], new_partition[u]);
                }

                total_migration_cost_from_start += migration_after;

                rebuild_q(q);

                std::cout << "Comm Cost: " << get_qap(g, p_manager.partition.get_vector(), oracle) << "\n";
                std::cout << "Migration Cost: Before=" << migration_before << " After=" << migration_after << " Improvement=" << (migration_before - migration_after) << "\n";
            }

            g.clear_dirty_status();
            g.clear_dirty_status_partition();
            g.clear_new_vertices();
            std::cout << "PARTITION_DONE\n" << std::flush;
        }

        void run_refine_fast(u32 num_iterations) {
            HEIPROMAP_PROFILE_SCOPE("solver", "refine-fast", "refine-fast");

            u64 num_blocks = p_manager.k;

            // Step 2: Resize structures & compute updated lmax
            if (p_manager.n < g.n) {
                p_manager.resize(g.n, 0);
            }
            p_manager.recalculate_weights_dyn(g);

            weight_t lmax = (num_blocks > 0)
                ? (weight_t) std::ceil((1.0 + config.imbalance) * ((f64) g.g_weight / (f64) num_blocks))
                : 0;
            for (partition_t i = 0; i < num_blocks; ++i) {
                p_manager.lmax[i] = lmax;
            }

            size_t num_dirty = g.dirty_list.size();
            size_t num_new = g.new_vertices.size();

            weight_t comm_cost_before = get_qap(g, p_manager.partition.get_vector(), oracle);
            f64 balance_before = (f64) p_manager.max_weight() / ((f64) g.g_weight / (f64) p_manager.k);

            // Step 3: Greedily assign new vertices (has access to updated lmax)
            if (!g.new_vertices.empty()) {
                total_migration_cost_from_start += DynGreedyAssign::assign(g, p_manager, oracle, lmax);
            }

            std::vector<partition_t> partition_before = p_manager.partition.get_vector();

            // Create temporary CSRGraph for rebalancer and refinement
            CSRGraph csr_g(g.n, g.m, g.g_weight);
            for (vertex_t u = 0; u < g.n; ++u) {
                csr_g.v_weights[u] = g.v_weights[u];
                csr_g.neighborhoods[u + 1] = csr_g.neighborhoods[u] + g.neighbors[u].size();
                for (size_t i = 0; i < g.neighbors[u].size(); ++i) {
                    csr_g.edges_v[csr_g.neighborhoods[u] + i] = g.neighbors[u][i].u;
                    csr_g.edges_w[csr_g.neighborhoods[u] + i] = g.neighbors[u][i].w;
                }
            }

            // Populate structures for rebalancing and refinement
            bv_manager.initialize(csr_g.n, num_blocks);
            b_conn.initialize(csr_g.n, csr_g.m, num_blocks);
            b_conn.reset_build();
            q.initialize(num_blocks);

            for (vertex_t u = 0; u < csr_g.n; ++u) {
                b_conn.begin_vertex(csr_g, u);
                partition_t u_id = p_manager[u];
                for (size_t i = csr_g.neighborhoods[u]; i < csr_g.neighborhoods[u + 1]; ++i) {
                    vertex_t v = csr_g.edges_v[i];
                    weight_t w = csr_g.edge_w(i);
                    partition_t v_id = p_manager[v];
                    b_conn.add_connection(u, v_id, w);
                    if (u_id != v_id) {
                        bv_manager.add(u, u_id);
                        if (u < v) q.add_edge(u_id, v_id, w);
                    }
                }
            }

            // Step 4: Rebalance
            rebalancer.initialize(csr_g.n, csr_g.m, num_blocks, config.seed);
            rebalancer.rebalance_last_layer<true, true>(csr_g, p_manager, bv_manager, q, oracle, b_conn);

            // Step 5: Refinement
            LabelPropagationConfiguration lp_config("Label Propagation");
            lp_config.enabled = true;
            lp_config.max_iteration = num_iterations;

            LabelPropagationRefinement<false> lp_refine;
            lp_refine.initialize(csr_g.n, csr_g.m, num_blocks, (u32) config.n_threads, (u32) config.seed, lp_config);

            lp_refine.refine(csr_g, oracle, bv_manager, p_manager, q, b_conn);

            weight_t step_migration = calculate_migration_cost(partition_before, p_manager.partition.get_vector());
            total_migration_cost_from_start += step_migration;
            previous_partition = partition_before;

            g.clear_dirty_status();
            g.clear_new_vertices();

            weight_t comm_cost_after = get_qap(g, p_manager.partition.get_vector(), oracle);
            f64 balance_after = (f64) p_manager.max_weight() / ((f64) g.g_weight / (f64) p_manager.k);

            std::cout << "Incremental Refinement done (" << num_dirty << " dirty, " << num_new << " new vertices):" << std::endl;
            std::cout << "  Comm Cost: Before=" << comm_cost_before << " After=" << comm_cost_after << " (Diff=" << comm_cost_before - comm_cost_after << ")" << std::endl;
            std::cout << "  Balance:   Before=" << balance_before << " After=" << balance_after << std::endl;
            std::cout << "  Migration: Step=" << step_migration << " Total=" << total_migration_cost_from_start << std::endl;
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
