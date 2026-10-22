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
#include "../utility/partition_aligner.h"

namespace HeiProMap {
    class DynHeiProMapSolver {
        DynGraph g;

        PartitionManager p_manager;
        BoundaryVertexManager bv_manager;
        BlockConn b_conn;
        QuotientGraph q_graph;

        PartitionManager prev_p_manager;

        DynConfiguration config;
        DistanceOracle oracle;

        Rebalancer rebalancer;

        // Vertices that exist in the graph but have not been placed into a block yet.
        std::vector<char> unplaced;

        weight_t total_migration_cost_from_start = 0;
        weight_t total_comm_cost_cost_from_start = 0;

        // Applies the buffered graph modifications and keeps p_manager in sync: it grows with the
        // graph and every block weight / vertex count follows the net weight change of its vertices.
        // New vertices stay in block 0 (their placeholder) until they are placed.
        void apply_pending() {
            vertex_t n_before = p_manager.n;
            std::vector<DynGraph::VertexChange> changes = g.apply_pending();
            if (changes.empty()) return;

            if (p_manager.n < g.n) {
                p_manager.resize(g.n, 0);
            }
            unplaced.resize(g.n, 0);

            for (const auto &c: changes) {
                bool was_unplaced = unplaced[c.v];
                bool unplaced_now = c.is_new || (was_unplaced && c.after > 0);
                partition_t block = p_manager[c.v];

                // Unplaced vertices are not counted in any block until DynGreedyAssign places them.
                weight_t counted_before = was_unplaced ? 0 : c.before;
                weight_t counted_after = unplaced_now ? 0 : c.after;
                p_manager.bweights[block] += counted_after - counted_before;
                if (c.is_new && !was_unplaced && c.v < n_before) p_manager.n_vertices[block] -= 1;

                unplaced[c.v] = unplaced_now;
            }
        }

        void exec_file(const std::string &file_path) {
            HEIPROMAP_PROFILE_SCOPE("solver", "exec_file", "exec_file");
            std::ifstream file(file_path);
            if (!file.is_open()) {
                std::cerr << "Cannot open file: " << file_path << std::endl;
                return;
            }

            std::string line;
            while (std::getline(file, line)) {
                if (line.empty() || line[0] == '#') continue;
                process_command(line);
            }
        }

    public:
        void process_command(const std::string &line) {
            std::istringstream iss(line);
            std::string cmd;
            iss >> cmd;
            if (cmd.empty()) return;

            // Graph modifications are only buffered here, every other command applies them first.
            if (cmd == "+v") {
                DynGraph::GraphOp op{DynGraph::GraphOp::ADD_V};
                iss >> op.u >> op.w;
                g.pending_ops.push_back(op);
            } else if (cmd == "+e") {
                DynGraph::GraphOp op{DynGraph::GraphOp::ADD_E};
                iss >> op.u >> op.v >> op.w;
                g.pending_ops.push_back(op);
            } else if (cmd == "-v") {
                DynGraph::GraphOp op{DynGraph::GraphOp::REMOVE_V};
                iss >> op.u;
                op.has_w = static_cast<bool>(iss >> op.w);
                g.pending_ops.push_back(op);
            } else if (cmd == "-e") {
                DynGraph::GraphOp op{DynGraph::GraphOp::REMOVE_E};
                iss >> op.u >> op.v;
                op.has_w = static_cast<bool>(iss >> op.w);
                g.pending_ops.push_back(op);
            } else if (cmd == "vw") {
                DynGraph::GraphOp op{DynGraph::GraphOp::SET_VW};
                iss >> op.u >> op.w;
                g.pending_ops.push_back(op);
            } else if (cmd == "ew") {
                DynGraph::GraphOp op{DynGraph::GraphOp::SET_EW};
                iss >> op.u >> op.v >> op.w;
                g.pending_ops.push_back(op);
            } else if (cmd == "exec") {
                std::string file_path;
                iss >> file_path;
                exec_file(file_path);
            } else if (cmd == "set-hierarchy") {
                apply_pending();
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
                q_graph.initialize(num_blocks);

                // Everything sits in block 0 again, so nothing is unplaced anymore.
                unplaced.assign(g.n, 0);
                g.clear_new_vertices();
                g.clear_removed_vertices();
            } else if (cmd == "partition") {
                apply_pending();
                std::string config_str;
                iss >> config_str;
                run_partition(config_str);
            } else if (cmd == "save-partition" || cmd == "save_partition") {
                apply_pending();
                std::string file_path;
                iss >> file_path;
                std::ofstream file(file_path);
                if (file.is_open()) {
                    // One line per external vertex id, -1 for ids that are not in the graph.
                    for (vertex_t ext = 0; ext < g.ext_to_int.size(); ++ext) {
                        vertex_t slot = g.ext_to_int[ext];
                        if (slot == DynGraph::INVALID_ID) file << -1 << "\n";
                        else file << p_manager[slot] << "\n";
                    }
                    std::cout << "Partition saved to " << file_path << std::endl;
                } else {
                    std::cerr << "Cannot open file: " << file_path << std::endl;
                }
            } else if (cmd == "load-partition" || cmd == "load_partition") {
                apply_pending();
                std::string file_path;
                iss >> file_path;
                std::ifstream file(file_path);
                if (file.is_open()) {
                    // One line per external vertex id, negative values mean "not in the graph".
                    s64 value;
                    vertex_t ext = 0;
                    u64 loaded = 0;
                    while (file >> value) {
                        vertex_t slot = g.find_slot(ext++);
                        if (value < 0 || slot == DynGraph::INVALID_ID) continue;
                        p_manager.partition[slot] = (partition_t) value;
                        loaded++;
                    }
                    p_manager.recalculate_weights_dyn(g);
                    rebuild_q(q_graph);
                    unplaced.assign(g.n, 0);
                    g.clear_new_vertices();
                    g.clear_removed_vertices();
                    std::cout << "Partition loaded from " << file_path << " (" << loaded << " vertices)" << std::endl;
                } else {
                    std::cerr << "Cannot open file: " << file_path << std::endl;
                }
            } else if (cmd == "refine-fast") {
                apply_pending();
                u32 iters = 2;
                iss >> iters;
                run_refine_fast(iters);
            } else if (cmd == "autorefine") {
                apply_pending();
                f64 threshold_percent = 10.0;
                std::string config_str = "fast";
                iss >> threshold_percent >> config_str;

                f64 current_dirty_percent = (g.n > 0) ? (100.0 * (f64) g.dirty_list_partition.size() / (f64) g.n) : 0.0;
                if (current_dirty_percent >= threshold_percent) {
                    run_partition(config_str);
                } else {
                    u32 iters = 5;
                    if (config_str == "eco") iters = 10;
                    else if (config_str == "strong") iters = 20;

                    run_refine_fast(iters);
                }
            } else if (cmd == "greedy-assign") {
                apply_pending();
                run_greedy_assign();
            } else if (cmd == "no-op") {
                // The application advances without refinement, so it still pays the current comm cost.
                apply_pending();
                total_comm_cost_cost_from_start += get_qap(g, p_manager.partition.get_vector(), oracle);
            } else if (cmd == "stats") {
                apply_pending();
                print_dyn_stats(g, p_manager, oracle, total_migration_cost_from_start, config.imbalance, total_comm_cost_cost_from_start);
            }
        }

        explicit DynHeiProMapSolver(const DynConfiguration &t_config) : config(t_config) {
            oracle.initialize(t_config.hierarchy, t_config.distance);

            u64 num_blocks = 1;
            for (auto h: config.hierarchy) num_blocks *= h;

            p_manager.initialize(0, num_blocks, 0);
            bv_manager.initialize(0, num_blocks);
            b_conn.initialize(0, 0, num_blocks);
            q_graph.initialize(num_blocks);

            g.dirty_callback = [this](vertex_t v) {
                if (v >= p_manager.n) {
                    p_manager.resize(v + 1, 0);
                }
            };
        }

    private:
        void snapshot_previous_p_manager() {
            HEIPROMAP_PROFILE_SCOPE("solver", "misc", "snapshot_previous_p_manager");
            prev_p_manager = p_manager;
        }

        // Unit cost a new vertex pays on placement: it is charged the largest hierarchy distance.
        weight_t max_distance() const {
            return config.distance.empty() ? 0 : *std::max_element(config.distance.begin(), config.distance.end());
        }

        weight_t new_vertices_cost() const {
            std::vector<char> is_new = g.new_vertex_mask();
            weight_t cost = 0;
            for (vertex_t v = 0; v < g.n; ++v) {
                if (is_new[v] && g.vertex_exists(v)) cost += g.v_weights[v] * max_distance();
            }
            return cost;
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

            // prev_p_manager is only filled by a snapshot, so it is empty before the first step
            const bool first_partition = (prev_p_manager.n == 0);
            snapshot_previous_p_manager();

            std::vector<partition_t> new_partition;
            heipromap_partition(g, config.hierarchy, config.distance, config.imbalance, config.seed, config.n_threads, config_str, new_partition);

            if (first_partition) {
                // First partition
                p_manager.initialize(g.n, p_manager.k, g.g_weight);
                p_manager.activate_all_blocks();
                p_manager.reset_weights();
                for (vertex_t u = 0; u < g.n; ++u) {
                    p_manager.set(u, g.v_weights[u], new_partition[u]);
                }

                rebuild_q(q_graph);
            } else {
                // Subsequent partition
                // prev_p_manager holds the partition from before this step (already grown to g.n).
                align_partitions(g, prev_p_manager.partition, new_partition, config.hierarchy);
                weight_t migration_after = calculate_migration_cost(g, oracle, prev_p_manager.partition, new_partition);

                p_manager.initialize(g.n, p_manager.k, g.g_weight);
                p_manager.activate_all_blocks();
                p_manager.reset_weights();
                for (vertex_t u = 0; u < g.n; ++u) {
                    p_manager.set(u, g.v_weights[u], new_partition[u]);
                }

                total_migration_cost_from_start += migration_after + new_vertices_cost();

                rebuild_q(q_graph);
            }

            total_comm_cost_cost_from_start += get_qap(g, p_manager.partition.get_vector(), oracle);

            g.clear_dirty_status();
            g.clear_dirty_status_partition();
            g.clear_new_vertices();
            unplaced.assign(g.n, 0);
            g.clear_removed_vertices();
        }

        // Updates lmax for the current total weight and greedily places all unplaced new vertices.
        void assign_new_vertices() {
            if (p_manager.n < g.n) {
                p_manager.resize(g.n, 0);
            }

            weight_t lmax = (p_manager.k > 0) ? (weight_t) std::ceil((1.0 + config.imbalance) * ((f64) g.g_weight / (f64) p_manager.k)) : 0;
            for (partition_t i = 0; i < p_manager.k; ++i) {
                p_manager.lmax[i] = lmax;
            }

            if (!g.new_vertices.empty()) {
                total_migration_cost_from_start += DynGreedyAssign::assign(g, p_manager, oracle, g.new_vertices, lmax, max_distance());
                std::fill(unplaced.begin(), unplaced.end(), 0);
            }
        }

        // Temporary CSRGraph for the rebalancer and the refinement.
        CSRGraph build_csr() const {
            HEIPROMAP_PROFILE_SCOPE("solver", "refine-fast", "create_csr");
            CSRGraph csr_g(g.n, g.m, g.g_weight);
            for (vertex_t u = 0; u < g.n; ++u) {
                csr_g.v_weights[u] = g.v_weights[u];
                csr_g.neighborhoods[u + 1] = csr_g.neighborhoods[u] + g.neighbors[u].size();
                for (size_t i = 0; i < g.neighbors[u].size(); ++i) {
                    csr_g.edges_v[csr_g.neighborhoods[u] + i] = g.neighbors[u][i].u;
                    csr_g.edges_w[csr_g.neighborhoods[u] + i] = g.neighbors[u][i].w;
                }
            }
            return csr_g;
        }

        // Populate bv_manager, b_conn and q_graph for rebalancing and refinement.
        void populate_structures(const CSRGraph &csr_g) {
            HEIPROMAP_PROFILE_SCOPE("solver", "refine-fast", "update_structures");
            bv_manager.initialize(csr_g.n, p_manager.k);
            b_conn.initialize(csr_g.n, csr_g.m, p_manager.k);
            b_conn.reset_build();
            q_graph.initialize(p_manager.k);

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
                        if (u < v) q_graph.add_edge(u_id, v_id, w);
                    }
                }
            }
        }

        void rebalance_graph(const CSRGraph &csr_g) {
            rebalancer.initialize(csr_g.n, csr_g.m, p_manager.k, config.seed);
            rebalancer.rebalance_last_layer<true, true>(csr_g, p_manager, bv_manager, q_graph, oracle, b_conn);
        }

        // Updates the bookkeeping that follows every step that changed the partition.
        void finish_step() {
            total_migration_cost_from_start += calculate_migration_cost(g, oracle, prev_p_manager.partition, p_manager.partition);

            g.clear_dirty_status();
            g.clear_new_vertices();
            g.clear_removed_vertices();

            total_comm_cost_cost_from_start += get_qap(g, p_manager.partition.get_vector(), oracle);
        }

        // Places the new vertices greedily and rebalances, without any refinement.
        void run_greedy_assign() {
            HEIPROMAP_PROFILE_SCOPE("solver", "greedy-assign", "greedy-assign");

            snapshot_previous_p_manager();
            assign_new_vertices();

            CSRGraph csr_g = build_csr();
            populate_structures(csr_g);
            rebalance_graph(csr_g);

            finish_step();
        }

        void run_refine_fast(u32 num_iterations) {
            HEIPROMAP_PROFILE_SCOPE("solver", "refine-fast", "refine-fast");

            snapshot_previous_p_manager();

            // Step 2 & 3: Update lmax and greedily assign new vertices
            assign_new_vertices();

            CSRGraph csr_g = build_csr();
            populate_structures(csr_g);

            // Step 4: Rebalance
            rebalance_graph(csr_g);

            // Step 5: Refinement
            LabelPropagationConfiguration lp_config("Label Propagation");
            lp_config.enabled = true;
            lp_config.max_iteration = num_iterations;

            LabelPropagationRefinement<false> lp_refine;
            lp_refine.initialize(csr_g.n, csr_g.m, p_manager.k, (u32) config.n_threads, (u32) config.seed, lp_config);

            lp_refine.refine(csr_g, oracle, bv_manager, p_manager, q_graph, b_conn);

            finish_step();
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
