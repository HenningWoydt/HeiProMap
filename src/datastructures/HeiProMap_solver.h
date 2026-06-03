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

#ifndef HEIPROMAP_SOLVER_H
#define HEIPROMAP_SOLVER_H

#include <cmath>

#include "boundary_vertex_manger.h"
#include "partition_manager.h"
#include "quotient_graph.h"
#include "block_conn.h"
#include "../definitions.h"
#include "../utility/macros.h"
#include "../utility/random_engine.h"
#include "../utility/utils.h"
#include "../coarsening/global_path_algorithm.h"
#include "../coarsening/size_constrained_lp.h"
#include "../coarsening/heavy_edge_matching.h"
#include "../rebalance/rebalancer.h"
#include "../partitioning/global_multisection.h"
#include "../refinement/flow_based_refinement.h"
#include "../HeiProMap_configuration.h"
#include "../utility/assert_state.h"
#include "../utility/qap.h"
#include "distance_oracle.h"
#include "../partitioning/kaffpa_partitioner.h"
#include "../partitioning/recursive_bisection.h"
#include "../utility/translation_table.h"
#include "HeiPa_solver.h"

namespace HeiProMap {
    class HeiProMapSolver {
        AlgorithmConfiguration ac;
        RandomEngine random_engine;
        f64 init_time = 0.0;

        // statistics
        weight_t initial_qap = 0;
        weight_t initial_max_block_weight = 0;

        std::vector<graph_t> graphs;

        PartitionManager p_manager;
        BoundaryVertexManager bv_manager;
        QuotientGraph q_graph;
        BlockConn block_conn;
        DistanceOracle d_oracle;

        // matching
        std::vector<Mapping> mappings;
        GlobalPathAlgorithmMatcher gpa_matcher;
        SizeConstrainedLP size_constrained_lp;
        HeavyEdgeMatching heavy_edge_matching;

        Rebalancer rebalancer;

        // refinement
        LabelPropagationRefinement lp_refine;
        QuotientGraphRefinement qg_refine;
        FlowBasedRefinement flow_based_refinement;

        f64 io_ms = 0.0;
        f64 misc_ms = 0.0;
        f64 coarsening_ms = 0.0;
        f64 contraction_ms = 0.0;
        f64 initial_partitioning_ms = 0.0;
        f64 uncontraction_ms = 0.0;
        f64 rebalance_ms = 0.0;
        f64 refinement_ms = 0.0;

        struct level_info {
            u32 level = -1;
            vertex_t n = -1;
            vertex_t m = -1;

            weight_t comm_cost = -1;
            weight_t max_b_weight = -1;
            weight_t lmax = -1;
            f64 imb = -1.0;
            partition_t empty_partitions = -1;
            partition_t oload_partitions = -1;
            weight_t sum_oload_weights = -1;

            f64 t_coarsening = 0.0;
            f64 t_contraction = 0.0;
            f64 t_uncontraction = 0.0;
            f64 t_rebalance = 0.0;
            f64 t_refinement = 0.0;
        };

        std::vector<level_info> level_infos;

        inline void print_level_row(const level_info &L) {
            std::cout
                    << std::setw(3) << L.level << " | "
                    << std::setw(8) << L.n << " | "
                    << std::setw(11) << L.m << " | "
                    << std::setw(10) << L.comm_cost << " | "
                    << std::setw(7) << L.lmax << " | "
                    << std::setw(7) << L.max_b_weight << " | "
                    << std::setw(8) << L.imb << " | "
                    << std::setw(6) << (u32) L.empty_partitions << " | "
                    << std::setw(6) << (u32) L.oload_partitions << " | "
                    << std::setw(8) << L.sum_oload_weights << " | "
                    << std::setw(10) << L.t_coarsening << " | "
                    << std::setw(10) << L.t_contraction << " | "
                    << std::setw(10) << L.t_uncontraction << " | "
                    << std::setw(10) << L.t_rebalance << " | "
                    << std::setw(10) << L.t_refinement
                    << "\n";
        }

        inline void print_all_levels() {
            std::cout
                    << std::setw(3) << "Lvl" << " | "
                    << std::setw(8) << "n" << " | "
                    << std::setw(11) << "m" << " | "
                    << std::setw(10) << "comm cost" << " | "
                    << std::setw(7) << "lmax" << " | "
                    << std::setw(7) << "maxW" << " | "
                    << std::setw(8) << "imb" << " | "
                    << std::setw(6) << "empty" << " | "
                    << std::setw(6) << "oload" << " | "
                    << std::setw(8) << "w_oload" << " | "
                    << std::setw(10) << "t_c" << " | "
                    << std::setw(10) << "t_con" << " | "
                    << std::setw(10) << "t_unc" << " | "
                    << std::setw(10) << "t_reb" << " | "
                    << std::setw(10) << "t_ref"
                    << "\n";

            std::cout << std::string(100, '-') << "\n";
            for (const auto &L: level_infos) {
                print_level_row(L);
            }
        }

    public:
        explicit HeiProMapSolver(const AlgorithmConfiguration &t_ac) {
            graphs.reserve(100);
            auto sp_io = get_time_point();
            graphs.emplace_back(t_ac.graph_in);
            io_ms += get_milli_seconds(sp_io, get_time_point());

            auto sp = get_time_point();
            ac = t_ac;
            random_engine = RandomEngine(ac.seed);

            // manager
            p_manager.initialize(graphs[0].n, ac.k, graphs[0].g_weight);
            bv_manager.initialize(graphs[0].n, ac.k);
            q_graph.initialize(ac.k);
            block_conn.initialize(graphs[0].n, graphs[0].m, ac.k);
            HEAVYASSERT(assert_state_pre_partitioning(graphs[0], p_manager, ac.k));

            // distance
            d_oracle.initialize(ac.hierarchy, ac.distance);

            // matching
            gpa_matcher.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine, ac.global_path_algorithm_config);
            size_constrained_lp.initialize(graphs[0].n, graphs[0].m, ac.k, random_engine.get_u64(), ac.size_constrained_lp_config);

            rebalancer.initialize(graphs[0].n, graphs[0].m, ac.k, random_engine.get_u64());

            // refinement
            if (ac.label_propagation_config.enabled) {
                lp_refine.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine.get_u64(), ac.label_propagation_config);
            }
            if (ac.quotient_graph_refinement_config.enabled) {
                qg_refine.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine.get_u64(), ac.quotient_graph_refinement_config);
            }
            if (ac.flow_based_refinement_config.enabled) {
                flow_based_refinement.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine.get_u64(), ac.flow_based_refinement_config);
            }

            auto ep = get_time_point();
            init_time += get_seconds(sp, ep);
        }

        explicit HeiProMapSolver(graph_t &&g, const AlgorithmConfiguration &t_ac) {
            graphs.reserve(100);
            graphs.emplace_back(std::move(g));

            ac = t_ac;
            random_engine = RandomEngine(ac.seed);

            // manager
            p_manager.initialize(graphs[0].n, ac.k, graphs[0].g_weight);
            bv_manager.initialize(graphs[0].n, ac.k);
            q_graph.initialize(ac.k);
            block_conn.initialize(graphs[0].n, graphs[0].m, ac.k);
            HEAVYASSERT(assert_state_pre_partitioning(graphs[0], p_manager, ac.k));

            // distance
            d_oracle.initialize(ac.hierarchy, ac.distance);

            // matching
            gpa_matcher.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine, ac.global_path_algorithm_config);
            size_constrained_lp.initialize(graphs[0].n, graphs[0].m, ac.k, random_engine.get_u64(), ac.size_constrained_lp_config);

            rebalancer.initialize(graphs[0].n, graphs[0].m, ac.k, random_engine.get_u64());

            // refinement
            if (ac.label_propagation_config.enabled) {
                lp_refine.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine.get_u64(), ac.label_propagation_config);
            }
            if (ac.quotient_graph_refinement_config.enabled) {
                qg_refine.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine.get_u64(), ac.quotient_graph_refinement_config);
            }
            if (ac.flow_based_refinement_config.enabled) {
                flow_based_refinement.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine.get_u64(), ac.flow_based_refinement_config);
            }
        }

        const PartitionManager &solve_subproblem() {
            internal_solve();
            return p_manager;
        }

        const graph_t &get_graph(u64 level) const { return graphs[level]; }

        const PartitionManager &get_p_manager() const { return p_manager; }

        std::vector<vertex_t> solve() {
            const auto sp = std::chrono::high_resolution_clock::now();

            if (ac.hm_level > 0) {
                const weight_t total_weight = graphs[0].g_weight;
                HEAVYASSERT(assert_graph(graphs[0]));
                TranslationTable<vertex_t> tt;
                tt.reserve(graphs[0].n, graphs[0].n);
                for (vertex_t u = 0; u < graphs[0].n; ++u) {
                    tt.add(u, u);
                }

                recursive_solve(graphs[0], p_manager, ac.hierarchy, ac.distance, 0, 0, tt, total_weight);
                p_manager.recalculate_weights(graphs[0]);
                HEAVYASSERT(assert_state_after_partitioning(graphs[0], p_manager, ac.k));

                if (ac.get("--config") == "super-strong") {
                    HEIPROMAP_PROFILE_SCOPE("final_refinement", "Solver", "refine");

                    bv_manager.reset();
                    q_graph.initialize(ac.k);
                    block_conn.initialize(graphs[0].n, graphs[0].m, ac.k);

                    for (vertex_t u = 0; u < graphs[0].n; ++u) {
                        partition_t u_id = p_manager[u];
                        for (size_t i = graphs[0].neighborhoods[u]; i < graphs[0].neighborhoods[u + 1]; ++i) {
                            vertex_t v = graphs[0].edges_v[i];
                            partition_t v_id = p_manager[v];
                            if (u_id != v_id) {
                                bv_manager.add(u, u_id);
                                if (u < v) {
                                    q_graph.add_edge(u_id, v_id, graphs[0].edges_w[i]);
                                }
                            }
                        }
                    }
                    block_conn.compute_from_scratch(graphs[0], p_manager);

                    AlignedArray<weight_t> lmax_constraints;
                    lmax_constraints.initialize(ac.k);
                    weight_t lmax = std::ceil((1.0 + ac.imbalance) * ((f64) graphs[0].g_weight / (f64) ac.k));
                    for (partition_t i = 0; i < ac.k; ++i) {
                        lmax_constraints[i] = lmax;
                    }
                    flow_based_refinement.refine(graphs[0], d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints, graphs[0].uniform_v_weights, graphs[0].uniform_e_weights);
                }
            } else {
                internal_solve();
            }

            weight_t qap = get_qap(graphs[0], p_manager, d_oracle);

            std::vector<partition_t> p(graphs[0].n);
            for (vertex_t u = 0; u < graphs[0].n; ++u) { p[u] = p_manager[u]; }
            write_partition(p, ac.mapping_out);

            const auto ep = std::chrono::high_resolution_clock::now();
            f64 duration = get_seconds(sp, ep);

            weight_t lmax = std::ceil((1.0 + ac.imbalance) * ((f64) graphs[0].g_weight / (f64) ac.k));

            std::cout << "Total time        : " << duration + init_time << std::endl;
            std::cout << "#Nodes            : " << graphs[0].n << std::endl;
            std::cout << "#Edges            : " << graphs[0].m << std::endl;
            std::cout << "k                 : " << ac.k << std::endl;
            std::cout << "Lmax              : " << lmax << std::endl;
            std::cout << "Init. QAP         : " << initial_qap << std::endl;
            std::cout << "Init. max block w : " << initial_max_block_weight << std::endl;
            std::cout << "Final QAP         : " << qap << std::endl;
            std::cout << "max block w       : " << max(p_manager.get_bweights()) << std::endl;

            size_t n_empty_partitions = 0;
            size_t n_overloaded_partitions = 0;
            weight_t sum_too_much = 0;
            for (partition_t id = 0; id < ac.k; ++id) {
                n_empty_partitions += p_manager.get_bweight(id) == 0;
                n_overloaded_partitions += p_manager.get_bweight(id) > lmax;
                sum_too_much += std::max((weight_t) 0, p_manager.get_bweight(id) - lmax);
            }
            std::cout << "#empty partitions : " << n_empty_partitions << std::endl;
            std::cout << "#oload partitions : " << n_overloaded_partitions << std::endl;
            std::cout << "Sum oload weights : " << sum_too_much << std::endl;

            if (ac.hm_level == 0) {
                std::cout << "IO                : " << io_ms << std::endl;
                std::cout << "Misc              : " << misc_ms << std::endl;
                std::cout << "Coarsening        : " << coarsening_ms << std::endl;
                std::cout << "Contraction       : " << contraction_ms << std::endl;
                std::cout << "Init. Part.       : " << initial_partitioning_ms << std::endl;
                std::cout << "Uncontraction     : " << uncontraction_ms << std::endl;
                std::cout << "Rebalance         : " << rebalance_ms << std::endl;
                std::cout << "Refinement        : " << refinement_ms << std::endl;
                std::cout << "ALL               : " << io_ms + misc_ms + coarsening_ms + contraction_ms + initial_partitioning_ms + uncontraction_ms + rebalance_ms + refinement_ms << std::endl;

                #if ENABLE_PROFILER
                print_all_levels();
                #endif
            }

            return p;
        }

    private:
        void internal_solve() {
            u64 level = 0;
            u64 mult = ac.initial_c;

            f64 level_imbalance = 0.0;
            [[maybe_unused]] weight_t level_lmax = 0;

            while (graphs.back().n > ac.k * mult) {
                #if ENABLE_PROFILER
                level_infos.emplace_back();
                level_infos[level].level = level;
                level_infos[level].n = graphs.back().n;
                level_infos[level].m = graphs.back().m;
                #endif

                level_imbalance = ac.imbalance;
                level_lmax = std::ceil((1.0 + level_imbalance) * ((f64) graphs[0].g_weight / (f64) ac.k));

                coarsening(level, level_imbalance);
                contraction(level);

                level += 1;
            }

            #if ENABLE_PROFILER
            level_infos.emplace_back();
            level_infos[level].level = level;
            level_infos[level].n = graphs.back().n;
            level_infos[level].m = graphs.back().m;
            #endif

            level_imbalance = ac.imbalance;
            level_lmax = std::ceil((1.0 + level_imbalance) * ((f64) graphs[0].g_weight / (f64) ac.k));

            partition(level, level_imbalance);
            refinement(level, level_imbalance);

            #if ENABLE_PROFILER
            level_infos[level].max_b_weight = p_manager.max_weight();
            level_infos[level].lmax = level_lmax;
            level_infos[level].imb = (f64) level_infos[level].max_b_weight / ((f64) graphs[0].g_weight / (f64) ac.k);
            level_infos[level].comm_cost = get_qap(graphs.back(), p_manager, d_oracle);
            level_infos[level].empty_partitions = p_manager.n_empty_blocks();
            level_infos[level].oload_partitions = p_manager.n_oload_blocks(level_lmax);
            level_infos[level].sum_oload_weights = p_manager.sum_oload_weight(level_lmax);
            #endif

            while (!mappings.empty()) {
                level -= 1;

                level_imbalance = ac.imbalance;
                level_lmax = std::ceil((1.0 + level_imbalance) * ((f64) graphs[0].g_weight / (f64) ac.k));
                uncontraction(level);

                rebalancing(level, level_imbalance);

                refinement(level, level_imbalance);

                #if ENABLE_PROFILER
                level_infos[level].max_b_weight = p_manager.max_weight();
                level_infos[level].lmax = level_lmax;
                level_infos[level].imb = (f64) level_infos[level].max_b_weight / ((f64) graphs[0].g_weight / (f64) ac.k);
                level_infos[level].comm_cost = get_qap(graphs.back(), p_manager, d_oracle);
                level_infos[level].empty_partitions = p_manager.n_empty_blocks();
                level_infos[level].oload_partitions = p_manager.n_oload_blocks(level_lmax);
                level_infos[level].sum_oload_weights = p_manager.sum_oload_weight(level_lmax);
                #endif
            }
        }

        void partition(u64 level, f64 level_imbalance) {
            auto sp = get_time_point();
            HEIPROMAP_PROFILE_SCOPE("partition", "partition", "partition");

            for (u64 iteration = 0; iteration < 1; ++iteration) {
                if (ac.partitioning_algorithm_id == PARTITIONING_ALG_MULTISECTION) {
                    if (ac.threads > 1) {
                        std::vector<PartitionManager> local_p_managers;
                        std::vector<weight_t> local_qaps;

                        local_p_managers.reserve(ac.threads);
                        local_qaps.reserve(ac.threads);

                        for (u64 thread_id = 0; thread_id < ac.threads; ++thread_id) {
                            local_p_managers.emplace_back();
                            local_p_managers.back().initialize(graphs.back().n, graphs.back().m, graphs.back().g_weight);
                            local_qaps.emplace_back(std::numeric_limits<weight_t>::max());
                        }

                        #pragma omp parallel for schedule(static) num_threads(ac.threads)
                        for (u64 thread_id = 0; thread_id < ac.threads; ++thread_id) {
                            GlobalMultisectionPartitioner partitioner;
                            partitioner.partition(graphs.back(), local_p_managers[thread_id], ac.hierarchy, ac.distance, level_imbalance, ac.global_multisection_config, thread_id);

                            local_qaps[thread_id] = get_qap(graphs.back(), local_p_managers[thread_id], d_oracle);
                        }

                        size_t best_idx = 0;
                        for (u64 thread_id = 0; thread_id < ac.threads; ++thread_id) {
                            if (local_qaps[thread_id] < local_qaps[best_idx]) {
                                best_idx = thread_id;
                            }
                        }

                        p_manager.copy_from(local_p_managers[best_idx]);
                    } else {
                        GlobalMultisectionPartitioner::partition(graphs.back(), p_manager, ac.hierarchy, ac.distance, level_imbalance, ac.global_multisection_config, 0);
                    }
                } else if (ac.partitioning_algorithm_id == PARTITIONING_ALG_GREEDY) {
                    if (ac.threads > 1) {
                        std::vector<PartitionManager> local_p_managers;
                        std::vector<weight_t> local_qaps;

                        local_p_managers.reserve(ac.threads);
                        local_qaps.reserve(ac.threads);

                        for (u64 thread_id = 0; thread_id < ac.threads; ++thread_id) {
                            local_p_managers.emplace_back();
                            local_p_managers.back().initialize(graphs.back().n, ac.k, graphs.back().g_weight);
                            local_qaps.emplace_back(std::numeric_limits<weight_t>::max());
                        }

                        #pragma omp parallel for schedule(static) num_threads(ac.threads)
                        for (u64 thread_id = 0; thread_id < ac.threads; ++thread_id) {
                            std::vector<partition_t> part(graphs.back().n);
                            greedy_partition(graphs.back(), d_oracle, level_imbalance, ac.seed + thread_id, part);
                            for (vertex_t u = 0; u < graphs.back().n; ++u) {
                                local_p_managers[thread_id].set(u, graphs.back().v_weights[u], part[u]);
                            }
                            local_qaps[thread_id] = get_qap(graphs.back(), local_p_managers[thread_id], d_oracle);
                        }

                        size_t best_idx = 0;
                        for (u64 thread_id = 1; thread_id < ac.threads; ++thread_id) {
                            if (local_qaps[thread_id] < local_qaps[best_idx]) {
                                best_idx = thread_id;
                            }
                        }

                        p_manager.copy_from(local_p_managers[best_idx]);
                    } else {
                        std::vector<partition_t> part(graphs.back().n);
                        greedy_partition(graphs.back(), d_oracle, level_imbalance, ac.seed, part);
                        p_manager.reset_weights();
                        for (vertex_t u = 0; u < graphs.back().n; ++u) {
                            p_manager.set(u, graphs.back().v_weights[u], part[u]);
                        }
                    }
                } else {
                    std::cerr << "Partitioning algorithm " << partitioning_algorithm_to_string(ac.partitioning_algorithm_id) << " with id " << ac.partitioning_algorithm_id << " not known!" << std::endl;
                    exit(EXIT_FAILURE);
                }

                // initialize boundary vertices and quotient graph
                HEIPROMAP_PROFILE_SCOPE("partition", "misc", "initialize_datastructures");
                p_manager.reset_weights();
                bv_manager.reset();
                q_graph.initialize(ac.k);
                block_conn.initialize(graphs.back().n, graphs.back().m, ac.k);
                block_conn.reset_build();

                for (vertex_t u = 0; u < graphs.back().n; ++u) {
                    block_conn.begin_vertex(graphs.back(), u);

                    const partition_t u_id = p_manager[u];
                    const weight_t u_w = graphs.back().v_weights[u];
                    p_manager.set(u, u_w, u_id);

                    for (size_t i = graphs.back().neighborhoods[u]; i < graphs.back().neighborhoods[u + 1]; ++i) {
                        const vertex_t v = graphs.back().edges_v[i];
                        const weight_t w = graphs.back().edges_w[i];
                        const partition_t v_id = p_manager[v];

                        // build block_conns directly here
                        block_conn.add_connection(u, v_id, w);

                        if (u_id != v_id) {
                            bv_manager.add(u, u_id); // boundary vertex
                            if (u < v) {
                                q_graph.add_edge(u_id, v_id, w); // quotient graph
                            }
                        }
                    }
                }

                initial_qap = get_qap(graphs.back(), p_manager, d_oracle);
                initial_max_block_weight = max(p_manager.get_bweights());
            }
            auto ep = get_time_point();
            initial_partitioning_ms += get_milli_seconds(sp, ep);

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
        }

        void coarsening(const u64 level, f64 level_imbalance) {
            auto sp = get_time_point();
            mappings.emplace_back();
            mappings.back().initialize(graphs.back().n);

            if (ac.coarsening_algorithm_id == COARSENING_ALG_GLOBAL_PATHS) {
                gpa_matcher.match(level, graphs.back(), p_manager, mappings.back(), level_imbalance);
            } else if (ac.coarsening_algorithm_id == COARSENING_ALG_SIZE_CONSTRAINED_LP) {
                size_constrained_lp.cluster(level, graphs.back(), p_manager, mappings.back(), level_imbalance, ac.threads);
            } else if (ac.coarsening_algorithm_id == COARSENING_ALG_HEAVY_EDGE) {
                heavy_edge_matching.match(graphs.back(), p_manager, mappings.back(), level_imbalance, random_engine.get_u64(), ac.heavy_edge_matching_config);
            } else {
                std::cerr << "Coarsening algorithm " << coarsening_algorithm_to_string(ac.coarsening_algorithm_id) << " with id " << ac.coarsening_algorithm_id << " not known!" << std::endl;
                exit(EXIT_FAILURE);
            }


            auto ep = get_time_point();
            coarsening_ms += get_milli_seconds(sp, ep);

            #if ENABLE_PROFILER
            level_infos[level].t_coarsening += get_milli_seconds(sp, ep);
            #endif
        }

        void contraction([[maybe_unused]] const u64 level) {
            auto sp = get_time_point();

            size_t prev_idx = graphs.size() - 1;
            graphs.emplace_back(); // coarse the graph

            if (ac.threads == 1) {
                graphs.back().initialize<false, false>(graphs[prev_idx], mappings.back());
            } else {
                graphs.back().parallel_initialize<false, false>(graphs[prev_idx], mappings.back(), ac.threads);
            }
            p_manager.contract(mappings.back());

            auto ep = get_time_point();
            contraction_ms += get_milli_seconds(sp, ep);

            #if ENABLE_PROFILER
            level_infos[level].t_contraction += get_milli_seconds(sp, ep);
            #endif

            HEAVYASSERT(assert_state_pre_partitioning(graphs.back(), p_manager, ac.k));
        }

        void uncontraction([[maybe_unused]] const u64 level) {
            auto sp = get_time_point();

            p_manager.uncontract(mappings.back());
            //
            {
                HEIPROMAP_PROFILE_SCOPE("uncontraction", "misc", "compute_from_scratch");
                const graph_t &g_uncontracted = graphs[graphs.size() - 2];

                bv_manager.reset();
                block_conn.initialize(g_uncontracted.n, g_uncontracted.m, ac.k);
                block_conn.reset_build();

                for (vertex_t u = 0; u < g_uncontracted.n; ++u) {
                    const partition_t u_id = p_manager[u];
                    size_t n_different = 0;

                    block_conn.begin_vertex(g_uncontracted, u);
                    weight_t own_weight = 0;

                    for (size_t i = g_uncontracted.neighborhoods[u]; i < g_uncontracted.neighborhoods[u + 1]; ++i) {
                        const vertex_t v = g_uncontracted.edges_v[i];
                        const weight_t w = g_uncontracted.edges_w[i];
                        const partition_t v_id = p_manager[v];

                        // rebuild block connections
                        if (u_id == v_id) {
                            own_weight += w;
                        } else {
                            block_conn.add_connection(u, v_id, w);
                        }

                        // rebuild boundary information
                        n_different += (u_id != v_id);
                    }
                    if (own_weight > 0) {
                        block_conn.add_connection(u, u_id, own_weight);
                    }
                    bv_manager.add_boundary_vertex_from_count(u, u_id, n_different);
                }
            }
            //
            {
                HEIPROMAP_PROFILE_SCOPE("uncontraction", "misc", "free_graph");
                graphs.pop_back(); // this is doing uncontraction
                mappings.pop_back();
            }

            auto ep = get_time_point();

            uncontraction_ms += get_milli_seconds(sp, ep);

            #if ENABLE_PROFILER
            level_infos[level].t_uncontraction += get_milli_seconds(sp, ep);
            #endif

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(),p_manager, bv_manager, q_graph, block_conn , ac.k));
        }

        void rebalancing(const u64 level, f64 level_imbalance) {
            auto sp = get_time_point();

            if (level == 0) {
                rebalancer.rebalance_last_layer(graphs.back(), p_manager, bv_manager, q_graph, d_oracle, block_conn, level_imbalance);
            } else {
                rebalancer.rebalance(graphs.back(), p_manager, bv_manager, q_graph, d_oracle, block_conn, level_imbalance);
            }

            auto ep = get_time_point();

            rebalance_ms += get_milli_seconds(sp, ep);

            #if ENABLE_PROFILER
            level_infos[level].t_rebalance += get_milli_seconds(sp, ep);
            #endif

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
        }

        void refinement(const u64 level, const f64 level_imbalance) {
            auto sp = get_time_point();

            AlignedArray<weight_t> lmax_constraints;
            lmax_constraints.initialize(ac.k);
            weight_t lmax = std::ceil((1.0 + level_imbalance) * ((f64) graphs.back().g_weight / (f64) ac.k));
            for (partition_t i = 0; i < ac.k; ++i) {
                lmax_constraints[i] = lmax;
            }

            if (ac.label_propagation_config.enabled) {
                lp_refine.refine(graphs.back(), d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints, graphs.back().uniform_v_weights, graphs.back().uniform_e_weights);
                HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
            }

            if (ac.quotient_graph_refinement_config.enabled) {
                qg_refine.refine(graphs.back(), d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints, graphs.back().uniform_v_weights, graphs.back().uniform_e_weights);
                HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
            }

            if (ac.flow_based_refinement_config.enabled) {
                flow_based_refinement.refine(graphs.back(), d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints, graphs.back().uniform_v_weights, graphs.back().uniform_e_weights);
                HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
            }

            auto ep = get_time_point();

            refinement_ms += get_milli_seconds(sp, ep);

            #if ENABLE_PROFILER
            level_infos[level].t_refinement += get_milli_seconds(sp, ep);
            #endif

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
        }

        void recursive_solve(graph_t &g, PartitionManager &p_manager, std::vector<partition_t> hierarchy,
                             std::vector<weight_t> distance, u64 current_level, u64 offset, const TranslationTable<vertex_t> &tt,
                             const weight_t total_weight) {
            partition_t k_of_subgraph = prod<partition_t>(hierarchy);
            f64 total_remaining_slack = ((1.0 + ac.imbalance) * (f64) k_of_subgraph * (f64) total_weight) / ((f64) ac.k * (f64) g.g_weight) - 1.0;
            total_remaining_slack = std::max(0.0, total_remaining_slack);

            if (current_level >= ac.hm_level) {
                AlgorithmConfiguration sub_ac = ac;
                sub_ac.hierarchy = hierarchy;
                sub_ac.distance = distance;
                sub_ac.k = k_of_subgraph;
                sub_ac.imbalance = total_remaining_slack;

                if (ac.get("--config") == "fast") {
                    sub_ac.set_fast();
                } else if (ac.get("--config") == "eco") {
                    sub_ac.set_eco();
                } else if (ac.get("--config") == "strong") {
                    sub_ac.set_strong();
                } else if (ac.get("--config") == "super-strong") {
                    sub_ac.set_super_strong();
                }

                std::vector<weight_t> v_weights(g.n);
                const weight_t *v_weights_ptr = g.v_weights.get_ptr();
                std::copy(v_weights_ptr, v_weights_ptr + g.n, v_weights.begin());
                HeiProMapSolver sub_solver(std::move(g), sub_ac);
                const PartitionManager &sub_p_manager = sub_solver.solve_subproblem();
                HEAVYASSERT(assert_state_after_partitioning(sub_solver.get_graph(0), sub_p_manager, sub_p_manager.k));

                for (vertex_t u = 0; u < sub_p_manager.n; ++u) {
                    p_manager.set(tt.get_o(u), v_weights[u], offset + sub_p_manager[u]);
                }
                HEAVYASSERT(assert_state_partial(sub_solver.get_graph(0), p_manager, tt, offset, k_of_subgraph));
                return;
            }

            partition_t k = hierarchy.back();
            hierarchy.pop_back();
            distance.pop_back();

            f64 per_level_epsilon = std::pow(1.0 + total_remaining_slack, 1.0 / (f64) (hierarchy.size() + 1)) - 1.0;
            per_level_epsilon = std::max(0.0, per_level_epsilon);

            AlignedArray<partition_t> partition;
            partition.initialize(g.n, 0);
            //
            {
                HEIPROMAP_PROFILE_SCOPE("adaptive_solver", "adaptive_solver", "partition");
                if (ac.global_multisection_config.mode == GLOBAL_MULTISECTION_KAFFPA_STRONG) {
                    kaffpa_partition(g, k, per_level_epsilon, KAFFPA_PARTITION_STRONG, ac.seed, partition, ac.global_multisection_config.kappa, ac.collect_dataset, ac.data_dir);
                } else if (ac.global_multisection_config.mode == GLOBAL_MULTISECTION_KAFFPA_ECO) {
                    kaffpa_partition(g, k, per_level_epsilon, KAFFPA_PARTITION_ECO, ac.seed, partition, ac.global_multisection_config.kappa, ac.collect_dataset, ac.data_dir);
                } else if (ac.global_multisection_config.mode == GLOBAL_MULTISECTION_KAFFPA_FAST) {
                    kaffpa_partition(g, k, per_level_epsilon, KAFFPA_PARTITION_FAST, ac.seed, partition, ac.global_multisection_config.kappa, ac.collect_dataset, ac.data_dir);
                } else if (ac.global_multisection_config.mode >= GLOBAL_MULTISECTION_HEIPA_FAST && ac.global_multisection_config.mode <= GLOBAL_MULTISECTION_HEIPA_SUPER_STRONG) {
                    heipa_partition(g, k, per_level_epsilon, ac.seed, partition, ac.global_multisection_config.mode, ac.global_multisection_config.kappa);
                } else if (ac.global_multisection_config.mode == GLOBAL_MULTISECTION_METIS_KWAY) {
                    kway_partition(g, k, per_level_epsilon, ac.seed, partition, ac.global_multisection_config.kappa);
                } else {
                    std::cerr << "Mode " << ac.global_multisection_config.mode << " not implemented" << std::endl;
                    abort();
                }
            }

            std::vector<vertex_t> new_ns(k, 0);
            std::vector<vertex_t> new_ms(k, 0);
            std::vector<weight_t> new_ws(k, 0);
            for (vertex_t u = 0; u < g.n; ++u) {
                partition_t u_id = partition[u];
                new_ns[u_id] += 1;
                new_ws[u_id] += g.v_weights[u];
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    vertex_t v = g.edges_v[i];
                    partition_t v_id = partition[v];
                    if (v_id == u_id) { new_ms[u_id] += 1; }
                }
            }

            partition_t k_per_subgraph = prod<partition_t>(hierarchy);

            for (partition_t i = 0; i < k; ++i) {
                HEIPROMAP_PROFILE_SCOPE("adaptive_solver", "adaptive_solver", "extrcact_graph");
                graph_t sub_g(new_ns[i], new_ms[i], new_ws[i]);
                TranslationTable<vertex_t> sub_tt;
                sub_tt.reserve(new_ns[i], p_manager.n);

                std::vector<vertex_t> new_us(k, 0);
                for (vertex_t old_u = 0; old_u < g.n; ++old_u) {
                    if (partition[old_u] == i) {
                        sub_tt.add(tt.get_o(old_u), new_us[i]);
                        new_us[i] += 1;
                    }
                }

                std::vector<vertex_t> degrees(sub_g.n, 0);
                for (vertex_t old_u = 0; old_u < g.n; ++old_u) {
                    if (partition[old_u] == i) {
                        vertex_t new_u = sub_tt.get_n(tt.get_o(old_u));
                        for (size_t j = g.neighborhoods[old_u]; j < g.neighborhoods[old_u + 1]; ++j) {
                            vertex_t old_v = g.edges_v[j];
                            if (partition[old_v] == i) {
                                degrees[new_u]++;
                            }
                        }
                    }
                }

                sub_g.neighborhoods[0] = 0;
                for (vertex_t j = 0; j < sub_g.n; ++j) {
                    sub_g.neighborhoods[j + 1] = sub_g.neighborhoods[j] + degrees[j];
                }

                std::vector<vertex_t> cursor(sub_g.n, 0);
                for (vertex_t old_u = 0; old_u < g.n; ++old_u) {
                    if (partition[old_u] == i) {
                        vertex_t new_u = sub_tt.get_n(tt.get_o(old_u));
                        sub_g.v_weights[new_u] = g.v_weights[old_u];

                        for (size_t j = g.neighborhoods[old_u]; j < g.neighborhoods[old_u + 1]; ++j) {
                            vertex_t old_v = g.edges_v[j];
                            if (partition[old_v] == i) {
                                vertex_t new_v = sub_tt.get_n(tt.get_o(old_v));
                                size_t pos = sub_g.neighborhoods[new_u] + cursor[new_u];
                                sub_g.edges_v[pos] = new_v;
                                sub_g.edges_w[pos] = g.edges_w[j];
                                cursor[new_u]++;
                            }
                        }
                    }
                }
                HEAVYASSERT(assert_graph(sub_g));
                recursive_solve(sub_g, p_manager, hierarchy, distance, current_level + 1, offset + i * k_per_subgraph, sub_tt, total_weight);
            }

            HEAVYASSERT(assert_state_partial(g, p_manager, tt, offset, k_of_subgraph));
        }
    };
}

#endif //HEIPROMAP_SOLVER_H
