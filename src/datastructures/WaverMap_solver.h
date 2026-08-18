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

#include <omp.h>
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
#include "../partitioning/greedy_partitioner.h"
#include "../refinement/flow_based_refinement.h"
#include "../WaverMap_configuration.h"
#include "../utility/assert_state.h"
#include "../utility/qap.h"
#include "distance_oracle.h"
#include "../partitioning/kaffpa_partitioner.h"
// #include "../partitioning/recursive_bisection.h"
#include "../utility/translation_table.h"
// #include "HeiPa_solver.h"

namespace HeiProMap {
    class WaverMapSolver {
        WaverMapConfiguration ac;
        RandomEngine random_engine;
        f64 init_time = 0.0;

        // statistics
        weight_t initial_qap = 0;
        weight_t initial_max_block_weight = 0;
        size_t initial_n_empty_partitions = 0;
        size_t initial_n_overloaded_partitions = 0;
        weight_t initial_sum_too_much = 0;

        std::vector<graph_t> graphs;
        graph_t topology_graph;

        PartitionManager p_manager;
        BoundaryVertexManager bv_manager;
        QuotientGraph q_graph;
        BlockConn block_conn;
        d_oracle_t d_oracle;

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
        NegativeCycleRefinement negative_cycle_refinement;

        f64 io_ms = 0.0;
        f64 misc_ms = 0.0;
        f64 coarsening_ms = 0.0;
        f64 contraction_ms = 0.0;
        f64 initial_partitioning_ms = 0.0;
        f64 uncontraction_ms = 0.0;
        f64 rebalance_ms = 0.0;
        f64 refinement_ms = 0.0;
        f64 lp_refine_ms = 0.0;
        f64 qg_refine_ms = 0.0;
        f64 negative_cycle_refine_ms = 0.0;
        f64 flow_refine_ms = 0.0;


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
        explicit WaverMapSolver(const WaverMapConfiguration &t_ac) {
            graphs.reserve(100);
            auto sp_io = get_time_point();
            graphs.emplace_back(t_ac.graph_in);
            io_ms += get_milli_seconds(sp_io, get_time_point());

            auto sp = get_time_point();
            ac = t_ac;
            random_engine = RandomEngine(ac.seed);

            // distance
            topology_graph = CSRGraph(ac.topology_in);
            ac.k = topology_graph.n;
            d_oracle.initialize(topology_graph, ac.threads);

            // manager
            p_manager.initialize(graphs[0].n, ac.k, graphs[0].g_weight);
            bv_manager.initialize(graphs[0].n, ac.k);
            q_graph.initialize(ac.k);
            block_conn.initialize(graphs[0].n, graphs[0].m, ac.k);
            HEAVYASSERT(assert_state_pre_partitioning(graphs[0], p_manager, ac.k));


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
            if (ac.negative_cycle_config.enabled) {
                negative_cycle_refinement.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine.get_u64(), ac.negative_cycle_config);
            }

            auto ep = get_time_point();
            init_time += get_seconds(sp, ep);
        }

        explicit WaverMapSolver(graph_t &&g, const WaverMapConfiguration &t_ac) {
            graphs.reserve(100);
            graphs.emplace_back(std::move(g));

            ac = t_ac;
            random_engine = RandomEngine(ac.seed);

            // distance
            topology_graph = CSRGraph(ac.topology_in);
            ac.k = topology_graph.n;
            d_oracle.initialize(topology_graph, ac.threads);

            // manager
            p_manager.initialize(graphs[0].n, ac.k, graphs[0].g_weight);
            bv_manager.initialize(graphs[0].n, ac.k);
            q_graph.initialize(ac.k);
            block_conn.initialize(graphs[0].n, graphs[0].m, ac.k);
            HEAVYASSERT(assert_state_pre_partitioning(graphs[0], p_manager, ac.k));


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
            if (ac.negative_cycle_config.enabled) {
                negative_cycle_refinement.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine.get_u64(), ac.negative_cycle_config);
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

            internal_solve();

            weight_t qap = get_qap(graphs[0], p_manager, d_oracle);

            std::vector<partition_t> p(graphs[0].n);
            for (vertex_t u = 0; u < graphs[0].n; ++u) { p[u] = p_manager[u]; }
            write_partition(p, ac.mapping_out);

            const auto ep = std::chrono::high_resolution_clock::now();
            f64 duration = get_seconds(sp, ep);

            weight_t lmax = std::ceil((1.0 + ac.imbalance) * ((f64) graphs[0].g_weight / (f64) ac.k));

            std::cout << "Graph                   : " << ac.graph_in << std::endl;
            std::cout << "Total time (s)          : " << duration + init_time << std::endl;
            std::cout << "#Nodes                  : " << graphs[0].n << std::endl;
            std::cout << "#Edges                  : " << graphs[0].m << std::endl;
            std::cout << "k                       : " << ac.k << std::endl;
            std::cout << "Hierarchy               : " << ac.hierarchy_string << std::endl;
            std::cout << "Distances               : " << ac.distance_string << std::endl;
            std::cout << "Lmax                    : " << lmax << std::endl;
            std::cout << "Threads                 : " << ac.threads << std::endl;
            std::cout << "--------------------------" << std::endl;
            std::cout << "Init. QAP               : " << initial_qap << std::endl;
            std::cout << "Init. max block w       : " << initial_max_block_weight << std::endl;
            std::cout << "Init. #empty partitions : " << initial_n_empty_partitions << std::endl;
            std::cout << "Init. #oload partitions : " << initial_n_overloaded_partitions << std::endl;
            std::cout << "Init. Sum oload weights : " << initial_sum_too_much << std::endl;
            std::cout << "--------------------------" << std::endl;
            std::cout << "Final QAP               : " << qap << std::endl;
            std::cout << "max block w             : " << max(p_manager.get_bweights()) << std::endl;

            size_t n_empty_partitions = 0;
            size_t n_overloaded_partitions = 0;
            weight_t sum_too_much = 0;
            for (partition_t id = 0; id < ac.k; ++id) {
                n_empty_partitions += p_manager.get_bweight(id) == 0;
                n_overloaded_partitions += p_manager.get_bweight(id) > lmax;
                sum_too_much += std::max((weight_t) 0, p_manager.get_bweight(id) - lmax);
            }
            std::cout << "#empty partitions       : " << n_empty_partitions << std::endl;
            std::cout << "#oload partitions       : " << n_overloaded_partitions << std::endl;
            std::cout << "Sum oload weights       : " << sum_too_much << std::endl;

            if (ac.hm_level == 0) {
                std::cout << "IO (ms)                 : " << io_ms << std::endl;
                std::cout << "Misc (ms)               : " << misc_ms << std::endl;
                std::cout << "Coarsening (ms)         : " << coarsening_ms << std::endl;
                std::cout << "Contraction (ms)        : " << contraction_ms << std::endl;
                std::cout << "Init. Part. (ms)        : " << initial_partitioning_ms << std::endl;
                std::cout << "Uncontraction (ms)      : " << uncontraction_ms << std::endl;
                std::cout << "Rebalance (ms)          : " << rebalance_ms << std::endl;
                std::cout << "Refinement (ms)         : " << refinement_ms << std::endl;
                std::cout << "  Label Propagation (ms): " << lp_refine_ms << std::endl;
                std::cout << "  Quotient Graph (ms)   : " << qg_refine_ms << std::endl;
                // std::cout << "  Negative Cycle (ms)   : " << negative_cycle_refine_ms << std::endl;
                std::cout << "  Flow (ms)             : " << flow_refine_ms << std::endl;
                std::cout << "ALL (ms)                : " << io_ms + misc_ms + coarsening_ms + contraction_ms + initial_partitioning_ms + uncontraction_ms + rebalance_ms + refinement_ms << std::endl;


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

                std::cout << "[Coarsening] Level " << level << " | Nodes: " << graphs.back().n << std::endl;
                coarsening(level, level_imbalance);
                contraction(level);

                if (graphs.back().n == graphs[graphs.size() - 2].n) {
                    graphs.pop_back();
                    mappings.pop_back();
                    break;
                }

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

            std::cout << "[Initial Partitioning] Level " << level << " | Nodes: " << graphs.back().n << std::endl;
            partition(level, level_imbalance);
            initial_n_empty_partitions = 0;
            initial_n_overloaded_partitions = 0;
            initial_sum_too_much = 0;
            for (partition_t id = 0; id < ac.k; ++id) {
                initial_n_empty_partitions += p_manager.get_bweight(id) == 0;
                initial_n_overloaded_partitions += p_manager.get_bweight(id) > level_lmax;
                initial_sum_too_much += std::max((weight_t) 0, p_manager.get_bweight(id) - level_lmax);
            }

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

            while (!mappings.empty()) {
                level -= 1;

                level_imbalance = ac.imbalance;
                level_lmax = std::ceil((1.0 + level_imbalance) * ((f64) graphs[0].g_weight / (f64) ac.k));
                
                std::cout << "[Uncoarsening] Level " << level << " | Nodes: " << graphs[graphs.size() - 2].n << std::endl;
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
                        greedy_partition(graphs.back(), d_oracle, level_imbalance, ac.seed + thread_id, local_p_managers[thread_id]);

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
                    greedy_partition(graphs.back(), d_oracle, level_imbalance, ac.seed, p_manager);
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

            graphs.back().contract(graphs[prev_idx], mappings.back(), ac.threads);
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

            HEIPROMAP_PROFILE_SCOPE("uncontraction", "misc", "compute_from_scratch");
            const graph_t &g_uncontracted = graphs[graphs.size() - 2];

            if (ac.threads > 1) {
                bv_manager.parallel_reset(ac.threads);
                block_conn.parallel_initialize_offsets(g_uncontracted, ac.threads);

                std::vector<std::vector<std::vector<vertex_t> > > thread_boundaries(ac.threads, std::vector<std::vector<vertex_t> >(ac.k));

                #pragma omp parallel num_threads(ac.threads)
                {
                    u64 tid = omp_get_thread_num();
                    vertex_t chunk = (g_uncontracted.n + ac.threads - 1) / ac.threads;
                    vertex_t start_u = tid * chunk;
                    vertex_t end_u = std::min(g_uncontracted.n, start_u + chunk);

                    for (vertex_t u = start_u; u < end_u; ++u) {
                        const partition_t u_id = p_manager[u];
                        size_t n_different = 0;
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

                        if (n_different > 0) {
                            bv_manager.set_boundary_edges_count(u, n_different);
                            thread_boundaries[tid][u_id].push_back(u);
                        }
                    }
                }

                bv_manager.parallel_import_boundary_vertices(thread_boundaries, ac.threads);
            } else {
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

            HEIPROMAP_PROFILE_SCOPE("uncontraction", "misc", "free_graph");
            graphs.pop_back(); // this is doing uncontraction
            mappings.pop_back();

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
                auto sp_local = get_time_point();
                lp_refine.refine(graphs.back(), d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
                lp_refine_ms += get_milli_seconds(sp_local, get_time_point());
                HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
            }

            if (ac.quotient_graph_refinement_config.enabled) {
                auto sp_local = get_time_point();
                qg_refine.refine(graphs.back(), d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
                qg_refine_ms += get_milli_seconds(sp_local, get_time_point());
                HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
            }

            // if (ac.negative_cycle_config.enabled) {
            //     auto sp_local = get_time_point();
            //     negative_cycle_refinement.refine(graphs.back(), d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints, graphs.back().uniform_v_weights, graphs.back().uniform_e_weights);
            //     negative_cycle_refine_ms += get_milli_seconds(sp_local, get_time_point());
            //     HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
            // }

            if (ac.flow_based_refinement_config.enabled) {
                auto sp_local = get_time_point();
                flow_based_refinement.refine(graphs.back(), d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
                flow_refine_ms += get_milli_seconds(sp_local, get_time_point());
                HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
            }

            auto ep = get_time_point();

            refinement_ms += get_milli_seconds(sp, ep);

            #if ENABLE_PROFILER
            level_infos[level].t_refinement += get_milli_seconds(sp, ep);
            #endif

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
        }

    };
}

#endif //HEIPROMAP_SOLVER_H
