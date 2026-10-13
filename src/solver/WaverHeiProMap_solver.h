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

#ifndef HEIPROMAP_WAVER_MAP_SOLVER_H
#define HEIPROMAP_WAVER_MAP_SOLVER_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <omp.h>

#include "../datastructures/block_conn.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/large_quotient_graph.h"
#include "../datastructures/partition_manager.h"
#include "../definitions.h"
#include "../configuration/WaverHeiProMap_configuration.h"
#include "../utility/macros.h"
#include "../utility/random_engine.h"
#include "../utility/assert_state.h"
#include "../utility/qap.h"
#include "../utility/utils.h"

#include "../coarsening/global_path_algorithm.h"
#include "../coarsening/size_constrained_lp.h"
#include "../coarsening/heavy_edge_matching.h"
#include "../rebalance/rebalancer.h"
#include "../init_partitioning/greedy_partitioner.h"
#include "../refinement/flow_based_refinement.h"
#include "../distance_oracles/grid_distance_oracle.h"

namespace HeiProMap {
    class WaverHeiProMapSolver {
        WaverHeiProMapConfiguration ac;
        RandomEngine random_engine;

        std::chrono::high_resolution_clock::time_point sp;

        std::vector<graph_t> graphs;

        p_manager_t p_manager;
        bv_manager_t bv_manager;
        LargeQuotientGraph q_graph;
        block_conn_t block_conn;
        GridDistanceOracle d_oracle;

        // matching
        std::vector<Mapping> mappings;
        GlobalPathAlgorithmMatcher gpa_matcher;
        SizeConstrainedLP size_constrained_lp;
        HeavyEdgeMatching heavy_edge_matching;

        Rebalancer rebalancer;

        // refinement
        LabelPropagationRefinement<false> lp_refine;
        QuotientGraphRefinement<false> qg_refine;
        FlowBasedRefinement<false> flow_based_refinement;

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
        f64 flow_refine_ms = 0.0;

        void activate_all_blocks() {
            for (partition_t id = 0; id < ac.k; ++id) {
                p_manager.set_hierarchy_level(id, 0);
            }
        }

        void update_lmax_array() {
            const CSRGraph &topo = d_oracle.current_topology();
            std::vector<weight_t> lmax_vec(ac.k);
            if (topo.uniform_v_weights) {
                weight_t lmax = std::ceil((1.0 + ac.imbalance) * ((f64) graphs[0].g_weight / (f64) ac.k));
                for (partition_t i = 0; i < ac.k; ++i) {
                    lmax_vec[i] = lmax;
                }
            } else {
                f64 topo_total_w = topo.g_weight > 0 ? (f64) topo.g_weight : (f64) ac.k;
                for (partition_t i = 0; i < ac.k; ++i) {
                    double fair_share = (double) topo.v_weight(i) / topo_total_w;
                    lmax_vec[i] = std::ceil((1.0 + ac.imbalance) * graphs[0].g_weight * fair_share);
                }
            }
            for (partition_t i = 0; i < ac.k; ++i) {
                p_manager.set_lmax(i, lmax_vec[i]);
            }
        }

        void reinit_structures_for_new_k() {
            p_manager.initialize(graphs[0].n, ac.k, graphs[0].g_weight);
            activate_all_blocks();
            update_lmax_array();
            bv_manager.initialize(graphs[0].n, ac.k);
            q_graph.initialize(ac.k);
            block_conn.initialize(graphs[0].n, graphs[0].m, ac.k);

            rebalancer.initialize(graphs[0].n, graphs[0].m, ac.k, random_engine.get_u64());
            if (ac.label_propagation_config.enabled) lp_refine.update_k(ac.k);
            if (ac.quotient_graph_refinement_config.enabled) qg_refine.update_k(ac.k);
            if (ac.flow_based_refinement_config.enabled) flow_based_refinement.update_k(ac.k);
        }

        void initialize_common(CSRGraph topology) {
            random_engine = RandomEngine(ac.seed);

            d_oracle.initialize(std::move(topology), ac.threads);
            ac.k = d_oracle.current_k();

            p_manager.initialize(graphs[0].n, ac.k, graphs[0].g_weight);
            activate_all_blocks();
            update_lmax_array();
            bv_manager.initialize(graphs[0].n, ac.k);
            q_graph.initialize(ac.k);
            block_conn.initialize(graphs[0].n, graphs[0].m, ac.k);
            HEAVYASSERT(assert_state_pre_partitioning(graphs[0], p_manager, ac.k));

            gpa_matcher.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine, ac.global_path_algorithm_config);
            size_constrained_lp.initialize(graphs[0].n, graphs[0].m, ac.k, random_engine.get_u64(), ac.size_constrained_lp_config);

            rebalancer.initialize(graphs[0].n, graphs[0].m, ac.k, random_engine.get_u64());

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

    public:
        explicit WaverHeiProMapSolver(const WaverHeiProMapConfiguration &t_ac) {
            sp = std::chrono::high_resolution_clock::now();

            ac = t_ac;

            auto sp_io = get_time_point();
            graphs.emplace_back(ac.graph_in);
            io_ms += get_milli_seconds(sp_io, get_time_point());

            if (ac.topology_in.empty()) {
                std::cerr << "Topology required for WaverMap!" << std::endl;
                exit(1);
            }

            initialize_common(CSRGraph(ac.topology_in));
        }

        explicit WaverHeiProMapSolver(graph_t &&g, graph_t &&topology, const WaverHeiProMapConfiguration &t_ac) {
            sp = std::chrono::high_resolution_clock::now();

            ac = t_ac;
            graphs.emplace_back(std::move(g));

            initialize_common(std::move(topology));
        }

        std::vector<vertex_t> solve() {
            internal_solve();

            std::vector<partition_t> p(graphs[0].n);
            for (vertex_t u = 0; u < graphs[0].n; ++u) { p[u] = p_manager[u]; }

            const auto ep = std::chrono::high_resolution_clock::now();
            f64 duration = get_seconds(sp, ep);

            weight_t qap = get_qap(graphs[0], p_manager, d_oracle);
            MessageCostBreakdown cost = get_message_cost(q_graph, ac.k, d_oracle,
                                                         ac.alpha_end(), ac.beta_end(),
                                                         ac.alpha_link, ac.beta_link);
            weight_t lmax = std::ceil((1.0 + ac.imbalance) * ((f64) graphs[0].g_weight / (f64) ac.k));

            size_t n_empty_partitions = 0;
            size_t n_overloaded_partitions = 0;
            weight_t sum_too_much = 0;
            for (partition_t id = 0; id < ac.k; ++id) {
                n_empty_partitions += p_manager.get_bweight(id) == 0;
                n_overloaded_partitions += p_manager.get_bweight(id) > lmax;
                sum_too_much += std::max((weight_t) 0, p_manager.get_bweight(id) - lmax);
            }

            std::cout << "Graph                 : " << ac.graph_in << std::endl;
            std::cout << "Total time (s)        : " << duration << std::endl;
            std::cout << "#Nodes                : " << graphs[0].n << std::endl;
            std::cout << "#Edges                : " << graphs[0].m << std::endl;
            std::cout << "k                     : " << ac.k << std::endl;
            std::cout << "Lmax                  : " << lmax << std::endl;
            std::cout << "Threads               : " << ac.threads << std::endl;
            std::cout << "-----------------------" << std::endl;
            std::cout << "Final QAP             : " << qap << std::endl;
            cost.print(std::cout);
            std::cout << "max block w           : " << max(p_manager.get_bweights()) << std::endl;
            std::cout << "#empty partitions     : " << n_empty_partitions << std::endl;
            std::cout << "#oload partitions     : " << n_overloaded_partitions << std::endl;
            std::cout << "Sum oload weights     : " << sum_too_much << std::endl;
            std::cout << "------- Time -------" << std::endl;
            std::cout << "IO                    : " << io_ms << std::endl;
            std::cout << "Coarsening            : " << coarsening_ms << std::endl;
            std::cout << "Contraction           : " << contraction_ms << std::endl;
            std::cout << "Init. Part.           : " << initial_partitioning_ms << std::endl;
            std::cout << "Uncontraction         : " << uncontraction_ms << std::endl;
            std::cout << "Rebalance             : " << rebalance_ms << std::endl;
            std::cout << "Refinement            : " << refinement_ms << std::endl;
            std::cout << "  Label Prop.         : " << lp_refine_ms << std::endl;
            std::cout << "  Quotient Graph      : " << qg_refine_ms << std::endl;
            std::cout << "  Flow                : " << flow_refine_ms << std::endl;
            std::cout << "Misc                  : " << misc_ms << std::endl;
            std::cout << "ALL                   : " << io_ms + misc_ms + coarsening_ms + contraction_ms + initial_partitioning_ms + uncontraction_ms + rebalance_ms + refinement_ms << std::endl;

            return p;
        }

    private:
        void internal_solve() {
            u64 level = 0;

            update_lmax_array();

            while (graphs.back().n > ac.k * ac.initial_c) {
                coarsening(level);
                contraction();

                if (d_oracle.can_coarsen()) {
                    d_oracle.coarsen();
                    ac.k = d_oracle.current_k();
                    reinit_structures_for_new_k();
                }

                if (graphs.back().n == graphs[graphs.size() - 2].n) {
                    graphs.pop_back();
                    mappings.pop_back();
                    break;
                }

                level += 1;
            }

            initial_partitioning();
            rebalance(level);
            refinement();

            while (!mappings.empty()) {
                level -= 1;

                // If oracle was coarsened at this graph level, uncoarsen it first
                if (d_oracle.depth() == mappings.size()) {
                    uncoarsen_topology();
                }

                uncoarsening();
                rebalance(level);
                refinement();
            }
        }

        void initial_partitioning() {
            auto p = get_time_point();
            HEIPROMAP_PROFILE_SCOPE("partition", "partition", "partition");

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
                    AlignedArray<vertex_t> vertices;
                    vertices.initialize(graphs.back().n);
                    greedy_partition(graphs.back(), d_oracle, ac.imbalance, ac.seed + thread_id, local_p_managers[thread_id], vertices);
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
                AlignedArray<vertex_t> vertices;
                vertices.initialize(graphs.back().n);
                greedy_partition(graphs.back(), d_oracle, ac.imbalance, ac.seed, p_manager, vertices);
            }

            HEIPROMAP_PROFILE_SCOPE("partition", "misc", "initialize_datastructures");
            p_manager.reset_weights();
            bv_manager.reset();
            q_graph.initialize(ac.k);
            block_conn.initialize(graphs.back().n, graphs.back().m, ac.k);
            block_conn.reset_build();

            for (vertex_t u = 0; u < graphs.back().n; ++u) {
                block_conn.begin_vertex(graphs.back(), u);

                const partition_t u_id = p_manager[u];
                const weight_t u_w = graphs.back().v_weight(u);
                p_manager.set(u, u_w, u_id);

                for (size_t i = graphs.back().neighborhoods[u]; i < graphs.back().neighborhoods[u + 1]; ++i) {
                    const vertex_t v = graphs.back().edges_v[i];
                    const weight_t w = graphs.back().edge_w(i);
                    const partition_t v_id = p_manager[v];

                    block_conn.add_connection(u, v_id, w);

                    if (u_id != v_id) {
                        bv_manager.add(u, u_id);
                        if (u < v) {
                            q_graph.add_edge(u_id, v_id, w);
                        }
                    }
                }
            }

            initial_partitioning_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
        }

        void coarsening(const u64 level) {
            auto p = get_time_point();
            mappings.emplace_back();
            mappings.back().initialize(graphs.back().n);

            if (ac.coarsening_algorithm_id == COARSENING_ALG_GLOBAL_PATHS) {
                gpa_matcher.match(level, graphs.back(), p_manager, mappings.back(), ac.imbalance);
            } else if (ac.coarsening_algorithm_id == COARSENING_ALG_SIZE_CONSTRAINED_LP) {
                size_constrained_lp.cluster(level, graphs.back(), p_manager, mappings.back(), ac.imbalance, ac.threads);
            } else if (ac.coarsening_algorithm_id == COARSENING_ALG_HEAVY_EDGE) {
                heavy_edge_matching.match(graphs.back(), p_manager, mappings.back(), ac.imbalance, random_engine.get_u64(), ac.heavy_edge_matching_config);
            } else {
                std::cerr << "Coarsening algorithm " << coarsening_algorithm_to_string(ac.coarsening_algorithm_id) << " with id " << ac.coarsening_algorithm_id << " not known!" << std::endl;
                exit(EXIT_FAILURE);
            }

            coarsening_ms += get_milli_seconds(p, get_time_point());
        }

        void contraction() {
            auto p = get_time_point();
            graphs.emplace_back();
            graphs.back().contract(graphs[graphs.size() - 2], mappings.back(), ac.threads);
            p_manager.contract(mappings.back());

            contraction_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_pre_partitioning(graphs.back(), p_manager, ac.k));
        }

        void uncoarsen_topology() {
            auto p = get_time_point();

            // Build coarse->fine map for partition expansion before popping oracle level
            partition_t coarse_k = d_oracle.current_k();
            partition_t fine_k = d_oracle.fine_k();
            const Mapping &topo_map = d_oracle.top_mapping();

            std::vector<std::vector<partition_t>> coarse_to_fine(coarse_k);
            for (partition_t f = 0; f < fine_k; ++f) {
                coarse_to_fine[topo_map.get(f)].push_back(f);
            }

            d_oracle.uncoarsen();
            ac.k = d_oracle.current_k();

            // Expand partition: assign each graph vertex from coarse block to least-loaded fine block
            std::vector<partition_t> old_partition(graphs.back().n);
            for (vertex_t u = 0; u < graphs.back().n; ++u) {
                old_partition[u] = p_manager[u];
            }

            reinit_structures_for_new_k();
            p_manager.reset_weights();

            for (vertex_t u = 0; u < graphs.back().n; ++u) {
                const auto &fine_blocks = coarse_to_fine[old_partition[u]];
                weight_t u_weight = graphs.back().v_weight(u);

                partition_t best = fine_blocks[0];
                weight_t min_w = p_manager.get_bweight(best);
                for (size_t i = 1; i < fine_blocks.size(); ++i) {
                    if (p_manager.get_bweight(fine_blocks[i]) < min_w) {
                        min_w = p_manager.get_bweight(fine_blocks[i]);
                        best = fine_blocks[i];
                    }
                }
                p_manager.set(u, u_weight, best);
            }

            q_graph.compute_from_scratch(graphs.back(), p_manager);

            uncontraction_ms += get_milli_seconds(p, get_time_point());
        }

        void uncoarsening() {
            auto p = get_time_point();

            p_manager.uncontract(mappings.back());

            HEIPROMAP_PROFILE_SCOPE("uncontraction", "misc", "compute_from_scratch");
            const graph_t &g_uncontracted = graphs[graphs.size() - 2];

            if (ac.threads > 1) {
                bv_manager.parallel_reset(ac.threads);
                block_conn.parallel_initialize_offsets(g_uncontracted, ac.threads);

                std::vector<std::vector<std::vector<vertex_t>>> thread_boundaries(ac.threads, std::vector<std::vector<vertex_t>>(ac.k));

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
                            const weight_t w = g_uncontracted.edge_w(i);
                            const partition_t v_id = p_manager[v];

                            if (u_id == v_id) {
                                own_weight += w;
                            } else {
                                block_conn.add_connection(u, v_id, w);
                            }

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
                        const weight_t w = g_uncontracted.edge_w(i);
                        const partition_t v_id = p_manager[v];

                        if (u_id == v_id) {
                            own_weight += w;
                        } else {
                            block_conn.add_connection(u, v_id, w);
                        }

                        n_different += (u_id != v_id);
                    }
                    if (own_weight > 0) {
                        block_conn.add_connection(u, u_id, own_weight);
                    }
                    bv_manager.add_boundary_vertex_from_count(u, u_id, n_different);
                }
            }

            graphs.pop_back();
            mappings.pop_back();

            uncontraction_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
        }

        void rebalance(const u64 level) {
            auto p = get_time_point();

            if (level == 0) {
                rebalancer.template rebalance_last_layer<true, true>(graphs.back(), p_manager, bv_manager, q_graph, d_oracle, block_conn);
            } else {
                rebalancer.template rebalance<true, true>(graphs.back(), p_manager, bv_manager, q_graph, d_oracle, block_conn);
            }

            rebalance_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
        }

        void refinement() {
            auto p = get_time_point();

            if (ac.label_propagation_config.enabled) {
                auto sp_local = get_time_point();
                lp_refine.refine(graphs.back(), d_oracle, bv_manager, p_manager, q_graph, block_conn);
                lp_refine_ms += get_milli_seconds(sp_local, get_time_point());
            }

            if (ac.quotient_graph_refinement_config.enabled) {
                auto sp_local = get_time_point();
                qg_refine.refine(graphs.back(), d_oracle, bv_manager, p_manager, q_graph, block_conn);
                qg_refine_ms += get_milli_seconds(sp_local, get_time_point());
            }

            if (ac.flow_based_refinement_config.enabled) {
                auto sp_local = get_time_point();
                flow_based_refinement.refine(graphs.back(), d_oracle, bv_manager, p_manager, q_graph, block_conn);
                flow_refine_ms += get_milli_seconds(sp_local, get_time_point());
            }

            refinement_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
        }
    };
}

#endif //HEIPROMAP_WAVER_MAP_SOLVER_H
