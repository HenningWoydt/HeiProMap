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

#ifndef HEIPROMAP_DEEP_HEIPROMAP_SOLVER_H
#define HEIPROMAP_DEEP_HEIPROMAP_SOLVER_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>
#include <omp.h>

#include "../datastructures/block_conn.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../distance_oracles/matrix_distance_oracle.h"
#include "../distance_oracles/binary_distance_oracle.h"
#include "../distance_oracles/division_distance_oracle.h"
#include "../distance_oracles/stored_division_distance_oracle.h"
#include "../distance_oracles/verify_distance_oracles.h"
#include "../datastructures/partition_manager.h"
#include "../datastructures/quotient_graph.h"
#include "../datastructures/large_quotient_graph.h"
#include "../datastructures/subgraph_extractor.h"

#include "../definitions.h"
#include "../configuration/DeepHeiProMap_configuration.h"
#include "../utility/macros.h"
#include "../utility/mapping.h"
#include "../utility/random_engine.h"
#include "../utility/assert_state.h"
#include "../utility/small_translation_table.h"
#include "../utility/profiler.h"
#include "../utility/qap.h"
#include "../utility/utils.h"

#include "../rebalance/rebalancer.h"
#include "../coarsening/heavy_edge_matching.h"
#include "../coarsening/global_path_algorithm.h"
#include "../coarsening/size_constrained_lp.h"
#include "../refinement/quotient_graph_refinement.h"
#include "../refinement/simple_quotient_graph_refinement.h"
#include "../refinement/flow_based_refinement.h"
#include "../refinement/simple_label_propagation_refinement.h"
#include "HeiPa_solver.h"

namespace HeiProMap {
    /**
     * Solver for Deep Process Mapping.
     */
    template<bool LARGE_K, typename DistanceOracleT, typename QuotientGraphT>
    class DeepHeiProMapSolver {
        DeepHeiProMapConfiguration ac;
        RandomEngine random_engine;

        // statistics
        std::chrono::high_resolution_clock::time_point sp;

        std::vector<graph_t> graphs;
        p_manager_t p_manager;
        DistanceOracleT d_oracle;
        bv_manager_t bv_manager;
        QuotientGraphT q_graph;
        block_conn_t block_conn;
        SubgraphExtractor subgraph_extractor;

        Rebalancer rebalancer;

        // balance
        weight_t lmax = 0;
        std::vector<weight_t> lmax_vec;
        std::vector<partition_t> k_rem;

        // matching
        std::vector<Mapping> mappings;
        GlobalPathAlgorithmMatcher gpa_matcher;
        HeavyEdgeMatching heavy_edge_matcher;
        SizeConstrainedLP size_constrained_lp_clustering;

        // refinement
        SimpleLabelPropagationRefinement<LARGE_K> lp_refine;
        QuotientGraphRefinement<LARGE_K> qg_refine;
        SimpleQuotientGraphRefinement<LARGE_K> simple_qg_refine;
        FlowBasedRefinement<LARGE_K> flow_based_refinement;

        std::vector<partition_t> inter_ids;
        std::vector<partition_t> inter_id_to_dense;
        std::vector<PartitionManager> thread_sub_pm;

        std::vector<std::vector<std::tuple<partition_t, partition_t, weight_t> > > thread_edges;

        f64 misc_ms = 0.0;
        f64 coarsening_ms = 0.0;
        f64 contraction_ms = 0.0;
        f64 initial_partitioning_ms = 0.0;
        f64 intermediate_partitioning_ms = 0.0;
        f64 recompute_datastructures_ms = 0.0;
        f64 uncontraction_ms = 0.0;
        f64 rebalance_ms = 0.0;
        f64 refinement_ms = 0.0;
        f64 lp_refine_ms = 0.0;
        f64 qg_refine_ms = 0.0;
        f64 flow_refine_ms = 0.0;

    public:
        explicit DeepHeiProMapSolver(const DeepHeiProMapConfiguration &t_ac) {
            sp = std::chrono::high_resolution_clock::now();

            ac = t_ac;
            random_engine = RandomEngine(ac.seed);

            graphs.emplace_back(ac.graph_in);

            // balance
            lmax = std::ceil((1.0 + ac.imbalance) * ((f64) graphs[0].g_weight / (f64) ac.k));
            lmax_vec.resize(ac.hierarchy.size());
            lmax_vec[0] = lmax;
            for (u64 i = 1; i < ac.hierarchy.size(); ++i) {
                lmax_vec[i] = lmax_vec[i - 1] * (weight_t) ac.hierarchy[i - 1];
            }

            partition_t temp_k = 1;
            k_rem.push_back(temp_k);
            for (u64 i = 0; i < ac.hierarchy.size(); ++i) {
                temp_k *= ac.hierarchy[ac.hierarchy.size() - 1 - i];
                k_rem.push_back(k_rem.back() * ac.hierarchy[i]);
            }

            // manager
            p_manager.initialize(graphs[0].n, ac.k, graphs[0].g_weight);
            p_manager.set_hierarchy_level(0, ac.hierarchy.size());
            p_manager.set_lmax(0, lmax_vec.back());

            subgraph_extractor.initialize(ac.k, ac.threads);
            inter_id_to_dense.assign(ac.k, std::numeric_limits<partition_t>::max());
            thread_sub_pm.resize(ac.threads);
            thread_edges.resize(ac.threads);

            rebalancer.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, ac.seed);
            bv_manager.initialize(graphs[0].n, ac.k);
            if (ac.enable_q_graph) { q_graph.initialize(ac.k); }
            if (ac.enable_block_conn) { block_conn.initialize(graphs[0].n, graphs[0].m, ac.k); }
            d_oracle.initialize(ac.hierarchy, ac.distance);
            HEAVYASSERT(verify_distance_oracles());

            HEAVYASSERT(assert_state_pre_partitioning(graphs[0], p_manager, ac.k, ac.threads));

            // matching
            gpa_matcher.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine, ac.global_path_algorithm_config);
            size_constrained_lp_clustering.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, ac.size_constrained_lp_clustering_configuration);

            // refinement
            if (ac.deep_label_propagation_refinement_config.enabled) {
                lp_refine.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine.get_u64(), ac.deep_label_propagation_refinement_config);
            }
            if (ac.deep_quotient_graph_refinement_config.enabled) {
                qg_refine.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine.get_u64(), ac.deep_quotient_graph_refinement_config);
            }
            if (ac.deep_simple_qg_refinement_config.enabled) {
                simple_qg_refine.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine.get_u64(), ac.deep_simple_qg_refinement_config);
            }
            if (ac.deep_flow_based_refinement_config.enabled) {
                flow_based_refinement.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine.get_u64(), ac.deep_flow_based_refinement_config);
            }
        }

        explicit DeepHeiProMapSolver(graph_t &&g, const DeepHeiProMapConfiguration &t_ac) {
            sp = std::chrono::high_resolution_clock::now();

            ac = t_ac;
            random_engine = RandomEngine(ac.seed);

            graphs.emplace_back(std::move(g));

            // balance
            lmax = std::ceil((1.0 + ac.imbalance) * ((f64) graphs[0].g_weight / (f64) ac.k));
            lmax_vec.resize(ac.hierarchy.size());
            lmax_vec[0] = lmax;
            for (u64 i = 1; i < ac.hierarchy.size(); ++i) {
                lmax_vec[i] = lmax_vec[i - 1] * (weight_t) ac.hierarchy[i - 1];
            }

            partition_t temp_k = 1;
            k_rem.push_back(temp_k);
            for (u64 i = 0; i < ac.hierarchy.size(); ++i) {
                temp_k *= ac.hierarchy[ac.hierarchy.size() - 1 - i];
                k_rem.push_back(k_rem.back() * ac.hierarchy[i]);
            }

            // manager
            p_manager.initialize(graphs[0].n, ac.k, graphs[0].g_weight);
            p_manager.set_hierarchy_level(0, ac.hierarchy.size());
            p_manager.set_lmax(0, lmax_vec.back());

            subgraph_extractor.initialize(ac.k, ac.threads);
            inter_id_to_dense.assign(ac.k, std::numeric_limits<partition_t>::max());
            thread_sub_pm.resize(ac.threads);
            thread_edges.resize(ac.threads);

            rebalancer.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, ac.seed);
            bv_manager.initialize(graphs[0].n, ac.k);
            if (ac.enable_q_graph) { q_graph.initialize(ac.k); }
            if (ac.enable_block_conn) { block_conn.initialize(graphs[0].n, graphs[0].m, ac.k); }
            d_oracle.initialize(ac.hierarchy, ac.distance);
            HEAVYASSERT(verify_distance_oracles());

            HEAVYASSERT(assert_state_pre_partitioning(graphs[0], p_manager, ac.k, ac.threads));

            // matching
            gpa_matcher.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine, ac.global_path_algorithm_config);
            size_constrained_lp_clustering.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, ac.size_constrained_lp_clustering_configuration);

            // refinement
            if (ac.deep_label_propagation_refinement_config.enabled) {
                lp_refine.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine.get_u64(), ac.deep_label_propagation_refinement_config);
            }
            if (ac.deep_quotient_graph_refinement_config.enabled) {
                qg_refine.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine.get_u64(), ac.deep_quotient_graph_refinement_config);
            }
            if (ac.deep_simple_qg_refinement_config.enabled) {
                simple_qg_refine.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine.get_u64(), ac.deep_simple_qg_refinement_config);
            }
            if (ac.deep_flow_based_refinement_config.enabled) {
                flow_based_refinement.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine.get_u64(), ac.deep_flow_based_refinement_config);
            }
        }

        std::vector<partition_t> solve() {
            internal_solve();

            std::vector<partition_t> p(graphs.back().n);
            for (vertex_t u = 0; u < graphs.back().n; ++u) { p[u] = p_manager[u]; }

            const auto ep = std::chrono::high_resolution_clock::now();
            f64 duration = get_seconds(sp, ep);

            weight_t qap = get_qap(graphs.back(), p_manager, d_oracle);
            size_t n_empty_partitions = 0;
            size_t n_overloaded_partitions = 0;
            weight_t sum_too_much = 0;
            for (partition_t id = 0; id < ac.k; ++id) {
                n_empty_partitions += p_manager.get_bweight(id) == 0;
                n_overloaded_partitions += p_manager.get_bweight(id) > p_manager.lmax[id];
                sum_too_much += std::max((weight_t) 0, p_manager.get_bweight(id) - p_manager.lmax[id]);
            }

            std::cout << "Graph                 : " << ac.graph_in << std::endl;
            std::cout << "Total time (s)        : " << duration << std::endl;
            std::cout << "#Nodes                : " << graphs[0].n << std::endl;
            std::cout << "#Edges                : " << graphs[0].m << std::endl;
            std::cout << "k                     : " << ac.k << std::endl;
            std::cout << "Hierarchy             : " << ac.hierarchy_string << std::endl;
            std::cout << "Distances             : " << ac.distance_string << std::endl;
            std::cout << "Distance Oracle       : ";
            if constexpr (std::is_same_v<DistanceOracleT, BinaryDistanceOracle>) std::cout << "binary";
            else if constexpr (std::is_same_v<DistanceOracleT, DivisionDistanceOracle>) std::cout << "division";
            else if constexpr (std::is_same_v<DistanceOracleT, StoredDivisionDistanceOracle>) std::cout << "stored-division";
            else std::cout << "matrix";
            std::cout << std::endl;
            std::cout << "Lmax                  : " << lmax << std::endl;
            std::cout << "Threads               : " << ac.threads << std::endl;
            std::cout << "Coarsening Mode       : " << (ac.hierarchical_coarsening ? "hierarchical" : "flat") << std::endl;
            std::cout << "-----------------------" << std::endl;
            std::cout << "Final QAP             : " << qap << std::endl;
            std::cout << "max block w           : " << max(p_manager.get_bweights()) << std::endl;
            std::cout << "#empty partitions     : " << n_empty_partitions << std::endl;
            std::cout << "#oload partitions     : " << n_overloaded_partitions << std::endl;
            std::cout << "Sum oload weights     : " << sum_too_much << std::endl;
            std::cout << "------- Time -------" << std::endl;
            std::cout << "Total solve time      : " << duration * 1000.0 << std::endl;
            std::cout << "Coarsening            : " << coarsening_ms << std::endl;
            std::cout << "Contraction           : " << contraction_ms << std::endl;
            std::cout << "Init. Part.           : " << initial_partitioning_ms << std::endl;
            std::cout << "Inter. Part.          : " << intermediate_partitioning_ms << std::endl;
            std::cout << "Recompute             : " << recompute_datastructures_ms << std::endl;
            std::cout << "Uncontraction         : " << uncontraction_ms << std::endl;
            std::cout << "Rebalance             : " << rebalance_ms << std::endl;
            std::cout << "Refinement            : " << refinement_ms << std::endl;
            std::cout << "  Label Prop.         : " << lp_refine_ms << std::endl;
            std::cout << "  Quotient Graph      : " << qg_refine_ms << std::endl;
            std::cout << "  Flow                : " << flow_refine_ms << std::endl;
            std::cout << "Misc                  : " << misc_ms << std::endl;
            std::cout << "ALL                   : " << coarsening_ms + contraction_ms + initial_partitioning_ms + intermediate_partitioning_ms + recompute_datastructures_ms + uncontraction_ms + rebalance_ms + refinement_ms + misc_ms << std::endl;

            return p;
        }

    private:
        void internal_solve() {
            u64 level = 0;
            [[maybe_unused]] u64 max_level = 0;

            if (ac.hierarchical_coarsening) {
                const size_t l = ac.hierarchy.size();
                for (u64 h_level = 0; h_level < l; ++h_level) {
                    while (graphs.back().n > k_rem[l - h_level] * ac.initial_C) {
                        coarsening(level, lmax_vec[h_level]);
                        contraction();

                        if (graphs.back().n == graphs[graphs.size() - 2].n) {
                            graphs.pop_back();
                            mappings.pop_back();
                            break;
                        }

                        level += 1;
                    }
                }
            } else {
                const partition_t z = ac.hierarchy.back();
                const weight_t lmax_z = lmax * (ac.k / z);
                const vertex_t target_n = ac.initial_C * z;

                while (graphs.back().n > target_n) {
                    coarsening(level, lmax_z);
                    contraction();

                    if (graphs.back().n == graphs[graphs.size() - 2].n) {
                        graphs.pop_back();
                        mappings.pop_back();
                        break;
                    }

                    level += 1;
                }
            }

            max_level = level > 0 ? level - 1 : 0;

            initial_partitioning();

            while (level > 0) {
                level -= 1;

                uncoarsening();

                intermediate_partitioning(level);

                recompute_datastructures();

                rebalance(level);

                refinement(level, max_level);
            }
        }

        void initial_partitioning() {
            auto p = get_time_point();
            HEIPROMAP_PROFILE_SCOPE("initial_partition", "misc", "initial_partitioning");

            partition_t block_k = ac.hierarchy.back();
            partition_t block_level = p_manager.get_hierarchy_level(0);
            partition_t k_spacing = k_rem[block_level - 1];
            weight_t child_lmax = lmax_vec[block_level - 1];

            f64 imb = ((f64) (child_lmax * block_k) / (f64) graphs.back().g_weight) - 1.0;
            if (imb <= 0.0) imb = 0.0001;

            graph_t g_copy = graphs.back();

            PartitionManager sub_pm(g_copy.n, block_k, g_copy.g_weight);
            UniformDistanceOracle temp_d_oracle(block_k);
            greedy_partition(g_copy, temp_d_oracle, imb, 42, sub_pm);

            for (vertex_t u = 0; u < graphs.back().n; ++u) {
                partition_t sub_p = sub_pm[u];
                if (sub_p == 0) continue;
                partition_t move_id = 0 + k_spacing * sub_p;
                weight_t w = graphs.back().v_weights[u];
                p_manager.move(u, w, 0, move_id);
            }

            for (partition_t i = 0; i < block_k; ++i) {
                partition_t move_id = 0 + k_spacing * i;
                p_manager.set_lmax(move_id, child_lmax);
                p_manager.set_hierarchy_level(move_id, block_level - 1);
            }

            initial_partitioning_ms += get_milli_seconds(p, get_time_point());
        }

        void intermediate_partitioning(const u64 level) {
            HEIPROMAP_PROFILE_SCOPE("intermediate_partitioning", "misc", "intermediate_partitioning");

            auto p = get_time_point();

            while (true) {
                HEIPROMAP_PROFILE_SCOPE("intermediate_partitioning", "misc", "collect");
                inter_ids.clear();
                for (partition_t id = 0; id < ac.k; ++id) {
                    if (!p_manager.is_active(id)) continue;
                    if (p_manager.get_hierarchy_level(id) == 0) continue;

                    vertex_t threshold = ac.intermediate_C * ac.hierarchy[p_manager.get_hierarchy_level(id) - 1];
                    if (p_manager.size(id) >= threshold || level == 0) {
                        inter_id_to_dense[id] = inter_ids.size();
                        inter_ids.push_back(id);
                    }
                }

                if (inter_ids.empty()) break;

                HEIPROMAP_PROFILE_SCOPE("intermediate_partitioning", "subgraph_extractor", "extract");
                subgraph_extractor.extract(graphs.back(), p_manager, inter_ids, inter_id_to_dense);

                HEIPROMAP_PROFILE_SCOPE("intermediate_partitioning", "misc", "reset");
                // Reset only used entries in inter_id_to_dense
                for (partition_t id: inter_ids) {
                    inter_id_to_dense[id] = std::numeric_limits<partition_t>::max();
                }

                HEIPROMAP_PROFILE_SCOPE("intermediate_partitioning", "greedy_partitioner", "sub_block_partitioning");
                #pragma omp parallel for num_threads(ac.threads) schedule(dynamic)
                for (size_t i = 0; i < inter_ids.size(); ++i) {
                    partition_t block_id = inter_ids[i];
                    partition_t block_level = p_manager.get_hierarchy_level(block_id);
                    partition_t block_k = ac.hierarchy[block_level - 1];
                    partition_t k_spacing = k_rem[block_level - 1];
                    weight_t child_lmax = (block_level >= 2) ? lmax_vec[block_level - 1] : lmax_vec[0];

                    CSRGraph &sub_g = subgraph_extractor.graphs[i];
                    auto &tt = subgraph_extractor.tts[i];

                    if (sub_g.n > 0) {
                        f64 imb = ((f64) (child_lmax * block_k) / (f64) sub_g.g_weight) - 1.0;
                        if (imb <= 0.0) imb = 0.0001;

                        u64 tid = omp_get_thread_num();
                        PartitionManager &sub_pm = thread_sub_pm[tid];
                        sub_pm.initialize(sub_g.n, block_k, sub_g.g_weight);
                        UniformDistanceOracle temp_d_oracle(block_k);
                        greedy_partition(sub_g, temp_d_oracle, imb, 42 + block_id, sub_pm);

                        for (vertex_t u = 0; u < sub_g.n; ++u) {
                            partition_t sub_p = sub_pm[u];
                            if (sub_p == 0) continue;
                            partition_t move_id = block_id + k_spacing * sub_p;
                            vertex_t orig_u = tt.get_o(u);
                            weight_t w = sub_g.v_weights[u];
                            p_manager.move(orig_u, w, block_id, move_id);
                        }
                    }

                    for (partition_t p_id = 0; p_id < block_k; ++p_id) {
                        partition_t move_id = block_id + k_spacing * p_id;
                        p_manager.set_lmax(move_id, child_lmax);
                        p_manager.set_hierarchy_level(move_id, block_level - 1);
                    }
                }
            }

            intermediate_partitioning_ms += get_milli_seconds(p, get_time_point());
        }

        void recompute_datastructures() {
            HEIPROMAP_PROFILE_SCOPE("recompute_datastructures", "datastructures", "all");
            auto p = get_time_point();

            const graph_t &g = graphs.back();

            // Phase 1: Recompute block connections
            if (ac.enable_block_conn) {
                block_conn.compute_from_scratch(g, p_manager, ac.threads);
            }

            // Phase 2: Recompute boundary vertices
            bv_manager.compute_from_scratch(g, p_manager, ac.threads);

            // Phase 3: Recompute quotient graph using boundary vertices from Phase 2
            if (ac.enable_q_graph) {
                q_graph.compute_from_scratch(g, p_manager, bv_manager, ac.threads);
            }

            recompute_datastructures_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k, ac.threads, ac.enable_q_graph, ac.enable_block_conn));
        }

        void coarsening(const u64 level, const weight_t max_v_weight) {
            auto p = get_time_point();
            mappings.emplace_back();
            mappings.back().initialize(graphs.back().n);

            if (ac.coarsening_algorithm_id == COARSENING_ALG_HEAVY_MATCHING) {
                heavy_edge_matcher.match(graphs.back(), p_manager, mappings.back(), ac.imbalance, random_engine.get_u64(), ac.parallel_heavy_edge_matching_configuration, max_v_weight);
            } else if (ac.coarsening_algorithm_id == COARSENING_ALG_GLOBAL_PATHS) {
                gpa_matcher.match(level, graphs.back(), p_manager, mappings.back(), ac.imbalance, max_v_weight);
            } else if (ac.coarsening_algorithm_id == COARSENING_ALG_SIZE_CONSTRAINED_LP) {
                size_constrained_lp_clustering.cluster(level, graphs.back(), p_manager, mappings.back(), ac.imbalance, ac.threads, max_v_weight);
            } else {
                std::cerr << "Coarsening Algorithm not recognized : " << ac.coarsening_algorithm_id << std::endl;
                abort();
            }

            coarsening_ms += get_milli_seconds(p, get_time_point());
        }

        void contraction() {
            auto p = get_time_point();
            graphs.emplace_back(); // coarse the graph
            graphs.back().contract(graphs[graphs.size() - 2], mappings.back(), ac.threads, ac.threads > 1 || ac.force_parallel_contraction, ac.use_kaminpar_contraction);
            p_manager.contract(mappings.back());

            contraction_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_pre_partitioning(graphs.back(), p_manager, ac.k, ac.threads));
        }

        void uncoarsening() {
            auto p = get_time_point();
            p_manager.uncontract(mappings.back());
            mappings.pop_back();
            graphs.pop_back();

            uncontraction_ms += get_milli_seconds(p, get_time_point());
        }

        void refinement([[maybe_unused]] const u64 level, [[maybe_unused]] const u64 max_level) {
            auto p = get_time_point();
            if (ac.deep_label_propagation_refinement_config.enabled) {
                auto sp_local = get_time_point();
                lp_refine.refine(graphs.back(), d_oracle, bv_manager, p_manager, q_graph, block_conn);
                lp_refine_ms += get_milli_seconds(sp_local, get_time_point());
            }

            if (ac.deep_quotient_graph_refinement_config.enabled) {
                auto sp_local = get_time_point();
                qg_refine.refine(graphs.back(), d_oracle, bv_manager, p_manager, q_graph, block_conn);
                qg_refine_ms += get_milli_seconds(sp_local, get_time_point());
                HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k, ac.threads, ac.enable_q_graph, ac.enable_block_conn));
            }

            if (ac.deep_simple_qg_refinement_config.enabled) {
                auto sp_local = get_time_point();
                simple_qg_refine.refine(graphs.back(), d_oracle, bv_manager, p_manager);
                qg_refine_ms += get_milli_seconds(sp_local, get_time_point());
            }

            if (ac.deep_flow_based_refinement_config.enabled) {
                auto sp_local = get_time_point();
                flow_based_refinement.refine(graphs.back(), d_oracle, bv_manager, p_manager, q_graph, block_conn);
                flow_refine_ms += get_milli_seconds(sp_local, get_time_point());
                HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k, ac.threads, ac.enable_q_graph, ac.enable_block_conn));
            }

            refinement_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k, ac.threads, ac.enable_q_graph, ac.enable_block_conn));
        }

        void rebalance(const u64 level) {
            auto p = get_time_point();
            if (!p_manager.is_overloaded()) {
                misc_ms += get_milli_seconds(p, get_time_point());
                return;
            }

            if (level == 0) {
                rebalancer.rebalance_last_layer(graphs.back(), p_manager, bv_manager, q_graph, d_oracle, block_conn, ac.imbalance, ac.enable_q_graph, ac.enable_block_conn);
            } else {
                rebalancer.rebalance(graphs.back(), p_manager, bv_manager, q_graph, d_oracle, block_conn, ac.imbalance, ac.enable_q_graph, ac.enable_block_conn);
            }

            rebalance_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k, ac.threads, ac.enable_q_graph, ac.enable_block_conn));
        }
    };
}

#endif //HEIPROMAP_DEEP_HEIPROMAP_SOLVER_H
