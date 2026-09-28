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
#include "../distance_oracles/distance_oracle.h"
#include "../distance_oracles/binary_distance_oracle.h"
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
#include "../refinement/flow_based_refinement.h"
#include "HeiPa_solver.h"

namespace HeiProMap {
    /**
     * Solver for Deep Process Mapping.
     */
    template<typename DistanceOracleT = DistanceOracle, typename QuotientGraphT = LargeQuotientGraph>
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
        LabelPropagationRefinement lp_refine;
        QuotientGraphRefinement qg_refine;
        FlowBasedRefinement flow_based_refinement;

        bool structs_up_to_date = true;

        f64 misc_ms = 0.0;
        f64 coarsening_ms = 0.0;
        f64 contraction_ms = 0.0;
        f64 initial_partitioning_ms = 0.0;
        f64 intermediate_partitioning_ms = 0.0;
        f64 uncontraction_ms = 0.0;
        f64 rebalance_ms = 0.0;
        f64 refinement_ms = 0.0;

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

            rebalancer.initialize(graphs[0].n, graphs[0].m, ac.k, ac.seed);
            bv_manager.initialize(graphs[0].n, ac.k);
            q_graph.initialize(ac.k);
            block_conn.initialize(graphs[0].n, graphs[0].m, ac.k);
            d_oracle.initialize(ac.hierarchy, ac.distance);

            HEAVYASSERT(assert_state_pre_partitioning(graphs[0], p_manager, ac.k));

            // matching
            gpa_matcher.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine, ac.global_path_algorithm_config);
            size_constrained_lp_clustering.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, ac.size_constrained_lp_clustering_configuration);

            // refinement
            LabelPropagationConfiguration label_propagation_config("label propy");
            lp_refine.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine.get_u64(), ac.deep_label_propagation_refinement_config);
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

            rebalancer.initialize(graphs[0].n, graphs[0].m, ac.k, ac.seed);
            bv_manager.initialize(graphs[0].n, ac.k);
            q_graph.initialize(ac.k);
            block_conn.initialize(graphs[0].n, graphs[0].m, ac.k);
            d_oracle.initialize(ac.hierarchy, ac.distance);

            HEAVYASSERT(assert_state_pre_partitioning(graphs[0], p_manager, ac.k));

            // matching
            gpa_matcher.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine, ac.global_path_algorithm_config);
            size_constrained_lp_clustering.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, ac.size_constrained_lp_clustering_configuration);

            // refinement
            if (ac.deep_label_propagation_refinement_config.enabled) {
                lp_refine.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine.get_u64(), ac.deep_label_propagation_refinement_config);
            }
        }

        std::vector<partition_t> solve() {
            internal_solve();

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));

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
            std::cout << "Distance Oracle       : " << (std::is_same_v<DistanceOracleT, BinaryDistanceOracle> ? "binary" : "matrix") << std::endl;
            std::cout << "Lmax                  : " << lmax << std::endl;
            std::cout << "Threads               : " << ac.threads << std::endl;
            std::cout << "-----------------------" << std::endl;
            std::cout << "Final QAP             : " << qap << std::endl;
            std::cout << "max block w           : " << max(p_manager.get_bweights()) << std::endl;
            std::cout << "#empty partitions     : " << n_empty_partitions << std::endl;
            std::cout << "#oload partitions     : " << n_overloaded_partitions << std::endl;
            std::cout << "------- Time -------" << std::endl;
            std::cout << "Total solve time      : " << duration * 1000.0 << std::endl;
            std::cout << "Coarsening            : " << coarsening_ms << std::endl;
            std::cout << "Contraction           : " << contraction_ms << std::endl;
            std::cout << "Init. Part.           : " << initial_partitioning_ms << std::endl;
            std::cout << "Inter. Part.          : " << intermediate_partitioning_ms << std::endl;
            std::cout << "Uncontraction         : " << uncontraction_ms << std::endl;
            std::cout << "Rebalance             : " << rebalance_ms << std::endl;
            std::cout << "Refinement            : " << refinement_ms << std::endl;
            std::cout << "Misc                  : " << misc_ms << std::endl;
            std::cout << "ALL                   : " << coarsening_ms + contraction_ms + initial_partitioning_ms + intermediate_partitioning_ms + uncontraction_ms + rebalance_ms + refinement_ms + misc_ms << std::endl;

            return p;
        }

    private:
        void internal_solve() {
            u64 level = 0;
            [[maybe_unused]] u64 max_level = 0;

            size_t l = ac.hierarchy.size();

            bool stagnated = false;
            for (u64 h_level = 0; h_level < l && !stagnated; ++h_level) {
                while (graphs.back().n > k_rem[l - h_level] * ac.initial_C) {
                    coarsening(level, lmax_vec[h_level]);
                    contraction();

                    if (graphs.back().n == graphs[graphs.size() - 2].n) {
                        graphs.pop_back();
                        mappings.pop_back();
                        stagnated = true;
                        break;
                    }

                    level += 1;
                }
            }

            max_level = level > 0 ? level - 1 : 0;

            initial_partitioning();

            rebalance(level);

            while (level > 0) {
                level -= 1;

                uncoarsening();

                intermediate_partitioning(level);

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

            PartitionManager sub_pm;
            sub_pm.initialize(g_copy.n, block_k, g_copy.g_weight);
            UniformDistanceOracle temp_d_oracle(block_k);
            greedy_partition(g_copy, temp_d_oracle, imb, 42, sub_pm);

            for (vertex_t u = 0; u < graphs.back().n; ++u) {
                partition_t sub_p = sub_pm[u];
                if (sub_p == 0) continue;
                partition_t move_id = 0 + k_spacing * sub_p;
                weight_t w = graphs.back().v_weights[u];
                p_manager.move(u, w, 0, move_id);
                block_conn.move(graphs.back(), u, 0, move_id);
            }

            for (partition_t i = 0; i < block_k; ++i) {
                partition_t move_id = 0 + k_spacing * i;
                p_manager.set_lmax(move_id, child_lmax);
                p_manager.set_hierarchy_level(move_id, block_level - 1);
            }

            bv_manager.compute_from_scratch(graphs.back(), p_manager);
            q_graph.compute_from_scratch(graphs.back(), p_manager);
            block_conn.compute_from_scratch(graphs.back(), p_manager);

            initial_partitioning_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
        }

        void intermediate_partitioning(const u64 level) {
            auto p = get_time_point();
            HEIPROMAP_PROFILE_SCOPE("intermediate_partitioning", "misc", "intermediate_partitioning");

            std::vector<partition_t> ids = {0};
            std::vector<u32> one_hot(ac.k);

            while (!ids.empty()) {
                ids.clear();
                std::fill(one_hot.begin(), one_hot.end(), 0);

                for (partition_t id = 0; id < ac.k; ++id) {
                    if (!p_manager.is_active(id)) continue;
                    if (p_manager.get_hierarchy_level(id) == 0) continue;

                    vertex_t threshold = ac.intermediate_C * ac.hierarchy[p_manager.get_hierarchy_level(id) - 1];
                    if (p_manager.size(id) >= threshold || level == 0) {
                        ids.push_back(id);
                        one_hot[id] = 1;
                    }
                }

                if (ids.empty()) break;

                structs_up_to_date = false;

                HEIPROMAP_PROFILE_SCOPE("intermediate_partitioning", "subgraph_extractor", "extract");
                subgraph_extractor.extract(graphs.back(), p_manager, one_hot);

                HEIPROMAP_PROFILE_SCOPE("intermediate_partitioning", "greedy_partitioner", "sub_block_partitioning");
                #pragma omp parallel for num_threads(ac.threads) schedule(dynamic)
                for (size_t i = 0; i < ids.size(); ++i) {
                    partition_t block_id = ids[i];
                    partition_t block_level = p_manager.get_hierarchy_level(block_id);
                    partition_t block_k = ac.hierarchy[block_level - 1];
                    partition_t k_spacing = k_rem[block_level - 1];
                    weight_t child_lmax = (block_level >= 2) ? lmax_vec[block_level - 1] : lmax_vec[0];

                    CSRGraph &sub_g = subgraph_extractor.graphs[block_id];
                    auto &tt = subgraph_extractor.tts[block_id];

                    if (sub_g.n > 0) {
                        f64 imb = ((f64) (child_lmax * block_k) / (f64) sub_g.g_weight) - 1.0;
                        if (imb <= 0.0) imb = 0.0001;

                        CSRGraph sub_g_copy = sub_g;

                        PartitionManager sub_pm;
                        sub_pm.initialize(sub_g_copy.n, block_k, sub_g_copy.g_weight);
                        UniformDistanceOracle temp_d_oracle(block_k);
                        greedy_partition(sub_g_copy, temp_d_oracle, imb, 42 + block_id, sub_pm);

                        for (vertex_t u = 0; u < sub_g.n; ++u) {
                            partition_t sub_p = sub_pm[u];
                            if (sub_p == 0) continue;
                            partition_t move_id = block_id + k_spacing * sub_p;
                            vertex_t orig_u = tt.get_o(u);
                            weight_t w = sub_g.v_weights[u];
                            p_manager.move(orig_u, w, block_id, move_id);
                        }
                    }

                    for (partition_t p = 0; p < block_k; ++p) {
                        partition_t move_id = block_id + k_spacing * p;
                        p_manager.set_lmax(move_id, child_lmax);
                        p_manager.set_hierarchy_level(move_id, block_level - 1);
                    }
                }

                HEIPROMAP_PROFILE_SCOPE("intermediate_partitioning", "datastructures", "compute_from_scratch");
                bv_manager.compute_from_scratch(graphs.back(), p_manager);
                q_graph.compute_from_scratch(graphs.back(), p_manager);
                block_conn.compute_from_scratch(graphs.back(), p_manager);
                structs_up_to_date = true;
            }

            if (!structs_up_to_date) {
                HEIPROMAP_PROFILE_SCOPE("intermediate_partitioning", "datastructures", "compute_from_scratch");
                bv_manager.compute_from_scratch(graphs.back(), p_manager);
                q_graph.compute_from_scratch(graphs.back(), p_manager);
                block_conn.compute_from_scratch(graphs.back(), p_manager);
                structs_up_to_date = true;
            }

            intermediate_partitioning_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
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
            graphs.back().contract(graphs[graphs.size() - 2], mappings.back(), ac.threads);
            p_manager.contract(mappings.back());

            contraction_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_pre_partitioning(graphs.back(), p_manager, ac.k));
        }

        void uncoarsening() {
            auto p = get_time_point();
            p_manager.uncontract(mappings.back());
            mappings.pop_back();
            graphs.pop_back();

            structs_up_to_date = false;
            uncontraction_ms += get_milli_seconds(p, get_time_point());
        }

        void refinement([[maybe_unused]] const u64 level, [[maybe_unused]] const u64 max_level) {
            auto p = get_time_point();
            if (ac.deep_label_propagation_refinement_config.enabled) {
                lp_refine.refine(graphs.back(), d_oracle, bv_manager, p_manager, q_graph, block_conn);
            }
            refinement_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
        }

        void rebalance(const u64 level) {
            auto p = get_time_point();
            if (!p_manager.is_overloaded()) {
                misc_ms += get_milli_seconds(p, get_time_point());
                return;
            }

            if (level == 0) {
                rebalancer.rebalance_last_layer(graphs.back(), p_manager, bv_manager, q_graph, d_oracle, block_conn, ac.imbalance);
            } else {
                rebalancer.rebalance(graphs.back(), p_manager, bv_manager, q_graph, d_oracle, block_conn, ac.imbalance);
            }

            rebalance_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
        }
    };
}

#endif //HEIPROMAP_DEEP_HEIPROMAP_SOLVER_H
