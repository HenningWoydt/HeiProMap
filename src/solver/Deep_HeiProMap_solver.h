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
#include "../datastructures/distance_oracle.h"
#include "../datastructures/binary_distance_oracle.h"
#include "../datastructures/partition_manager.h"
#include "../datastructures/quotient_graph.h"
#include "../datastructures/large_quotient_graph.h"
#include "../datastructures/subgraph_extractor.h"

#include "../definitions.h"
#include "../configuration/Deep_HeiProMap_configuration.h"
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

    public:
        explicit DeepHeiProMapSolver(const DeepHeiProMapConfiguration &t_ac) {
            sp = std::chrono::high_resolution_clock::now();

            ac = t_ac;
            random_engine = RandomEngine(ac.seed);

            auto sp_temp = std::chrono::high_resolution_clock::now();
            graphs.emplace_back(ac.graph_in);
            auto ep_temp = std::chrono::high_resolution_clock::now();
            std::cout << "load_graph level " << 0 << " in " << get_milli_seconds(sp_temp, ep_temp) << std::endl;

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
            rebalancer.initialize(graphs[0].n, graphs[0].m, ac.k, ac.seed);
            bv_manager.initialize(graphs[0].n, ac.k);
            q_graph.initialize(ac.k);
            block_conn.initialize(graphs[0].n, graphs[0].m, ac.k);
            d_oracle.initialize(ac.hierarchy, ac.distance);

            // matching
            gpa_matcher.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine, ac.global_path_algorithm_config);
            size_constrained_lp_clustering.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, ac.size_constrained_lp_clustering_configuration);
        }

        std::vector<partition_t> solve() {
            internal_solve();

            weight_t qap = get_qap(graphs.back(), p_manager, d_oracle);

            std::vector<partition_t> p(graphs.back().n);
            for (vertex_t u = 0; u < graphs.back().n; ++u) { p[u] = p_manager[u]; }
            write_partition(p, ac.mapping_out);

            const auto ep = std::chrono::high_resolution_clock::now();
            f64 duration = get_seconds(sp, ep);

            std::cout << "Graph                   : " << ac.graph_in << std::endl;
            std::cout << "Total time (s)          : " << duration << std::endl;
            std::cout << "#Nodes                  : " << graphs[0].n << std::endl;
            std::cout << "#Edges                  : " << graphs[0].m << std::endl;
            std::cout << "k                       : " << ac.k << std::endl;
            std::cout << "Hierarchy               : " << ac.hierarchy_string << std::endl;
            std::cout << "Distances               : " << ac.distance_string << std::endl;
            std::cout << "Distance Oracle         : " << (std::is_same_v<DistanceOracleT, BinaryDistanceOracle> ? "binary" : "matrix") << std::endl;
            std::cout << "Lmax                    : " << lmax << std::endl;
            std::cout << "Threads                 : " << ac.threads << std::endl;
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

            return p;
        }

    private:
        void internal_solve() {
            u64 level = 0;
            [[maybe_unused]] u64 max_level = 0;

            auto sp_time = get_time_point();
            auto ep_time = get_time_point();

            size_t l = ac.hierarchy.size();
            std::cout << "[DeepHeiProMap] Starting coarsening to " << k_rem[1] * ac.initial_C << " vertices" << std::endl;

            print(k_rem);
            print(lmax_vec);
            print(ac.hierarchy);

            while (graphs.back().n > k_rem[1] * ac.initial_C) {
                std::cout << "[DeepHeiProMap] -> Coarsening level " << level << " (n=" << graphs.back().n << ", m=" << graphs.back().m << ", max_w=" << lmax_vec[l - 1] << ")..." << std::flush;
                sp_time = get_time_point();
                coarsening(level, lmax_vec[l - 1]);
                ep_time = get_time_point();
                std::cout << " done in " << get_milli_seconds(sp_time, ep_time) << " ms" << std::endl;

                std::cout << "[DeepHeiProMap] -> Contracting level " << level << "..." << std::flush;
                sp_time = get_time_point();
                contraction();
                ep_time = get_time_point();
                std::cout << " done in " << get_milli_seconds(sp_time, ep_time) << " ms (coarse n=" << graphs.back().n << ", m=" << graphs.back().m << ")" << std::endl;

                if (graphs.back().n == graphs[graphs.size() - 2].n) {
                    std::cout << "[DeepHeiProMap] Coarsening stagnated (no contraction possible). Popping duplicate level and terminating coarsening." << std::endl;
                    graphs.pop_back();
                    mappings.pop_back();
                    break;
                }

                level += 1;
            }

            max_level = level > 0 ? level - 1 : 0;
            std::cout << "[DeepHeiProMap] Coarsening finished at level " << level
                    << " (coarsest graph: n=" << graphs.back().n << ", m=" << graphs.back().m << ")" << std::endl;

            std::cout << "[DeepHeiProMap] Running initial_partitioning (level " << level << ")..." << std::flush;
            sp_time = get_time_point();
            initial_partitioning();
            ep_time = get_time_point();
            std::cout << " done in " << get_milli_seconds(sp_time, ep_time) << " ms" << std::endl;

            std::cout << "[DeepHeiProMap] Running rebalance (level " << level << ")..." << std::flush;
            sp_time = get_time_point();
            rebalance(level);
            ep_time = get_time_point();
            std::cout << " done in " << get_milli_seconds(sp_time, ep_time) << " ms" << std::endl;

            while (level > 0) {
                level -= 1;
                std::cout << "[DeepHeiProMap] Uncoarsening to level " << level << " (n=" << graphs[level].n << ")..." << std::flush;

                sp_time = get_time_point();
                uncoarsening();
                ep_time = get_time_point();
                std::cout << " done in " << get_milli_seconds(sp_time, ep_time) << " ms" << std::endl;

                std::cout << "[DeepHeiProMap] Intermediate partitioning level " << level << "..." << std::flush;
                sp_time = get_time_point();
                intermediate_partitioning(level);
                ep_time = get_time_point();
                std::cout << " done in " << get_milli_seconds(sp_time, ep_time) << " ms" << std::endl;

                std::cout << "[DeepHeiProMap] Rebalancing level " << level << "..." << std::flush;
                sp_time = get_time_point();
                rebalance(level);
                ep_time = get_time_point();
                std::cout << " done in " << get_milli_seconds(sp_time, ep_time) << " ms" << std::endl;

                std::cout << "[DeepHeiProMap] Refinement level " << level << "..." << std::flush;
                sp_time = get_time_point();
                refinement(level, max_level);
                ep_time = get_time_point();
                std::cout << " done in " << get_milli_seconds(sp_time, ep_time) << " ms" << std::endl;
            }
            std::cout << "[DeepHeiProMap] internal_solve completed!" << std::endl;
        }

        void initial_partitioning() {
            ScopedTimer _t("initial_partition", "misc", "initial_partitioning");

            partition_t k_init = ac.hierarchy.back();
            std::cout << "\n  [initial_partitioning] Partitioning coarsest graph (n=" << graphs.back().n
                      << ", m=" << graphs.back().m << ") into k=" << k_init << " blocks using HeiPaSolver..." << std::flush;

            HeiPaConfiguration h_ac;
            h_ac.k = k_init;
            h_ac.imbalance = ac.imbalance;
            h_ac.threads = ac.threads;
            h_ac.seed = random_engine.get_u64();

            if (ac.config_string == "fast") {
                h_ac.set_fast();
            } else if (ac.config_string == "eco") {
                h_ac.set_eco();
            } else if (ac.config_string == "strong") {
                h_ac.set_strong();
            } else if (ac.config_string == "experimental") {
                h_ac.set_experimental();
            } else {
                h_ac.set_fast();
            }
            h_ac.k = k_init;
            h_ac.imbalance = ac.imbalance;
            h_ac.threads = ac.threads;

            graph_t g_copy = graphs.back();
            HeiPaSolver heipa_solver(std::move(g_copy), h_ac);
            const PartitionManager &sub_pm = heipa_solver.solve_subproblem();

            p_manager.reset_weights();
            for (vertex_t u = 0; u < graphs.back().n; ++u) {
                p_manager.set(u, graphs.back().v_weights[u], sub_pm[u]);
            }

            bv_manager.compute_from_scratch(graphs.back(), p_manager);
            q_graph.compute_from_scratch(graphs.back(), p_manager);
            block_conn.compute_from_scratch(graphs.back(), p_manager);

            std::cout << " done." << std::flush;
        }

        void intermediate_partitioning([[maybe_unused]] const u64 level) {
            ScopedTimer _t("intermediate_partitioning", "misc", "compute_bv_manager_from_scratch");
            std::cout << "\n  [intermediate_partitioning] Computing bv_manager from scratch..." << std::flush;
            bv_manager.compute_from_scratch(graphs.back(), p_manager);
            std::cout << " done." << std::flush;
        }

        void coarsening(const u64 level, const weight_t max_v_weight) {
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
        }

        void contraction() {
            graphs.emplace_back(); // coarse the graph
            graphs.back().contract(graphs[graphs.size() - 2], mappings.back(), ac.threads);
            p_manager.contract(mappings.back());
        }

        void uncoarsening() {
            p_manager.uncontract(mappings.back());
            mappings.pop_back();

            ScopedTimer _t("uncontraction", "misc", "compute_from_scratch");
            graphs.pop_back();
            bv_manager.compute_from_scratch(graphs.back(), p_manager);
        }

        void refinement([[maybe_unused]] const u64 level, [[maybe_unused]] const u64 max_level) {
        }

        void rebalance(const u64 level) {
            if (level == 0) {
                rebalancer.rebalance_last_layer(graphs.back(), p_manager, bv_manager, q_graph, d_oracle, block_conn, ac.imbalance);
            } else {
                rebalancer.rebalance(graphs.back(), p_manager, bv_manager, q_graph, d_oracle, block_conn, ac.imbalance);
            }
        }
    };
}

#endif //HEIPROMAP_DEEP_HEIPROMAP_SOLVER_H
