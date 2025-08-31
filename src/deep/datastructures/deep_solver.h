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

#ifndef HEIPROMAP_DEEP_SOLVER_H
#define HEIPROMAP_DEEP_SOLVER_H

#include <cmath>
#include <omp.h>

#include "deep_boundary_vertex_manager.h"
#include "deep_partition_manager.h"
#include "deep_quotient_graph.h"
#include "../../commons/definitions.h"
#include "../../commons/macros.h"
#include "../../commons/matching.h"
#include "../../commons/random_engine.h"
#include "../../commons/small_statistic_collector.h"
#include "../../commons/utils.h"
#include "../../serial/coarsening/greedy_edge_matcher.h"
#include "../coarsening/parallel_global_path_algorithm.h"
#include "../coarsening/suitor_algorithm.h"
#include "../coarsening/heavy_edge_matching.h"
#include "../partition/kaffpa_kway_partitioner.h"
#include "../rebalance/deep_rebalancer.h"
#include "../../serial/refinement/flow_based_refinement.h"
#include "../refinement/deep_flow_based_refinement.h"
#include "../refinement/deep_lightning_refinement.h"
#include "../refinement/deep_quotient_graph_refinement.h"
#include "../../serial/utility/algorithm_configuration.h"
#include "../../serial/utility/assert_state.h"
#include "../../serial/utility/qap.h"
#include "../utility/deep_algorithm_configuration.h"
#include "../utility/deep_assert_state.h"

namespace HeiProMap {
    /**
     * Solver for serial Process Mapping.
     */
    class DeepSolver {
        DeepAlgorithmConfiguration ac;
        RandomEngine random_engine;

        // statistics
        SmallStatisticCollector small_stat_collect;
        s64 initial_qap = 0;
        weight_t initial_max_block_weight = 0;
        std::chrono::high_resolution_clock::time_point sp;

        std::vector<deep_graph_t> graphs;
        deep_p_manager_t p_manager;
        deep_d_oracle_t d_oracle;
        deep_bv_manager_t bv_manager;
        deep_q_graph_t q_graph;

        KaffpaKWayPartitioner partitioner;
        DeepRebalancer deep_rebalancer;

        // balance
        weight_t lmax = 0;
        std::vector<weight_t> lmax_vec;
        std::vector<partition_t> k_rem;

        // matching
        std::vector<Matching> matches;
        // GlobalPathAlgorithmMatcher gpa_matcher;
        ParallelGlobalPathAlgorithmMatcher parallel_gpa_matcher;
        SuitorMatcher suitor_matcher;
        ParallelHeavyEdgeMatching heavy_edge_matcher;

        // refinement
        std::vector<std::pair<ISerialDeepRefiner *, ISerialDeepRefinerConfiguration *> > refinements;

        DeepQuotientGraphRefinement deep_quotient_graph_refinement;
        DeepFlowBasedRefinement deep_flow_based_refinement;

    public:
        explicit DeepSolver(const DeepAlgorithmConfiguration &t_ac) {
            sp = std::chrono::high_resolution_clock::now();

            ac = t_ac;
            random_engine = RandomEngine(ac.seed);

            const auto sp_graph_io = std::chrono::high_resolution_clock::now();
            graphs.emplace_back(ac.graph_in);
            const auto ep_graph_io = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("graph_io", get_seconds(sp_graph_io, ep_graph_io));

            const auto sp_io = std::chrono::high_resolution_clock::now();

            // balance
            lmax = std::ceil((1.0 + ac.imbalance) * ((f64) graphs[0].weight() / (f64) ac.k));
            lmax_vec.resize(ac.hierarchy.size());
            lmax_vec[0] = lmax;
            for (u64 i = 1; i < ac.hierarchy.size(); ++i) {
                lmax_vec[i] = lmax_vec[i - 1] * ac.hierarchy[i - 1];
            }

            partition_t temp_k = 1;
            k_rem.push_back(temp_k);
            for (u64 i = 0; i < ac.hierarchy.size(); ++i) {
                temp_k *= ac.hierarchy[ac.hierarchy.size() - 1 - i];
                k_rem.push_back(k_rem.back() * ac.hierarchy[i]);
            }

            // manager
            p_manager.initialize(graphs[0].get_n(), ac.k, lmax_vec.back());
            p_manager.set_hierarchy_level(0, ac.hierarchy.size());
            p_manager.set_lmax(0, lmax_vec.back());
            forall_gu(graphs.back(), u)
                {
                    p_manager.set(u, graphs.back().weight(u), 0);
                }
            endfor

            partitioner.initialize(graphs.back().get_n(), ac.k, ac.threads, random_engine, ac.kaffpa_kway_partitioner_config);

            deep_rebalancer.initialize(ac.threads, random_engine);
            bv_manager.initialize(graphs[0].get_n(), ac.k);
            q_graph.initialize(ac.hierarchy, ac.k);
            // distance
            d_oracle.initialize(ac.hierarchy, ac.distance, ac.threads);

            // matching
            parallel_gpa_matcher.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax_vec.back(), ac.threads, random_engine, ac.global_path_algorithm_config);
            suitor_matcher.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax_vec.back(), ac.threads, random_engine, SuitorMatcherConfiguration());
            heavy_edge_matcher.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax_vec.back(), ac.threads, random_engine, ac.parallel_heavy_edge_matching_configuration);

            refinements.emplace_back(&deep_quotient_graph_refinement, &ac.deep_quotient_graph_refinement_config);
            refinements.emplace_back(&deep_flow_based_refinement, &ac.deep_flow_based_refinement_config);

            for (auto &[refiner, config]: refinements) {
                if (config->enabled) {
                    refiner->initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, ac.imbalance, ac.threads, ac.hierarchy, ac.distance, random_engine, *config);
                }
            }

            const auto ep_io = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("io", get_seconds(sp_io, ep_io));
        }

        std::vector<partition_t> solve() {
            internal_solve();

            weight_t qap = get_qap(graphs.back(), p_manager, d_oracle, ac.threads);

            std::vector<partition_t> p(graphs.back().get_n());
            for (vertex_t u = 0; u < graphs.back().get_n(); ++u) { p[u] = p_manager[u]; }
            write_partition(p, ac.mapping_out);

            small_stat_collect.print();

            const auto ep = std::chrono::high_resolution_clock::now();
            f64 duration = get_seconds(sp, ep);

            std::cout << "Total time        : " << duration << std::endl;
            std::cout << "#Nodes            : " << graphs.back().get_n() << std::endl;
            std::cout << "#Edges            : " << graphs.back().get_m() << std::endl;
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

            return p;
        }

    private:
        void internal_solve() {
            u64 level = 0;
            u64 max_level = 0;

            u64 threshold = 16;

            while (graphs.back().get_n() > ac.hierarchy.back() * threshold) {
                auto sp = std::chrono::high_resolution_clock::now();
                matching(level);
                if (matches.back().size() == 0) {
                    matches.pop_back();
                    break;
                }

                coarsening(level);
                auto ep = std::chrono::high_resolution_clock::now();
                std::cout << level << " " << graphs.back().get_n() << " " << matches.back().size() << " " << get_seconds(sp, ep) << " " << get_memory_usage_gb() << std::endl;

                level += 1;
            }

            max_level = level - 1;

            partition();
            rebalance(level);
            std::cout << level << " First partition " << graphs.back().get_n() << " " << get_qap(graphs.back(), p_manager, d_oracle, ac.threads) << " " << get_memory_usage_gb() << std::endl;
            print(get_qap_per_layer(graphs.back(), p_manager, d_oracle, ac.hierarchy.size()));

            while (level > 0) {
                level -= 1;

                std::cout << level << " A " << graphs.back().get_n() << " " << get_memory_usage_gb() << std::endl;

                uncoarsening(level);

                std::cout << level << " B " << graphs.back().get_n() << " " << get_memory_usage_gb() << std::endl;

                partition_subgraphs(level, threshold);

                std::cout << level << " C " << graphs.back().get_n() << " " << get_memory_usage_gb() << std::endl;

                rebalance(level);

                std::cout << level << " D " << graphs.back().get_n() << " " << get_memory_usage_gb() << std::endl;

                refinement(level, max_level);

                std::cout << level << " E " << graphs.back().get_n() << " " << get_memory_usage_gb() << std::endl;
            }
        }

        void partition() {
            const auto sp_partition = std::chrono::high_resolution_clock::now();

            partition_t id = 0;
            u64 thread_id = 0;

            partitioner.determine_all_blocks(graphs.back(), p_manager);
            partitioner.partition(graphs.back(), p_manager, thread_id, id, k_rem[p_manager.get_hierarchy_level(id) - 1], ac.hierarchy[p_manager.get_hierarchy_level(id) - 1], lmax_vec[p_manager.get_hierarchy_level(id) - 1], p_manager.get_hierarchy_level(id));
            q_graph.compute_from_scratch(graphs.back(), p_manager);
            bv_manager.compute_from_scratch(graphs.back(), p_manager);

            initial_qap = get_qap(graphs.back(), p_manager, d_oracle);
            initial_max_block_weight = max(p_manager.get_bweights());

            const auto ep_partition = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("partition", get_seconds(sp_partition, ep_partition));

            HEAVYASSERT(deep_assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, ac.k));
        }

        void partition_subgraphs(const u64 level, const u64 threshold) {
            const auto sp_partition = std::chrono::high_resolution_clock::now();

            std::vector<partition_t> ids;
            for (partition_t id = 0; id < ac.k; ++id) {
                if (p_manager.is_active(id) && p_manager.get_hierarchy_level(id) > 0 && (level == 0 || p_manager.size(id) > threshold * ac.hierarchy[p_manager.get_hierarchy_level(id) - 1])) {
                    ids.push_back(id);
                }
            }

            if (!ids.empty()) {
                partitioner.determine_all_blocks(graphs.back(), p_manager);
#pragma omp parallel for num_threads(ac.threads) schedule(dynamic)
                for (size_t i = 0; i < ids.size(); ++i) {
                    partition_t id = ids[i];
                    partition_t thread_id = omp_get_thread_num();
                    partitioner.partition(graphs.back(), p_manager, thread_id, id, k_rem[p_manager.get_hierarchy_level(id) - 1], ac.hierarchy[p_manager.get_hierarchy_level(id) - 1], lmax_vec[p_manager.get_hierarchy_level(id) - 1], p_manager.get_hierarchy_level(id));
                }
                q_graph.compute_from_scratch(graphs.back(), p_manager);
                bv_manager.compute_from_scratch(graphs.back(), p_manager);
            }

            const auto ep_partition = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("partition", get_seconds(sp_partition, ep_partition));

            HEAVYASSERT(deep_assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, ac.k));
        }

        void matching(const u64 level) {
            const auto sp_match = std::chrono::high_resolution_clock::now();

            matches.emplace_back();
            matches.back().initialize(graphs.back().get_n());

            // parallel_gpa_matcher.match(level, graphs.back(), p_manager, matches.back());
            std::cout << "matching level: " << level << std::endl;
            // suitor_matcher.match(level, graphs.back(), p_manager, matches.back());
            heavy_edge_matcher.match(level, graphs.back(), p_manager, matches.back());

            const auto ep_match = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("matching", get_seconds(sp_match, ep_match));

            HEAVYASSERT(deep_assert_matching(graphs.back(), matches.back()));
        }

        void coarsening(const u64 level) {
            const auto sp_coarse = std::chrono::high_resolution_clock::now();

            graphs.emplace_back(); // coarse the graph
            graphs.back().parallel_initialize(graphs[graphs.size() - 2], matches.back(), ac.threads);
            p_manager.contract(matches.back());

            const auto ep_coarse = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("coarsening", get_seconds(sp_coarse, ep_coarse));

            HEAVYASSERT(deep_assert_state_pre_partitioning(graphs.back(), p_manager, ac.k));
        }

        void uncoarsening(const u64 level) {
            const auto sp_uncoarse = std::chrono::high_resolution_clock::now();

            p_manager.uncontract(matches.back());
            matches.pop_back();
            graphs.pop_back();
            bv_manager.compute_from_scratch(graphs.back(), p_manager);

            const auto ep_uncoarse = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("uncoarsening", get_seconds(sp_uncoarse, ep_uncoarse));

            HEAVYASSERT(deep_assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, ac.k));
        }

        void refinement(const u64 level, const u64 max_level) {
            const auto sp_refinement = std::chrono::high_resolution_clock::now();

            SMALL_METRICS(s64 qap_before = get_qap(graphs.back(), p_manager, d_oracle, ac.threads);)

            u64 refinement_max_iterations = 1;
            for (u64 refinement_i = 0; refinement_i < refinement_max_iterations; ++refinement_i) {
                for (auto [refiner, config]: refinements) {
                    if (config->enabled) {
                        const auto sp = std::chrono::high_resolution_clock::now();

                        HEAVYASSERT(deep_assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, ac.k));
                        refiner->refine(level, max_level, graphs.back(), d_oracle, bv_manager, p_manager, q_graph);

                        const auto ep = std::chrono::high_resolution_clock::now();
                        HEAVYASSERT(deep_assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, ac.k));
                        SMALL_METRICS(s64 qap_after = get_qap(graphs.back(), p_manager, d_oracle, ac.threads);)
                        s64 qap_delta = 0;
                        SMALL_METRICS(qap_delta = qap_before - qap_after;)

                        small_stat_collect.add_refinement(config->name, get_seconds(sp, ep), qap_delta);

                        SMALL_METRICS(qap_before = qap_after;)
                    }
                }
            }

            const auto ep_refinement = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("refinement", get_seconds(sp_refinement, ep_refinement));

            HEAVYASSERT(deep_assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, ac.k));
        }

        void rebalance(const u64 level) {
            const auto sp_rebalance = std::chrono::high_resolution_clock::now();

            if (level == 0) {
                deep_rebalancer.rebalance_last_layer(graphs.back(), p_manager, bv_manager, q_graph, d_oracle, ac.k);
            } else {
                deep_rebalancer.rebalance(graphs.back(), p_manager, bv_manager, q_graph, d_oracle, ac.k);
            }

            const auto ep_rebalance = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("rebalance", get_seconds(sp_rebalance, ep_rebalance));

            HEAVYASSERT(deep_assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, ac.k));
        }
    };
}

#endif //HEIPROMAP_DEEP_SOLVER_H
