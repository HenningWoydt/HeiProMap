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
#include "deep_distance_oracle.h"
#include "deep_partition_manager.h"
#include "deep_quotient_graph.h"
#include "../../../commons/definitions.h"
#include "../../../commons/macros.h"
#include "../../../commons/matching.h"
#include "../../../commons/random_engine.h"
#include "../../../commons/small_statistic_collector.h"
#include "../../../commons/statistic_collector.h"
#include "../../../commons/utils.h"
#include "../../coarsening/global_path_algorithm.h"
#include "../../coarsening/greedy_edge_matcher.h"
#include "../../coarsening/heavy_edge_matcher.h"
#include "../../coarsening/random_edge_matcher.h"
#include "../../partitioning/global_multisection.h"
#include "../../partitioning/greedy_kway_partitioner.h"
#include "../../partitioning/kaffpa_kway_partitioner.h"
#include "../../partitioning/kaffpa_partitioner.h"
#include "../../refinement/flow_based_refinement.h"
#include "../../refinement/hierarchy_aware_multi_try_multi_way_fm_refinement.h"
#include "../../refinement/three_vertex_label_propagation_refinement.h"
#include "../../refinement/two_vertex_label_propagation_refinement.h"
#include "../../refinement/zero_gain_perturbator.h"
#include "../../refinement/deep/deep_flow_based_refinement.h"
#include "../../refinement/deep/deep_lightning_refinement.h"
#include "../../refinement/deep/deep_quotient_graph_refinement.h"
#include "../../rebalance/deep/deep_rebalancer.h"
#include "../../utility/algorithm_configuration.h"
#include "../../utility/assert_state.h"
#include "../../utility/qap.h"
#include "../../utility/deep/deep_assert_state.h"

namespace HeiProMap {
    /**
     * Solver for serial Process Mapping.
     */
    class DeepSolver {
        AlgorithmConfiguration ac;
        RandomEngine random_engine;

        // statistics
        SmallStatisticCollector small_stat_collect;
        s64 initial_qap                   = 0;
        weight_t initial_max_block_weight = 0;
        std::chrono::high_resolution_clock::time_point sp;

        std::vector<graph_t> graphs;
        DeepPartitionManager p_manager;
        DeepDistanceOracle d_oracle;
        DeepBoundaryVertexManager bv_manager;
        DeepQuotientGraph q_graph;

        KaffpaKWayPartitioner partitioner;

        // balance
        weight_t lmax = 0;
        std::vector<weight_t> lmax_vec;
        std::vector<partition_t> k_rem;

        // matching
        std::vector<Matching> matches;
        GlobalPathAlgorithmMatcher gpa_matcher;

        // refinement
        std::vector<std::pair<ISerialDeepRefiner*, ISerialDeepRefinerConfiguration*>> refinements;

        DeepQuotientGraphRefinementConfiguration deep_quotient_graph_refinement_config = DeepQuotientGraphRefinementConfiguration("Deep Quotient Graph Refinement");
        DeepQuotientGraphRefinement deep_quotient_graph_refinement;

        DeepFlowBasedRefinementConfiguration deep_flow_based_refinement_config = DeepFlowBasedRefinementConfiguration("Deep Flow Based Refinement");
        DeepFlowBasedRefinement deep_flow_based_refinement;

        DeepLightningRefinementConfiguration deep_lightning_refinement_config = DeepLightningRefinementConfiguration("Deep Lighting Refinement");
        DeepLightningRefinement deep_lightning_refinement;

    public:
        explicit DeepSolver(const AlgorithmConfiguration& t_ac) {
            sp = std::chrono::high_resolution_clock::now();

            ac            = t_ac;
            random_engine = RandomEngine(ac.seed);

            const auto sp_graph_io = std::chrono::high_resolution_clock::now();
            graphs.emplace_back(ac.graph_in);
            const auto ep_graph_io = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("graph_io", get_seconds(sp_graph_io, ep_graph_io));

            const auto sp_io = std::chrono::high_resolution_clock::now();

            // balance
            lmax = std::ceil((1.0 + ac.imbalance) * ((f64)graphs[0].weight() / (f64)ac.k));
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

            bv_manager.initialize(graphs[0].get_n(), ac.k);
            q_graph.initialize(ac.k);

            // distance
            d_oracle.initialize(ac.hierarchy, ac.distance, ac.threads);

            // matching
            gpa_matcher.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax_vec.back(), ac.threads, random_engine, ac.global_path_algorithm_config);

            // refinement
            deep_quotient_graph_refinement_config.enabled       = true;
            deep_quotient_graph_refinement_config.max_iteration = 5;
            deep_quotient_graph_refinement_config.alpha         = 1000.0;

            deep_flow_based_refinement_config.enabled                    = true;
            deep_flow_based_refinement_config.min_level                  = 0;
            deep_flow_based_refinement_config.max_level                  = 100;
            deep_flow_based_refinement_config.max_global_iteration       = 2;
            deep_flow_based_refinement_config.max_local_iteration        = 5;
            deep_flow_based_refinement_config.alpha                      = 2.0;
            deep_flow_based_refinement_config.alpha_upper_bound          = 16.0;
            deep_flow_based_refinement_config.alpha_modifier             = 2.0;
            deep_flow_based_refinement_config.use_closed_vertex_set      = true;
            deep_flow_based_refinement_config.closed_vertex_sets_repeats = 100;

            deep_lightning_refinement_config.enabled = false;

            refinements.emplace_back(&deep_quotient_graph_refinement, &deep_quotient_graph_refinement_config);
            refinements.emplace_back(&deep_flow_based_refinement, &deep_flow_based_refinement_config);
            refinements.emplace_back(&deep_lightning_refinement, &deep_lightning_refinement_config);

            for (auto& [refiner, config] : refinements) {
                if (config->enabled) {
                    refiner->initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, ac.imbalance, ac.hierarchy, ac.distance, random_engine, *config);
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
            f64 duration  = get_seconds(sp, ep);

            std::cout << "Total time        : " << duration << std::endl;
            std::cout << "#Nodes            : " << graphs.back().get_n() << std::endl;
            std::cout << "#Edges            : " << graphs.back().get_m() << std::endl;
            std::cout << "k                 : " << ac.k << std::endl;
            std::cout << "Lmax              : " << lmax << std::endl;
            std::cout << "Init. QAP         : " << initial_qap << std::endl;
            std::cout << "Init. max block w : " << initial_max_block_weight << std::endl;
            std::cout << "Final QAP         : " << qap << std::endl;
            std::cout << "max block w       : " << max(p_manager.get_bweights()) << std::endl;

            size_t n_empty_partitions      = 0;
            size_t n_overloaded_partitions = 0;
            weight_t sum_too_much          = 0;
            for (partition_t id = 0; id < ac.k; ++id) {
                n_empty_partitions += p_manager.get_bweight(id) == 0;
                n_overloaded_partitions += p_manager.get_bweight(id) > lmax;
                sum_too_much += std::max((weight_t)0, p_manager.get_bweight(id) - lmax);
            }
            std::cout << "#empty partitions : " << n_empty_partitions << std::endl;
            std::cout << "#oload partitions : " << n_overloaded_partitions << std::endl;
            std::cout << "Sum oload weights : " << sum_too_much << std::endl;

            return p;
        }

    private:
        void internal_solve() {
            DeepRebalancer deep_rebalancer;

            u64 level     = 0;
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
                std::cout << level << " " << graphs.back().get_n() << " " << get_seconds(sp, ep) << std::endl;

                level += 1;
            }

            max_level = level - 1;

            partition();
            deep_rebalancer.rebalance(graphs.back(), p_manager, bv_manager, q_graph, ac.k);
            std::cout << level << " " << graphs.back().get_n() << " " << get_qap(graphs.back(), p_manager, d_oracle, ac.threads) << std::endl;
            print(get_qap_per_layer(graphs.back(), p_manager, d_oracle, ac.hierarchy.size()));

            while (level > 0) {
                auto sp = std::chrono::high_resolution_clock::now();
                level -= 1;
                uncoarsening(level);

                partition_subgraphs(level, threshold);
                deep_rebalancer.rebalance(graphs.back(), p_manager, bv_manager, q_graph, ac.k);

                refinement(level, max_level);
                auto ep = std::chrono::high_resolution_clock::now();

                std::cout << level << " " << graphs.back().get_n() << " " << get_qap(graphs.back(), p_manager, d_oracle, ac.threads) << " " << get_seconds(sp, ep) << std::endl;
                print(get_qap_per_layer(graphs.back(), p_manager, d_oracle, ac.hierarchy.size(), ac.threads));
            }
        }

        void partition() {
            const auto sp_partition = std::chrono::high_resolution_clock::now();

            partition_t id = 0;
            u64 thread_id = 0;
            partitioner.partition(graphs.back(), p_manager, bv_manager, q_graph, thread_id, id, k_rem[p_manager.get_hierarchy_level(id) - 1], ac.hierarchy[p_manager.get_hierarchy_level(id) - 1], lmax_vec[p_manager.get_hierarchy_level(id) - 1], p_manager.get_hierarchy_level(id));

            initial_qap              = get_qap(graphs.back(), p_manager, d_oracle);
            initial_max_block_weight = max(p_manager.get_bweights());

            const auto ep_partition = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("partition", get_seconds(sp_partition, ep_partition));

            HEAVYASSERT(deep_assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, ac.k));
        }

        void partition_subgraphs(const u64 level, const u64 threshold) {
            const auto sp_partition = std::chrono::high_resolution_clock::now();

            std::vector<partition_t> ids;
            for (partition_t id = 0; id < ac.k; ++id) {
                if (p_manager.get_hierarchy_level(id) != ac.k && p_manager.get_hierarchy_level(id) > 0 && (level == 0 || p_manager.size(id) > threshold * ac.hierarchy[p_manager.get_hierarchy_level(id) - 1])) {
                    ids.push_back(id);
                }
            }

            std::cout << ids.size() << std::endl;
            partitioner.determine_all_blocks(graphs.back(), p_manager);
#pragma omp parallel for num_threads(ac.threads)
            for (size_t i = 0; i < ids.size(); ++i) {
                partition_t id = ids[i];
                partition_t thread_id = omp_get_thread_num();
                partitioner.partition(graphs.back(), p_manager, bv_manager, q_graph, thread_id, id, k_rem[p_manager.get_hierarchy_level(id) - 1], ac.hierarchy[p_manager.get_hierarchy_level(id) - 1], lmax_vec[p_manager.get_hierarchy_level(id) - 1], p_manager.get_hierarchy_level(id));
            }

            const auto ep_partition = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("partition", get_seconds(sp_partition, ep_partition));

            HEAVYASSERT(deep_assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, ac.k));
        }

        void matching(const u64 level) {
            const auto sp_match = std::chrono::high_resolution_clock::now();

            matches.emplace_back();
            matches.back().initialize(graphs.back().get_n());

            gpa_matcher.match(level, graphs.back(), p_manager, matches.back());

            const auto ep_match = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("matching", get_seconds(sp_match, ep_match));
        }

        void coarsening(const u64 level) {
            const auto sp_coarse = std::chrono::high_resolution_clock::now();

            graphs.emplace_back(); // coarse the graph
            graphs.back().initialize(graphs[graphs.size() - 2], matches.back());
            p_manager.contract(matches.back());

            const auto ep_coarse = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("coarsening", get_seconds(sp_coarse, ep_coarse));

            HEAVYASSERT(deep_assert_state_pre_partitioning(graphs.back(), p_manager, ac.k));
        }

        void uncoarsening(const u64 level) {
            const auto sp_uncoarse = std::chrono::high_resolution_clock::now();

            p_manager.uncontract(matches.back());
            bv_manager.compute_from_scratch(graphs[graphs.size() - 2], p_manager);
            graphs.pop_back(); // this is doing uncontraction
            matches.pop_back();

            const auto ep_uncoarse = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("uncoarsening", get_seconds(sp_uncoarse, ep_uncoarse));

            HEAVYASSERT(deep_assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, ac.k));
        }

        void refinement(const u64 level, const u64 max_level) {
            const auto sp_refinement = std::chrono::high_resolution_clock::now();

            SMALL_METRICS(s64 qap_before = get_qap(graphs.back(), p_manager, d_oracle, ac.threads);)

            u64 refinement_max_iterations = 1;
            for (u64 refinement_i = 0; refinement_i < refinement_max_iterations; ++refinement_i) {
                for (auto [refiner, config] : refinements) {
                    if (config->enabled) {
                        const auto sp = std::chrono::high_resolution_clock::now();

                        refiner->refine(level, max_level, graphs.back(), d_oracle, bv_manager, p_manager, q_graph);
                        HEAVYASSERT(deep_assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, ac.k));

                        const auto ep = std::chrono::high_resolution_clock::now();
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
    };
}

#endif //HEIPROMAP_DEEP_SOLVER_H
