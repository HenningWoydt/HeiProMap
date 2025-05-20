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

#include "deep_boundary_vertex_manager.h"
#include "deep_distance_oracle.h"
#include "deep_partition_manager.h"
#include "deep_quotient_graph.h"
#include "../../commons/definitions.h"
#include "../../commons/macros.h"
#include "../../commons/matching.h"
#include "../../commons/random_engine.h"
#include "../../commons/small_statistic_collector.h"
#include "../../commons/statistic_collector.h"
#include "../../commons/utils.h"
#include "../coarsening/global_path_algorithm.h"
#include "../coarsening/greedy_edge_matcher.h"
#include "../coarsening/heavy_edge_matcher.h"
#include "../coarsening/random_edge_matcher.h"
#include "../partitioning/global_multisection.h"
#include "../partitioning/greedy_kway_partitioner.h"
#include "../partitioning/kaffpa_kway_partitioner.h"
#include "../partitioning/kaffpa_partitioner.h"
#include "../refinement/flow_based_refinement.h"
#include "../refinement/hierarchy_aware_multi_try_multi_way_fm_refinement.h"
#include "../refinement/three_vertex_label_propagation_refinement.h"
#include "../refinement/two_vertex_label_propagation_refinement.h"
#include "../refinement/zero_gain_perturbator.h"
#include "../utility/algorithm_configuration.h"
#include "../utility/assert_state.h"
#include "../utility/qap.h"

namespace HeiProMap {
    /**
     * Solver for serial Process Mapping.
     */
    class DeepSolver {
        AlgorithmConfiguration ac;
        RandomEngine           random_engine;

        // statistics
        StatisticCollector                             stat_collect;
        SmallStatisticCollector                        small_stat_collect;
        s64                                            initial_qap              = 0;
        weight_t                                       initial_max_block_weight = 0;
        std::chrono::high_resolution_clock::time_point sp;

        std::vector<graph_t> graphs;
        DeepPartitionManager p_manager;
        DeepDistanceOracle d_oracle;
        DeepBoundaryVertexManager bv_manager;
        DeepQuotientGraph q_graph;

        // KaffpaKWayPartitioner partitioner;
        GreedyKWayPartitioner partitioner;

        // balance
        weight_t                 lmax = 0;
        std::vector<weight_t>    lmax_vec;
        std::vector<partition_t> k_rem;

        // matching
        std::vector<Matching>      matches;
        GlobalPathAlgorithmMatcher gpa_matcher;

        // refinement
        std::vector<std::pair<ISerialRefiner *, ISerialRefinerConfiguration *>> refinements;

    public:
        explicit DeepSolver(const AlgorithmConfiguration &t_ac) {
            sp = std::chrono::high_resolution_clock::now();

            ac            = t_ac;
            random_engine = RandomEngine(ac.seed);

            const auto sp_graph_io = std::chrono::high_resolution_clock::now();
            graphs.emplace_back(ac.graph_in);
            const auto ep_graph_io = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("graph_io", get_seconds(sp_graph_io, ep_graph_io));

            const auto sp_io = std::chrono::high_resolution_clock::now();

            // balance
            lmax = std::ceil((1.0 + ac.imbalance) * ((f64) graphs[0].weight() / (f64) ac.k));
            lmax_vec.resize(ac.hierarchy.size());
            partition_t temp_k = 1;
            k_rem.push_back(temp_k);
            for (u64 i = 0; i < ac.hierarchy.size(); ++i) {
                temp_k *= ac.hierarchy[ac.hierarchy.size() - 1 - i];
                k_rem.push_back(k_rem.back() * ac.hierarchy[i]);
                lmax_vec[ac.hierarchy.size() - 1 - i] = std::ceil((1.0 + ac.imbalance) * ((f64) graphs[0].weight() / (f64) temp_k));
            }

            // manager
            p_manager.initialize(graphs[0].get_n(), ac.k, lmax_vec.back());
            p_manager.set_hierarchy_level(0, (s32) ac.hierarchy.size() - 1);
            p_manager.set_lmax(0, lmax_vec.back());
            forall_gu(graphs.back(), u)
                {
                    p_manager.set(u, graphs.back().weight(u), 0);
                }
            endfor

            bv_manager.initialize(graphs[0].get_n(), ac.k);
            q_graph.initialize(ac.k);

            // distance
            d_oracle.initialize(ac.hierarchy, ac.distance);

            // matching
            gpa_matcher.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax_vec.back(), random_engine, ac.global_path_algorithm_config, stat_collect);

            // refinement

            for (auto &[refiner, config]: refinements) {
                if (config->enabled) {
                    refiner->initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, ac.imbalance, lmax, ac.hierarchy, ac.distance, random_engine, *config, stat_collect);
                }
            }

            const auto ep_io = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("io", get_seconds(sp_io, ep_io));
            METRICS(stat_collect.set_io(get_seconds(sp_graph_io, ep_graph_io), get_seconds(sp_io, ep_io));)
        }

        std::vector<vertex_t> solve() {
            internal_solve();

            weight_t qap = get_qap(graphs.back(), p_manager, d_oracle);

            METRICS(stat_collect.set_final(qap, p_manager.get_bweights(), lmax);)
            METRICS(stat_collect.finalize();)
            METRICS(std::ofstream stat_file(ac.statistics_out);)
            METRICS(stat_file << stat_collect.to_JSON();)

            std::vector<partition_t> p(graphs.back().get_n());
            for (vertex_t            u = 0; u < graphs.back().get_n(); ++u) { p[u] = p_manager[u]; }
            write_partition(p, ac.mapping_out);

            small_stat_collect.print();

            const auto ep       = std::chrono::high_resolution_clock::now();
            f64        duration = get_seconds(sp, ep);

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
            for(partition_t id = 0; id < ac.k; ++id){
                n_empty_partitions += p_manager.get_bweight(id) == 0;
                n_overloaded_partitions += p_manager.get_bweight(id) > lmax;
                sum_too_much += std::max(0, p_manager.get_bweight(id) - lmax);
            }
            std::cout << "#empty partitions : " << n_empty_partitions << std::endl;
            std::cout << "#oload partitions : " << n_overloaded_partitions << std::endl;
            std::cout << "Sum oload weights : " << sum_too_much << std::endl;

            return p;
        }

    private:
        void internal_solve() {
            print(lmax_vec);

            u64 level     = 0;
            u64 max_level = 0;

            while (graphs.back().get_n() > ac.hierarchy.back() * 64) {
                matching(level);
                if (matches.back().size() == 0) {
                    std::cout << "No matching found!" << std::endl;
                    graphs.back().write_graph("temp.graph");
                    matches.pop_back();
                    break;
                }

                coarsening(level);

                level += 1;
            }

            max_level = level - 1;

            partition();

            while (level > 0) {
                level -= 1;
                uncoarsening(level);

                partition_subgraphs(level);

                refinement(level, max_level);
            }

            METRICS(stat_collect.add_matching_method_stats(gpa_matcher.get_stats());)
            for (auto [refiner, config]: refinements) {
                if (config->enabled) {
                    METRICS(stat_collect.add_refinement_method_stats(config->name, refiner->get_stats());)
                }
            }
        }

        void partition() {
            const auto sp_partition = std::chrono::high_resolution_clock::now();

            std::cout << "partition" << std::endl;

            partition_t id = 0;
            partitioner.partition(graphs.back(), p_manager, bv_manager, q_graph, id, k_rem[p_manager.get_hierarchy_level(id)], ac.hierarchy[p_manager.get_hierarchy_level(id)], lmax_vec[p_manager.get_hierarchy_level(id)], p_manager.get_hierarchy_level(id), random_engine, ac.greedy_kway_partitioner_config, stat_collect);

            initial_qap              = get_qap(graphs.back(), p_manager, d_oracle);
            initial_max_block_weight = max(p_manager.get_bweights());

            const auto ep_partition = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("partition", get_seconds(sp_partition, ep_partition));

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, ac.k));

            METRICS(stat_collect.set_partition_time(get_seconds(sp_partition, ep_partition));)
            METRICS(stat_collect.set_partition_stats(get_qap(graphs.back(), p_manager, d_oracle), p_manager.get_bweights(), lmax);)
        }

        void partition_subgraphs(const u64 level) {
            const auto sp_partition = std::chrono::high_resolution_clock::now();

            std::cout << "partition subgraphs level " << level << std::endl;

            std::vector<partition_t> ids;
            for (partition_t         id = 0; id < ac.k; ++id) {
                if (p_manager.get_hierarchy_level(id) >= 0 && (level == 0 || p_manager.size(id) > 64 * ac.hierarchy[p_manager.get_hierarchy_level(id)])) {
                    ids.push_back(id);
                }
            }

            for (partition_t id: ids) {
                partitioner.partition(graphs.back(), p_manager, bv_manager, q_graph, id, k_rem[p_manager.get_hierarchy_level(id)], ac.hierarchy[p_manager.get_hierarchy_level(id)], lmax_vec[p_manager.get_hierarchy_level(id)], p_manager.get_hierarchy_level(id), random_engine, ac.greedy_kway_partitioner_config, stat_collect);
            }

            const auto ep_partition = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("partition", get_seconds(sp_partition, ep_partition));

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, ac.k));
        }

        void matching(const u64 level) {
            const auto sp_match = std::chrono::high_resolution_clock::now();

            matches.emplace_back();
            matches.back().initialize(graphs.back().get_n());

            gpa_matcher.match(level, graphs.back(), p_manager, matches.back());

            const auto ep_match = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("matching", get_seconds(sp_match, ep_match));

            METRICS(stat_collect.set_matching_time(get_seconds(sp_match, ep_match), level);)
            METRICS(stat_collect.set_matching_stats(level, matches.back().size());)
        }

        void coarsening(const u64 level) {
            const auto sp_coarse = std::chrono::high_resolution_clock::now();

            graphs.emplace_back(); // coarse the graph
            graphs.back().initialize(graphs[graphs.size() - 2], matches.back());
            p_manager.contract(matches.back());

            const auto ep_coarse = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("coarsening", get_seconds(sp_coarse, ep_coarse));

            HEAVYASSERT(assert_state_pre_partitioning(graphs.back(), p_manager, ac.k));
            METRICS(stat_collect.set_coarsening_time(get_seconds(sp_coarse, ep_coarse), level);)
            METRICS(stat_collect.set_coarsening_stats(graphs.back().get_n(), level);)
        }

        void uncoarsening(const u64 level) {
            const auto sp_uncoarse = std::chrono::high_resolution_clock::now();

            std::cout << "uncoarsening" << std::endl;

            p_manager.uncontract(matches.back());
            bv_manager.compute_from_scratch(graphs[graphs.size() - 2], p_manager);
            graphs.pop_back(); // this is doing uncontraction
            matches.pop_back();

            const auto ep_uncoarse = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("uncoarsening", get_seconds(sp_uncoarse, ep_uncoarse));

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, ac.k));
            METRICS(stat_collect.set_uncoarsening_time(get_seconds(sp_uncoarse, ep_uncoarse), level);)
            METRICS(stat_collect.set_uncoarsening_stats(level, graphs.back().get_n());)
        }

        void refinement(const u64 level, const u64 max_level) {}
    };
}

#endif //HEIPROMAP_DEEP_SOLVER_H
