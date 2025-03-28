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
#include "../../commons/definitions.h"
#include "../../commons/macros.h"
#include "../../commons/random_engine.h"
#include "../../commons/small_statistic_collector.h"
#include "../../commons/statistic_collector.h"
#include "../../commons/utils.h"
#include "../coarsening/global_path_algorithm.h"
#include "../coarsening/greedy_edge_matcher.h"
#include "../coarsening/heavy_edge_matcher.h"
#include "../coarsening/matching.h"
#include "../coarsening/random_edge_matcher.h"
#include "../partitioning/global_multisection.h"
#include "../partitioning/kaffpa_partitioner.h"
#include "../refinement/flow_based_refinement.h"
#include "../refinement/label_propagation_refinement_Faraj20.h"
#include "../refinement/quotient_graph_refinement_Faraj20.h"
#include "../refinement/three_vertex_label_propagation_refinement.h"
#include "../refinement/two_vertex_label_propagation_refinement.h"
#include "../refinement/pertubation.h"
#include "../utility/algorithm_configuration.h"
#include "../utility/assert_state.h"
#include "../utility/qap.h"

namespace HeiProMap {
    /**
     * Solver for serial Process Mapping.
     */
    class Solver {
        AlgorithmConfiguration ac;
        RandomEngine random_engine;

        // statistics
        StatisticCollector stat_collect;
        SmallStatisticCollector small_stat_collect;
        s64 initial_qap = 0;
        weight_t initial_max_block_weight = 0;
        std::chrono::high_resolution_clock::time_point sp;

        std::vector<graph_t> graphs;
        PartitionManager p_manager;
        BoundaryVertexManager bv_manager;
        QuotientGraph q_graph;
        DistanceOracle d_oracle;

        // balance
        weight_t lmax = 0;

        // matching
        std::vector<Matching> matches;
        GreedyEdgeMatcher ge_matcher;
        HeavyEdgeMatcher he_matcher;
        RandomEdgeMatcher rnd_matcher;
        GlobalPathAlgorithmMatcher gpa_matcher;

        // refinement
        LabelPropagationRefinementFaraj20 lp_refine_faraj20;
        QuotientGraphRefinementFaraj20 qg_refine_faraj20;
        KWayFMRefinementFaraj20 k_way_refine_faraj20;
        MultiTryFMRefinementFaraj20 multi_try_fm_refinement_faraj20;

        LabelPropagationRefinement lp_refine;
        TwoVertexLabelPropagationRefinement two_vertex_lp_refine;
        ThreeVertexLabelPropagationRefinement three_vertex_lp_refine;
        QuotientGraphRefinement qg_refine;
        KWayFMRefinement k_way_refine;
        MultiTryFMRefinement multi_try_fm_refinement;
        FlowBasedRefinement flow_based_refinement;

        HierarchyAwareMultiWayFMRefinement hierarchy_aware_fm_refinement;

        Pertubation pertubation;

        std::vector<std::pair<ISerialRefiner*, ISerialRefinerConfiguration*>> refinements;

    public:
        explicit Solver(const AlgorithmConfiguration& t_ac) {
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

            // manager
            p_manager.initialize(graphs[0].get_n(), ac.k, lmax);
            bv_manager.initialize(graphs[0].get_n(), ac.k);
            q_graph.initialize(ac.k);
            HEAVYASSERT(assert_state_pre_partitioning(graphs[0]));

            // distance
            d_oracle.initialize(ac.hierarchy, ac.distance);

            // matching
            ge_matcher.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, random_engine, ac.greedy_edge_matcher_config, stat_collect);
            he_matcher.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, random_engine, ac.heavy_edge_matcher_config, stat_collect);
            rnd_matcher.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, random_engine, ac.random_edge_matcher_config, stat_collect);
            gpa_matcher.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, random_engine, ac.global_path_algorithm_config, stat_collect);

            // refinement
            refinements.emplace_back(&qg_refine_faraj20, &ac.quotient_graph_refinement_faraj20_config);
            refinements.emplace_back(&lp_refine_faraj20, &ac.label_propagation_faraj20_config);
            refinements.emplace_back(&k_way_refine_faraj20, &ac.k_way_fm_refinement_faraj20_config);
            refinements.emplace_back(&multi_try_fm_refinement_faraj20, &ac.multi_try_fm_refinement_faraj20_config);

            refinements.emplace_back(&qg_refine, &ac.quotient_graph_refinement_config);
            refinements.emplace_back(&lp_refine, &ac.label_propagation_config);
            refinements.emplace_back(&k_way_refine, &ac.k_way_fm_refinement_config);
            refinements.emplace_back(&multi_try_fm_refinement, &ac.multi_try_fm_refinement_config);
            refinements.emplace_back(&two_vertex_lp_refine, &ac.two_vertex_label_propagation_config);
            refinements.emplace_back(&three_vertex_lp_refine, &ac.three_vertex_label_propagation_config);
            refinements.emplace_back(&flow_based_refinement, &ac.flow_based_refinement_config);

            refinements.emplace_back(&hierarchy_aware_fm_refinement, &ac.hierarchy_aware_multi_way_fm_config);

            for (auto& [refiner, config] : refinements) {
                if (config->enabled) {
                    refiner->initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, ac.imbalance, lmax, ac.hierarchy, ac.distance, random_engine, *config, stat_collect);
                }
            }

            pertubation.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, ac.hierarchy, ac.distance, random_engine, stat_collect);

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
            for (vertex_t u = 0; u < graphs.back().get_n(); ++u) { p[u] = p_manager[u]; }
            write_partition(p, ac.mapping_out);

            small_stat_collect.print();

            const auto ep = std::chrono::high_resolution_clock::now();
            f64 duration  = get_seconds(sp, ep);

            std::cout << "Total time        : " << duration << std::endl;
            std::cout << "#Nodes            : " << graphs.back().get_n() << std::endl;
            std::cout << "#Edges            : " << graphs.back().get_m() << std::endl;
            std::cout << "Lmax              : " << lmax << std::endl;
            std::cout << "Init. QAP         : " << initial_qap << std::endl;
            std::cout << "Init. max block w : " << initial_max_block_weight << std::endl;
            std::cout << "Final QAP         : " << qap << std::endl;
            std::cout << "max block w       : " << max(p_manager.get_bweights()) << std::endl;

            return p;
        }

    private:
        void internal_solve() {
            u64 level     = 0;
            u64 max_level = 0;

            while (graphs.back().get_n() > ac.k * 16) {
                matching(level);
                if (matches.back().size() == 0) {
                    matches.pop_back();
                    break;
                }

                coarsening(level);

                level += 1;
            }

            max_level = level - 1;
            partition();

            ASSERT(max(p_manager.get_bweights()) <= lmax);
            if(p_manager.is_overloaded()){
                print(p_manager.get_bweights());
                std::cout << max(p_manager.get_bweights()) << std::endl;
            }

            while (level > 0) {
                level -= 1;
                uncoarsening(level);
                refinement(level, max_level);
            }

            METRICS(stat_collect.add_matching_method_stats(gpa_matcher.get_stats());)
            for (auto [refiner, config] : refinements) {
                if (config->enabled) {
                    METRICS(stat_collect.add_refinement_method_stats(config->name, refiner->get_stats());)
                }
            }
        }

        void partition() {
            const auto sp_partition = std::chrono::high_resolution_clock::now();

            if (ac.partitioning_algorithm_id == PARTITIONING_ALG_KAFFPA) {
                KaffpaPartitioner partitioner;
                partitioner.partition(graphs.back(), p_manager, ac.hierarchy, ac.distance, ac.imbalance, random_engine, ac.kaffpa_partitioner_config, stat_collect);
            } else if (ac.partitioning_algorithm_id == PARTITIONING_ALG_MULTISECTION) {
                GlobalMultisectionPartitioner partitioner;
                partitioner.partition(graphs.back(), p_manager, ac.hierarchy, ac.distance, ac.imbalance, random_engine, ac.global_multisection_config, stat_collect);
            } else {
                std::cerr << "Partitioning algorithm " << partitioning_algorithm_to_string(ac.partitioning_algorithm_id) << " with id " << ac.partitioning_algorithm_id << " not known!" << std::endl;
                exit(EXIT_FAILURE);
            }

            // initialize boundary vertices and quotient graph
            forall_gu(graphs.back(), u)
                {
                    const partition_t u_id = p_manager[u];

                    forall_guivw(graphs.back(), u, i, v, w)
                        {
                            const partition_t v_id = p_manager[v];

                            if (u_id != v_id) {
                                bv_manager.add(u, u_id); // boundary vertex
                                q_graph.add_edge(u_id, v_id, w); // quotient graph
                            }
                        }
                    endfor
                }
            endfor

            initial_qap = get_qap(graphs.back(), p_manager, d_oracle);
            initial_max_block_weight = max(p_manager.get_bweights());

            const auto ep_partition = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("partition", get_seconds(sp_partition, ep_partition));

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, ac.k));

            METRICS(stat_collect.set_partition_time(get_seconds(sp_partition, ep_partition));)
            METRICS(stat_collect.set_partition_stats(get_qap(graphs.back(), p_manager, d_oracle), p_manager.get_bweights(), lmax);)
        }

        void matching(const u64 level) {
            const auto sp_match = std::chrono::high_resolution_clock::now();

            matches.emplace_back();
            matches.back().initialize(graphs.back().get_n());

            if (ac.coarsening_algorithm_id == COARSENING_ALG_GREEDY_MATCHING) {
                ge_matcher.match(level, graphs.back(), matches.back());
            } else if (ac.coarsening_algorithm_id == COARSENING_ALG_HEAVY_MATCHING) {
                he_matcher.match(level, graphs.back(), matches.back());
            } else if (ac.coarsening_algorithm_id == COARSENING_ALG_RANDOM_MATCHING) {
                rnd_matcher.match(level, graphs.back(), matches.back());
            } else if (ac.coarsening_algorithm_id == COARSENING_ALG_GLOBAL_PATHS) {
                gpa_matcher.match(level, graphs.back(), matches.back());
            } else {
                std::cerr << "Coarsening algorithm " << coarsening_algorithm_to_string(ac.coarsening_algorithm_id) << " with id " << ac.coarsening_algorithm_id << " not known!" << std::endl;
                exit(EXIT_FAILURE);
            }

            const auto ep_match = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("matching", get_seconds(sp_match, ep_match));

            METRICS(stat_collect.set_matching_time(get_seconds(sp_match, ep_match), level);)
            METRICS(stat_collect.set_matching_stats(level, matches.back().size());)
        }

        void coarsening(const u64 level) {
            const auto sp_coarse = std::chrono::high_resolution_clock::now();

            graphs.emplace_back(); // coarse the graph
            graphs.back().initialize(graphs[graphs.size() - 2], matches.back());

            const auto ep_coarse = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("coarsening", get_seconds(sp_coarse, ep_coarse));

            HEAVYASSERT(assert_state_pre_partitioning(graphs.back()));
            METRICS(stat_collect.set_coarsening_time(get_seconds(sp_coarse, ep_coarse), level);)
            METRICS(stat_collect.set_coarsening_stats(graphs.back().get_n(), level);)
        }

        void uncoarsening(const u64 level) {
            const auto sp_uncoarse = std::chrono::high_resolution_clock::now();

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

        void refinement(const u64 level, const u64 max_level) {
            const auto sp_refinement = std::chrono::high_resolution_clock::now();

            SMALL_METRICS(s64 qap_before = get_qap(graphs.back(), p_manager, d_oracle);)
            // s64 qap_1 = get_qap(graphs.back(), p_manager, d_oracle);
            // pertubation.pertubate(level, max_level, graphs.back(), d_oracle, bv_manager, p_manager, q_graph);
            // s64 qap_2 = get_qap(graphs.back(), p_manager, d_oracle);
            // std::cout << qap_2 - qap_1 << " " << qap_2 << std::endl;

            for (auto [refiner, config] : refinements) {
                if (config->enabled) {
                    const auto sp = std::chrono::high_resolution_clock::now();

                    refiner->refine(level, max_level, graphs.back(), d_oracle, bv_manager, p_manager, q_graph);

                    const auto ep = std::chrono::high_resolution_clock::now();
                    SMALL_METRICS(s64 qap_after = get_qap(graphs.back(), p_manager, d_oracle);)
                    s64 qap_delta = 0;
                    SMALL_METRICS(qap_delta = qap_before - qap_after;)

                    small_stat_collect.add_refinement(config->name, get_seconds(sp, ep), qap_delta);

                    SMALL_METRICS(qap_before = qap_after;)
                }
            }

            const auto ep_refinement = std::chrono::high_resolution_clock::now();
            small_stat_collect.add("refinement", get_seconds(sp_refinement, ep_refinement));

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, ac.k));
            METRICS(stat_collect.set_refinement_time(get_seconds(sp_refinement, ep_refinement), level);)
            METRICS(stat_collect.set_refinement_stats(level, get_qap(graphs.back(), p_manager, d_oracle));)
        }
    };
}

#endif //HEIPROMAP_SOLVER_H
