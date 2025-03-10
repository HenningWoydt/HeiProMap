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

#include "active_vertex_manager.h"
#include "boundary_vertex_manger.h"
#include "graph.h"
#include "partition_manager.h"
#include "quotient_graph.h"
#include "../../definitions.h"
#include "../../macros.h"
#include "../../commons/random_engine.h"
#include "../../commons/statistic_collector.h"
#include "../coarsening/global_path_algorithm.h"
#include "../coarsening/greedy_edge_matcher.h"
#include "../coarsening/heavy_edge_matcher.h"
#include "../partitioning/global_multisection.h"
#include "../partitioning/kaffpa_partitioner.h"
#include "../rebalance/simple_rebalancer.h"
#include "../refinement/k_way_fm_refinement_boost.h"
#include "../refinement/label_propagation_refinement_Faraj20.h"
#include "../refinement/quotient_graph_refinement_boost.h"
#include "../refinement/quotient_graph_refinement_Faraj20.h"
#include "../refinement/quotient_graph_refinement_simple_pq.h"
#include "../refinement/two_vertex_label_propagation_refinement.h"
#include "../utility/AlgorithmConfiguration.h"
#include "../utility/assert_state.h"
#include "../utility/qap.h"
#include "../../commons/utils.h"

namespace HeiProMap {
    /**
     * Solver for serial Process Mapping.
     */
    class Solver {
        AlgorithmConfiguration ac;
        RandomEngine random_engine;

        // statistics
        StatisticCollector stat_collect;

        std::vector<graph_t> graphs;
        ActiveVertexManager av_manager;
        PartitionManager p_manager;
        BoundaryVertexManager bv_manager;
        QuotientGraph q_graph;
        DistanceOracle d_oracle;

        // balance
        weight_t lmax = 0;

        // matching
        std::vector<EdgeUV*> matches;
        std::vector<size_t> matches_size;
        GreedyEdgeMatcher ge_matcher;
        HeavyEdgeMatcher he_matcher;
        GlobalPathAlgorithmMatcher gpa_matcher;

        // refinement
        LabelPropagationRefinementFaraj20 lp_refine_faraj20;
        LabelPropagationRefinement lp_refine;
        TwoVertexLabelPropagationRefinement two_vertex_lp_refine;
        QuotientGraphRefinementFaraj20 qg_refine_faraj20;
        QuotientGraphRefinement qg_refine;
        KWayFMRefinementFaraj20 k_way_refine_faraj20;
        KWayFMRefinement k_way_refine;
        MultiTryFMRefinementFaraj20 multi_try_fm_refinement_faraj20;
        MultiTryFMRefinement multi_try_fm_refinement;
        HierarchyAwareCycleRefinement hierarchy_aware_cycle_refinement;

    public:
        explicit Solver(const AlgorithmConfiguration& t_ac) {
            ac            = t_ac;
            random_engine = RandomEngine(ac.seed);

            METRICS_TIME(sp_graph_io);
            graphs.emplace_back(ac.graph_in);
            METRICS_TIME(ep_graph_io);

            METRICS_TIME(sp_io);
            // balance
            lmax = std::ceil((1.0 + ac.imbalance) * ((f64)graphs[0].get_weight() / (f64)ac.k));

            // manager
            av_manager.initialize(graphs[0].get_n());
            p_manager.initialize(graphs[0].get_n(), ac.k, lmax);
            bv_manager.initialize(graphs[0].get_n(), ac.k);
            q_graph.initialize(ac.k);
            HEAVYASSERT(assert_state_pre_partitioning(graphs[0], av_manager));

            // distance
            d_oracle.initialize(ac.hierarchy, ac.distance);

            // matching
            ge_matcher.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, random_engine, ac.greedy_edge_matcher_config, stat_collect);
            he_matcher.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, random_engine, ac.heavy_edge_matcher_config, stat_collect);
            gpa_matcher.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, random_engine, ac.global_path_algorithm_config, stat_collect);

            // refinement
            lp_refine_faraj20.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, ac.hierarchy, ac.distance, random_engine, ac.label_propagation_faraj20_config, stat_collect);
            lp_refine.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, ac.hierarchy, ac.distance, random_engine, ac.label_propagation_config, stat_collect);
            two_vertex_lp_refine.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, ac.hierarchy, ac.distance, random_engine, ac.two_vertex_label_propagation_config, stat_collect);
            qg_refine_faraj20.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, ac.hierarchy, ac.distance, random_engine, ac.quotient_graph_refinement_faraj20_config, stat_collect);
            qg_refine.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, ac.hierarchy, ac.distance, random_engine, ac.quotient_graph_refinement_config, stat_collect);
            k_way_refine_faraj20.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, ac.hierarchy, ac.distance, random_engine, ac.k_way_fm_refinement_faraj20_config, stat_collect);
            k_way_refine.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, ac.hierarchy, ac.distance, random_engine, ac.k_way_fm_refinement_config, stat_collect);
            multi_try_fm_refinement_faraj20.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, ac.hierarchy, ac.distance, random_engine, ac.multi_try_fm_refinement_faraj20_config, stat_collect);
            multi_try_fm_refinement.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, ac.hierarchy, ac.distance, random_engine, ac.multi_try_fm_refinement_config, stat_collect);
            hierarchy_aware_cycle_refinement.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, ac.hierarchy, ac.distance, random_engine, ac.hierarchy_aware_cycles_config, stat_collect);

            METRICS_TIME(ep_io);
            METRICS(stat_collect.set_io(get_seconds(sp_graph_io, ep_graph_io), get_seconds(sp_io, ep_io)));
        }

        std::vector<vertex_t> solve() {
            const auto sp = std::chrono::high_resolution_clock::now();

            internal_solve();

            weight_t qap = get_qap(graphs.back(), av_manager, p_manager, d_oracle);

            METRICS(stat_collect.set_final(qap, p_manager.get_bweights(), lmax));
            METRICS(stat_collect.finalize());
            METRICS(std::ofstream stat_file(ac.statistics_out));
            METRICS(stat_file << stat_collect.to_JSON());

            std::vector<partition_t> p(graphs.back().get_n());
            for (vertex_t u = 0; u < graphs.back().get_n(); ++u) { p[u] = p_manager[u]; }
            write_partition(p, ac.mapping_out);

            const auto ep = std::chrono::high_resolution_clock::now();
            std::cout << "Total time : " << get_seconds(sp, ep) << std::endl;
            std::cout << "Final QAP  : " << qap << std::endl;
            std::cout << "Lmax       : " << lmax << std::endl;
            std::cout << "max block w: " << max(p_manager.get_bweights()) << std::endl;

            return p;
        }

    private:
        void internal_solve() {
            u64 level = 0;

            while (av_manager.size() > ac.k * 64) {
                matching(level);
                if (matches_size.back() == 0) { break; }

                coarsening(level);

                level += 1;
            }

            partition();

            while (level > 0) {
                level -= 1;
                uncoarsening(level);
                refinement(level);
            }
        }

        void partition() {
            METRICS_TIME(sp_partition)

            if (ac.partitioning_algorithm_id == PARTITIONING_ALG_KAFFPA) {
                KaffpaPartitioner partitioner;
                partitioner.partition(graphs.back(), av_manager, p_manager, ac.hierarchy, ac.distance, ac.imbalance, random_engine, ac.kaffpa_partitioner_config, stat_collect);
            } else if (ac.partitioning_algorithm_id == PARTITIONING_ALG_MULTISECTION) {
                GlobalMultisectionPartitioner partitioner;
                partitioner.partition(graphs.back(), av_manager, p_manager, ac.hierarchy, ac.distance, ac.imbalance, random_engine, ac.global_multisection_config, stat_collect);
            } else {
                std::cerr << "Partitioning algorithm " << partitioning_algorithm_to_string(ac.partitioning_algorithm_id) << " with id " << ac.partitioning_algorithm_id << " not known!" << std::endl;
                exit(EXIT_FAILURE);
            }

            if (ac.rebalancing_algorithm_id == REBALANCING_ALG_SIMPLE) {
                SimpleRebalancer simple_rebalancer(graphs[0].get_n(), ac.k, lmax);
                simple_rebalancer.rebalance(ac.simple_rebalancer_configuration, graphs.back(), av_manager, bv_manager, p_manager, d_oracle, q_graph);
            } else {
                std::cerr << "Rebalancing algorithm " << rebalancing_algorithm_to_string(ac.rebalancing_algorithm_id) << " with id " << ac.rebalancing_algorithm_id << " not known!" << std::endl;
                exit(EXIT_FAILURE);
            }

            // initialize boundary vertices and quotient graph
            forall_av_iu(av_manager, j, u)
                {
                    for (size_t i = 0; i < graphs.back().size(u); ++i) {
                        const vertex_t v       = graphs.back().neighbor(u, i);
                        const weight_t w       = graphs.back().get_weight(u, i);
                        const partition_t u_id = p_manager[u];
                        const partition_t v_id = p_manager[v];

                        if (u_id != v_id) {
                            bv_manager.add(u, u_id); // boundary vertex
                            q_graph.add_edge(u_id, v_id, w); // quotient graph
                        }
                    }
                }
            endfor

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), av_manager, p_manager, bv_manager, q_graph, ac.k));

            METRICS_TIME(ep_partition);
            METRICS(stat_collect.set_partition_time(get_seconds(sp_partition, ep_partition));)
            METRICS(stat_collect.set_partition_stats(get_qap(graphs.back(), av_manager, p_manager, d_oracle), p_manager.get_bweights(), lmax);)
        }

        void matching(const u64 level) {
            METRICS_TIME(sp_match)

            size_t t_n_64       = round_up_64(av_manager.size() / 2);
            EdgeUV* matches_arr = (EdgeUV*)aligned_alloc(64, t_n_64 * sizeof(EdgeUV));
            matches.push_back(matches_arr);
            matches_size.push_back(0);

            if (ac.coarsening_algorithm_id == COARSENING_ALG_GREEDY_MATCHING) {
                ge_matcher.match(level, graphs.back(), av_manager, matches.back(), matches_size.back());
            } else if (ac.coarsening_algorithm_id == COARSENING_ALG_HEAVY_MATCHING) {
                he_matcher.match(level, graphs.back(), av_manager, matches.back(), matches_size.back());
            } else if (ac.coarsening_algorithm_id == COARSENING_ALG_GLOBAL_PATHS) {
                gpa_matcher.match(level, graphs.back(), av_manager, matches.back(), matches_size.back());
            } else {
                std::cerr << "Coarsening algorithm " << coarsening_algorithm_to_string(ac.coarsening_algorithm_id) << " with id " << ac.coarsening_algorithm_id << " not known!" << std::endl;
                exit(EXIT_FAILURE);
            }

            METRICS_TIME(ep_match)
            METRICS(stat_collect.set_matching_time(get_seconds(sp_match, ep_match), level);)
            METRICS(stat_collect.set_matching_stats(level, matches_size.back());)
        }

        void coarsening(const u64 level) {
            METRICS_TIME(sp_coarse)

            graphs.emplace_back(graphs.back(), matches.back(), matches_size.back()); // coarse the graph
            av_manager.contract(matches.back(), matches_size.back());

            HEAVYASSERT(assert_state_pre_partitioning(graphs.back(), av_manager));

            METRICS_TIME(ep_coarse)
            METRICS(stat_collect.set_coarsening_time(get_seconds(sp_coarse, ep_coarse), level);)
            METRICS(stat_collect.set_coarsening_stats(av_manager.size(), level);)
        }

        void uncoarsening(const u64 level) {
            METRICS_TIME(sp_uncoarse)

            p_manager.uncontract(matches.back(), matches_size.back());
            av_manager.uncontract(matches.back(), matches_size.back());
            bv_manager.uncontract(graphs[graphs.size() - 2], av_manager, p_manager);
            graphs.pop_back(); // this is doing uncontraction

            free(matches.back());
            matches.pop_back(); // throw away the matching, not needed anymore
            matches_size.pop_back();

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), av_manager, p_manager, bv_manager, q_graph, ac.k));

            METRICS_TIME(ep_uncoarse)
            METRICS(stat_collect.set_uncoarsening_time(get_seconds(sp_uncoarse, ep_uncoarse), level);)
            METRICS(stat_collect.set_uncoarsening_stats(level, av_manager.size());)
        }

        void refinement(const u64 level) {
            METRICS_TIME(sp_refinement)

            if (ac.do_refinement_quotient_graph_faraj20) {
                qg_refine_faraj20.refine(level, graphs.back(), av_manager, d_oracle, bv_manager, p_manager, q_graph);
            }

            if (ac.do_refinement_quotient_graph) {
                qg_refine.refine(level, graphs.back(), av_manager, d_oracle, bv_manager, p_manager, q_graph);
            }

            if (ac.do_refinement_label_propagation_faraj20) {
                lp_refine_faraj20.refine(level, graphs.back(), av_manager, d_oracle, bv_manager, p_manager, q_graph);
            }

            if (ac.do_refinement_label_propagation) {
                lp_refine.refine(level, graphs.back(), av_manager, d_oracle, bv_manager, p_manager, q_graph);
            }

            if (ac.do_refinement_two_vertex_label_propagation_enable) {
                two_vertex_lp_refine.refine(level, graphs.back(), av_manager, d_oracle, bv_manager, p_manager, q_graph);
            }

            if (ac.do_refinement_k_way_fm_faraj20) {
                k_way_refine_faraj20.refine(level, graphs.back(), av_manager, d_oracle, bv_manager, p_manager, q_graph);
            }

            if (ac.do_refinement_k_way_fm) {
                k_way_refine.refine(level, graphs.back(), av_manager, d_oracle, bv_manager, p_manager, q_graph);
            }

            if (ac.do_refinement_multi_try_fm_faraj20) {
                multi_try_fm_refinement_faraj20.refine(level, graphs.back(), av_manager, d_oracle, bv_manager, p_manager, q_graph);
            }

            if (ac.do_refinement_multi_try_fm) {
                multi_try_fm_refinement.refine(level, graphs.back(), av_manager, d_oracle, bv_manager, p_manager, q_graph);
            }

            if (ac.do_refinement_hierarchy_aware_cycles_enable) {
                hierarchy_aware_cycle_refinement.refine(level, graphs.back(), av_manager, d_oracle, bv_manager, p_manager, q_graph);
            }

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), av_manager, p_manager, bv_manager, q_graph, ac.k));

            METRICS_TIME(ep_refinement)
            METRICS(stat_collect.set_refinement_time(get_seconds(sp_refinement, ep_refinement), level);)
            METRICS(stat_collect.set_refinement_stats(level, get_qap(graphs.back(), av_manager, p_manager, d_oracle));)
        }
    };
}

#endif //HEIPROMAP_SOLVER_H
