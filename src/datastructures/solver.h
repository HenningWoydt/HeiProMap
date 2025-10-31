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
#include "../definitions.h"
#include "../utility/macros.h"
#include "../utility/random_engine.h"
#include "../utility/utils.h"
#include "../coarsening/global_path_algorithm.h"
#include "../coarsening/greedy_edge_matcher.h"
#include "../coarsening/heavy_edge_matcher.h"
#include "../coarsening/random_edge_matcher.h"
#include "../rebalance/rebalancer.h"
#include "../partitioning/global_multisection.h"
#include "../partitioning/kaffpa_partitioner.h"
#include "../refinement/flow_based_refinement.h"
#include "../refinement/hierarchy_aware_multi_try_multi_way_fm_refinement.h"
#include "../refinement/three_vertex_label_propagation_refinement.h"
#include "../refinement/two_vertex_label_propagation_refinement.h"
#include "../utility/algorithm_configuration.h"
#include "../utility/assert_state.h"
#include "../utility/qap.h"

namespace HeiProMap {
    class Solver {
        AlgorithmConfiguration ac;
        RandomEngine random_engine;

        // statistics
        s64 initial_qap = 0;
        weight_t initial_max_block_weight = 0;

        std::vector<graph_t> graphs;
        PartitionManager p_manager;
        BoundaryVertexManager bv_manager;
        QuotientGraph q_graph;
        DistanceOracle d_oracle;

        // balance
        weight_t lmax = 0;

        // matching
        std::vector<Mapping> mappings;
        GreedyEdgeMatcher ge_matcher;
        HeavyEdgeMatcher he_matcher;
        RandomEdgeMatcher rnd_matcher;
        GlobalPathAlgorithmMatcher gpa_matcher;

        Rebalancer rebalancer;

        // refinement
        LabelPropagationRefinement lp_refine;
        TwoVertexLabelPropagationRefinement two_vertex_lp_refine;
        ThreeVertexLabelPropagationRefinement three_vertex_lp_refine;
        QuotientGraphRefinement qg_refine;
        KWayFMRefinement k_way_refine;
        MultiTryFMRefinement multi_try_fm_refinement;
        FlowBasedRefinement flow_based_refinement;
        // ILPRefinement ilp_refinement;

        HierarchyAwareMultiWayFMRefinement hierarchy_aware_fm_refinement;
        HierarchyAwareMultiTryMultiWayFMRefinement hierarchy_aware_multi_try_multi_way_fm_refinement;
        // HierarchyAwareILPRefinement hierarchy_aware_ilp_refinement;
        HierarchyAwareFlowBasedRefinement hierarchy_aware_flow_based_refinement;
        HierarchyAwareQuotientGraphRefinement hierarchy_aware_quotient_graph_refinement;

        WaveRefinement wave_refinement;
        LightningRefinement lightning_refinement;

        std::vector<std::pair<ISerialRefiner *, ISerialRefinerConfiguration *> > refinements;

        std::vector<std::pair<ISerialRefiner *, ISerialRefinerConfiguration *> > hierarchy_refinements;

    public:
        explicit Solver(const AlgorithmConfiguration &t_ac) {
            ac = t_ac;
            random_engine = RandomEngine(ac.seed);

            graphs.emplace_back(ac.graph_in);

            // balance
            lmax = std::ceil((1.0 + ac.imbalance) * ((f64) graphs[0].weight() / (f64) ac.k));

            // manager
            p_manager.initialize(graphs[0].get_n(), ac.k, lmax);
            bv_manager.initialize(graphs[0].get_n(), ac.k);
            q_graph.initialize(ac.k);
            HEAVYASSERT(assert_state_pre_partitioning(graphs[0], p_manager, ac.k));

            // distance
            d_oracle.initialize(ac.hierarchy, ac.distance);

            // matching
            ge_matcher.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, random_engine, ac.greedy_edge_matcher_config);
            he_matcher.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, random_engine, ac.heavy_edge_matcher_config);
            rnd_matcher.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, random_engine, ac.random_edge_matcher_config);
            gpa_matcher.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, random_engine, ac.global_path_algorithm_config);

            rebalancer.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, ac.hierarchy, ac.distance, random_engine);

            // refinement
            // refinements.emplace_back(&lightning_refinement, &ac.lightning_refinement_configuration);
            refinements.emplace_back(&lp_refine, &ac.label_propagation_config);
            refinements.emplace_back(&k_way_refine, &ac.k_way_fm_refinement_config);
            refinements.emplace_back(&qg_refine, &ac.quotient_graph_refinement_config);
            refinements.emplace_back(&flow_based_refinement, &ac.flow_based_refinement_config);
            refinements.emplace_back(&two_vertex_lp_refine, &ac.two_vertex_label_propagation_config);
            refinements.emplace_back(&three_vertex_lp_refine, &ac.three_vertex_label_propagation_config);
            refinements.emplace_back(&multi_try_fm_refinement, &ac.multi_try_fm_refinement_config);

            hierarchy_refinements.emplace_back(&hierarchy_aware_fm_refinement, &ac.hierarchy_aware_multi_way_fm_config);
            hierarchy_refinements.emplace_back(&hierarchy_aware_multi_try_multi_way_fm_refinement, &ac.hierarchy_aware_multi_try_multi_way_fm_config);
            hierarchy_refinements.emplace_back(&hierarchy_aware_flow_based_refinement, &ac.hierarchy_aware_flow_based_refinement_configuration);
            hierarchy_refinements.emplace_back(&hierarchy_aware_quotient_graph_refinement, &ac.hierarchy_aware_quotient_graph_refinement_configuration);
            // hierarchy_refinements.emplace_back(&ilp_refinement, &ac.ilp_refinement_configuration);
            // hierarchy_refinements.emplace_back(&hierarchy_aware_ilp_refinement, &ac.hierarchy_aware_ilp_refinement_configuration);

            refinements.emplace_back(&wave_refinement, &ac.wave_refinement_configuration);
            refinements.emplace_back(&lightning_refinement, &ac.lightning_refinement_configuration);

            for (auto &[refiner, config]: refinements) {
                if (config->enabled) {
                    refiner->initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, ac.imbalance, lmax, ac.hierarchy, ac.distance, random_engine, *config);
                }
            }

            for (auto &[refiner, config]: hierarchy_refinements) {
                if (config->enabled) {
                    refiner->initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, ac.imbalance, lmax, ac.hierarchy, ac.distance, random_engine, *config);
                }
            }
        }

        std::vector<vertex_t> solve() {
            const auto sp = std::chrono::high_resolution_clock::now();
            internal_solve();

            weight_t qap = get_qap(graphs.back(), p_manager, d_oracle);

            std::vector<partition_t> p(graphs.back().get_n());
            for (vertex_t u = 0; u < graphs.back().get_n(); ++u) { p[u] = p_manager[u]; }
            write_partition(p, ac.mapping_out);

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
            for (u64 v_cycle = 0; v_cycle < ac.n_v_cycle; ++v_cycle) {
                u64 level = 0;
                u64 max_level = 0;

                u64 mult = 64;
                if (v_cycle > 0) { mult = 16; }

                while (graphs.back().get_n() > ac.k * mult) {
                    matching(level, v_cycle);
                    if (mappings.back().get_coarse_n() >= 0.9 * graphs.back().get_n()) {
                        mappings.pop_back();
                        break;
                    }
                    coarsening(level);

                    level += 1;
                }

                max_level = level - 1;
                partition(v_cycle);

                while (level > 0) {
                    level -= 1;
                    uncoarsening(level);
                    rebalancing(level);
                    refinement(level, max_level);
                }
            }
        }

        void partition(u64 v_cycle) {
            for (u64 iteration = 0; iteration < 1; ++iteration) {
                if (v_cycle == 0) {
                    if (ac.partitioning_algorithm_id == PARTITIONING_ALG_KAFFPA) {
                        KaffpaPartitioner partitioner;
                        partitioner.partition(graphs.back(), p_manager, ac.hierarchy, ac.distance, ac.imbalance, random_engine, ac.kaffpa_partitioner_config);
                    } else if (ac.partitioning_algorithm_id == PARTITIONING_ALG_MULTISECTION) {
                        GlobalMultisectionPartitioner partitioner;
                        partitioner.partition(graphs.back(), p_manager, ac.hierarchy, ac.distance, ac.imbalance, random_engine, ac.global_multisection_config);
                    } else {
                        std::cerr << "Partitioning algorithm " << partitioning_algorithm_to_string(ac.partitioning_algorithm_id) << " with id " << ac.partitioning_algorithm_id << " not known!" << std::endl;
                        exit(EXIT_FAILURE);
                    }
                }

                // initialize boundary vertices and quotient graph
                p_manager.reset_weights();
                bv_manager.reset();
                q_graph.initialize(ac.k);
                forall_gu(graphs.back(), u)
                    {
                        const partition_t u_id = p_manager[u];
                        const weight_t u_w = graphs.back().weight(u);
                        p_manager.set(u, u_w, u_id);

                        forall_guivw(graphs.back(), u, i, v, w) {
                                const partition_t v_id = p_manager[v];

                                if (u_id != v_id) {
                                    bv_manager.add(u, u_id); // boundary vertex
                                    if (u < v) {
                                        q_graph.add_edge(u_id, v_id, w); // quotient graph
                                    }
                                }
                            }
                        endfor
                    }
                endfor

                initial_qap = get_qap(graphs.back(), p_manager, d_oracle);
                initial_max_block_weight = max(p_manager.get_bweights());
            }

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, ac.k));
        }

        void matching(const u64 level, u64 v_cycle) {
            mappings.emplace_back();
            mappings.back().initialize(graphs.back().get_n());

            if (v_cycle == 0) {
                if (ac.coarsening_algorithm_id == COARSENING_ALG_GREEDY_MATCHING) {
                    ge_matcher.match(level, graphs.back(), p_manager, mappings.back());
                } else if (ac.coarsening_algorithm_id == COARSENING_ALG_HEAVY_MATCHING) {
                    he_matcher.match(level, graphs.back(), p_manager, mappings.back());
                } else if (ac.coarsening_algorithm_id == COARSENING_ALG_RANDOM_MATCHING) {
                    rnd_matcher.match(level, graphs.back(), p_manager, mappings.back());
                } else if (ac.coarsening_algorithm_id == COARSENING_ALG_GLOBAL_PATHS) {
                    gpa_matcher.match(level, graphs.back(), p_manager, mappings.back());
                } else {
                    std::cerr << "Coarsening algorithm " << coarsening_algorithm_to_string(ac.coarsening_algorithm_id) << " with id " << ac.coarsening_algorithm_id << " not known!" << std::endl;
                    exit(EXIT_FAILURE);
                }
            } else {
                rnd_matcher.match(level, graphs.back(), p_manager, mappings.back());
            }
        }

        void coarsening([[maybe_unused]] const u64 level) {
            HEAVYASSERT(assert_state_pre_partitioning(graphs.back(), p_manager, ac.k));

            graphs.emplace_back(); // coarse the graph
            graphs.back().initialize(graphs[graphs.size() - 2], mappings.back());
            p_manager.contract(mappings.back());

            HEAVYASSERT(assert_state_pre_partitioning(graphs.back(), p_manager, ac.k));
        }

        void uncoarsening([[maybe_unused]] const u64 level) {
            p_manager.uncontract(mappings.back());
            bv_manager.compute_from_scratch(graphs[graphs.size() - 2], p_manager);

            ScopedTimer _t("uncontraction", "misc", "free_graph");
            graphs.pop_back(); // this is doing uncontraction
            mappings.pop_back();
            _t.stop();

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, ac.k));
        }

        void rebalancing([[maybe_unused]] const u64 level) {
            if (level == 0) {
                rebalancer.rebalance_last_layer(graphs.back(), p_manager, bv_manager, q_graph, d_oracle);
            } else {
                rebalancer.rebalance(graphs.back(), p_manager, bv_manager, q_graph, d_oracle);
            }

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, ac.k));
        }

        void refinement(const u64 level, const u64 max_level) {
            u64 refinement_max_iterations = 1;
            for (u64 refinement_i = 0; refinement_i < refinement_max_iterations; ++refinement_i) {
                for (auto [refiner, config]: refinements) {
                    if (config->enabled) {
                        refiner->refine(level, max_level, graphs.back(), d_oracle, bv_manager, p_manager, q_graph);
                        HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, ac.k));
                    }
                }
            }

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, ac.k));
        }
    };
}

#endif //HEIPROMAP_SOLVER_H
