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
#include "block_conn.h"
#include "../definitions.h"
#include "../utility/macros.h"
#include "../utility/random_engine.h"
#include "../utility/utils.h"
#include "../coarsening/global_path_algorithm.h"
#include "../coarsening/greedy_edge_matcher.h"
#include "../coarsening/heavy_edge_matcher.h"
#include "../coarsening/random_edge_matcher.h"
#include "../coarsening/size_constrained_lp.h"
#include "../rebalance/rebalancer.h"
#include "../partitioning/global_multisection.h"
#include "../refinement/flow_based_refinement.h"
#include "../refinement/three_vertex_label_propagation_refinement.h"
#include "../refinement/two_vertex_label_propagation_refinement.h"
#include "../refinement/pertubator.h"
#include "../utility/algorithm_configuration.h"
#include "../utility/assert_state.h"
#include "../utility/qap.h"

namespace HeiProMap {
    class Solver {
        AlgorithmConfiguration ac;
        RandomEngine random_engine;
        f64 init_time = 0.0;

        // statistics
        s64 initial_qap = 0;
        weight_t initial_max_block_weight = 0;

        std::vector<graph_t> graphs;

        u64 max_n_partitions = 1;
        std::vector<PartitionManager> p_managers;
        std::vector<BoundaryVertexManager> bv_managers;
        std::vector<QuotientGraph> q_graphs;
        std::vector<BlockConn> block_conns;
        DistanceOracle d_oracle;

        // matching
        std::vector<Mapping> mappings;
        GreedyEdgeMatcher ge_matcher;
        HeavyEdgeMatcher he_matcher;
        RandomEdgeMatcher rnd_matcher;
        GlobalPathAlgorithmMatcher gpa_matcher;
        SizeConstrainedLP size_constrained_lp;

        Rebalancer rebalancer;

        // refinement
        LabelPropagationRefinement lp_refine;
        // TwoVertexLabelPropagationRefinement two_vertex_lp_refine;
        // ThreeVertexLabelPropagationRefinement three_vertex_lp_refine;
        QuotientGraphRefinement qg_refine;
        KWayFMRefinement k_way_refine;
        MultiTryFMRefinement multi_try_fm_refinement;
        FlowBasedRefinement flow_based_refinement;
        // ILPRefinement ilp_refinement;

        // WaveRefinement wave_refinement;
        // LightningRefinement lightning_refinement;

        std::vector<std::pair<ISerialRefiner *, ISerialRefinerConfiguration *> > refinements;

        f64 io_ms = 0.0;
        f64 misc_ms = 0.0;
        f64 coarsening_ms = 0.0;
        f64 contraction_ms = 0.0;
        f64 initial_partitioning_ms = 0.0;
        f64 uncontraction_ms = 0.0;
        f64 rebalance_ms = 0.0;
        f64 refinement_ms = 0.0;

        struct level_info {
            u32 level = -1;
            vertex_t n = -1;
            vertex_t m = -1;

            weight_t comm_cost = -1;
            weight_t max_b_weight = -1;
            weight_t lmax = -1;
            f64 imb = -1.0;
            partition_t empty_partitions = -1;
            partition_t oload_partitions = -1;
            weight_t sum_oload_weights = -1;

            f64 t_coarsening = 0.0;
            f64 t_contraction = 0.0;
            f64 t_uncontraction = 0.0;
            f64 t_rebalance = 0.0;
            f64 t_refinement = 0.0;
        };

        std::vector<std::vector<level_info> > level_infos;

        inline void print_level_row(const level_info &L) {
            std::cout
                    << std::setw(3) << L.level << " | "
                    << std::setw(8) << L.n << " | "
                    << std::setw(11) << L.m << " | "
                    << std::setw(10) << L.comm_cost << " | "
                    << std::setw(7) << L.lmax << " | "
                    << std::setw(7) << L.max_b_weight << " | "
                    << std::setw(8) << L.imb << " | "
                    << std::setw(6) << (u32) L.empty_partitions << " | "
                    << std::setw(6) << (u32) L.oload_partitions << " | "
                    << std::setw(8) << L.sum_oload_weights << " | "
                    << std::setw(10) << L.t_coarsening << " | "
                    << std::setw(10) << L.t_contraction << " | "
                    << std::setw(10) << L.t_uncontraction << " | "
                    << std::setw(10) << L.t_rebalance << " | "
                    << std::setw(10) << L.t_refinement
                    << "\n";
        }

        inline void print_all_levels(const std::vector<level_info> &infos) {
            std::cout
                    << std::setw(3) << "Lvl" << " | "
                    << std::setw(8) << "n" << " | "
                    << std::setw(11) << "m" << " | "
                    << std::setw(10) << "comm cost" << " | "
                    << std::setw(7) << "lmax" << " | "
                    << std::setw(7) << "maxW" << " | "
                    << std::setw(8) << "imb" << " | "
                    << std::setw(6) << "empty" << " | "
                    << std::setw(6) << "oload" << " | "
                    << std::setw(8) << "w_oload" << " | "
                    << std::setw(10) << "t_c" << " | "
                    << std::setw(10) << "t_con" << " | "
                    << std::setw(10) << "t_unc" << " | "
                    << std::setw(10) << "t_reb" << " | "
                    << std::setw(10) << "t_ref"
                    << "\n";

            std::cout << std::string(100, '-') << "\n";
            for (const auto &L: infos) {
                print_level_row(L);
            }
        }

    public:
        explicit Solver(const AlgorithmConfiguration &t_ac) {
            auto sp_io = get_time_point();
            graphs.emplace_back(t_ac.graph_in);
            io_ms += get_milli_seconds(sp_io, get_time_point());

            auto sp = get_time_point();
            ac = t_ac;
            random_engine = RandomEngine(ac.seed);

            // manager
            max_n_partitions = ac.n_max_partitions;
            p_managers.resize(max_n_partitions);
            bv_managers.resize(max_n_partitions);
            q_graphs.resize(max_n_partitions);
            block_conns.resize(max_n_partitions);
            for (u64 i = 0; i < max_n_partitions; ++i) {
                p_managers[i].initialize(graphs[0].n, ac.k, graphs[0].g_weight);
                bv_managers[i].initialize(graphs[0].n, ac.k);
                q_graphs[i].initialize(ac.k);
                block_conns[i].initialize(graphs[0].n, ac.k);
            }
            HEAVYASSERT(assert_state_pre_partitioning(graphs[0], p_managers[0], ac.k));

            // distance
            d_oracle.initialize(ac.hierarchy, ac.distance);

            // matching
            ge_matcher.initialize(graphs[0].n, graphs[0].m, ac.k, random_engine, ac.greedy_edge_matcher_config);
            he_matcher.initialize(graphs[0].n, graphs[0].m, ac.k, random_engine, ac.heavy_edge_matcher_config);
            rnd_matcher.initialize(graphs[0].n, graphs[0].m, ac.k, random_engine, ac.random_edge_matcher_config);
            gpa_matcher.initialize(graphs[0].n, graphs[0].m, ac.k, random_engine, ac.global_path_algorithm_config);
            size_constrained_lp.initialize(graphs[0].n, graphs[0].m, ac.k, random_engine,
                                           ac.size_constrained_lp_config);

            rebalancer.initialize(graphs[0].n, graphs[0].m, ac.k, random_engine.get_u64());

            // refinement
            // refinements.emplace_back(&lightning_refinement, &ac.lightning_refinement_configuration);
            refinements.emplace_back(&lp_refine, &ac.label_propagation_config);
            refinements.emplace_back(&k_way_refine, &ac.k_way_fm_refinement_config);
            refinements.emplace_back(&qg_refine, &ac.quotient_graph_refinement_config);
            refinements.emplace_back(&flow_based_refinement, &ac.flow_based_refinement_config);
            // refinements.emplace_back(&two_vertex_lp_refine, &ac.two_vertex_label_propagation_config);
            // refinements.emplace_back(&three_vertex_lp_refine, &ac.three_vertex_label_propagation_config);
            refinements.emplace_back(&multi_try_fm_refinement, &ac.multi_try_fm_refinement_config);

            // refinements.emplace_back(&wave_refinement, &ac.wave_refinement_configuration);
            // refinements.emplace_back(&lightning_refinement, &ac.lightning_refinement_configuration);

            for (auto &[refiner, config]: refinements) {
                if (config->enabled) {
                    refiner->initialize(graphs[0].n, graphs[0].m, ac.k, *config);
                }
            }
            lp_refine.initialize(graphs[0].n, graphs[0].m, ac.k, ac.label_propagation_config);

            auto ep = get_time_point();
            init_time += get_seconds(sp, ep);
        }

        std::vector<vertex_t> solve() {
            const auto sp = std::chrono::high_resolution_clock::now();
            internal_solve();

            weight_t qap = get_qap(graphs.back(), p_managers[0], d_oracle);

            std::vector<partition_t> p(graphs.back().n);
            for (vertex_t u = 0; u < graphs.back().n; ++u) { p[u] = p_managers[0][u]; }
            write_partition(p, ac.mapping_out);

            const auto ep = std::chrono::high_resolution_clock::now();
            f64 duration = get_seconds(sp, ep);

            weight_t lmax = std::ceil((1.0 + ac.imbalance) * ((f64) graphs[0].g_weight / (f64) ac.k));

            std::cout << "Total time        : " << duration + init_time << std::endl;
            std::cout << "#Nodes            : " << graphs.back().n << std::endl;
            std::cout << "#Edges            : " << graphs.back().m << std::endl;
            std::cout << "k                 : " << ac.k << std::endl;
            std::cout << "Lmax              : " << lmax << std::endl;
            std::cout << "Init. QAP         : " << initial_qap << std::endl;
            std::cout << "Init. max block w : " << initial_max_block_weight << std::endl;
            std::cout << "Final QAP         : " << qap << std::endl;
            std::cout << "max block w       : " << max(p_managers[0].get_bweights()) << std::endl;

            size_t n_empty_partitions = 0;
            size_t n_overloaded_partitions = 0;
            weight_t sum_too_much = 0;
            for (partition_t id = 0; id < ac.k; ++id) {
                n_empty_partitions += p_managers[0].get_bweight(id) == 0;
                n_overloaded_partitions += p_managers[0].get_bweight(id) > lmax;
                sum_too_much += std::max((weight_t) 0, p_managers[0].get_bweight(id) - lmax);
            }
            std::cout << "#empty partitions : " << n_empty_partitions << std::endl;
            std::cout << "#oload partitions : " << n_overloaded_partitions << std::endl;
            std::cout << "Sum oload weights : " << sum_too_much << std::endl;
            std::cout << "IO                : " << io_ms << std::endl;
            std::cout << "Misc              : " << misc_ms << std::endl;
            std::cout << "Coarsening        : " << coarsening_ms << std::endl;
            std::cout << "Contraction       : " << contraction_ms << std::endl;
            std::cout << "Init. Part.       : " << initial_partitioning_ms << std::endl;
            std::cout << "Uncontraction     : " << uncontraction_ms << std::endl;
            std::cout << "Rebalance         : " << rebalance_ms << std::endl;
            std::cout << "Refinement        : " << refinement_ms << std::endl;
            std::cout << "ALL               : " << io_ms + misc_ms + coarsening_ms + contraction_ms +
                    initial_partitioning_ms + uncontraction_ms + rebalance_ms + refinement_ms << std::endl;

#if ENABLE_PROFILER
            for (size_t i = 0; i < level_infos.size(); ++i) {
                print_all_levels(level_infos[i]);
            }
#endif

            return p;
        }

    private:
        void internal_solve() {
            u64 level = 0;
            u64 mult = ac.initial_c;

            run_cycle(0, level, mult);

            u64 n_v_cycle = ac.n_v_cycle;
            for (u64 i_v_cycle = 0; i_v_cycle < n_v_cycle; ++i_v_cycle) {
                run_cycle(i_v_cycle + 1, level, mult);
            }
        }

        void run_cycle(const u64 v_cycle, u64 &level, const u64 mult) {
            const bool is_initial = (v_cycle == 0);

            f64 level_imbalance = 0.0;
            f64 per_level_imb_add = is_initial ? 1.0 / 400.0 : 0.0;
            weight_t level_lmax = 0;

            level_infos.emplace_back();

            while (graphs.back().n > ac.k * mult && (is_initial || graphs.size() <= ac.v_cycle_max_depth)) {
#if ENABLE_PROFILER
                level_infos[v_cycle].emplace_back();
                level_infos[v_cycle][level].level = level;
                level_infos[v_cycle][level].n = graphs.back().n;
                level_infos[v_cycle][level].m = graphs.back().m;
#endif

                level_imbalance = ac.imbalance + (f64) level * per_level_imb_add;
                level_lmax = std::ceil((1.0 + level_imbalance) * ((f64) graphs[0].g_weight / (f64) ac.k));

                coarsening(v_cycle, level, level_imbalance, !is_initial);
                contraction(v_cycle, level);

                level += 1;
            }

#if ENABLE_PROFILER
            level_infos[v_cycle].emplace_back();
            level_infos[v_cycle][level].level = level;
            level_infos[v_cycle][level].n = graphs.back().n;
            level_infos[v_cycle][level].m = graphs.back().m;
#endif

            level_imbalance = ac.imbalance + (f64) level * per_level_imb_add;
            level_lmax = std::ceil((1.0 + level_imbalance) * ((f64) graphs[0].g_weight / (f64) ac.k));

            if (is_initial) {
                partition(level, level_imbalance);
            } else {
                block_conns[0].compute_from_scratch(graphs.back(), p_managers[0]);
            }

#if ENABLE_PROFILER
            level_infos[v_cycle][level].max_b_weight = p_managers[0].max_weight();
            level_infos[v_cycle][level].lmax = level_lmax;
            level_infos[v_cycle][level].imb = (f64) level_infos[v_cycle][level].max_b_weight / (
                                                  (f64) graphs[0].g_weight / (f64) ac.k);
            level_infos[v_cycle][level].comm_cost = get_qap(graphs.back(), p_managers[0], d_oracle);
            level_infos[v_cycle][level].empty_partitions = p_managers[0].n_empty_blocks();
            level_infos[v_cycle][level].oload_partitions = p_managers[0].n_oload_blocks(level_lmax);
            level_infos[v_cycle][level].sum_oload_weights = p_managers[0].sum_oload_weight(level_lmax);
#endif

            while (!mappings.empty()) {
                level -= 1;

                level_imbalance = ac.imbalance + (f64) level * per_level_imb_add;
                level_lmax = std::ceil((1.0 + level_imbalance) * ((f64) graphs[0].g_weight / (f64) ac.k));
                uncoarsening(v_cycle, level);

                u64 n_partitions = is_initial ? 1 : std::min(1 + level * level, max_n_partitions);
                lp_refine.min_improvement = -1;
                for (u64 i = 1; i < n_partitions; ++i) {
                    p_managers[i].copy_from(p_managers[0]);
                    bv_managers[i].copy_from(bv_managers[0]);
                    q_graphs[i].copy_from(q_graphs[0]);
                    block_conns[i].copy_from(block_conns[0]);

                    // perturbate(graphs.back(), d_oracle, bv_managers[i], p_managers[i], q_graphs[i], block_conns[i], level_imbalance);
                    lp_refine.refine(graphs.back(), d_oracle, bv_managers[i], p_managers[i], q_graphs[i], block_conns[i], level_imbalance);
                }
                lp_refine.min_improvement = 0;

                for (u64 i = 0; i < n_partitions; ++i) {
                    std::cout << "pre: " << v_cycle << " " << level << " " << i << " " << get_qap(
                        graphs.back(), p_managers[i], d_oracle) << std::endl;
                }

                rebalancing(v_cycle, level, level_imbalance, n_partitions);
                refinement(v_cycle, level, level_imbalance, n_partitions);

                for (u64 i = 0; i < n_partitions; ++i) {
                    std::cout << "aft: " << v_cycle << " " << level << " " << i << " " << get_qap(
                        graphs.back(), p_managers[i], d_oracle) << std::endl;
                }

                weight_t best_qap = get_qap(graphs.back(), p_managers[0], d_oracle);
                u64 best_id = 0;
                for (u64 i = 1; i < n_partitions; ++i) {
                    weight_t qap = get_qap(graphs.back(), p_managers[i], d_oracle);
                    if (qap < best_qap) {
                        best_qap = qap;
                        best_id = i;
                    }
                }

                if (best_id != 0) {
                    std::swap(p_managers[0], p_managers[best_id]);
                    std::swap(bv_managers[0], bv_managers[best_id]);
                    std::swap(q_graphs[0], q_graphs[best_id]);
                    std::swap(block_conns[0], block_conns[best_id]);
                }

#if ENABLE_PROFILER
                level_infos[v_cycle][level].max_b_weight = p_managers[0].max_weight();
                level_infos[v_cycle][level].lmax = level_lmax;
                level_infos[v_cycle][level].imb =                        (f64) level_infos[v_cycle][level].max_b_weight / ((f64) graphs[0].g_weight / (f64) ac.k);
                level_infos[v_cycle][level].comm_cost = get_qap(graphs.back(), p_managers[0], d_oracle);
                level_infos[v_cycle][level].empty_partitions = p_managers[0].n_empty_blocks();
                level_infos[v_cycle][level].oload_partitions = p_managers[0].n_oload_blocks(level_lmax);
                level_infos[v_cycle][level].sum_oload_weights = p_managers[0].sum_oload_weight(level_lmax);
#endif
            }
        }

        void partition(u64 level, f64 level_imbalance) {
            auto sp = get_time_point();

            for (u64 iteration = 0; iteration < 1; ++iteration) {
                if (ac.partitioning_algorithm_id == PARTITIONING_ALG_MULTISECTION) {
                    GlobalMultisectionPartitioner partitioner;
                    partitioner.partition(graphs.back(), p_managers[0], ac.hierarchy, ac.distance, level_imbalance, random_engine, ac.global_multisection_config);
                } else {
                    std::cerr << "Partitioning algorithm " << partitioning_algorithm_to_string(ac.partitioning_algorithm_id) << " with id " << ac.partitioning_algorithm_id << " not known!" << std::endl;
                    exit(EXIT_FAILURE);
                }

                // initialize boundary vertices and quotient graph
                ScopedTimer _t_allocate("partition", "misc", "initialize_datastructures");
                p_managers[0].reset_weights();
                bv_managers[0].reset();
                q_graphs[0].initialize(ac.k);
                forall_gu(graphs.back(), u)
                    {
                        const partition_t u_id = p_managers[0][u];
                        const weight_t u_w = graphs.back().v_weights[u];
                        p_managers[0].set(u, u_w, u_id);

                        forall_guivw(graphs.back(), u, i, v, w)
                            {
                                const partition_t v_id = p_managers[0][v];

                                if (u_id != v_id) {
                                    bv_managers[0].add(u, u_id); // boundary vertex
                                    if (u < v) {
                                        q_graphs[0].add_edge(u_id, v_id, w); // quotient graph
                                    }
                                }
                            }
                        endfor
                    }
                endfor
                block_conns[0].compute_from_scratch(graphs.back(), p_managers[0]);
                _t_allocate.stop();

                initial_qap = get_qap(graphs.back(), p_managers[0], d_oracle);
                initial_max_block_weight = max(p_managers[0].get_bweights());
            }
            auto ep = get_time_point();
            initial_partitioning_ms += get_milli_seconds(sp, ep);

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_managers[0], bv_managers[0], q_graphs[0], block_conns[0], ac.k));
        }

        void coarsening(const u64 v_cycle, const u64 level, f64 level_imbalance, bool random = false) {
            auto sp = get_time_point();
            mappings.emplace_back();
            mappings.back().initialize(graphs.back().n);

            if (random) {
                rnd_matcher.match(level, graphs.back(), p_managers[0], mappings.back(), level_imbalance);
            } else if (ac.coarsening_algorithm_id == COARSENING_ALG_GREEDY_MATCHING) {
                ge_matcher.match(level, graphs.back(), p_managers[0], mappings.back(), level_imbalance);
            } else if (ac.coarsening_algorithm_id == COARSENING_ALG_HEAVY_MATCHING) {
                he_matcher.match(level, graphs.back(), p_managers[0], mappings.back(), level_imbalance);
            } else if (ac.coarsening_algorithm_id == COARSENING_ALG_RANDOM_MATCHING) {
                rnd_matcher.match(level, graphs.back(), p_managers[0], mappings.back(), level_imbalance);
            } else if (ac.coarsening_algorithm_id == COARSENING_ALG_GLOBAL_PATHS) {
                gpa_matcher.match(level, graphs.back(), p_managers[0], mappings.back(), level_imbalance);
            } else if (ac.coarsening_algorithm_id == COARSENING_ALG_SIZE_CONSTRAINED_LP) {
                size_constrained_lp.cluster(level, graphs.back(), p_managers[0], mappings.back(), level_imbalance);
            } else {
                std::cerr << "Coarsening algorithm " << coarsening_algorithm_to_string(ac.coarsening_algorithm_id) <<
                        " with id " << ac.coarsening_algorithm_id << " not known!" << std::endl;
                exit(EXIT_FAILURE);
            }


            auto ep = get_time_point();
            coarsening_ms += get_milli_seconds(sp, ep);

#if ENABLE_PROFILER
            level_infos[v_cycle][level].t_coarsening += get_milli_seconds(sp, ep);
#endif
        }

        void contraction(const u64 v_cycle, [[maybe_unused]] const u64 level) {
            auto sp = get_time_point();

            graphs.emplace_back(); // coarse the graph
            graphs.back().initialize(graphs[graphs.size() - 2], mappings.back());
            p_managers[0].contract(mappings.back());

            auto ep = get_time_point();
            contraction_ms += get_milli_seconds(sp, ep);

#if ENABLE_PROFILER
            level_infos[v_cycle][level].t_contraction += get_milli_seconds(sp, ep);
#endif

            HEAVYASSERT(assert_state_pre_partitioning(graphs.back(), p_managers[0], ac.k));
        }

        void uncoarsening(const u64 v_cycle, [[maybe_unused]] const u64 level) {
            auto sp = get_time_point();

            p_managers[0].uncontract(mappings.back());
            bv_managers[0].compute_from_scratch(graphs[graphs.size() - 2], p_managers[0]);
            block_conns[0].compute_from_scratch(graphs[graphs.size() - 2], p_managers[0]);

            ScopedTimer _t("uncontraction", "misc", "free_graph");
            graphs.pop_back(); // this is doing uncontraction
            mappings.pop_back();
            _t.stop();

            auto ep = get_time_point();

            uncontraction_ms += get_milli_seconds(sp, ep);

#if ENABLE_PROFILER
            level_infos[v_cycle][level].t_uncontraction += get_milli_seconds(sp, ep);
#endif

            HEAVYASSERT(
                assert_state_after_partitioning(graphs.back(),p_managers[0], bv_managers[0], q_graphs[0], block_conns[0]
                    , ac.k));
        }

        void rebalancing(const u64 v_cycle, const u64 level, f64 level_imbalance, const u64 n_partitions) {
            auto sp = get_time_point();

            for (u64 i = 0; i < n_partitions; ++i) {
                if (level == 0) {
                    rebalancer.rebalance_last_layer(graphs.back(), p_managers[i], bv_managers[i], q_graphs[i], d_oracle, block_conns[i], level_imbalance);
                } else {
                    rebalancer.rebalance(graphs.back(), p_managers[i], bv_managers[i], q_graphs[i], d_oracle, block_conns[i], level_imbalance);
                }
            }

            auto ep = get_time_point();

            rebalance_ms += get_milli_seconds(sp, ep);

#if ENABLE_PROFILER
            level_infos[v_cycle][level].t_rebalance += get_milli_seconds(sp, ep);
#endif

            HEAVYASSERT(                assert_state_after_partitioning(graphs.back(), p_managers[0], bv_managers[0], q_graphs[0], block_conns[0                ], ac.k));
        }

        void refinement(const u64 v_cycle, const u64 level, const f64 level_imbalance, const u64 n_partitions) {
            auto sp = get_time_point();

            for (u64 i = 0; i < n_partitions; ++i) {
                u64 refinement_max_iterations = ac.n_refinement_iterations;
                for (u64 refinement_i = 0; refinement_i < refinement_max_iterations; ++refinement_i) {
                    for (auto [refiner, config]: refinements) {
                        if (config->enabled) {
                            refiner->refine(graphs.back(), d_oracle, bv_managers[i], p_managers[i], q_graphs[i], block_conns[i], level_imbalance);
                            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_managers[0], bv_managers[0], q_graphs[0], block_conns[0], ac.k));
                        }
                    }
                }
            }

            auto ep = get_time_point();

            refinement_ms += get_milli_seconds(sp, ep);

#if ENABLE_PROFILER
            level_infos[v_cycle][level].t_refinement += get_milli_seconds(sp, ep);
#endif

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_managers[0], bv_managers[0], q_graphs[0], block_conns[0 ], ac.k));
        }
    };
}

#endif //HEIPROMAP_SOLVER_H
