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

#ifndef HEIPROMAP_WAVER_MAP_SOLVER_H
#define HEIPROMAP_WAVER_MAP_SOLVER_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <thread>
#include <iostream>
#include <numeric>
#include <vector>
#include <omp.h>

#include "../datastructures/block_conn.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/large_quotient_graph.h"
#include "../datastructures/partition_manager.h"

#include "../definitions.h"
#include "../configuration/WaverHeiProMap_configuration.h"
#include "../utility/macros.h"
#include "../utility/random_engine.h"
#include "../utility/assert_state.h"
#include "../utility/qap.h"
#include "../utility/spectral_tracker.h"
#include "../utility/kaminpar_partitioner.h"
#include "../utility/waver_report.h"
#include "../utility/utils.h"

#include "../coarsening/global_path_algorithm.h"
#include "../coarsening/size_constrained_lp.h"
#include "../coarsening/heavy_edge_matching.h"
#include "../rebalance/rebalancer.h"
#include "../rebalance/empty_block_diffusion.h"
#include "../uncoarsening/topology_uncoarsener.h"
#include "../init_partitioning/greedy_partitioner.h"
#include "../init_partitioning/waver_greedy_partitioner.h"
#include "../refinement/simple_label_propagation_refinement.h"
#include "../refinement/waver_label_propagation_refinement.h"
#include "../refinement/waver_message_refinement.h"
#include "../refinement/waver_long_message_refinement.h"
#include "../refinement/waver_sa_refinement.h"
#include "../refinement/waver_swap_refinement.h"
#include "../distance_oracles/grid_distance_oracle.h"

namespace HeiProMap {
    class WaverHeiProMapSolver {
        // ---- configuration and run state ----
        WaverHeiProMapConfiguration ac;
        RandomEngine random_engine;
        std::chrono::high_resolution_clock::time_point sp;
        WaverTimings timings;

        // ---- graph hierarchy and partition state ----
        std::vector<graph_t> graphs;
        std::vector<Mapping> mappings; // mappings[i]: graphs[i] -> graphs[i + 1]
        p_manager_t p_manager;
        bv_manager_t bv_manager;
        LargeQuotientGraph q_graph;
        block_conn_t block_conn;
        GridDistanceOracle d_oracle;

        // ---- spectral solution (follows the graph hierarchy) ----
        SpectralTracker spectral;

        // ---- input reduction: set when the input graph had more vertices than PEs ----
        std::vector<partition_t> kaminpar_part; // KaMinPar block of every original vertex
        vertex_t orig_n = 0;                    // number of original vertices

        // ---- coarsening ----
        GlobalPathAlgorithmMatcher gpa_matcher;
        SizeConstrainedLP size_constrained_lp;
        HeavyEdgeMatching heavy_edge_matching;

        // ---- rebalancing and refinement ----
        Rebalancer rebalancer;
        WaverLabelPropagationRefinement<false> lp_refine;

        void reinit_structures_for_new_k() {
            p_manager.initialize(graphs[0].n, ac.k, graphs[0].g_weight);
            p_manager.activate_all();
            p_manager.update_limits(d_oracle.current_topology(), graphs[0].g_weight, graphs[0].n,
                                    ac.imbalance, ac.min_imbalance);
            bv_manager.initialize(graphs[0].n, ac.k);

            rebalancer.initialize(graphs[0].n, graphs[0].m, ac.k, random_engine.get_u64());
        }

        void initialize_common(CSRGraph topology) {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "initialize_common");
            random_engine = RandomEngine(ac.seed);

            d_oracle.initialize(std::move(topology), ac.threads);
            ac.k = d_oracle.current_k();

            p_manager.initialize(graphs[0].n, ac.k, graphs[0].g_weight);
            p_manager.activate_all();
            p_manager.update_limits(d_oracle.current_topology(), graphs[0].g_weight, graphs[0].n,
                                    ac.imbalance, ac.min_imbalance);
            bv_manager.initialize(graphs[0].n, ac.k);
            HEAVYASSERT(assert_state_pre_partitioning(graphs[0], p_manager, ac.k));

            gpa_matcher.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine, ac.global_path_algorithm_config);
            size_constrained_lp.initialize(graphs[0].n, graphs[0].m, ac.k, random_engine.get_u64(), ac.size_constrained_lp_config);

            rebalancer.initialize(graphs[0].n, graphs[0].m, ac.k, random_engine.get_u64());

            // Note: the non-LARGE_K label propagation allocates a k x k table (10^10 entries at k = 99856),
            // so it is not set up when refinement is disabled.
            if (ac.waver_label_propagation_config.enabled && !ac.spectral_only) {
                lp_refine.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine.get_u64(), ac.waver_label_propagation_config);
            }
        }

        // If g has more vertices than there are PEs, partition it into k blocks with KaMinPar and replace it by
        // its quotient graph (one unit-weight vertex per block, summed edge weights), so n == k afterwards.
        void reduce_to_quotient(const partition_t k) {
            CSRGraph &g = graphs[0];
            if (g.n <= k) { return; }
            HEIPROMAP_PROFILE_SCOPE("partition", "WaverSolver", "reduce_to_quotient");

            orig_n = g.n;
            QuotientResult r = partition_into_quotient(g, k, ac.kaminpar_imbalance, ac.kaminpar_threads,
                                                       ac.threads, ac.seed, ac.kaminpar_preset);
            kaminpar_part = std::move(r.part);
            graphs[0] = std::move(r.quotient);
        }

    public:
        explicit WaverHeiProMapSolver(graph_t &&g, graph_t &&topology, const WaverHeiProMapConfiguration &t_ac) {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "constructor");
            sp = std::chrono::high_resolution_clock::now();

            ac = t_ac;
            graphs.emplace_back(std::move(g));

            reduce_to_quotient(topology.n);
            initialize_common(std::move(topology));
        }

        std::vector<vertex_t> solve() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "solve");
            internal_solve();

            std::vector<partition_t> p(graphs[0].n);
            for (vertex_t u = 0; u < graphs[0].n; ++u) { p[u] = p_manager[u]; }

            const auto ep = std::chrono::high_resolution_clock::now();
            const f64 duration = get_seconds(sp, ep);

            q_graph.compute_from_scratch(graphs[0], p_manager);
            block_conn.compute_from_scratch(graphs[0], p_manager);
            HEAVYASSERT(assert_state_after_partitioning(graphs[0], p_manager, bv_manager, q_graph, block_conn, ac.k));
            const MessageCostBreakdown cost = get_message_cost(q_graph, ac.k, d_oracle,
                                                               ac.alpha_end(), ac.beta_end(),
                                                               ac.alpha_link, ac.beta_link);
            const f64 max_link_load = get_max_link_load(q_graph, ac.k, d_oracle.current_W(),
                                                        ac.alpha_link, ac.beta_link);
            const weight_t lmax = std::ceil((1.0 + ac.imbalance) * ((f64) graphs[0].g_weight / (f64) ac.k));

            print_waver_summary(std::cout, ac, graphs[0], duration, lmax, cost, max_link_load, p_manager, timings);

            if (!kaminpar_part.empty()) {
                // Propagate the block -> PE mapping back to the original vertices.
                std::vector<partition_t> orig_p(orig_n);
                for (vertex_t u = 0; u < orig_n; ++u) { orig_p[u] = p[kaminpar_part[u]]; }
                return orig_p;
            }
            return p;
        }


        void recompute_datastructures() {
            const auto &g_cur = graphs.back();
            bv_manager.reset();
            for (vertex_t u = 0; u < g_cur.n; ++u) {
                const partition_t u_id = p_manager[u];
                size_t n_different = 0;
                for (size_t i = g_cur.neighborhoods[u]; i < g_cur.neighborhoods[u + 1]; ++i) {
                    n_different += (u_id != p_manager[g_cur.edges_v[i]]);
                }
                bv_manager.add_boundary_vertex_from_count(u, u_id, n_different);
            }
            q_graph.compute_from_scratch(g_cur, p_manager);
            block_conn.compute_from_scratch(g_cur, p_manager);
        }

        // Multilevel mode: if the tracked spectral placement is cheaper than the current multilevel partition,
        // adopt it (rebalanced and refined). The spectral tracker itself is never overwritten.
        void sync_with_spectral_solution(const u64 level) {
            const graph_t &g = graphs.back();
            if (!ac.spectral_solution || ac.spectral_only || spectral.size() != g.n) { return; }
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "sync_with_spectral_solution");

            std::vector<partition_t> sp_part;
            spectral.place(g, d_oracle.current_topology(), d_oracle.current_W(), sp_part);
            const f64 ml = SpectralTracker::message_cost(g, p_manager.partition, d_oracle,
                                                         ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
            const f64 sp = SpectralTracker::message_cost(g, sp_part, d_oracle,
                                                         ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
            if (sp >= ml) { return; }

            std::vector<partition_t> ml_part(g.n);
            for (vertex_t u = 0; u < g.n; ++u) ml_part[u] = p_manager[u];
            load_partition(sp_part);
            recompute_datastructures();
            rebalance(level);
            EmptyBlockDiffusion::diffuse(graphs.back(), p_manager, static_cast<size_t>(d_oracle.current_W()));
            recompute_datastructures();
            refinement();
            if (current_message_cost() >= ml) {
                load_partition(ml_part);
                recompute_datastructures();
            }
        }

        void load_partition(const std::vector<partition_t> &part) {
            for (vertex_t u = 0; u < graphs.back().n; ++u) p_manager.partition[u] = part[u];
            p_manager.recalculate_weights(graphs.back());
        }

        f64 current_message_cost() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "current_message_cost");
            recompute_datastructures();
            return get_message_cost(q_graph, ac.k, d_oracle, ac.alpha_end(), ac.beta_end(),
                                    ac.alpha_link, ac.beta_link).total;
        }

        bool try_uncoarsen_topology() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "try_uncoarsen_topology");
            if (d_oracle.depth() == 0) return false;

            const partition_t coarse_k = d_oracle.current_k();
            const partition_t fine_k_val = d_oracle.fine_k();
            const Mapping &topo_map = d_oracle.top_mapping();

            std::vector<partition_t> children_count(coarse_k, 0);
            for (partition_t f = 0; f < fine_k_val; ++f) {
                children_count[topo_map.get(f)]++;
            }

            // While the graph is coarse, keep at least topo_split_factor vertices per finer cell on average,
            // otherwise the finer topology is mostly empty/overloaded cells.
            if (!mappings.empty() && graphs.back().n < static_cast<u64>(ac.topo_split_factor) * fine_k_val) { return false; }
            for (partition_t c = 0; c < coarse_k; ++c) {
                if (p_manager.get_bweight(c) < static_cast<weight_t>(children_count[c])) { return false; }
            }

            uncoarsen_topology();
            p_manager.update_limits(d_oracle.current_topology(), graphs[0].g_weight, graphs[0].n,
                                    ac.imbalance, ac.min_imbalance);
            return true;
        }

    private:
        void internal_solve() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "internal_solve");
            u64 level = 0;

            // Pre-coarsen the topology oracle as far as possible.
            while (d_oracle.can_coarsen()) {
                d_oracle.coarsen();
                ac.k = d_oracle.current_k();
                reinit_structures_for_new_k();
            }

            p_manager.update_limits(d_oracle.current_topology(), graphs[0].g_weight, graphs[0].n,
                                    ac.imbalance, ac.min_imbalance);

            // Coarsen the graph.
            while (graphs.back().n > d_oracle.current_k() * ac.initial_c) {
                coarsening(level);
                contraction();

                if (graphs.back().n == graphs[graphs.size() - 2].n) {
                    graphs.pop_back();
                    mappings.pop_back();
                    break;
                }

                level += 1;
            }

            if (ac.spectral_solution) spectral.init(graphs.back(), ac.seed);

            // Initial partitioning at coarsest level.
            initial_partitioning();
            recompute_datastructures();
            rebalance(level);
            refinement();
            sync_with_spectral_solution(level);

            // Graph uncoarsening loop.
            while (!mappings.empty()) {
                level -= 1;

                // Try to uncoarsen topology levels before graph uncoarsening.
                while (d_oracle.depth() > 0) {
                    if (!try_uncoarsen_topology()) break;

                    rebalance(level);
                    EmptyBlockDiffusion::diffuse(graphs.back(), p_manager, static_cast<size_t>(d_oracle.current_W()));
                    recompute_datastructures();
                    refinement();
                }

                if (ac.spectral_solution) spectral.prolong(mappings.back(), graphs[graphs.size() - 2]);
                uncoarsening();

                rebalance(level);
                EmptyBlockDiffusion::diffuse(graphs.back(), p_manager, static_cast<size_t>(d_oracle.current_W()));
                recompute_datastructures();
                refinement();
                sync_with_spectral_solution(level);
            }

            // Uncoarsen remaining topology levels after graph is fully uncoarsened.
            while (d_oracle.depth() > 0) {
                if (!try_uncoarsen_topology()) {
                    rebalance(0);
                    EmptyBlockDiffusion::diffuse(graphs.back(), p_manager, static_cast<size_t>(d_oracle.current_W()));
                    recompute_datastructures();
                    refinement();
                    sync_with_spectral_solution(0);
                    break;
                }

                rebalance(0);
                EmptyBlockDiffusion::diffuse(graphs.back(), p_manager, static_cast<size_t>(d_oracle.current_W()));
                recompute_datastructures();
                refinement();
                sync_with_spectral_solution(0);
            }

            if (ac.spectral_only && ac.spectral_solution && spectral.size() == graphs.back().n) {
                adopt_spectral_only();
            }
        }

        // Spectral-only mode: replace the multilevel partition by the spectral placement (rebalanced, not refined).
        void adopt_spectral_only() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "adopt_spectral_only");
            const graph_t &g = graphs.back();
            std::vector<partition_t> sp_part;
            spectral.place(g, d_oracle.current_topology(), d_oracle.current_W(), sp_part);
            load_partition(sp_part);
            recompute_datastructures();
            rebalance(0);
            EmptyBlockDiffusion::diffuse(graphs.back(), p_manager, static_cast<size_t>(d_oracle.current_W()));
            recompute_datastructures();
            if (ac.waver_message_refinement_config.enabled) {
                WaverMessageRefinement::refine(graphs.back(), p_manager, bv_manager, q_graph, d_oracle, block_conn,
                                               ac.waver_message_refinement_config,
                                               ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
            }
            if (ac.swap_radius > 0) {
                WaverSwapRefinement::swap_refine(graphs.back(), p_manager, d_oracle,
                                                 ac.swap_radius, ac.swap_extended,
                                                 ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
            }
            if (ac.long_msg_padding > 0) {
                WaverLongMessageRefinement::refine(graphs.back(), p_manager, d_oracle,
                                                   ac.long_msg_padding,
                                                   ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
            }
        }

        void initial_partitioning() {
            auto p = get_time_point();
            HEIPROMAP_PROFILE_SCOPE("partition", "partition", "partition");

            const u64 kappa = ac.initial_kappa;

            if (ac.threads > 1) {
                std::vector<PartitionManager> local_best_pm(ac.threads);
                std::vector<f64> local_best_cost(ac.threads, std::numeric_limits<f64>::max());

                for (u64 t = 0; t < ac.threads; ++t) {
                    local_best_pm[t].initialize(graphs.back().n, ac.k, graphs.back().g_weight);
                }

                #pragma omp parallel num_threads(ac.threads)
                {
                    u64 tid = omp_get_thread_num();
                    PartitionManager cur_pm;
                    cur_pm.initialize(graphs.back().n, ac.k, graphs.back().g_weight);
                    AlignedArray<vertex_t> vertices;
                    vertices.initialize(graphs.back().n);

                    for (u64 trial = 0; trial < kappa; ++trial) {
                        // tid=0, trial=0 → seed=0 (natural vertex order, no shuffle)
                        const u64 trial_seed = (tid == 0 && trial == 0) ? 0 : (1 + ac.seed + tid * kappa + trial);
                        waver_greedy_partition(graphs.back(), d_oracle, ac.imbalance, trial_seed, cur_pm, vertices, ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
                        f64 cur_cost = get_waver_cost(graphs.back(), cur_pm, d_oracle, ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);

                        if (trial == 0 || cur_cost < local_best_cost[tid]) {
                            local_best_cost[tid] = cur_cost;
                            local_best_pm[tid].copy_from(cur_pm);
                        }
                    }
                }

                size_t best_idx = 0;
                for (u64 t = 1; t < ac.threads; ++t) {
                    if (local_best_cost[t] < local_best_cost[best_idx]) {
                        best_idx = t;
                    }
                }

                p_manager.copy_from(local_best_pm[best_idx]);
            } else {
                PartitionManager best_pm;
                best_pm.initialize(p_manager.n, ac.k, graphs.back().g_weight);
                f64 best_cost = std::numeric_limits<f64>::max();
                AlignedArray<vertex_t> vertices;
                vertices.initialize(graphs.back().n);

                for (u64 trial = 0; trial < kappa; ++trial) {
                    // trial=0 → seed=0 (natural vertex order, no shuffle)
                    const u64 trial_seed = (trial == 0) ? 0 : (1 + ac.seed + trial);
                    waver_greedy_partition(graphs.back(), d_oracle, ac.imbalance, trial_seed, p_manager, vertices, ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
                    f64 cur_cost = get_waver_cost(graphs.back(), p_manager, d_oracle, ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);

                    if (trial == 0 || cur_cost < best_cost) {
                        best_cost = cur_cost;
                        best_pm.copy_from(p_manager);
                    }
                }

                p_manager.copy_from(best_pm);
            }

            HEIPROMAP_PROFILE_SCOPE("partition", "misc", "initialize_datastructures");
            p_manager.reset_weights();
            bv_manager.reset();

            for (vertex_t u = 0; u < graphs.back().n; ++u) {
                const partition_t u_id = p_manager[u];
                const weight_t u_w = graphs.back().v_weight(u);
                p_manager.set(u, u_w, u_id);

                for (size_t i = graphs.back().neighborhoods[u]; i < graphs.back().neighborhoods[u + 1]; ++i) {
                    const vertex_t v = graphs.back().edges_v[i];
                    const partition_t v_id = p_manager[v];

                    if (u_id != v_id) {
                        bv_manager.add(u, u_id);
                    }
                }
            }


            q_graph.compute_from_scratch(graphs.back(), p_manager);
            block_conn.compute_from_scratch(graphs.back(), p_manager);
            timings.initial_partitioning_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
        }

        void coarsening(const u64 level) {
            auto p = get_time_point();
            mappings.emplace_back();
            mappings.back().initialize(graphs.back().n);

            if (ac.coarsening_algorithm_id == COARSENING_ALG_GLOBAL_PATHS) {
                gpa_matcher.match(level, graphs.back(), p_manager, mappings.back(), ac.imbalance);
            } else if (ac.coarsening_algorithm_id == COARSENING_ALG_SIZE_CONSTRAINED_LP) {
                size_constrained_lp.cluster(level, graphs.back(), p_manager, mappings.back(), ac.imbalance, ac.threads);
            } else if (ac.coarsening_algorithm_id == COARSENING_ALG_HEAVY_EDGE) {
                heavy_edge_matching.match(graphs.back(), p_manager, mappings.back(), ac.imbalance, random_engine.get_u64(), ac.heavy_edge_matching_config);
            } else {
                std::cerr << "Coarsening algorithm " << coarsening_algorithm_to_string(ac.coarsening_algorithm_id) << " with id " << ac.coarsening_algorithm_id << " not known!" << std::endl;
                exit(EXIT_FAILURE);
            }

            timings.coarsening_ms += get_milli_seconds(p, get_time_point());
        }

        void contraction() {
            auto p = get_time_point();
            graphs.emplace_back();
            graphs.back().contract(graphs[graphs.size() - 2], mappings.back(), ac.threads);
            p_manager.contract(mappings.back());

            timings.contraction_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_pre_partitioning(graphs.back(), p_manager, ac.k));
        }

        void uncoarsen_topology() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "uncoarsen_topology");
            auto p = get_time_point();

            // The child cells must be read before the oracle level is popped.
            const auto coarse_to_fine = TopologyUncoarsener::children(d_oracle);

            d_oracle.uncoarsen();
            ac.k = d_oracle.current_k();

            std::vector<partition_t> old_partition(graphs.back().n);
            for (vertex_t u = 0; u < graphs.back().n; ++u) {
                old_partition[u] = p_manager[u];
            }

            reinit_structures_for_new_k();
            p_manager.reset_weights();

            const TopologyUncoarsener::CostParams costs{ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link};
            TopologyUncoarsener(graphs.back(), p_manager, d_oracle, costs, coarse_to_fine).run(old_partition);

            timings.uncontraction_ms += get_milli_seconds(p, get_time_point());
        }

        void uncoarsening() {
            auto p = get_time_point();

            p_manager.uncontract(mappings.back());

            HEIPROMAP_PROFILE_SCOPE("uncontraction", "misc", "compute_from_scratch");
            const graph_t &g_uncontracted = graphs[graphs.size() - 2];

            if (ac.threads > 1) {
                bv_manager.parallel_reset(ac.threads);

                std::vector<std::vector<std::vector<vertex_t> > > thread_boundaries(ac.threads, std::vector<std::vector<vertex_t> >(ac.k));

                #pragma omp parallel num_threads(ac.threads)
                {
                    u64 tid = omp_get_thread_num();
                    vertex_t chunk = (g_uncontracted.n + ac.threads - 1) / ac.threads;
                    vertex_t start_u = tid * chunk;
                    vertex_t end_u = std::min(g_uncontracted.n, start_u + chunk);

                    for (vertex_t u = start_u; u < end_u; ++u) {
                        const partition_t u_id = p_manager[u];
                        size_t n_different = 0;

                        for (size_t i = g_uncontracted.neighborhoods[u]; i < g_uncontracted.neighborhoods[u + 1]; ++i) {
                            const vertex_t v = g_uncontracted.edges_v[i];
                            const partition_t v_id = p_manager[v];
                            n_different += (u_id != v_id);
                        }

                        if (n_different > 0) {
                            bv_manager.set_boundary_edges_count(u, n_different);
                            thread_boundaries[tid][u_id].push_back(u);
                        }
                    }
                }

                bv_manager.parallel_import_boundary_vertices(thread_boundaries, ac.threads);
            } else {
                bv_manager.reset();

                for (vertex_t u = 0; u < g_uncontracted.n; ++u) {
                    const partition_t u_id = p_manager[u];
                    size_t n_different = 0;

                    for (size_t i = g_uncontracted.neighborhoods[u]; i < g_uncontracted.neighborhoods[u + 1]; ++i) {
                        const vertex_t v = g_uncontracted.edges_v[i];
                        const partition_t v_id = p_manager[v];
                        n_different += (u_id != v_id);
                    }
                    bv_manager.add_boundary_vertex_from_count(u, u_id, n_different);
                }
            }

            graphs.pop_back();
            mappings.pop_back();


            q_graph.compute_from_scratch(graphs.back(), p_manager);
            block_conn.compute_from_scratch(graphs.back(), p_manager);
            timings.uncontraction_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
        }

        void rebalance(const u64 level) {
            auto p = get_time_point();

            if (level == 0) {
                rebalancer.rebalance_last_layer<true, true>(graphs.back(), p_manager, bv_manager, q_graph, d_oracle, block_conn);
            } else {
                rebalancer.rebalance<true, true>(graphs.back(), p_manager, bv_manager, q_graph, d_oracle, block_conn);
            }

            timings.rebalance_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
        }

        void refinement() {
            if (ac.spectral_only) { return; }
            auto p = get_time_point();

            if (ac.waver_label_propagation_config.enabled) {
                auto sp_local = get_time_point();
                lp_refine.refine(graphs.back(), d_oracle, bv_manager, p_manager, q_graph, block_conn);
                timings.lp_refine_ms += get_milli_seconds(sp_local, get_time_point());
            }

            if (ac.waver_message_refinement_config.enabled) {
                auto sp_local = get_time_point();
                WaverMessageRefinement::refine(graphs.back(), p_manager, bv_manager, q_graph, d_oracle, block_conn,
                                               ac.waver_message_refinement_config,
                                               ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
                timings.msg_refine_ms += get_milli_seconds(sp_local, get_time_point());
            }

            if (ac.swap_radius > 0) {
                WaverSwapRefinement::swap_refine(graphs.back(), p_manager, d_oracle,
                                                 ac.swap_radius, ac.swap_extended,
                                                 ac.alpha_end(), ac.beta_end(),
                                                 ac.alpha_link, ac.beta_link);
            }

            if (ac.long_msg_padding > 0) {
                WaverLongMessageRefinement::refine(graphs.back(), p_manager, d_oracle,
                                                   ac.long_msg_padding,
                                                   ac.alpha_end(), ac.beta_end(),
                                                   ac.alpha_link, ac.beta_link);
            }

            if (ac.sa_iterations > 0 && d_oracle.depth() == 0) {
                WaverSARefinement::refine(graphs.back(), p_manager, d_oracle,
                                          ac.sa_iterations, ac.sa_radius,
                                          ac.seed + 99999,
                                          ac.alpha_end(), ac.beta_end(),
                                          ac.alpha_link, ac.beta_link);
            }

            // WaverSwapRefinement::vertex_move_refine(graphs.back(), p_manager, d_oracle,
            //                                         ac.vmove_radius, ac.vmove_extended,
            //                                         ac.alpha_end(), ac.beta_end(),
            //                                         ac.alpha_link, ac.beta_link);

            timings.refinement_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
        }
    };
}

#endif //HEIPROMAP_WAVER_MAP_SOLVER_H
