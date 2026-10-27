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

#ifndef HEIPROMAP_WAVERHEIPROMAP_SOLVER_H
#define HEIPROMAP_WAVERHEIPROMAP_SOLVER_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <thread>
#include <iostream>
#include <numeric>
#include <vector>

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
#include "../utility/engine_load.h"
#include "../utility/spectral_tracker.h"
#include "../utility/kaminpar_partitioner.h"
#include "../utility/spectral_prepartitioner.h"
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
#include "../refinement/waver_engine_refinement.h"
#include "../refinement/waver_sa_refinement.h"
#include "../refinement/waver_swap_refinement.h"
#include "../distance_oracles/grid_distance_oracle.h"

namespace HeiProMap {
    // Maps a graph onto a PE grid, minimizing the message cost.
    //
    // Pipeline:
    //   1. prepartition   n > k: partition the input graph into k blocks and contract them into the quotient graph with
    //                     one vertex per PE. Everything up to 6. works on a graph with n == k.
    //   2. coarsen        the topology and the graph
    //   3. initial        coarsest level: greedy partition (multilevel, both) or spectral placement (spectral), refined
    //   4. uncoarsen      level by level, legalize and refine on every level. The mapping is a bijection (one vertex
    //                     per PE) on the finest level.
    //   5. compare        spectral and both: the spectral placement of the finest level competes with the result;
    //                     the cheaper one is kept
    //   6. V-cycle        optional, on the finished bijection
    //   7. undo prepartition  projects the mapping onto the input vertices, which may differ in weight by --imbalance,
    //                     then legalizes and refines that level like any other
    class WaverHeiProMapSolver {
        // The pipeline always maps a bijection (one quotient vertex per PE): no imbalance, and every PE gets
        // exactly one vertex. ac.imbalance is the imbalance of the input graph's partition into PEs (KaMinPar).
        static constexpr f64 kImbalance = 0.0;
        static constexpr f64 kMinImbalance = 1.0;

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

        // The input graph while its quotient is solved, and the contraction between them (see prepartition)
        graph_t input_graph;
        Mapping input_mapping;
        vertex_t max_n = 0; // size of the input graph, the largest graph: the state is allocated for it
        vertex_t max_m = 0;
        f64 imbalance = kImbalance; // PE weight limit of the current graph: only the input graph is not a bijection

        // ---- coarsening ----
        GlobalPathAlgorithmMatcher gpa_matcher;
        SizeConstrainedLP size_constrained_lp;
        HeavyEdgeMatching heavy_edge_matching;

        // ---- rebalancing and refinement ----
        Rebalancer rebalancer;
        WaverLabelPropagationRefinement lp_refine;
        bool lp_active = false; // lp_refine is only set up (and run) in the multilevel and both modes

    public:
        explicit WaverHeiProMapSolver(graph_t &&g, graph_t &&t_topology, const WaverHeiProMapConfiguration &t_ac) {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "constructor");
            sp = std::chrono::high_resolution_clock::now();

            ac = t_ac;
            graphs.reserve(100);
            max_n = g.n;
            max_m = g.m;
            input_graph = std::move(g);
            random_engine = RandomEngine(ac.seed);

            spectral.set_smoothing(ac.spectral_smoothing);
            d_oracle.initialize(std::move(t_topology), ac.threads);
            ac.k = d_oracle.current_k();

            // The helpers are sized for the input graph, the largest graph, and the finest topology; they are used on
            // the smaller graphs and coarser topologies as well.
            gpa_matcher.initialize(max_n, max_m, ac.k, ac.threads, random_engine, ac.global_path_algorithm_config);
            size_constrained_lp.initialize(max_n, max_m, ac.k, random_engine.get_u64(), ac.size_constrained_lp_config);

            // Note: the non-LARGE_K label propagation allocates a k x k table (10^10 entries at k = 99856),
            // so it is not set up when it is not used.
            lp_active = ac.waver_label_propagation_config.enabled && ac.do_multilevel;
            if (lp_active) {
                lp_refine.initialize(random_engine.get_u64(), ac.waver_label_propagation_config);
            }
        }

        std::vector<vertex_t> solve() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "solve");
            std::vector<partition_t> result = internal_solve();

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
            const EndpointLoad endpoint = get_max_endpoint_load(q_graph, ac.k, ac.alpha_end(), ac.beta_end());
            const weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) graphs[0].g_weight / (f64) ac.k));

            const EngineLoadStats engine_model = get_engine_load(graphs[0], p_manager, ac.k);

            print_waver_summary(std::cout, ac, graphs[0], duration, lmax, cost, max_link_load, endpoint, engine_model, p_manager, timings);

            return result;
        }

    private:
        // Runs the whole pipeline and returns the PE of every vertex of the input graph.
        std::vector<partition_t> internal_solve() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "internal_solve");

            prepartition();

            // From here on n == k: coarsen the topology as far as it goes, then the graph until it has about
            // initial_c vertices per (coarse) PE.
            while (d_oracle.can_coarsen()) {
                d_oracle.coarsen();
                ac.k = d_oracle.current_k();
            }
            enter_topology_level(); // the coarsest one

            for (u64 level = 0; graphs.back().n > d_oracle.current_k() * ac.initial_c; ++level) {
                coarsening(level);
                contraction();

                if (graphs.back().n == graphs[graphs.size() - 2].n) {
                    graphs.pop_back();
                    mappings.pop_back();
                    break;
                }
            }

            initial_partitioning();

            while (uncoarsen_one_level()) { settle_level(); }

            compare_with_spectral_placement();

            if (ac.vcycle) { vcycle(); }

            undo_prepartition();
            return snapshot();
        }

        // ------------------------------------------------------------------------------------------------
        // setup
        // ------------------------------------------------------------------------------------------------

        // Resets the partition state for the current topology level: its number of cells and its weight limits. The
        // arrays are sized for the input graph, the largest graph, so they fit every level. The weights and limits are
        // those of the graph with one vertex per PE (graphs[0]).
        void enter_topology_level() {
            const graph_t &g = graphs[0];
            p_manager.initialize(max_n, ac.k, g.g_weight);
            p_manager.activate_all();
            p_manager.update_limits(d_oracle.current_topology(), g.g_weight, g.n, kImbalance, kMinImbalance);
            bv_manager.initialize(max_n, ac.k);
            rebalancer.initialize(max_n, max_m, ac.k, random_engine.get_u64());
        }

        // If the input graph has more vertices than there are PEs, partition it into k blocks (KaMinPar or spectral)
        // and contract them into the quotient graph (one unit-weight vertex per block, summed edge weights), so that
        // the graph has one vertex per PE (n == k) from here on. The input graph and the contraction are kept for
        // undo_prepartition.
        void prepartition() {
            const partition_t k = ac.k; // the topology is not coarsened yet
            if (input_graph.n <= k) {
                graphs.emplace_back(std::move(input_graph));
                return;
            }
            HEIPROMAP_PROFILE_SCOPE("partition", "WaverSolver", "prepartition");

            QuotientResult r;
            if (ac.prepartitioner == "spectral3d") {
                const size_t W = d_oracle.current_W();
                if (W == 0 || k % W != 0) {
                    std::cout << "Error: --prepartitioner spectral3d needs a full grid topology (k = " << k
                            << ", inferred width " << W << ")" << std::endl;
                    exit(EXIT_FAILURE);
                }
                r.part = SpectralPrepartitioner::partition(input_graph, W, k / W, ac.seed, ac.threads,
                                                          ac.spectral_smoothing);
                r.quotient = contract_into_quotient(input_graph, r.part, k, ac.threads);
            } else {
                r = partition_into_quotient(input_graph, k, ac.imbalance, ac.kaminpar_threads,
                                            ac.threads, ac.seed, ac.kaminpar_preset);
            }

            input_mapping.initialize(input_graph.n);
            input_mapping.set_coarse_n(k);
            for (vertex_t u = 0; u < input_graph.n; ++u) { input_mapping.set(u, r.part[u]); }
            graphs.emplace_back(std::move(r.quotient));
        }

        // ------------------------------------------------------------------------------------------------
        // coarsening
        // ------------------------------------------------------------------------------------------------

        void coarsening(const u64 level) {
            auto p = get_time_point();
            mappings.emplace_back();
            mappings.back().initialize(graphs.back().n);

            if (ac.coarsening_algorithm_id == COARSENING_ALG_GLOBAL_PATHS) {
                gpa_matcher.match(level, graphs.back(), p_manager, mappings.back(), kImbalance);
            } else if (ac.coarsening_algorithm_id == COARSENING_ALG_SIZE_CONSTRAINED_LP) {
                size_constrained_lp.cluster(level, graphs.back(), p_manager, mappings.back(), kImbalance, ac.threads);
            } else if (ac.coarsening_algorithm_id == COARSENING_ALG_HEAVY_EDGE) {
                heavy_edge_matching.match(graphs.back(), p_manager, mappings.back(), kImbalance, random_engine.get_u64(), ac.heavy_edge_matching_config);
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

        // ------------------------------------------------------------------------------------------------
        // initial partition and uncoarsening
        // ------------------------------------------------------------------------------------------------

        std::vector<partition_t> spectral_placement() {
            std::vector<partition_t> part;
            spectral.place(graphs.back(), d_oracle.current_topology(), d_oracle.current_W(), part);
            return part;
        }

        // Partition of the coarsest level: greedy in the multilevel and both mode, the spectral placement in the
        // spectral mode. Then the same steps as on every other level.
        void initial_partitioning() {
            if (ac.do_spectral) { spectral.init(graphs.back(), ac.seed); }

            if (ac.do_multilevel) {
                greedy_partitioning();
            } else {
                load_partition(spectral_placement());
            }
            recompute_datastructures();
            rebalance(mappings.size());
            refinement();
        }

        // Uncoarsens one thing: a topology level if the coarse graph allows it, otherwise a graph level. Never
        // both at once. Returns false once nothing is left to uncoarsen down to the quotient graph.
        bool uncoarsen_one_level() {
            if (can_uncoarsen_topology()) {
                uncoarsen_topology();
                return true;
            }
            if (mappings.empty()) { return false; }

            if (ac.do_spectral) { spectral.prolong(mappings.back(), graphs[graphs.size() - 2]); }
            uncoarsening();
            return true;
        }

        // Restores a valid, refined partition after a level change.
        void settle_level() {
            legalize(mappings.size());
            refinement();
        }

        // Spectral and both mode, once on the finest level (the graph with one vertex per PE): the spectral placement
        // of that level is an independent candidate that ignores the partition found while uncoarsening. If it is
        // cheaper than that partition, it is adopted (rebalanced and refined); if the refined placement is still not
        // cheaper, the partition from uncoarsening is kept. The spectral tracker itself is never overwritten.
        void compare_with_spectral_placement() {
            const graph_t &g = graphs.back();
            if (!ac.do_spectral || spectral.size() != g.n) { return; }
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "compare_with_spectral_placement");

            const std::vector<partition_t> spectral_part = spectral_placement();
            const f64 multilevel_cost = placement_cost(p_manager.partition);
            if (placement_cost(spectral_part) >= multilevel_cost) { return; }

            const std::vector<partition_t> multilevel_part = snapshot();
            load_partition(spectral_part);
            recompute_datastructures();
            legalize(mappings.size());
            refinement();
            if (current_message_cost() >= multilevel_cost) {
                load_partition(multilevel_part);
                recompute_datastructures();
            }
        }

        // Greedy partition of the coarsest level: the best of initial_kappa attempts (the first one in natural vertex
        // order, the others shuffled).
        void greedy_partitioning() {
            auto p = get_time_point();
            HEIPROMAP_PROFILE_SCOPE("partition", "partition", "partition");

            const graph_t &g = graphs.back();

            // p_manager.n: the arrays keep the size of the finest graph. The working copy carries the PE limits the
            // greedy partitioner respects.
            PartitionManager cur_pm;
            cur_pm.initialize(p_manager.n, ac.k, g.g_weight);
            cur_pm.copy_from(p_manager);
            PartitionManager best_pm;
            best_pm.initialize(p_manager.n, ac.k, g.g_weight);
            f64 best_cost = std::numeric_limits<f64>::max();
            AlignedArray<vertex_t> vertices;
            vertices.initialize(g.n);

            for (u64 trial = 0; trial < ac.initial_kappa; ++trial) {
                const u64 trial_seed = trial == 0 ? 0 : 1 + ac.seed + trial;
                waver_greedy_partition(g, d_oracle, kImbalance, trial_seed, cur_pm, vertices, ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
                const f64 cur_cost = get_waver_cost(g, cur_pm, d_oracle, ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);

                if (trial == 0 || cur_cost < best_cost) {
                    best_cost = cur_cost;
                    best_pm.copy_from(cur_pm);
                }
            }

            p_manager.copy_from(best_pm);
            p_manager.recalculate_weights(g);

            timings.initial_partitioning_ms += get_milli_seconds(p, get_time_point());
        }

        // True if the topology can be uncoarsened one level now: there is a level left, every coarse cell holds at
        // least as many vertices (counted, not weighed) as it has child cells, and (while the graph is coarse) the
        // graph is not too small for the finer topology.
        bool can_uncoarsen_topology() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "can_uncoarsen_topology");
            if (d_oracle.depth() == 0) return false;

            const partition_t coarse_k = d_oracle.current_k();
            const partition_t fine_k_val = d_oracle.fine_k();
            const Mapping &topo_map = d_oracle.top_mapping();

            std::vector<size_t> children_count(coarse_k, 0);
            for (partition_t f = 0; f < fine_k_val; ++f) {
                children_count[topo_map.get(f)]++;
            }

            // While the graph is coarse, keep at least topo_split_factor vertices per finer cell on average,
            // otherwise the finer topology is mostly empty/overloaded cells.
            if (!mappings.empty() && graphs.back().n < ac.topo_split_factor * fine_k_val) { return false; }
            // Counted from the partition: p_manager.size() is not kept up to date by the swap refinements.
            std::vector<size_t> vertex_count(coarse_k, 0);
            for (vertex_t u = 0; u < graphs.back().n; ++u) { ++vertex_count[p_manager[u]]; }
            for (partition_t c = 0; c < coarse_k; ++c) {
                if (vertex_count[c] < children_count[c]) { return false; }
            }

            return true;
        }

        void uncoarsen_topology() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "uncoarsen_topology");
            auto p = get_time_point();

            // The child cells must be read before the oracle level is popped.
            const auto coarse_to_fine = TopologyUncoarsener::children(d_oracle);

            d_oracle.uncoarsen();
            ac.k = d_oracle.current_k();

            const std::vector<partition_t> old_partition = snapshot();

            enter_topology_level();
            p_manager.reset_weights();

            const TopologyUncoarsener::CostParams costs{ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link};
            TopologyUncoarsener(graphs.back(), p_manager, d_oracle, costs, coarse_to_fine).run(old_partition);

            timings.uncontraction_ms += get_milli_seconds(p, get_time_point());
        }

        void uncoarsening() {
            auto p = get_time_point();

            const graph_t &finer = graphs[graphs.size() - 2];
            p_manager.uncontract(mappings.back());

            HEIPROMAP_PROFILE_SCOPE("uncontraction", "misc", "compute_from_scratch");
            bv_manager.reset(finer.n);
            for (vertex_t u = 0; u < finer.n; ++u) {
                const partition_t u_id = p_manager[u];
                size_t n_different = 0;

                for (size_t i = finer.neighborhoods[u]; i < finer.neighborhoods[u + 1]; ++i) {
                    n_different += (u_id != p_manager[finer.edges_v[i]]);
                }
                bv_manager.add_boundary_vertex_from_count(u, u_id, n_different);
            }

            graphs.pop_back();
            mappings.pop_back();

            q_graph.compute_from_scratch(graphs.back(), p_manager);
            block_conn.compute_from_scratch(graphs.back(), p_manager);
            timings.uncontraction_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
        }

        // ------------------------------------------------------------------------------------------------
        // rebalancing and refinement
        // ------------------------------------------------------------------------------------------------

        // Rebuilds boundary vertices, quotient graph and block connectivity from the current partition.
        void recompute_datastructures() {
            const auto &g_cur = graphs.back();
            bv_manager.reset(g_cur.n);
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

        // Makes the partition valid for the grid: rebalance, fill empty PEs, rebuild the derived state.
        // level is the number of graph levels still to uncoarsen.
        void legalize(const u64 level) {
            rebalance(level);
            EmptyBlockDiffusion::diffuse(graphs.back(), p_manager, static_cast<size_t>(d_oracle.current_W()));
            recompute_datastructures();
        }

        // All refinements, in order. Each is switched on by its configuration. Label propagation refines the
        // partition of a coarse level and only exists in the multilevel and both modes; simulated annealing swaps
        // whole PEs and therefore only runs on the finest level.
        void refinement() {
            auto p = get_time_point();

            if (lp_active) {
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

            // On the individual input vertices of the final topology (more vertices than PEs), where a PE can lose a
            // row. Before the refinements that swap whole PEs: they relabel p_manager only, so the quotient graph,
            // block connectivity and boundary are stale afterwards.
            if (mappings.empty() && d_oracle.depth() == 0 && graphs.back().n > ac.k) {
                auto sp_local = get_time_point();
                if (ac.waver_engine_refinement_config.pair_penalty > 0.0) {
                    WaverEngineRefinement::refine_pairs(graphs.back(), p_manager, bv_manager, q_graph, block_conn,
                                                        ac.waver_engine_refinement_config);
                }
                if (ac.waver_engine_refinement_config.enabled) {
                    WaverEngineRefinement::refine(graphs.back(), p_manager, bv_manager, q_graph, block_conn,
                                                  ac.waver_engine_refinement_config);
                }
                timings.engine_refine_ms += get_milli_seconds(sp_local, get_time_point());
            }

            if (ac.waver_swap_refinement_config.enabled) {
                auto sp_local = get_time_point();
                WaverSwapRefinement::swap_refine(graphs.back(), p_manager, d_oracle,
                                                 ac.waver_swap_refinement_config,
                                                 ac.alpha_end(), ac.beta_end(),
                                                 ac.alpha_link, ac.beta_link);
                timings.swap_refine_ms += get_milli_seconds(sp_local, get_time_point());
            }

            if (ac.waver_long_message_refinement_config.enabled) {
                auto sp_local = get_time_point();
                WaverLongMessageRefinement::refine(graphs.back(), p_manager, d_oracle,
                                                   ac.waver_long_message_refinement_config,
                                                   ac.alpha_end(), ac.beta_end(),
                                                   ac.alpha_link, ac.beta_link);
                timings.long_msg_refine_ms += get_milli_seconds(sp_local, get_time_point());
            }

            if (ac.waver_sa_refinement_config.enabled && d_oracle.depth() == 0 && mappings.empty()) {
                auto sp_local = get_time_point();
                WaverSARefinement::refine(graphs.back(), p_manager, d_oracle,
                                          ac.waver_sa_refinement_config, ac.seed,
                                          ac.alpha_end(), ac.beta_end(),
                                          ac.alpha_link, ac.beta_link);
                timings.sa_refine_ms += get_milli_seconds(sp_local, get_time_point());
            }

            timings.refinement_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
        }

        // ------------------------------------------------------------------------------------------------
        // finished bijection
        // ------------------------------------------------------------------------------------------------

        // V-cycle on a finished bijection. Topology and graph are coarsened the same way: the vertices that sit in the
        // cells merged into a coarse cell are contracted into one coarse vertex, so every level has one vertex per
        // cell (n' == k'). Uncoarsening then restores graph and topology together, level by level, and legalizes and
        // refines every level like the main uncoarsening, starting at the coarsest. The previous mapping is restored
        // if the cycle does not lower the message cost.
        void vcycle() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "vcycle");
            const std::vector<partition_t> saved = snapshot();
            const f64 cost_before = current_message_cost();

            const u64 levels = coarsen_mapped_bijection();

            recompute_datastructures();
            settle_level();
            for (u64 l = 0; l < levels; ++l) {
                uncoarsening();
                uncoarsen_topology();
                recompute_datastructures();
                settle_level();
            }

            if (current_message_cost() >= cost_before) {
                load_partition(saved);
                recompute_datastructures();
            }
        }

        // Coarsens topology and graph together by the current placement: coarse vertex c is the union of the vertices
        // placed in the cells merged into coarse cell c, so the coarse graph has exactly one vertex per coarse cell.
        // Returns the number of levels.
        u64 coarsen_mapped_bijection() {
            u64 levels = 0;
            while (d_oracle.can_coarsen()) {
                d_oracle.coarsen();
                const Mapping &topo_map = d_oracle.top_mapping();
                const vertex_t fine_n = graphs.back().n;

                mappings.emplace_back();
                mappings.back().initialize(fine_n);
                mappings.back().set_coarse_n(d_oracle.current_k());
                for (vertex_t u = 0; u < fine_n; ++u) { mappings.back().set(u, topo_map.get(p_manager[u])); }

                graphs.emplace_back();
                graphs.back().contract(graphs[graphs.size() - 2], mappings.back(), ac.threads);

                ac.k = d_oracle.current_k();
                enter_topology_level();
                p_manager.reset_weights();
                for (vertex_t c = 0; c < graphs.back().n; ++c) { p_manager.set(c, graphs.back().v_weight(c), c); }
                ++levels;
            }
            return levels;
        }

        // The inverse of prepartition, one more uncoarsening level: the mapping of the quotient graph is projected
        // onto the individual input vertices, then legalized and refined like every other level. Here the PEs may
        // deviate from the average weight by --imbalance, and the weights are those of the input graph.
        void undo_prepartition() {
            if (input_mapping.get_old_n() == 0) { return; }
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "undo_prepartition");

            graphs.insert(graphs.begin(), std::move(input_graph));
            mappings.push_back(std::move(input_mapping));
            uncoarsening();

            imbalance = ac.imbalance;
            const graph_t &input = graphs.back();
            p_manager.update_limits(d_oracle.current_topology(), input.g_weight, input.n, imbalance, kMinImbalance);
            p_manager.recalculate_weights(input);
            settle_level();
        }

        // ------------------------------------------------------------------------------------------------
        // partition helpers
        // ------------------------------------------------------------------------------------------------

        std::vector<partition_t> snapshot() const {
            std::vector<partition_t> part(graphs.back().n);
            for (vertex_t u = 0; u < part.size(); ++u) { part[u] = p_manager[u]; }
            return part;
        }

        void load_partition(const std::vector<partition_t> &part) {
            for (vertex_t u = 0; u < graphs.back().n; ++u) p_manager.partition[u] = part[u];
            p_manager.recalculate_weights(graphs.back());
        }

        template<typename PartT>
        f64 placement_cost(const PartT &part) {
            return SpectralTracker::message_cost(graphs.back(), part, d_oracle,
                                                 ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
        }

        f64 current_message_cost() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "current_message_cost");
            recompute_datastructures();
            return get_message_cost(q_graph, ac.k, d_oracle, ac.alpha_end(), ac.beta_end(),
                                    ac.alpha_link, ac.beta_link).total;
        }
    };
}

#endif //HEIPROMAP_WAVERHEIPROMAP_SOLVER_H
