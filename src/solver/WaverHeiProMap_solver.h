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
#include "../refinement/waver_sa_refinement.h"
#include "../refinement/waver_swap_refinement.h"
#include "../distance_oracles/grid_distance_oracle.h"

namespace HeiProMap {
    // Maps a graph onto a PE grid, minimizing the message cost.
    //
    // Pipeline (the hierarchy is graphs[0] = input graph, graphs[1] = its quotient with one vertex per PE, then the
    // coarser graphs; mappings[i] contracts graphs[i] into graphs[i + 1]):
    //   1. coarsen     n > k: partition the input graph into k blocks and contract them (the first level);
    //                  then the topology and the quotient graph
    //   2. initial     coarsest level: greedy partition (multilevel, both) or spectral placement (spectral), refined
    //   3. uncoarsen   level by level, legalize and refine on every level; the spectral placement of a level replaces
    //                  the current partition whenever it is cheaper (spectral, both). Down to the quotient graph the
    //                  mapping is a bijection (one vertex per PE).
    //   4. V-cycle     optional, on the finished bijection of the quotient graph
    //   5. input level the last uncoarsening step projects onto the input vertices, which may differ in weight by
    //                  --imbalance
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
        graph_t topology; // PE topology; handed to the distance oracle by internal_solve
        std::vector<graph_t> graphs;
        std::vector<Mapping> mappings; // mappings[i]: graphs[i] -> graphs[i + 1]
        p_manager_t p_manager;
        bv_manager_t bv_manager;
        LargeQuotientGraph q_graph;
        block_conn_t block_conn;
        GridDistanceOracle d_oracle;

        // ---- spectral solution (follows the graph hierarchy) ----
        SpectralTracker spectral;

        // ---- what the last level (the individual input vertices) changed ----
        InputLevelStats input_level;
        bool reduced = false; // true if graphs[0] is the input graph and graphs[1] its quotient
        partition_t finest_k = 0; // number of PEs
        size_t finest_width = 0; // width of the PE grid

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
            graphs.emplace_back(std::move(g));
            topology = std::move(t_topology);
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
            const weight_t lmax = std::ceil((1.0 + level_imbalance()) * ((f64) graphs[0].g_weight / (f64) ac.k));

            print_waver_summary(std::cout, ac, graphs[0], duration, lmax, cost, max_link_load, p_manager, timings, input_level);

            return result;
        }

    private:
        // Runs the whole pipeline and returns the PE of every vertex of the input graph.
        std::vector<partition_t> internal_solve() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "internal_solve");

            initialize_common(std::move(topology));

            coarsen_topology();
            coarsen_graph();

            initial_partitioning();
            while (uncoarsen_one_level()) { settle_level(); }

            if (ac.vcycle && is_bijection()) { vcycle(); }
            if (reduced) { refine_input_level(); }
            return snapshot();
        }

        // ------------------------------------------------------------------------------------------------
        // setup
        // ------------------------------------------------------------------------------------------------

        // The graph with one vertex per PE: the largest graph the state is sized for until the input level.
        const graph_t &solver_graph() const { return graphs[reduced ? 1 : 0]; }

        // Number of graph levels that are uncoarsened by the main loop (the first level, the reduction of the
        // input graph to its quotient, is the last one and is handled separately).
        size_t base_levels() const { return reduced ? 1 : 0; }
        bool graph_levels_left() const { return mappings.size() > base_levels(); }

        // The input level may deviate from the average PE weight by the user's imbalance, all others are bijections.
        f64 level_imbalance() const { return mappings.empty() && reduced ? ac.imbalance : kImbalance; }

        // (Re)initializes the partition and boundary state for the current topology level.
        void init_state_for_k() {
            const graph_t &g = solver_graph();
            p_manager.initialize(g.n, ac.k, g.g_weight);
            p_manager.activate_all();
            p_manager.update_limits(d_oracle.current_topology(), g.g_weight, g.n, kImbalance, kMinImbalance);
            bv_manager.initialize(g.n, ac.k);
        }

        void reinit_structures_for_new_k() {
            init_state_for_k();
            rebalancer.initialize(solver_graph().n, solver_graph().m, ac.k, random_engine.get_u64());
        }

        void initialize_common(CSRGraph topology) {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "initialize_common");
            random_engine = RandomEngine(ac.seed);

            d_oracle.initialize(std::move(topology), ac.threads);
            ac.k = d_oracle.current_k();
            finest_k = ac.k;
            finest_width = d_oracle.current_W();
        }

        // Sets up the coarsening and refinement helpers and the partition state. They are sized by the graph with one
        // vertex per PE, so this happens once the input graph has been reduced to it.
        void initialize_structures() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "initialize_structures");
            const graph_t &g = solver_graph();

            gpa_matcher.initialize(g.n, g.m, finest_k, ac.threads, random_engine, ac.global_path_algorithm_config);
            size_constrained_lp.initialize(g.n, g.m, finest_k, random_engine.get_u64(), ac.size_constrained_lp_config);

            // Note: the non-LARGE_K label propagation allocates a k x k table (10^10 entries at k = 99856),
            // so it is not set up when it is not used.
            lp_active = ac.waver_label_propagation_config.enabled && ac.do_multilevel;
            if (lp_active) {
                lp_refine.initialize(random_engine.get_u64(), ac.waver_label_propagation_config);
            }

            reinit_structures_for_new_k(); // for the coarsest topology
            HEAVYASSERT(assert_state_pre_partitioning(g, p_manager, ac.k));
        }

        // First coarsening step: if the input graph has more vertices than there are PEs, partition it into k blocks
        // (KaMinPar or spectral) and contract them into the quotient graph (one unit-weight vertex per block, summed
        // edge weights), so that the graph has one vertex per PE from here on.
        void reduce_to_quotient() {
            const partition_t k = finest_k;
            if (graphs[0].n <= k) { return; }
            HEIPROMAP_PROFILE_SCOPE("partition", "WaverSolver", "reduce_to_quotient");

            QuotientResult r;
            if (ac.prepartitioner == "spectral3d") {
                const size_t W = finest_width;
                if (W == 0 || k % W != 0) {
                    std::cout << "Error: --prepartitioner spectral3d needs a full grid topology (k = " << k
                            << ", inferred width " << W << ")" << std::endl;
                    exit(EXIT_FAILURE);
                }
                r.part = SpectralPrepartitioner::partition(graphs[0], W, k / W, ac.seed, ac.threads);
                r.quotient = contract_into_quotient(graphs[0], r.part, k, ac.threads);
            } else {
                r = partition_into_quotient(graphs[0], k, ac.imbalance, ac.kaminpar_threads,
                                            ac.threads, ac.seed, ac.kaminpar_preset);
            }

            mappings.emplace_back();
            mappings.back().initialize(graphs[0].n);
            mappings.back().set_coarse_n(k);
            for (vertex_t u = 0; u < graphs[0].n; ++u) { mappings.back().set(u, r.part[u]); }
            graphs.emplace_back(std::move(r.quotient));
            reduced = true;
        }

        // ------------------------------------------------------------------------------------------------
        // coarsening
        // ------------------------------------------------------------------------------------------------

        // Coarsens the topology oracle as far as it goes.
        void coarsen_topology() {
            while (d_oracle.can_coarsen()) {
                d_oracle.coarsen();
                ac.k = d_oracle.current_k();
            }
        }

        // Coarsens the graph: first the input graph is reduced to its quotient with one vertex per PE (if it has
        // more vertices), then the quotient is coarsened until it has about initial_c vertices per (coarse) PE.
        void coarsen_graph() {
            reduce_to_quotient();
            initialize_structures();

            u64 level = 0;
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
        }

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
            sync_with_spectral_solution();
        }

        // Uncoarsens one thing: a topology level if the coarse graph allows it, otherwise a graph level. Never
        // both at once. Returns false once nothing is left to uncoarsen down to the quotient graph.
        bool uncoarsen_one_level() {
            if (can_uncoarsen_topology()) {
                uncoarsen_topology();
                return true;
            }
            if (!graph_levels_left()) { return false; }

            if (ac.do_spectral && !(reduced && mappings.size() == 1)) { spectral.prolong(mappings.back(), graphs[graphs.size() - 2]); }
            uncoarsening();
            return true;
        }

        // Restores a valid, refined partition after a level change.
        void settle_level() {
            legalize(mappings.size());
            refinement();
            sync_with_spectral_solution();
        }

        // Both mode: if the spectral placement of this level is cheaper than the current multilevel partition,
        // adopt it (rebalanced and refined); if the refined placement is still not cheaper, keep the multilevel
        // partition. The spectral tracker itself is never overwritten.
        void sync_with_spectral_solution() {
            const graph_t &g = graphs.back();
            if (!ac.do_spectral || spectral.size() != g.n) { return; }
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "sync_with_spectral_solution");

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
        // least as many vertices as it has children, and (while the graph is coarse) the graph is not too small for
        // the finer topology.
        bool can_uncoarsen_topology() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "can_uncoarsen_topology");
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
            if (graph_levels_left() && graphs.back().n < static_cast<u64>(ac.topo_split_factor) * fine_k_val) { return false; }
            for (partition_t c = 0; c < coarse_k; ++c) {
                if (p_manager.get_bweight(c) < static_cast<weight_t>(children_count[c])) { return false; }
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

            reinit_structures_for_new_k();
            p_manager.reset_weights();

            const TopologyUncoarsener::CostParams costs{ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link};
            TopologyUncoarsener(graphs.back(), p_manager, d_oracle, costs, coarse_to_fine).run(old_partition);

            timings.uncontraction_ms += get_milli_seconds(p, get_time_point());
        }

        // The state of the quotient levels is sized by the quotient graph; the input level needs arrays for all of its
        // vertices. The partition of the current (coarse) level is carried over.
        void grow_state_to(const graph_t &g) {
            const std::vector<partition_t> coarse = snapshot();
            p_manager.initialize(g.n, ac.k, g.g_weight);
            p_manager.activate_all();
            p_manager.update_limits(d_oracle.current_topology(), g.g_weight, g.n, ac.imbalance, kMinImbalance);
            for (vertex_t v = 0; v < coarse.size(); ++v) { p_manager.partition[v] = coarse[v]; }
            bv_manager.initialize(g.n, ac.k);
            rebalancer.initialize(g.n, g.m, ac.k, random_engine.get_u64());
        }

        void uncoarsening() {
            auto p = get_time_point();

            const graph_t &finer = graphs[graphs.size() - 2];
            const bool to_input = reduced && mappings.size() == 1;
            if (to_input) { grow_state_to(finer); }
            p_manager.uncontract(mappings.back());
            if (to_input) { p_manager.recalculate_weights(finer); }

            HEIPROMAP_PROFILE_SCOPE("uncontraction", "misc", "compute_from_scratch");
            bv_manager.reset();
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

        // True if the finished mapping is a bijection on the finest level: one vertex per PE, nothing else.
        bool is_bijection() const {
            if (graphs.size() != base_levels() + 1 || d_oracle.depth() != 0 || graphs.back().n != ac.k) { return false; }
            for (partition_t id = 0; id < ac.k; ++id) {
                if (p_manager.get_bweight(id) != 1) { return false; }
            }
            return true;
        }

        // One swap refinement pass on the current level of the V-cycle.
        void vcycle_swap_pass() {
            if (!ac.waver_vcycle_swap_refinement_config.enabled) { return; }
            WaverSwapRefinement::swap_refine(graphs.back(), p_manager, d_oracle,
                                             ac.waver_vcycle_swap_refinement_config,
                                             ac.alpha_end(), ac.beta_end(),
                                             ac.alpha_link, ac.beta_link);
        }

        // V-cycle on a finished bijection. The topology is coarsened (blocks of cells) and the vertices mapped to the merged
        // cells are contracted the same way, so every coarse vertex sits in its own coarse cell. Uncoarsening then
        // restores graph and topology level by level, with a swap refinement pass on every level, starting at the
        // coarsest. The previous mapping is restored if the cycle does not lower the message cost.
        void vcycle() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "vcycle");
            const std::vector<partition_t> saved = snapshot();
            const f64 cost_before = current_message_cost();

            const u64 levels = coarsen_mapped_bijection();

            recompute_datastructures();
            vcycle_swap_pass();
            for (u64 l = 0; l < levels; ++l) {
                uncoarsening();
                uncoarsen_topology();
                recompute_datastructures();
                rebalance(1);
                EmptyBlockDiffusion::diffuse(graphs.back(), p_manager, static_cast<size_t>(d_oracle.current_W()));
                vcycle_swap_pass();
            }

            if (!is_bijection() || current_message_cost() >= cost_before) {
                load_partition(saved);
                recompute_datastructures();
            }
        }

        // Coarsens topology and mapped vertices together; returns the number of levels.
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
                reinit_structures_for_new_k();
                p_manager.reset_weights();
                for (vertex_t c = 0; c < graphs.back().n; ++c) { p_manager.set(c, graphs.back().v_weight(c), c); }
                ++levels;
            }
            return levels;
        }

        struct LevelMeasure {
            f64 cost = 0.0;
            f64 max_link = 0.0;
            weight_t max_rows = 0;
            weight_t min_rows = 0;
        };

        // Message cost, maximum link load and PE weights of the current level.
        LevelMeasure measure_level(const bool recompute) {
            if (recompute) { recompute_datastructures(); }
            LevelMeasure m;
            m.cost = get_message_cost(q_graph, ac.k, d_oracle, ac.alpha_end(), ac.beta_end(),
                                      ac.alpha_link, ac.beta_link).total;
            m.max_link = get_max_link_load(q_graph, ac.k, d_oracle.current_W(), ac.alpha_link, ac.beta_link);
            m.max_rows = m.min_rows = p_manager.get_bweight(0);
            for (partition_t id = 1; id < ac.k; ++id) {
                m.max_rows = std::max(m.max_rows, p_manager.get_bweight(id));
                m.min_rows = std::min(m.min_rows, p_manager.get_bweight(id));
            }
            return m;
        }

        // The last uncoarsening level: the mapping of the quotient graph is projected onto the individual input
        // vertices, then legalized and refined like every other level. Here the PEs may deviate from the average
        // weight by --imbalance.
        void refine_input_level() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "refine_input_level");
            const auto p = get_time_point();

            uncoarsening();
            const LevelMeasure before = measure_level(false);
            settle_level();
            const LevelMeasure after = measure_level(true);

            input_level.ran = true;
            input_level.cost_before = before.cost;
            input_level.cost_after = after.cost;
            input_level.max_link_before = before.max_link;
            input_level.max_link_after = after.max_link;
            input_level.max_rows_before = before.max_rows;
            input_level.max_rows_after = after.max_rows;
            input_level.min_rows_after = after.min_rows;
            input_level.seconds = get_seconds(p, get_time_point());
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
