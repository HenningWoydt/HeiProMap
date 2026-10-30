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

#ifndef HEIPROMAP_WAVERBASE_SOLVER_H
#define HEIPROMAP_WAVERBASE_SOLVER_H

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <thread>
#include <iostream>
#include <numeric>
#include <type_traits>
#include <vector>

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
#include "../refinement/waver_message_refinement.h"
#include "../refinement/waver_long_message_refinement.h"
#include "../refinement/waver_sa_refinement.h"
#include "../refinement/waver_swap_refinement.h"
#include "../distance_oracles/grid_distance_oracle.h"

namespace HeiProMap {
    // Shared pipeline of the Waver solvers: maps a graph onto a PE grid, minimizing the message cost. The solvers differ
    // only in how the coarsest level is mapped and in what runs besides the common refinements (see the hooks).
    //
    // Pipeline:
    //   1. prepartition   n > k: partition the input graph into k blocks and contract them into the quotient graph with
    //                     one vertex per PE. Everything up to 6. works on a graph with n == k.
    //   2. coarsen        the topology and the graph
    //   3. initial        coarsest level: initial_partitioning() of the solver, then refined
    //   4. uncoarsen      level by level, legalize and refine on every level. The mapping is a bijection (one vertex
    //                     per PE) on the finest level.
    //   5. after          after_uncoarsening() of the solver
    //   6. V-cycle        optional, on the finished bijection
    //   7. undo prepartition  projects the mapping onto the input vertices, which may differ in weight by --imbalance,
    //                     then legalizes and refines that level like any other
    class WaverSolverBase {
    protected:
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
        // Not maintained: the rebalancer takes them as arguments but does not read them with the tags used here.
        // solve() builds q_graph once for the report.
        LargeQuotientGraph q_graph;
        block_conn_t block_conn;
        GridDistanceOracle d_oracle;
        // The boundary, quotient graph and block connectivity match the partition. Cleared whenever the partition
        // changes without updating them (topology uncoarsening, loading a partition, the refinements that relabel PEs).
        bool structures_fresh = false;

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

    public:
        WaverSolverBase(graph_t &&g, graph_t &&t_topology, const WaverHeiProMapConfiguration &t_ac) {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "constructor");
            sp = std::chrono::high_resolution_clock::now();

            ac = t_ac;
            graphs.reserve(100);
            max_n = g.n;
            max_m = g.m;
            input_graph = std::move(g);
            random_engine = RandomEngine(ac.seed);

            d_oracle.initialize(std::move(t_topology), ac.threads);
            ac.k = d_oracle.current_k();

            // The helpers are sized for the input graph, the largest graph, and the finest topology; they are used on
            // the smaller graphs and coarser topologies as well.
            gpa_matcher.initialize(max_n, max_m, ac.k, ac.threads, random_engine, ac.global_path_algorithm_config);
            size_constrained_lp.initialize(max_n, max_m, ac.k, random_engine.get_u64(), ac.size_constrained_lp_config);
        }

        virtual ~WaverSolverBase() = default;

        std::vector<vertex_t> solve() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "solve");
            std::vector<partition_t> result = internal_solve();

            const auto ep = std::chrono::high_resolution_clock::now();
            const f64 duration = get_seconds(sp, ep);
            // Everything the phase timers do not cover: constructor (oracle, helper allocation), structure rebuilds, etc.
            timings.misc_ms = std::max(0.0, duration * 1000.0 - timings.total_ms());

            // The report only needs the quotient graph; the boundary is not kept up to date after the last refinement.
            q_graph.compute_from_scratch(graphs[0], p_manager);
            HEAVYASSERT(assert_bweights(graphs[0], p_manager, ac.k));
            const MessageCostBreakdown cost = get_message_cost(q_graph, ac.k, d_oracle, ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
            const f64 max_link_load = get_max_link_load(q_graph, ac.k, d_oracle.current_W(), ac.alpha_link, ac.beta_link);
            const EndpointLoad endpoint = get_max_endpoint_load(q_graph, ac.k, ac.alpha_end(), ac.beta_end());
            const weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) graphs[0].g_weight / (f64) ac.k));

            print_waver_summary(std::cout, ac, graphs[0], duration, lmax, cost, max_link_load, endpoint, p_manager, timings);

            return result;
        }

    protected:
        // ------------------------------------------------------------------------------------------------
        // hooks of the solvers
        // ------------------------------------------------------------------------------------------------

        // Maps the coarsest level (graphs.back()) and settles it: set the partition, then finish_initial_partition().
        virtual void initial_partitioning() = 0;

        // Called before a graph level is uncoarsened: the mapping maps the finer graph to the current one.
        virtual void on_graph_uncoarsen(const Mapping &, const graph_t &) {
        }

        // Called once on the finest level of the bijection, after the last level is settled.
        virtual void after_uncoarsening() {
        }

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

            after_uncoarsening();

            if (ac.vcycle) { vcycle(); }

            undo_prepartition();
            corner_refinement();
            return snapshot();
        }

        // ------------------------------------------------------------------------------------------------
        // setup
        // ------------------------------------------------------------------------------------------------

        // Resets the partition state for the current topology level: its number of cells and its weight limits. The
        // arrays are sized for the input graph, the largest graph, so they fit every level. The weights and limits are
        // those of the graph with one vertex per PE (graphs[0]).
        void enter_topology_level() {
            structures_fresh = false;
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
            const auto t_start = get_time_point();
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
                                                           ac.spectral_smoothing, ac.prepartition_rotations, ac.prepartition_dims);
                r.quotient = contract_into_quotient(input_graph, r.part, k, ac.threads);
            } else {
                r = partition_into_quotient(input_graph, k, ac.imbalance, ac.kaminpar_threads,
                                            ac.threads, ac.seed, ac.kaminpar_preset);
            }

            input_mapping.initialize(input_graph.n);
            input_mapping.set_coarse_n(k);
            for (vertex_t u = 0; u < input_graph.n; ++u) { input_mapping.set(u, r.part[u]); }
            graphs.emplace_back(std::move(r.quotient));
            timings.prepartition_ms += get_milli_seconds(t_start, get_time_point());
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

        // The same steps as on every other level, after the solver has set the partition of the coarsest level.
        void finish_initial_partition() {
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

            on_graph_uncoarsen(mappings.back(), graphs[graphs.size() - 2]);
            uncoarsening();
            return true;
        }

        // Restores a valid, refined partition after a level change.
        void settle_level() {
            legalize(mappings.size());
            refinement();
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

        // recalculate_weights: the weights of the finer graph are not those of the coarser one (the prepartition level),
        // so the PE weights are recomputed before the derived state is checked.
        void uncoarsening(const bool recalculate_weights = false) {
            auto p = get_time_point();

            const graph_t &finer = graphs[graphs.size() - 2];
            p_manager.uncontract(mappings.back());
            if (recalculate_weights) { p_manager.recalculate_weights(finer); }

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

            structures_fresh = true;
            timings.uncontraction_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(state_is_consistent());
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
            structures_fresh = true;
        }

        void rebalance(const u64 level) {
            auto p = get_time_point();

            if (level == 0) {
                rebalancer.rebalance_last_layer<false, false>(graphs.back(), p_manager, bv_manager, q_graph, d_oracle, block_conn);
            } else {
                rebalancer.rebalance<false, false>(graphs.back(), p_manager, bv_manager, q_graph, d_oracle, block_conn);
            }

            timings.rebalance_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(state_is_consistent());
        }

        // Makes the partition valid for the grid: rebalance, fill empty PEs, rebuild the derived state.
        // level is the number of graph levels still to uncoarsen.
        void legalize(const u64 level) {
            // The rebalancer reads the boundary (and keeps all structures up to date while it moves vertices), so they
            // have to match the partition first: they do not after a topology level or a refinement.
            if (!structures_fresh) { recompute_datastructures(); }
            rebalance(level);
            const size_t moved = EmptyBlockDiffusion::diffuse(graphs.back(), p_manager, static_cast<size_t>(d_oracle.current_W()));
            // The diffusion does not keep them up to date.
            if (moved > 0) { recompute_datastructures(); }
        }

        // All refinements, in order. Each is switched on by its configuration. Simulated annealing swaps whole PEs.
        void refinement() {
            auto p = get_time_point();

            if (ac.waver_message_refinement_config.enabled) {
                auto sp_local = get_time_point();
                WaverMessageRefinement::refine(graphs.back(), p_manager, bv_manager, d_oracle,
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

            if (ac.waver_sa_refinement_config.enabled) {
                auto sp_local = get_time_point();
                WaverSARefinement::refine(graphs.back(), p_manager, d_oracle,
                                          ac.waver_sa_refinement_config, ac.seed,
                                          ac.alpha_end(), ac.beta_end(),
                                          ac.alpha_link, ac.beta_link);
                timings.sa_refine_ms += get_milli_seconds(sp_local, get_time_point());
            }

            structures_fresh = false; // the swap refinements relabel PEs in p_manager only
            timings.refinement_ms += get_milli_seconds(p, get_time_point());
            // The boundary, quotient graph and block connectivity are stale here by design, so only the partition
            // and its weights can be checked.
            HEAVYASSERT(assert_bweights(graphs.back(), p_manager, ac.k));
        }

        // The finished mapping of the input graph: single vertices move or swap between neighboring PEs, so that the
        // busiest engines get less work (see WaverCornerRefinement).
        void corner_refinement() {
            if (!ac.waver_corner_refinement_config.enabled) { return; }
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverSolver", "corner_refinement");
            const auto t_start = get_time_point();
            const graph_t &g = graphs.back();
            std::vector<partition_t> part = snapshot();
            const weight_t lmax = static_cast<weight_t>(std::ceil((1.0 + imbalance) * static_cast<f64>(g.g_weight) / static_cast<f64>(ac.k)));
            WaverCornerRefinement::refine(g, part, ac.k, lmax, ac.waver_corner_refinement_config);
            load_partition(part);
            timings.corner_refine_ms += get_milli_seconds(t_start, get_time_point());
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
            uncoarsening(true);

            imbalance = ac.imbalance;
            const graph_t &input = graphs.back();
            p_manager.update_limits(d_oracle.current_topology(), input.g_weight, input.n, imbalance, kMinImbalance);
            settle_level();
        }

        // The boundary matches the partition.
        // p_manager.size() is not kept up to date by the swap refinements, so the vertices per PE are not checked.
        bool state_is_consistent() const {
            return assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, ac.k, 1, false);
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
            structures_fresh = false;
            for (vertex_t u = 0; u < graphs.back().n; ++u) p_manager.partition[u] = part[u];
            p_manager.recalculate_weights(graphs.back());
        }

        template<typename PartT>
        f64 placement_cost(const PartT &part) {
            return SpectralTracker::message_cost(graphs.back(), part, d_oracle,
                                                 ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
        }

        // Computed from the graph and the partition: needs none of the derived structures, so they stay as they are.
        f64 current_message_cost() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "current_message_cost");
            return placement_cost(p_manager.partition);
        }
    };
}

#endif //HEIPROMAP_WAVERBASE_SOLVER_H
