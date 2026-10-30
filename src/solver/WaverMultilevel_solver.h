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

#ifndef HEIPROMAP_WAVERMULTILEVEL_SOLVER_H
#define HEIPROMAP_WAVERMULTILEVEL_SOLVER_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
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
#include "../utility/waver_report.h"
#include "../utility/utils.h"

#include "../coarsening/global_path_algorithm.h"
#include "../coarsening/size_constrained_lp.h"
#include "../coarsening/heavy_edge_matching.h"
#include "../rebalance/rebalancer.h"
#include "../rebalance/empty_block_diffusion.h"
#include "../utility/grid_utils.h"
#include "../uncoarsening/topology_uncoarsener.h"
#include "../refinement/waver_message_refinement.h"
#include "../refinement/waver_long_message_refinement.h"
#include "../refinement/waver_merge_refinement.h"
#include "../refinement/waver_busy_refinement.h"
#include "../refinement/waver_cycle_refinement.h"
#include "../refinement/waver_crowd_refinement.h"
#include "../refinement/waver_shift_refinement.h"
#include "../refinement/waver_smooth_refinement.h"
#include "../refinement/waver_sa_refinement.h"
#include "../refinement/waver_swap_refinement.h"
#include "../distance_oracles/grid_distance_oracle.h"
#include "../init_partitioning/waver_greedy_partitioner.h"

namespace HeiProMap {
    // Multilevel mapping of the input graph directly onto the PE grid, without a prepartition: there is no quotient
    // graph, so the input vertices (not blocks) are coarsened, mapped and uncoarsened, and --imbalance bounds the PE
    // weights on every level.
    //
    // Pipeline:
    //   1. coarsen        the topology as far as it goes, then the graph until it has about initial_c vertices per
    //                     (coarse) PE
    //   2. initial        coarsest level: the cheapest of initial_kappa greedy partitions and the spectral placement,
    //                     then rebalanced and refined
    //   3. uncoarsen      level by level (a topology level when the coarse graph allows it, otherwise a graph level),
    //                     legalize and refine on every level
    class WaverMultilevelSolver {
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
        // The boundary matches the partition. Cleared whenever the partition changes without updating it (topology
        // uncoarsening, the refinements that relabel PEs).
        bool structures_fresh = false;

        // cell_fine_counts[d][c]: the number of finest PEs merged into cell c of the topology at depth d. The balance a
        // coarse cell may have is that of its finest PEs together, so its limits are sums over these.
        std::vector<std::vector<partition_t> > cell_fine_counts;

        vertex_t max_n = 0; // size of the input graph, the largest graph: the state is allocated for it
        vertex_t max_m = 0;

        // ---- coarsening ----
        GlobalPathAlgorithmMatcher gpa_matcher;
        SizeConstrainedLP size_constrained_lp;
        HeavyEdgeMatching heavy_edge_matching;

        // ---- rebalancing ----
        Rebalancer rebalancer;

        // ---- initial placement ----
        SpectralTracker spectral; // the spectral placement of the coarsest level, a candidate besides the greedy ones

    public:
        WaverMultilevelSolver(graph_t &&g, graph_t &&t_topology, const WaverHeiProMapConfiguration &t_ac) {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "constructor");
            sp = std::chrono::high_resolution_clock::now();

            ac = t_ac;
            graphs.reserve(100);
            max_n = g.n;
            max_m = g.m;
            graphs.emplace_back(std::move(g));
            random_engine = RandomEngine(ac.seed);

            SpectralSmoothing smoothing = ac.spectral_smoothing;
            smoothing.threads = static_cast<int>(ac.threads);
            spectral.set_smoothing(smoothing);

            d_oracle.initialize(std::move(t_topology), ac.threads);
            ac.k = d_oracle.current_k();
            cell_fine_counts.emplace_back(ac.k, 1);

            // The helpers are sized for the input graph, the largest graph, and the finest topology; they are used on
            // the smaller graphs and coarser topologies as well.
            gpa_matcher.initialize(max_n, max_m, ac.k, ac.threads, random_engine, ac.global_path_algorithm_config);
            size_constrained_lp.initialize(max_n, max_m, ac.k, random_engine.get_u64(), ac.size_constrained_lp_config);
        }

        // Runs the whole pipeline and returns the PE of every vertex of the input graph.
        std::vector<vertex_t> solve() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "solve");

            // Coarsen the topology as far as it goes, then the graph until it has about initial_c vertices per
            // (coarse) PE.
            while (d_oracle.can_coarsen()) {
                d_oracle.coarsen();
                track_coarsened_cells();
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

            if (!ac.dump_hierarchy.empty()) { dump_hierarchy(); }
            greedy_partitioning();
            print_level("greedy partition");
            spectral_placement();
            print_level("initial partition (greedy / spectral)");
            finish_initial_partition();
            print_level("coarsest level refined");

            while (uncoarsen_one_level()) { settle_level(); }

            if (d_oracle.depth() != 0) {
                std::cerr << "Error: uncoarsening stopped with the topology at depth " << d_oracle.depth()
                        << " (" << ac.k << " cells instead of " << cell_fine_counts.front().size() << " PEs)" << std::endl;
                const Mapping &topo_map = d_oracle.top_mapping();
                std::vector<size_t> children(ac.k, 0), count(ac.k, 0);
                for (partition_t f = 0; f < d_oracle.fine_k(); ++f) { ++children[topo_map.get(f)]; }
                for (vertex_t u = 0; u < graphs.back().n; ++u) { ++count[p_manager[u]]; }
                size_t bad = 0, worst = 0;
                for (partition_t c = 0; c < ac.k; ++c) {
                    if (count[c] < children[c]) { ++bad; if (bad == 1 || children[c] - count[c] > worst) { worst = children[c] - count[c]; } }
                }
                for (partition_t c = 0, shown = 0; c < ac.k && shown < 2; ++c) {
                    if (count[c] >= children[c]) { continue; }
                    ++shown;
                    partition_t nb[4];
                    const size_t nn = grid_neighbors(c, d_oracle.current_W(), ac.k, nb);
                    std::cerr << "  cell " << c << ": vertices " << count[c] << ", bweight " << p_manager.bweights[c]
                            << ", lmin " << p_manager.get_lmin(c) << ", lmax " << p_manager.get_lmax(c)
                            << ", children " << children[c] << ", counted fine " << cell_fine_counts.back()[c]
                            << ", W " << d_oracle.current_W() << ", neighbors:";
                    for (size_t i = 0; i < nn; ++i) { std::cerr << " " << nb[i] << "(" << count[nb[i]] << ")"; }
                    std::cerr << std::endl;
                }
                std::cerr << "  cells with fewer vertices than child cells: " << bad << ", largest shortfall " << worst
                        << ", graph vertices " << graphs.back().n << ", lmin of cell 0 " << p_manager.get_lmin(0) << std::endl;
                exit(EXIT_FAILURE);
            }

            if (ac.smash_rounds > 0) { smash(); }
            if (ac.vcycle) { vcycle(); }

            std::vector<partition_t> result = snapshot();

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
            const weight_t lmax = std::ceil((1.0 + ac.imbalance) * ((f64) graphs[0].g_weight / (f64) ac.k));

            print_waver_summary(std::cout, ac, graphs[0], duration, lmax, cost, max_link_load, endpoint, p_manager, timings);

            return result;
        }

    private:
        // ------------------------------------------------------------------------------------------------
        // setup
        // ------------------------------------------------------------------------------------------------

        // Resets the partition state for the current topology level: its number of cells and its weight limits. The
        // arrays are sized for the input graph, the largest graph, so they fit every level. The weights and limits are
        // those of the input graph (graphs[0]).
        void enter_topology_level() {
            structures_fresh = false;
            const graph_t &g = graphs[0];
            p_manager.initialize(max_n, ac.k, g.g_weight);
            p_manager.activate_all();
            set_cell_limits(g.g_weight);
            bv_manager.initialize(max_n, ac.k);
            rebalancer.initialize(max_n, max_m, ac.k, random_engine.get_u64());
        }

        // Records how many finest PEs each cell of the topology just coarsened by the oracle stands for.
        void track_coarsened_cells() {
            const Mapping &topo_map = d_oracle.top_mapping();
            const std::vector<partition_t> &fine = cell_fine_counts.back();
            std::vector<partition_t> coarse(d_oracle.current_k(), 0);
            for (partition_t f = 0; f < fine.size(); ++f) { coarse[topo_map.get(f)] += fine[f]; }
            cell_fine_counts.push_back(std::move(coarse));
        }

        // The allowed balance of a cell is that of the finest PEs it stands for, added up: at most the finest weight
        // limit each, and at least --min-fill of the average weight each (but one unit). Without the lower limit the
        // refinement drains cells to save messages, and a drained cell has too few vertices to be split into its child
        // cells, so the topology could only be uncoarsened on the finest graph.
        void set_cell_limits(const weight_t total_weight) {
            const std::vector<partition_t> &counts = cell_fine_counts.back();
            const f64 finest_k = static_cast<f64>(cell_fine_counts.front().size());
            const weight_t finest_lmax = std::ceil((1.0 + ac.imbalance) * (static_cast<f64>(total_weight) / finest_k));
            const f64 finest_avg = static_cast<f64>(total_weight) / finest_k;
            for (partition_t c = 0; c < ac.k; ++c) {
                const weight_t floor_weight = std::max<weight_t>(1, std::floor(ac.min_fill * finest_avg));
                p_manager.set_lmax(c, counts[c] * finest_lmax);
                p_manager.set_lmin(c, counts[c] * std::min(floor_weight, finest_lmax));
            }
        }

        // ------------------------------------------------------------------------------------------------
        // coarsening
        // ------------------------------------------------------------------------------------------------

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

        // ------------------------------------------------------------------------------------------------
        // initial partition
        // ------------------------------------------------------------------------------------------------

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
                waver_greedy_partition(g, d_oracle, ac.imbalance, trial_seed, cur_pm, vertices, ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
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

        // The spectral placement of the coarsest level (Laplacian eigenvectors cut by recursive coordinate bisection onto
        // the coarse grid) replaces the greedy partition if it has a lower message cost.
        bool spectral_placement() {
            auto p = get_time_point();
            HEIPROMAP_PROFILE_SCOPE("partition", "partition", "spectral_placement");

            const graph_t &g = graphs.back();
            spectral.init(g, ac.seed);
            std::vector<partition_t> part;
            spectral.place(g, d_oracle.current_topology(), d_oracle.current_W(), part);
            // The turn of the two smoothest vectors with the lowest message cost (see place_rotated).
            f64 turned_cost = SpectralTracker::message_cost(g, part, d_oracle, ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
            for (u64 step = 1; step < ac.spectral_rotations; ++step) {
                const f64 theta = 0.5 * 3.14159265358979323846 * static_cast<f64>(step) / static_cast<f64>(ac.spectral_rotations);
                std::vector<partition_t> candidate;
                spectral.place_rotated(g, d_oracle.current_topology(), d_oracle.current_W(), theta, candidate);
                const f64 c = SpectralTracker::message_cost(g, candidate, d_oracle, ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
                if (c < turned_cost) {
                    turned_cost = c;
                    part = std::move(candidate);
                }
            }

            const f64 greedy_cost = SpectralTracker::message_cost(g, p_manager.partition, d_oracle, ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
            const f64 spectral_cost = SpectralTracker::message_cost(g, part, d_oracle, ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
            const bool replace = spectral_cost < greedy_cost;
            if (ac.level_stats) {
                std::vector<partition_t> greedy(g.n);
                for (vertex_t u = 0; u < g.n; ++u) { greedy[u] = p_manager.partition[u]; p_manager.partition[u] = part[u]; }
                p_manager.recalculate_weights(g);
                print_level("spectral placement");
                for (vertex_t u = 0; u < g.n; ++u) { p_manager.partition[u] = greedy[u]; }
                p_manager.recalculate_weights(g);
            }
            if (replace) {
                for (vertex_t u = 0; u < g.n; ++u) { p_manager.partition[u] = part[u]; }
                p_manager.recalculate_weights(g);
            }

            timings.spectral_ms += get_milli_seconds(p, get_time_point());
            return replace;
        }

        // The same steps as on every other level, after the partition of the coarsest level is set.
        void finish_initial_partition() {
            recompute_datastructures();
            rebalance(mappings.size());
            refinement();
            fill_underfull_cells();
        }

        // ------------------------------------------------------------------------------------------------
        // uncoarsening
        // ------------------------------------------------------------------------------------------------

        // Uncoarsens one thing: a topology level if the coarse graph allows it, otherwise a graph level. Never
        // both at once. Returns false once nothing is left to uncoarsen down to the input graph.
        bool uncoarsen_one_level() {
            if (can_uncoarsen_topology()) {
                uncoarsen_topology();
                return true;
            }
            if (mappings.empty()) { return false; }

            uncoarsening();
            return true;
        }

        // Mean number of neighbor blocks per PE (blocks with an edge between them), for the diagnostic print.
        f64 mean_block_degree() const {
            const graph_t &g = graphs.back();
            std::vector<std::vector<vertex_t> > blocks(ac.k);
            for (vertex_t u = 0; u < g.n; ++u) { blocks[p_manager[u]].push_back(u); }
            std::vector<partition_t> seen(ac.k, NO_ID);
            f64 total = 0.0;
            for (partition_t x = 0; x < ac.k; ++x) {
                for (const vertex_t u: blocks[x]) {
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const partition_t y = p_manager[g.edges_v[j]];
                        if (y != x && seen[y] != x) {
                            seen[y] = x;
                            total += 1.0;
                        }
                    }
                }
            }
            return total / static_cast<f64>(ac.k);
        }

        // Metrics of the mapping on the current level, printed with --level-stats (and the old one-line trace with
        // --busy-verbose): sizes, balance, message count and distance, the cost, and the exact engine load of the PEs.
        void print_level(const char *what) {
            if (ac.level_stats) {
                print_level_stats(what);
                return;
            }
            if (!ac.waver_busy_refinement_config.verbose) { return; }
            std::cout << "[level] " << what << ": graph n=" << graphs.back().n << " cells k=" << ac.k
                      << " mean neighbor cells " << mean_block_degree() << std::endl;
        }

        // One line per refinement step with --level-stats: what the step did to the mapping, for the analysis of the levels.
        void print_step(const char *what) {
            if (ac.level_stats) { print_level_stats(what, true); }
        }

        void print_level_stats(const char *what, const bool brief = false) {
            const graph_t &g = graphs.back();
            const partition_t k = ac.k;
            WaverRegionLoad tracker;
            tracker.initialize(g, p_manager, k, WaverRegionLoad::Weights());
            std::vector<std::vector<vertex_t> > blocks(k);
            std::vector<f64> weight(k, 0.0);
            for (vertex_t u = 0; u < g.n; ++u) {
                blocks[p_manager[u]].push_back(u);
                weight[p_manager[u]] += static_cast<f64>(g.v_weight(u));
            }
            const int64_t W = static_cast<int64_t>(d_oracle.current_W());
            std::vector<partition_t> seen(k, NO_ID);
            std::vector<f64> degrees(k, 0.0), loads(k, 0.0);
            f64 messages = 0.0, hops = 0.0, cut = 0.0, cut_hops = 0.0, far = 0.0, longest = 0.0;
            size_t empty = 0;
            f64 sr = 0, rr = 0, sv = 0, rv = 0;
            for (partition_t x = 0; x < k; ++x) {
                if (blocks[x].empty()) { ++empty; }
                size_t degree = 0;
                for (const vertex_t u: blocks[x]) {
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const partition_t y = p_manager[g.edges_v[j]];
                        if (y == x) { continue; }
                        const int64_t dx = static_cast<int64_t>(x) % W - static_cast<int64_t>(y) % W;
                        const int64_t dy = static_cast<int64_t>(x) / W - static_cast<int64_t>(y) / W;
                        const f64 d = static_cast<f64>((dx < 0 ? -dx : dx) + (dy < 0 ? -dy : dy));
                        cut += 0.5 * static_cast<f64>(g.edge_w(j));
                        cut_hops += 0.5 * static_cast<f64>(g.edge_w(j)) * d;
                        if (seen[y] != x) {
                            seen[y] = x;
                            ++degree;
                            messages += 0.5;
                            hops += 0.5 * d;
                            if (d > 1.0) { far += 0.5; }
                            longest = std::max(longest, d);
                        }
                    }
                }
                degrees[x] = static_cast<f64>(degree);
                loads[x] = tracker.load(x, degree);
                sr += static_cast<f64>(tracker.send_regions(x));
                rr += static_cast<f64>(tracker.receive_regions(x));
                sv += static_cast<f64>(tracker.send_values(x));
                rv += static_cast<f64>(tracker.receive_values(x));
            }
            const auto stat = [&](std::vector<f64> v, f64 &mean, f64 &p99, f64 &max) {
                mean = 0.0;
                for (const f64 x: v) { mean += x; }
                mean /= static_cast<f64>(v.size());
                std::sort(v.begin(), v.end());
                p99 = v[std::min(v.size() - 1, static_cast<size_t>(0.99 * static_cast<f64>(v.size())))];
                max = v.back();
            };
            f64 dm, dp, dx_, lm, lp, lx, wm, wp, wx;
            stat(degrees, dm, dp, dx_);
            stat(loads, lm, lp, lx);
            stat(weight, wm, wp, wx);
            const f64 kk = static_cast<f64>(k);
            const f64 cost = SpectralTracker::message_cost(g, p_manager.partition, d_oracle, ac.alpha_end(), ac.beta_end(),
                                                           ac.alpha_link, ac.beta_link);
            if (brief) {
                std::cout << "[step] " << std::setw(12) << what << ": weight max/mean=" << wx / wm << " degree=" << dm
                          << " messages=" << messages << " hops=" << (messages > 0 ? hops / messages : 0.0)
                          << " longest=" << longest << " cost=" << cost << " sendreg=" << sr / kk << " recvreg=" << rr / kk
                          << " sendval=" << sv / kk << " recvval=" << rv / kk << " load mean/max=" << lm << "/" << lx
                          << std::endl;
                return;
            }
            std::cout << "[level] " << what << ": n=" << g.n << " k=" << k << " empty=" << empty
                      << " | weight max/mean=" << wx / wm
                      << " | degree mean/p99/max=" << dm << "/" << dp << "/" << dx_
                      << " | messages=" << messages << " mean hops=" << (messages > 0 ? hops / messages : 0.0)
                      << " far(>1)=" << (messages > 0 ? far / messages : 0.0) << " longest=" << longest
                      << " | cut=" << cut << " cut/vertex=" << cut / static_cast<f64>(g.n)
                      << " cut hops=" << cut_hops << " | cost=" << cost
                      << " | per cell sendreg=" << sr / kk << " recvreg=" << rr / kk << " sendval=" << sv / kk
                      << " recvval=" << rv / kk << " | load mean/p99/max=" << lm << "/" << lp << "/" << lx << std::endl;
        }

        // Restores a valid, refined partition after a level change.
        void settle_level() {
            print_level("after the level change");
            legalize(mappings.size());
            refinement();
            print_level("settled");
            dump_partition("k" + std::to_string(ac.k) + "_n" + std::to_string(graphs.back().n) + "_settled");
            // The refinements keep the lower limits only as a limit on what they take from a cell; they can still leave
            // one short, which the next topology level cannot be split from.
            fill_underfull_cells();
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
            cell_fine_counts.pop_back();
            ac.k = d_oracle.current_k();

            const std::vector<partition_t> old_partition = snapshot();

            enter_topology_level();
            p_manager.reset_weights();

            const auto t_split = get_time_point();
            const std::vector<partition_t> assignment = local_spectral_split(old_partition, coarse_to_fine);
            timings.spectral_ms += get_milli_seconds(t_split, get_time_point());

            const TopologyUncoarsener::CostParams costs{ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link};
            TopologyUncoarsener(graphs.back(), p_manager, d_oracle, costs, coarse_to_fine).run_from_assignment(assignment);

            if (ac.global_split) { try_global_split(); }
            dump_partition("split" + std::to_string(ac.k) + "_n" + std::to_string(graphs.back().n) + "_raw");
            timings.uncontraction_ms += get_milli_seconds(p, get_time_point());
        }

        // Places the whole current graph spectrally onto the current (just split) topology; it replaces the partition if
        // its message cost is lower than that of the local splits.
        void try_global_split() {
            auto p = get_time_point();
            const graph_t &g = graphs.back();
            spectral.init(g, ac.seed);
            std::vector<partition_t> part;
            spectral.place(g, d_oracle.current_topology(), d_oracle.current_W(), part);
            const f64 local_cost = SpectralTracker::message_cost(g, p_manager.partition, d_oracle, ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
            const f64 global_cost = SpectralTracker::message_cost(g, part, d_oracle, ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
            if (ac.waver_busy_refinement_config.verbose) {
                std::cout << "[level] global split: local cost " << local_cost << ", global cost " << global_cost << std::endl;
            }
            if (global_cost < local_cost) {
                for (vertex_t u = 0; u < g.n; ++u) { p_manager.partition[u] = part[u]; }
                p_manager.recalculate_weights(g);
            }
            timings.spectral_ms += get_milli_seconds(p, get_time_point());
        }

        // Places the vertices verts (nv of them) onto the rectangle of fine cells with the top left cell
        // (row0_cell, col0_cell) and h rows and w columns, every cell about the same weight: a spectral placement of the
        // subgraph the vertices induce, the mirrored or transposed variant that puts the vertices with edges out of the
        // set closest to their outside neighbors, and (with --terminal-propagation) a pull towards those neighbors.
        // inside(v) tells whether a neighbor v belongs to the vertices; the position of an outside neighbor is its cell in
        // `assignment` if that is set, otherwise outside_centre(v), and outside_group(v) names the group it stands for.
        // The cell of every vertex is returned in result.
        template<typename InsideFn, typename CentreFn, typename GroupFn>
        void place_rectangle(const vertex_t *verts, const size_t nv, const size_t row0_cell, const size_t col0_cell,
                             const size_t h, const size_t w, const size_t fine_w, const InsideFn &inside,
                             const CentreFn &outside_centre, const GroupFn &outside_group,
                             const std::vector<partition_t> &assignment, std::vector<vertex_t> &local_id,
                             const u64 seed, std::vector<partition_t> &result) {
            const graph_t &g = graphs.back();
            SpectralSmoothing smoothing = ac.spectral_smoothing;
            smoothing.threads = 1;
            const size_t row0 = row0_cell, col0 = col0_cell;
            result.assign(nv, NO_ID);
        // The subgraph induced by the cell, and the edges that leave it.
        for (size_t i = 0; i < nv; ++i) { local_id[verts[i]] = static_cast<vertex_t>(i); }
        std::vector<weight_t> vw(nv), ew;
        std::vector<size_t> nb(nv + 1, 0);
        std::vector<vertex_t> ev;
        std::vector<vertex_t> out_vertex; // local vertex of every edge leaving the cell
        std::vector<f64> out_row, out_col, out_w;
        std::vector<vertex_t> out_neighbor; // the vertex outside the cell of every edge leaving it
        for (size_t i = 0; i < nv; ++i) {
            const vertex_t u = verts[i];
            vw[i] = g.v_weight(u);
            for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                const vertex_t v = g.edges_v[j];
                if (inside(v)) {
                    ev.push_back(local_id[v]);
                    ew.push_back(g.edge_w(j));
                } else {
                    out_vertex.push_back(static_cast<vertex_t>(i));
                    out_row.push_back(outside_centre(v).first);
                    out_col.push_back(outside_centre(v).second);
                    out_w.push_back(static_cast<f64>(g.edge_w(j)));
                    out_neighbor.push_back(v);
                }
            }
            nb[i + 1] = ev.size();
        }
        graph_t local;
        local.initialize(static_cast<vertex_t>(nv), static_cast<vertex_t>(ev.size()), vw.data(), nb.data(), ew.data(), ev.data());

        SpectralTracker tracker;
        tracker.set_smoothing(smoothing);
        tracker.init(local, seed);
        graph_t sub_grid; // only its number of cells matters: every child cell has the same capacity
        sub_grid.n = static_cast<vertex_t>(h * w);
        std::vector<partition_t> cell;
        tracker.place(local, sub_grid, w, cell);

        // The best of the mirrored and transposed variants of the placement.
        const size_t variants = h == w ? 8 : 4;
        const auto position = [&](const partition_t p, const size_t variant) {
            size_t r = p / w, q = p % w;
            if (variant & 4) { std::swap(r, q); }
            if (variant & 1) { r = h - 1 - r; }
            if (variant & 2) { q = w - 1 - q; }
            return std::pair<size_t, size_t>(row0 + r, col0 + q);
        };
        size_t best_variant = 0;
        f64 best_cost = std::numeric_limits<f64>::max();
        for (size_t variant = 0; variant < variants; ++variant) {
            f64 cost = 0.0;
            for (size_t e = 0; e < out_vertex.size(); ++e) {
                const auto pos = position(cell[out_vertex[e]], variant);
                cost += out_w[e] * (std::abs(static_cast<f64>(pos.first) - out_row[e])
                                    + std::abs(static_cast<f64>(pos.second) - out_col[e]));
            }
            if (cost < best_cost) {
                best_cost = cost;
                best_variant = variant;
            }
        }
        // Terminal propagation: the vertices are pulled towards their neighbors, inside the cell and outside it,
        // starting from the spectral placement; the neighbors in cells that are split already count at their
        // cell, the others at the centre of their coarse cell. The result is spread over the child cells again
        // by rank, and it replaces the spectral placement if it has a lower cost.
        std::vector<partition_t> chosen(nv);
        for (size_t i = 0; i < nv; ++i) {
            const auto pos = position(cell[i], best_variant);
            chosen[i] = static_cast<partition_t>(pos.first * fine_w + pos.second);
        }
        if (ac.terminal_propagation) {
            std::vector<f64> term_row(out_vertex.size()), term_col(out_vertex.size());
            for (size_t e = 0; e < out_vertex.size(); ++e) {
                const vertex_t v = out_neighbor[e];
                if (assignment[v] != NO_ID) {
                    term_row[e] = static_cast<f64>(assignment[v] / fine_w);
                    term_col[e] = static_cast<f64>(assignment[v] % fine_w);
                } else {
                    term_row[e] = out_row[e];
                    term_col[e] = out_col[e];
                }
            }
            // The message cost of a placement: every pair of cells that exchange values is one message (alpha_end
            // and the link cost of its distance), every cut edge adds its volume. A neighbor outside the cell
            // counts as its cell if that is split already, otherwise as its coarse cell.
            const auto placement_cost = [&](const std::vector<partition_t> &where) {
                std::vector<std::pair<unsigned long long, f64> > pairs;
                f64 total = 0.0;
                const auto edge = [&](const unsigned long long id_a, const f64 ra, const f64 qa, const unsigned long long id_b, const f64 rb,
                                      const f64 qb, const f64 weight) {
                    if (id_a == id_b) { return; }
                    const f64 d = std::abs(ra - rb) + std::abs(qa - qb);
                    pairs.emplace_back(std::min(id_a, id_b) * 4294967296ULL + std::max(id_a, id_b), ac.alpha_end() + ac.alpha_link * d);
                    total += (ac.beta_end() + ac.beta_link * d) * weight;
                };
                for (size_t i = 0; i < nv; ++i) {
                    const f64 r = static_cast<f64>(where[i] / fine_w), q = static_cast<f64>(where[i] % fine_w);
                    for (size_t j = nb[i]; j < nb[i + 1]; ++j) {
                        const partition_t o = where[ev[j]];
                        edge(where[i], r, q, o, static_cast<f64>(o / fine_w), static_cast<f64>(o % fine_w), 0.5 * static_cast<f64>(ew[j]));
                    }
                }
                for (size_t e = 0; e < out_vertex.size(); ++e) {
                    const partition_t o = where[out_vertex[e]];
                    const vertex_t v = out_neighbor[e];
                    const unsigned long long other = assignment[v] != NO_ID ? assignment[v] : (1ULL << 31) + outside_group(v);
                    edge(o, static_cast<f64>(o / fine_w), static_cast<f64>(o % fine_w), other, term_row[e], term_col[e], out_w[e]);
                }
                std::sort(pairs.begin(), pairs.end());
                for (size_t t = 0; t < pairs.size(); ++t) {
                    if (t == 0 || pairs[t].first != pairs[t - 1].first) { total += pairs[t].second; }
                }
                return total;
            };

            std::vector<f64> row(nv), col(nv), row0(nv), col0(nv), degree(nv, 0.0);
            for (size_t i = 0; i < nv; ++i) {
                row0[i] = row[i] = static_cast<f64>(chosen[i] / fine_w);
                col0[i] = col[i] = static_cast<f64>(chosen[i] % fine_w);
                for (size_t j = nb[i]; j < nb[i + 1]; ++j) { degree[i] += static_cast<f64>(ew[j]); }
            }
            for (size_t e = 0; e < out_vertex.size(); ++e) { degree[out_vertex[e]] += out_w[e]; }
            std::vector<std::vector<size_t> > out_of(nv);
            for (size_t e = 0; e < out_vertex.size(); ++e) { out_of[out_vertex[e]].push_back(e); }
            for (int sweep = 0; sweep < 20; ++sweep) {
                for (size_t i = 0; i < nv; ++i) {
                    // the anchor to the spectral position keeps the cell from collapsing onto its terminals
                    const f64 anchor = std::max(degree[i], 1.0);
                    f64 sr = anchor * row0[i], sc = anchor * col0[i], sw = anchor;
                    for (size_t j = nb[i]; j < nb[i + 1]; ++j) {
                        sr += static_cast<f64>(ew[j]) * row[ev[j]];
                        sc += static_cast<f64>(ew[j]) * col[ev[j]];
                        sw += static_cast<f64>(ew[j]);
                    }
                    for (const size_t e: out_of[i]) {
                        sr += out_w[e] * term_row[e];
                        sc += out_w[e] * term_col[e];
                        sw += out_w[e];
                    }
                    row[i] = sr / sw;
                    col[i] = sc / sw;
                }
            }
            // spread by rank: the child rows by the row coordinate, then the cells of a row by the column
            std::vector<size_t> by_row(nv);
            for (size_t i = 0; i < nv; ++i) { by_row[i] = i; }
            std::stable_sort(by_row.begin(), by_row.end(), [&](const size_t x, const size_t y) { return row[x] < row[y]; });
            std::vector<partition_t> moved(nv);
            for (size_t r = 0; r < h; ++r) {
                const size_t lo = r * nv / h, hi = (r + 1) * nv / h;
                std::vector<size_t> line(by_row.begin() + lo, by_row.begin() + hi);
                std::stable_sort(line.begin(), line.end(), [&](const size_t x, const size_t y) { return col[x] < col[y]; });
                for (size_t t = 0; t < line.size(); ++t) {
                    const size_t q = t * w / line.size();
                    moved[line[t]] = static_cast<partition_t>((row0_cell + r) * fine_w + col0_cell + q);
                }
            }
            if (placement_cost(moved) < placement_cost(chosen)) { chosen = std::move(moved); }
        }
        for (size_t i = 0; i < nv; ++i) { result[i] = chosen[i]; }

        }

        // Splits every coarse cell into its child cells: the vertices of the cell, and only those, are placed onto the
        // child cells by a spectral placement of the subgraph they induce (Laplacian eigenvectors cut by recursive
        // coordinate bisection onto the sub-grid of the children). The sign and the axes of such a placement are
        // arbitrary, so the mirrored and, for a square sub-grid, transposed variants are tried and the one that puts the
        // vertices with edges out of the cell closest to the neighboring cells is kept. Returns the fine cell of every
        // vertex; the oracle has already been uncoarsened.
        std::vector<partition_t> local_spectral_split(const std::vector<partition_t> &old_partition,
                                                      const std::vector<std::vector<partition_t> > &coarse_to_fine) {
            HEIPROMAP_PROFILE_SCOPE("spectral", "WaverSolver", "local_spectral_split");
            const graph_t &g = graphs.back();
            const size_t fine_w = d_oracle.current_W();
            const partition_t coarse_k = static_cast<partition_t>(coarse_to_fine.size());

            // The vertices of each coarse cell, and the centre of each cell in fine grid coordinates.
            std::vector<size_t> first(coarse_k + 1, 0);
            for (vertex_t u = 0; u < g.n; ++u) { ++first[old_partition[u] + 1]; }
            for (partition_t c = 0; c < coarse_k; ++c) { first[c + 1] += first[c]; }
            std::vector<vertex_t> members(g.n);
            {
                std::vector<size_t> cursor(first.begin(), first.end() - 1);
                for (vertex_t u = 0; u < g.n; ++u) { members[cursor[old_partition[u]]++] = u; }
            }
            std::vector<f64> centre_row(coarse_k, 0.0), centre_col(coarse_k, 0.0);
            for (partition_t c = 0; c < coarse_k; ++c) {
                for (const partition_t f: coarse_to_fine[c]) {
                    centre_row[c] += static_cast<f64>(f / fine_w);
                    centre_col[c] += static_cast<f64>(f % fine_w);
                }
                centre_row[c] /= static_cast<f64>(coarse_to_fine[c].size());
                centre_col[c] /= static_cast<f64>(coarse_to_fine[c].size());
            }

            SpectralSmoothing smoothing = ac.spectral_smoothing;
            smoothing.threads = 1;
            std::vector<partition_t> assignment(g.n, NO_ID);
            std::vector<vertex_t> local_id(g.n);

            for (partition_t c = 0; c < coarse_k; ++c) {
                const std::vector<partition_t> &children = coarse_to_fine[c];
                const size_t nv = first[c + 1] - first[c];
                if (nv == 0) { continue; }
                if (children.size() == 1) {
                    for (size_t i = first[c]; i < first[c + 1]; ++i) { assignment[members[i]] = children[0]; }
                    continue;
                }

                // The children form a rectangle of the fine grid.
                size_t row0 = children[0] / fine_w, row1 = row0, col0 = children[0] % fine_w, col1 = col0;
                for (const partition_t f: children) {
                    row0 = std::min(row0, f / fine_w);
                    row1 = std::max(row1, f / fine_w);
                    col0 = std::min(col0, f % fine_w);
                    col1 = std::max(col1, f % fine_w);
                }
                const size_t h = row1 - row0 + 1, w = col1 - col0 + 1;
                const size_t row0_cell = row0, col0_cell = col0;

                // Too few vertices (or an irregular set of children) to place: cut the cell into contiguous chunks.
                if (h * w != children.size() || nv < 8) {
                    for (size_t i = 0; i < nv; ++i) {
                        assignment[members[first[c] + i]] = children[i * children.size() / nv];
                    }
                    continue;
                }

                std::vector<partition_t> placed;
                place_rectangle(&members[first[c]], nv, row0_cell, col0_cell, h, w, fine_w,
                                [&](const vertex_t v) { return old_partition[v] == c; },
                                [&](const vertex_t v) { return std::pair<f64, f64>(centre_row[old_partition[v]], centre_col[old_partition[v]]); },
                                [&](const vertex_t v) { return static_cast<unsigned long long>(old_partition[v]); },
                                assignment, local_id, ac.seed + c, placed);
                for (size_t i = 0; i < nv; ++i) { assignment[members[first[c] + i]] = placed[i]; }
            }
            return assignment;
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

            structures_fresh = true;
            timings.uncontraction_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(state_is_consistent());
        }

        // ------------------------------------------------------------------------------------------------
        // rebalancing and refinement
        // ------------------------------------------------------------------------------------------------

        // Rebuilds the boundary vertices from the current partition.
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
            if (ac.diffusion_rebalance) {
                diffusion_rebalance();
                print_step("flow");
            }
            rebalance(level);
            print_step("rebalance");
            const size_t moved = EmptyBlockDiffusion::diffuse(graphs.back(), p_manager, static_cast<size_t>(d_oracle.current_W()), ac.diffusion_equalize);
            print_step("diffuse");
            // The diffusion does not keep them up to date.
            const size_t filled = fill_underfull_cells();
            print_step("fill");
            if (moved > 0 || filled > 0) { recompute_datastructures(); }
        }

        // Writes the cell of every input vertex in the current state to <dump prefix>_<tag>.txt.
        void dump_partition(const std::string &tag) const {
            if (ac.dump_hierarchy.empty()) { return; }
            std::vector<vertex_t> cluster(mappings.empty() ? graphs.back().n : mappings.front().get_old_n());
            for (size_t i = 0; i < cluster.size(); ++i) { cluster[i] = static_cast<vertex_t>(i); }
            for (const Mapping &m: mappings) { for (vertex_t &c: cluster) { c = m.get(c); } }
            std::ofstream out(ac.dump_hierarchy + "_p_" + tag + ".txt");
            for (const vertex_t c: cluster) { out << p_manager[c] << "\n"; }
        }

        // Writes, for every graph level, the cluster (vertex of that level) of every input vertex, one file per level.
        void dump_hierarchy() const {
            std::vector<vertex_t> cluster(mappings.empty() ? 0 : mappings.front().get_old_n());
            for (size_t i = 0; i < cluster.size(); ++i) { cluster[i] = static_cast<vertex_t>(i); }
            for (size_t level = 0; level < mappings.size(); ++level) {
                for (vertex_t &c: cluster) { c = mappings[level].get(c); }
                std::ofstream out(ac.dump_hierarchy + "_" + std::to_string(level + 1) + ".txt");
                for (const vertex_t c: cluster) { out << c << "\n"; }
            }
        }

        // Diffusion rebalancing: the cells outside their weight band get a flow towards the cells that have room, found as
        // the minimum-norm flow on the grid (a conjugate gradient solve of the grid Laplacian), and every unit of flow
        // across a grid edge is made by moving vertices of the sending cell into the receiving one, those with the best
        // gain first. Each vertex moves one hop, so no message gets long; the cells are handled from the highest
        // potential down, so a cell has received what it passes on. Does nothing if every cell is within its limits.
        void diffusion_rebalance() {
            HEIPROMAP_PROFILE_SCOPE("rebalance", "WaverSolver", "diffusion_rebalance");
            const graph_t &g = graphs.back();
            const size_t W = d_oracle.current_W();
            const partition_t k = ac.k;
            p_manager.recalculate_weights(g);
            bool over = false;
            for (partition_t x = 0; x < k && !over; ++x) { over = p_manager.bweights[x] > p_manager.get_lmax(x); }
            if (!over) { return; }

            // the right-hand side: the weight outside the band around the fair share
            std::vector<f64> b(k);
            f64 total_b = 0.0;
            const std::vector<partition_t> &counts = cell_fine_counts.back();
            f64 total_weight = 0.0, total_count = 0.0;
            for (partition_t x = 0; x < k; ++x) {
                total_weight += static_cast<f64>(p_manager.bweights[x]);
                total_count += static_cast<f64>(counts[x]);
            }
            for (partition_t x = 0; x < k; ++x) {
                const f64 fair = total_weight * static_cast<f64>(counts[x]) / total_count;
                const f64 w = static_cast<f64>(p_manager.bweights[x]);
                const f64 target = std::clamp(w, fair * (1.0 - 0.5 * ac.imbalance), fair * (1.0 + 0.5 * ac.imbalance));
                b[x] = w - target;
                total_b += b[x];
            }
            for (partition_t x = 0; x < k; ++x) { b[x] -= total_b / static_cast<f64>(k); }

            // conjugate gradient on L phi = b
            const auto multiply = [&](const std::vector<f64> &v, std::vector<f64> &out) {
                partition_t nbrs[4];
                for (partition_t x = 0; x < k; ++x) {
                    const size_t n = grid_neighbors(x, W, k, nbrs);
                    f64 sum = 0.0;
                    for (size_t i = 0; i < n; ++i) { sum += v[x] - v[nbrs[i]]; }
                    out[x] = sum;
                }
            };
            std::vector<f64> phi(k, 0.0), res = b, dir = b, mv(k);
            f64 rr = 0.0;
            for (const f64 v: res) { rr += v * v; }
            const f64 rr0 = rr;
            for (int it = 0; it < 4000 && rr > 1e-12 * rr0 && rr0 > 0.0; ++it) {
                multiply(dir, mv);
                f64 dm = 0.0;
                for (partition_t x = 0; x < k; ++x) { dm += dir[x] * mv[x]; }
                if (dm <= 0.0) { break; }
                const f64 step = rr / dm;
                f64 rr_new = 0.0;
                for (partition_t x = 0; x < k; ++x) {
                    phi[x] += step * dir[x];
                    res[x] -= step * mv[x];
                    rr_new += res[x] * res[x];
                }
                const f64 beta = rr_new / rr;
                for (partition_t x = 0; x < k; ++x) { dir[x] = res[x] + beta * dir[x]; }
                rr = rr_new;
            }

            std::vector<std::vector<vertex_t> > blocks(k);
            for (vertex_t u = 0; u < g.n; ++u) { blocks[p_manager[u]].push_back(u); }
            std::vector<partition_t> order(k);
            for (partition_t x = 0; x < k; ++x) { order[x] = x; }
            std::sort(order.begin(), order.end(), [&](const partition_t a, const partition_t c) { return phi[a] > phi[c]; });
            std::vector<std::pair<int64_t, vertex_t> > ranked;
            size_t moved = 0;
            for (const partition_t x: order) {
                partition_t nbrs[4];
                const size_t n = grid_neighbors(x, W, k, nbrs);
                for (size_t i = 0; i < n; ++i) {
                    const partition_t y = nbrs[i];
                    f64 flow = phi[x] - phi[y];
                    if (flow <= 0.5) { continue; }
                    // the vertices of x by their gain of moving to y
                    ranked.clear();
                    for (const vertex_t u: blocks[x]) {
                        if (p_manager[u] != x) { continue; }
                        int64_t gain = 0;
                        for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                            const partition_t c = p_manager[g.edges_v[j]];
                            gain += (c == y) - (c == x);
                        }
                        ranked.emplace_back(-gain, u);
                    }
                    std::sort(ranked.begin(), ranked.end());
                    for (const auto &[neg_gain, u]: ranked) {
                        if (flow <= 0.5) { break; }
                        const weight_t w = g.v_weight(u);
                        if (static_cast<f64>(w) > 2.0 * flow) { continue; } // would overshoot more than it helps
                        if (p_manager.bweights[x] - w < p_manager.get_lmin(x)) { break; }
                        p_manager.move_serial(u, w, x, y);
                        blocks[y].push_back(u);
                        flow -= static_cast<f64>(w);
                        ++moved;
                    }
                }
            }
            if (ac.waver_busy_refinement_config.verbose) { std::cout << "[diffusion] moved " << moved << " vertices" << std::endl; }
            recompute_datastructures();
        }

        // Brings every cell up to its lower limit (one unit of weight per finest PE it stands for). Without it, a coarse
        // cell can end up with fewer vertices than child cells and the topology cannot be uncoarsened. The vertices come
        // from the nearest cell on the grid that has weight above its own lower limit and travel along the shortest
        // path, one vertex per step, so that every cell on the way stays within its limits. Returns the number of
        // vertex moves; only the partition and the block weights are updated.
        size_t fill_underfull_cells() {
            const partition_t k = ac.k;
            const graph_t &g = graphs.back();
            // The swap refinements relabel PEs without updating the block weights, so take them from the partition.
            p_manager.recalculate_weights(g);
            bool any = false;
            for (partition_t e = 0; e < k && !any; ++e) { any = p_manager.bweights[e] < p_manager.get_lmin(e); }
            if (!any) { return 0; }

            HEIPROMAP_PROFILE_SCOPE("rebalance", "WaverMultilevel", "fill_underfull_cells");
            std::vector<std::vector<vertex_t> > blk_v(k);
            for (vertex_t u = 0; u < g.n; ++u) { blk_v[p_manager.partition[u]].push_back(u); }

            // Moves one vertex of cell from to cell to that keeps from at or above its lower limit and to within its
            // upper limit. Returns false if there is none.
            const auto move_one = [&](const partition_t from, const partition_t to) {
                std::vector<vertex_t> &list = blk_v[from];
                for (size_t i = list.size(); i-- > 0;) {
                    const vertex_t u = list[i];
                    if (p_manager.partition[u] != from) {
                        list[i] = list.back();
                        list.pop_back();
                        continue;
                    }
                    const weight_t uw = g.v_weight(u);
                    if (p_manager.bweights[from] - uw < p_manager.get_lmin(from)) { continue; }
                    if (p_manager.bweights[to] + uw > p_manager.get_lmax(to)) { continue; }
                    list[i] = list.back();
                    list.pop_back();
                    p_manager.bweights[from] -= uw;
                    p_manager.bweights[to] += uw;
                    p_manager.partition[u] = to;
                    blk_v[to].push_back(u);
                    return true;
                }
                return false;
            };

            std::vector<partition_t> parent(k, NO_ID);
            std::vector<partition_t> visited_in(k, NO_ID); // the deficient cell whose search visited the cell last
            std::vector<partition_t> queue;
            size_t moved = 0;
            for (partition_t e = 0; e < k; ++e) {
                while (p_manager.bweights[e] < p_manager.get_lmin(e)) {
                    // Breadth-first over the grid to the nearest cell with weight above its lower limit.
                    queue.assign(1, e);
                    visited_in[e] = e;
                    partition_t donor = NO_ID;
                    for (size_t head = 0; head < queue.size() && donor == NO_ID; ++head) {
                        const partition_t c = queue[head];
                        partition_t nbrs[4];
                        const size_t n_nbrs = grid_neighbors(c, d_oracle.current_W(), k, nbrs);
                        for (size_t i = 0; i < n_nbrs; ++i) {
                            const partition_t v = nbrs[i];
                            if (visited_in[v] == e) { continue; }
                            visited_in[v] = e;
                            parent[v] = c;
                            if (p_manager.bweights[v] > p_manager.get_lmin(v)) { donor = v; break; }
                            queue.push_back(v);
                        }
                    }
                    if (donor == NO_ID) { break; }

                    // Shift one vertex along the path, starting at the donor so that each cell gives only after it
                    // has received.
                    bool shifted = true;
                    for (partition_t c = donor; c != e && shifted; c = parent[c]) {
                        shifted = move_one(c, parent[c]);
                        moved += shifted;
                    }
                    if (!shifted) { break; }
                }
            }
            return moved;
        }

        // Warm-started relaxation: every vertex starts at the centre of its cell, the positions are smoothed over the graph
        // (each vertex moves towards the weighted mean of its neighbors, held back by an anchor to its start), and the
        // vertices are then given to cells again, every vertex to the nearest cell around its position that still has room,
        // nearest pairs first. The new mapping is kept if its message cost is lower; repeats while that holds.
        void relax_refine() {
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverSolver", "relax_refine");
            const graph_t &g = graphs.back();
            const int64_t W = static_cast<int64_t>(d_oracle.current_W()), k = ac.k, H = k / W;
            struct Offer {
                float distance;
                vertex_t vertex;
                partition_t cell;
            };
            size_t kept = 0;
            for (u64 round = 0; round < ac.relax_rounds; ++round) {
                std::vector<float> r0(g.n), c0(g.n), r(g.n), c(g.n), wsum(g.n, 0.f);
                for (vertex_t u = 0; u < g.n; ++u) {
                    r0[u] = r[u] = static_cast<float>(p_manager[u] / W);
                    c0[u] = c[u] = static_cast<float>(p_manager[u] % W);
                }
                // only neighbors in the same or an adjacent cell pull: the map is folded, and the mean of far positions says nothing
                const auto near = [&](const vertex_t u, const vertex_t v) {
                    const int64_t dr = static_cast<int64_t>(p_manager[u] / W) - static_cast<int64_t>(p_manager[v] / W);
                    const int64_t dc = static_cast<int64_t>(p_manager[u] % W) - static_cast<int64_t>(p_manager[v] % W);
                    return dr >= -1 && dr <= 1 && dc >= -1 && dc <= 1;
                };
                for (vertex_t u = 0; u < g.n; ++u) {
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        if (near(u, g.edges_v[j])) { wsum[u] += static_cast<float>(g.edge_w(j)); }
                    }
                }
                std::vector<float> nr(g.n), nc(g.n);
                for (u64 sweep = 0; sweep < ac.relax_sweeps; ++sweep) {
                    for (vertex_t u = 0; u < g.n; ++u) {
                        const float anchor = static_cast<float>(ac.relax_anchor) * wsum[u];
                        float sr = anchor * r0[u], sc = anchor * c0[u];
                        for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                            if (!near(u, g.edges_v[j])) { continue; }
                            const float w = static_cast<float>(g.edge_w(j));
                            sr += w * r[g.edges_v[j]];
                            sc += w * c[g.edges_v[j]];
                        }
                        const float total = anchor + wsum[u];
                        nr[u] = total > 0.f ? sr / total : r0[u];
                        nc[u] = total > 0.f ? sc / total : c0[u];
                    }
                    r.swap(nr);
                    c.swap(nc);
                }
                // nearest pairs first
                std::vector<Offer> offers;
                offers.reserve(static_cast<size_t>(g.n) * 9);
                for (vertex_t u = 0; u < g.n; ++u) {
                    const int64_t cr = std::clamp<int64_t>(std::llround(r[u]), 0, H - 1), cc = std::clamp<int64_t>(std::llround(c[u]), 0, W - 1);
                    for (int64_t dr = -1; dr <= 1; ++dr) {
                        for (int64_t dc = -1; dc <= 1; ++dc) {
                            const int64_t qr = cr + dr, qc = cc + dc;
                            if (qr < 0 || qr >= H || qc < 0 || qc >= W) { continue; }
                            const float er = r[u] - static_cast<float>(qr), ec = c[u] - static_cast<float>(qc);
                            offers.push_back({er * er + ec * ec, u, static_cast<partition_t>(qr * W + qc)});
                        }
                    }
                }
                std::sort(offers.begin(), offers.end(), [](const Offer &a, const Offer &b) { return a.distance < b.distance; });
                std::vector<partition_t> part(g.n, NO_ID);
                std::vector<weight_t> load(k, 0);
                for (const Offer &o: offers) {
                    if (part[o.vertex] != NO_ID) { continue; }
                    const weight_t w = g.v_weight(o.vertex);
                    if (load[o.cell] + w > p_manager.get_lmax(o.cell)) { continue; }
                    part[o.vertex] = o.cell;
                    load[o.cell] += w;
                }
                for (vertex_t u = 0; u < g.n; ++u) {
                    if (part[u] != NO_ID) { continue; }
                    // the nearest cell with room, searched in growing squares around the position
                    const weight_t w = g.v_weight(u);
                    const int64_t cr = std::clamp<int64_t>(std::llround(r[u]), 0, H - 1), cc = std::clamp<int64_t>(std::llround(c[u]), 0, W - 1);
                    f64 best = 1e300;
                    partition_t cell = NO_ID;
                    for (int64_t radius = 2; radius <= std::max(H, W) && cell == NO_ID; ++radius) {
                        for (int64_t qr = std::max<int64_t>(0, cr - radius); qr <= std::min(H - 1, cr + radius); ++qr) {
                            for (int64_t qc = std::max<int64_t>(0, cc - radius); qc <= std::min(W - 1, cc + radius); ++qc) {
                                const partition_t x = static_cast<partition_t>(qr * W + qc);
                                if (load[x] + w > p_manager.get_lmax(x)) { continue; }
                                const f64 er = r[u] - static_cast<f64>(qr), ec = c[u] - static_cast<f64>(qc);
                                if (er * er + ec * ec < best) {
                                    best = er * er + ec * ec;
                                    cell = x;
                                }
                            }
                        }
                    }
                    if (cell == NO_ID) { cell = p_manager[u]; }
                    part[u] = cell;
                    load[cell] += w;
                }
                const f64 before = SpectralTracker::message_cost(g, p_manager.partition, d_oracle, ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
                const f64 after = SpectralTracker::message_cost(g, part, d_oracle, ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
                if (ac.waver_busy_refinement_config.verbose) {
                    std::cout << "[relax] round " << round << ": cost " << before << " -> " << after << std::endl;
                }
                if (!(after < before)) { break; }
                for (vertex_t u = 0; u < g.n; ++u) { p_manager.partition[u] = part[u]; }
                p_manager.recalculate_weights(g);
                ++kept;
            }
            if (kept > 0) {
                fill_underfull_cells();
                recompute_datastructures();
            }
        }

        // Window re-embedding: every window of ac.window_size x ac.window_size cells is placed again from scratch (the
        // vertices of its cells are spread over its cells by place_rectangle, with the actual cells of the vertices
        // outside as terminals), and the new placement is kept if the message cost of the edges at the window drops and
        // every cell stays within its weight limits. The windows of a round do not overlap; later rounds shift them by
        // half a window.
        void window_refine() {
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverSolver", "window_refine");
            const graph_t &g = graphs.back();
            const size_t W = d_oracle.current_W(), k = ac.k, H = k / W;
            const size_t s = std::max<size_t>(2, ac.window_size);
            std::vector<partition_t> assignment(g.n);
            std::vector<vertex_t> local_id(g.n);
            const auto dist = [&](const partition_t a, const partition_t b) {
                const f64 dr = std::abs(static_cast<f64>(a / W) - static_cast<f64>(b / W));
                const f64 dc = std::abs(static_cast<f64>(a % W) - static_cast<f64>(b % W));
                return dr + dc;
            };
            size_t tried = 0, kept = 0;
            for (u64 round = 0; round < ac.window_rounds; ++round) {
                for (vertex_t u = 0; u < g.n; ++u) { assignment[u] = p_manager[u]; }
                std::vector<std::vector<vertex_t> > blocks(k);
                for (vertex_t u = 0; u < g.n; ++u) { blocks[assignment[u]].push_back(u); }
                const size_t shift = (round % 2 == 1) ? s / 2 : 0;
                size_t round_kept = 0;
                for (size_t r0 = 0; r0 < H + shift; r0 += s) {
                    for (size_t c0 = 0; c0 < W + shift; c0 += s) {
                        const size_t row0 = r0 < shift ? 0 : r0 - shift, col0 = c0 < shift ? 0 : c0 - shift;
                        const size_t row1 = std::min(H, r0 + s - shift), col1 = std::min(W, c0 + s - shift);
                        if (row1 <= row0 || col1 <= col0) { continue; }
                        const size_t h = row1 - row0, w = col1 - col0;
                        if (h * w < 2) { continue; }
                        std::vector<vertex_t> verts;
                        for (size_t r = row0; r < row1; ++r) {
                            for (size_t c = col0; c < col1; ++c) {
                                const auto &b = blocks[r * W + c];
                                verts.insert(verts.end(), b.begin(), b.end());
                            }
                        }
                        const size_t nv = verts.size();
                        if (nv < 8 * h * w / 2) { continue; }
                        const auto in_window = [&](const partition_t cell) {
                            const size_t r = cell / W, c = cell % W;
                            return r >= row0 && r < row1 && c >= col0 && c < col1;
                        };
                        std::vector<partition_t> placed;
                        place_rectangle(verts.data(), nv, row0, col0, h, w, W,
                                        [&](const vertex_t v) { return in_window(assignment[v]); },
                                        [&](const vertex_t v) { return std::pair<f64, f64>(static_cast<f64>(assignment[v] / W), static_cast<f64>(assignment[v] % W)); },
                                        [&](const vertex_t v) { return static_cast<unsigned long long>(assignment[v]); },
                                        assignment, local_id, ac.seed + round * 7919 + r0 * W + c0, placed);
                        ++tried;
                        // the cost of the edges at the window, before and after
                        const auto cost = [&](const bool after) {
                            std::vector<std::pair<unsigned long long, f64> > pairs;
                            f64 total = 0.0;
                            for (size_t i = 0; i < nv; ++i) {
                                const vertex_t u = verts[i];
                                const partition_t cu = after ? placed[i] : assignment[u];
                                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                                    const vertex_t v = g.edges_v[j];
                                    const bool inside = in_window(assignment[v]);
                                    const partition_t cv = (after && inside) ? placed[local_id[v]] : assignment[v];
                                    if (cu == cv) { continue; }
                                    const f64 d = dist(cu, cv);
                                    total += (ac.beta_end() + ac.beta_link * d) * static_cast<f64>(g.edge_w(j)) * (inside ? 0.5 : 1.0);
                                    pairs.emplace_back(std::min<unsigned long long>(cu, cv) * 4294967296ULL + std::max<unsigned long long>(cu, cv),
                                                       ac.alpha_end() + ac.alpha_link * d);
                                }
                            }
                            std::sort(pairs.begin(), pairs.end());
                            for (size_t t = 0; t < pairs.size(); ++t) {
                                if (t == 0 || pairs[t].first != pairs[t - 1].first) { total += pairs[t].second; }
                            }
                            return total;
                        };
                        if (!(cost(true) < cost(false) - 1e-6)) { continue; }
                        // the weight limits
                        std::vector<weight_t> weight_after(k, 0);
                        std::vector<partition_t> cells;
                        for (size_t r = row0; r < row1; ++r) { for (size_t c = col0; c < col1; ++c) { cells.push_back(static_cast<partition_t>(r * W + c)); } }
                        for (size_t i = 0; i < nv; ++i) { weight_after[placed[i]] += g.v_weight(verts[i]); }
                        bool fits = true;
                        for (const partition_t x: cells) {
                            if (weight_after[x] > p_manager.get_lmax(x) || weight_after[x] < p_manager.get_lmin(x)) { fits = false; }
                        }
                        if (!fits) { continue; }
                        for (size_t i = 0; i < nv; ++i) {
                            const vertex_t u = verts[i];
                            if (placed[i] != p_manager[u]) { p_manager.move_serial(u, g.v_weight(u), p_manager[u], placed[i]); }
                        }
                        for (const partition_t x: cells) { blocks[x].clear(); }
                        for (size_t i = 0; i < nv; ++i) { blocks[placed[i]].push_back(verts[i]); assignment[verts[i]] = placed[i]; }
                        ++kept;
                        ++round_kept;
                    }
                }
                if (round_kept == 0) { break; }
            }
            if (ac.waver_busy_refinement_config.verbose) {
                std::cout << "[window] kept " << kept << " of " << tried << " windows" << std::endl;
            }
            recompute_datastructures();
        }

        // All refinements, in order. Each is switched on by its configuration. Simulated annealing swaps whole PEs.
        void refinement() {
            auto p = get_time_point();
            print_step("start");
            if (ac.relax_rounds > 0 && d_oracle.depth() == 0 && graphs.back().n >= ac.relax_min_n) {
                relax_refine();
                print_step("relax");
            }
            if (ac.window_rounds > 0 && d_oracle.depth() == 0 && graphs.back().n >= ac.window_min_n &&
                (ac.window_max_n == 0 || graphs.back().n <= ac.window_max_n)) {
                window_refine();
                print_step("window");
            }

            if (ac.waver_message_refinement_config.enabled) {
                auto sp_local = get_time_point();
                WaverMessageRefinement::refine(graphs.back(), p_manager, bv_manager, d_oracle,
                                               ac.waver_message_refinement_config,
                                               ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
                timings.msg_refine_ms += get_milli_seconds(sp_local, get_time_point());
                print_step("msg");
            }

            if (ac.waver_merge_refinement_config.enabled) {
                auto sp_local = get_time_point();
                // After the message refinement: scatter, then chains, then scatter again (chains free room again).
                // With --merge-scatter 1 only the scatter pass runs.
                auto chain_cfg = ac.waver_merge_refinement_config;
                if (ac.intense_max_n > 0 && d_oracle.depth() == 0 && graphs.back().n <= ac.intense_max_n) {
                    chain_cfg.max_depth = ac.intense_depth;
                    chain_cfg.max_rounds = ac.intense_rounds;
                    chain_cfg.partial = ac.intense_partial;
                    chain_cfg.split_parts = ac.intense_split;
                }
                auto scatter_cfg = chain_cfg;
                scatter_cfg.scatter = true;
                const bool chains = !chain_cfg.scatter;
                WaverMergeRefinement::refine(graphs.back(), p_manager, bv_manager, d_oracle, scatter_cfg,
                                             ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
                if (chains) {
                    WaverMergeRefinement::refine(graphs.back(), p_manager, bv_manager, d_oracle, chain_cfg,
                                                 ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
                    WaverMergeRefinement::refine(graphs.back(), p_manager, bv_manager, d_oracle, scatter_cfg,
                                                 ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
                }
                timings.merge_refine_ms += get_milli_seconds(sp_local, get_time_point());
                print_step("merge");
            }

            if (ac.waver_smooth_refinement_config.enabled) {
                auto sp_local = get_time_point();
                WaverSmoothRefinement::refine(graphs.back(), p_manager, d_oracle, ac.waver_smooth_refinement_config);
                recompute_datastructures(); // the smoothing moved vertices without telling the boundary structures
                timings.shift_refine_ms += get_milli_seconds(sp_local, get_time_point());
                print_step("shift");
            }

            if (ac.waver_cycle_refinement_config.enabled) {
                auto sp_local = get_time_point();
                WaverCycleRefinement::refine(graphs.back(), p_manager, bv_manager, d_oracle,
                                             ac.waver_cycle_refinement_config,
                                             ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
                timings.cycle_refine_ms += get_milli_seconds(sp_local, get_time_point());
                print_step("cycle");
            }

            if (ac.waver_busy_refinement_config.enabled) {
                auto sp_local = get_time_point();
                WaverBusyRefinement::refine(graphs.back(), p_manager, bv_manager, d_oracle,
                                            ac.waver_busy_refinement_config,
                                            ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
                timings.busy_refine_ms += get_milli_seconds(sp_local, get_time_point());
                print_step("busy");
            }

            if (ac.waver_crowd_refinement_config.enabled) {
                auto sp_local = get_time_point();
                WaverCrowdRefinement::refine(graphs.back(), p_manager, bv_manager, d_oracle,
                                             ac.waver_crowd_refinement_config,
                                             ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
                timings.crowd_refine_ms += get_milli_seconds(sp_local, get_time_point());
                print_step("crowd");
            }

            if (ac.waver_swap_refinement_config.enabled) {
                auto sp_local = get_time_point();
                WaverSwapRefinement::swap_refine(graphs.back(), p_manager, d_oracle,
                                                 ac.waver_swap_refinement_config,
                                                 ac.alpha_end(), ac.beta_end(),
                                                 ac.alpha_link, ac.beta_link);
                timings.swap_refine_ms += get_milli_seconds(sp_local, get_time_point());
                print_step("swap");
            }

            if (ac.waver_shift_refinement_config.enabled) {
                auto sp_local = get_time_point();
                WaverShiftRefinement::refine(graphs.back(), p_manager, d_oracle, ac.waver_shift_refinement_config,
                                             ac.alpha_link, ac.beta_link);
                timings.shift_refine_ms += get_milli_seconds(sp_local, get_time_point());
                print_step("shift");
            }

            if (ac.waver_long_message_refinement_config.enabled) {
                auto sp_local = get_time_point();
                WaverLongMessageRefinement::refine(graphs.back(), p_manager, d_oracle,
                                                   ac.waver_long_message_refinement_config,
                                                   ac.alpha_end(), ac.beta_end(),
                                                   ac.alpha_link, ac.beta_link);
                timings.long_msg_refine_ms += get_milli_seconds(sp_local, get_time_point());
                print_step("long_msg");
            }

            if (ac.waver_sa_refinement_config.enabled) {
                auto sp_local = get_time_point();
                WaverSARefinement::refine(graphs.back(), p_manager, d_oracle,
                                          ac.waver_sa_refinement_config, ac.seed,
                                          ac.alpha_end(), ac.beta_end(),
                                          ac.alpha_link, ac.beta_link);
                timings.sa_refine_ms += get_milli_seconds(sp_local, get_time_point());
                print_step("sa");
            }

            structures_fresh = false; // the swap refinements relabel PEs in p_manager only
            timings.refinement_ms += get_milli_seconds(p, get_time_point());
            // The boundary is stale here by design, so only the partition and its weights can be checked.
            HEAVYASSERT(assert_bweights(graphs.back(), p_manager, ac.k));
        }

        // ------------------------------------------------------------------------------------------------
        // smash and repair
        // ------------------------------------------------------------------------------------------------

        // The exact engine load of every PE (see WaverRegionLoad) and its largest value, for the partition in p_manager.
        std::vector<f64> exact_loads() {
            const graph_t &g = graphs.back();
            const partition_t k = ac.k;
            WaverRegionLoad tracker;
            tracker.initialize(g, p_manager, k, WaverRegionLoad::Weights());
            std::vector<std::vector<vertex_t> > blocks(k);
            for (vertex_t u = 0; u < g.n; ++u) { blocks[p_manager[u]].push_back(u); }
            std::vector<partition_t> seen(k, NO_ID);
            std::vector<f64> loads(k, 0.0);
            for (partition_t x = 0; x < k; ++x) {
                size_t degree = 0;
                for (const vertex_t u: blocks[x]) {
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const partition_t y = p_manager[g.edges_v[j]];
                        if (y != x && seen[y] != x) {
                            seen[y] = x;
                            ++degree;
                        }
                    }
                }
                loads[x] = tracker.load(x, degree);
            }
            return loads;
        }

        // Smash and repair, on the finished mapping: the busiest PEs are fragmented for good. Around each of them a
        // window of PEs (--smash-radius on every side) is emptied and all the vertices in it are placed again from scratch
        // onto the window, spectrally with the pull of their neighbors outside it, so that its content is spread over the
        // window in a new way. The windows of all the busiest PEs are smashed together; then the mapping is repaired like
        // after any uncoarsening step (rebalance, all refinements). The result is kept if the busiest PE has less load
        // (and the sum of the 8th powers of the loads does not rise), otherwise the old mapping is restored. Repeats for
        // --smash-rounds, each round taking the busiest PEs that no earlier round has smashed.
        void smash() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "smash");
            auto p = get_time_point();
            const graph_t &g = graphs.back();
            const partition_t k = ac.k;
            const size_t W = d_oracle.current_W();
            const size_t H = (static_cast<size_t>(k) + W - 1) / W;
            const long radius = static_cast<long>(ac.smash_radius);

            const auto score = [&](const std::vector<f64> &loads) {
                f64 peak = 0.0, power = 0.0;
                for (const f64 l: loads) {
                    peak = std::max(peak, l);
                    power += std::pow(l / 1000.0, 8.0);
                }
                return std::pair<f64, f64>(peak, power);
            };

            std::vector<char> tried(k, 0);
            std::vector<vertex_t> local_id(g.n);
            for (u64 round = 0; round < ac.smash_rounds; ++round) {
                const std::vector<partition_t> saved = snapshot();
                load_partition(saved);
                const std::vector<f64> before = exact_loads();
                const auto [peak_before, power_before] = score(before);

                // the busiest PEs that were not smashed yet, with disjoint windows
                std::vector<partition_t> order(k);
                for (partition_t x = 0; x < k; ++x) { order[x] = x; }
                std::stable_sort(order.begin(), order.end(), [&](const partition_t a, const partition_t b) { return before[a] > before[b]; });
                std::vector<char> used(k, 0);
                struct Window {
                    size_t row0, col0, h, w;
                };
                std::vector<Window> windows;
                for (const partition_t x: order) {
                    if (windows.size() >= ac.smash_count) { break; }
                    if (tried[x]) { continue; }
                    const long row = static_cast<long>(x / W), col = static_cast<long>(x % W);
                    const long r0 = std::max(0L, row - radius), r1 = std::min(static_cast<long>(H) - 1, row + radius);
                    const long c0 = std::max(0L, col - radius), c1 = std::min(static_cast<long>(W) - 1, col + radius);
                    bool free_window = true;
                    for (long r = r0; r <= r1 && free_window; ++r) {
                        for (long c = c0; c <= c1; ++c) {
                            const long cell = r * static_cast<long>(W) + c;
                            if (cell >= static_cast<long>(k) || used[cell]) {
                                free_window = false;
                                break;
                            }
                        }
                    }
                    tried[x] = 1;
                    if (!free_window) { continue; }
                    for (long r = r0; r <= r1; ++r) {
                        for (long c = c0; c <= c1; ++c) { used[r * static_cast<long>(W) + c] = 1; }
                    }
                    windows.push_back({static_cast<size_t>(r0), static_cast<size_t>(c0), static_cast<size_t>(r1 - r0 + 1), static_cast<size_t>(c1 - c0 + 1)});
                }
                if (windows.empty()) { break; }

                std::vector<std::vector<vertex_t> > members(windows.size());
                std::vector<int> window_of(k, -1);
                for (size_t i = 0; i < windows.size(); ++i) {
                    for (size_t r = 0; r < windows[i].h; ++r) {
                        for (size_t c = 0; c < windows[i].w; ++c) { window_of[(windows[i].row0 + r) * W + windows[i].col0 + c] = static_cast<int>(i); }
                    }
                }
                for (vertex_t u = 0; u < g.n; ++u) {
                    const int wi = window_of[saved[u]];
                    if (wi >= 0) { members[wi].push_back(u); }
                }

                std::vector<partition_t> current = saved;
                const std::vector<partition_t> outside = saved; // the cells of the neighbors outside a window
                for (size_t i = 0; i < windows.size(); ++i) {
                    if (members[i].size() < 8) { continue; }
                    std::vector<partition_t> placed;
                    const int wi = static_cast<int>(i);
                    place_rectangle(members[i].data(), members[i].size(), windows[i].row0, windows[i].col0, windows[i].h,
                                    windows[i].w, W,
                                    [&](const vertex_t v) { return window_of[saved[v]] == wi; },
                                    [&](const vertex_t v) { return std::pair<f64, f64>(static_cast<f64>(saved[v] / W), static_cast<f64>(saved[v] % W)); },
                                    [&](const vertex_t v) { return static_cast<unsigned long long>(saved[v]); },
                                    outside, local_id, ac.seed + 7919 * (round + 1) + i, placed);
                    for (size_t j = 0; j < members[i].size(); ++j) { current[members[i][j]] = placed[j]; }
                }

                load_partition(current);
                // a window that does not respect the weight limits is not smashed after all
                for (size_t i = 0; i < windows.size(); ++i) {
                    bool ok = true;
                    for (size_t r = 0; r < windows[i].h && ok; ++r) {
                        for (size_t c = 0; c < windows[i].w; ++c) {
                            const partition_t cell = static_cast<partition_t>((windows[i].row0 + r) * W + windows[i].col0 + c);
                            if (p_manager.bweights[cell] > p_manager.get_lmax(cell) || p_manager.bweights[cell] < p_manager.get_lmin(cell)) {
                                ok = false;
                                break;
                            }
                        }
                    }
                    if (!ok) {
                        for (const vertex_t u: members[i]) { current[u] = saved[u]; }
                    }
                }
                load_partition(current);
                const std::vector<f64> smashed = exact_loads();
                recompute_datastructures();
                settle_level(); // the repair
                const std::vector<f64> repaired = exact_loads();
                const auto [peak_after, power_after] = score(repaired);
                const bool better = peak_after < peak_before || (peak_after == peak_before && power_after < power_before);
                if (ac.waver_busy_refinement_config.verbose || ac.waver_merge_refinement_config.verbose) {
                    std::cout << "[smash] round " << round << ": " << windows.size() << " windows, busiest PE " << peak_before
                              << " -> (smashed) " << score(smashed).first << " -> (repaired) " << peak_after
                              << (better ? ", kept" : ", restored") << std::endl;
                }
                if (!better) { load_partition(saved); }
            }
            timings.smash_ms += get_milli_seconds(p, get_time_point());
        }

        // ------------------------------------------------------------------------------------------------
        // V-cycle
        // ------------------------------------------------------------------------------------------------

        // V-cycle on the finished mapping: the graph is coarsened again, but only vertices placed on the same PE are
        // merged, so the mapping carries over to every level unchanged. There is no initial partitioning; the
        // coarsest level starts from the current mapping, which is refined there, and every finer level is legalized
        // and refined like in the main uncoarsening. The topology stays at its finest level throughout. The old mapping
        // is restored if the cycle does not lower the message cost.
        void vcycle() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "vcycle");
            auto p = get_time_point();
            const std::vector<partition_t> saved = snapshot();
            const f64 cost_before = current_message_cost();

            // the coarsening respects the partition: only vertices of one PE are matched
            recompute_datastructures();
            for (u64 level = 0; graphs.back().n > d_oracle.current_k() * ac.vcycle_c; ++level) {
                coarsening(level);
                contraction();
                if (graphs.back().n == graphs[graphs.size() - 2].n) {
                    graphs.pop_back();
                    mappings.pop_back();
                    break;
                }
            }

            finish_initial_partition();
            while (!mappings.empty()) {
                uncoarsening();
                settle_level();
            }

            const f64 cost_after = current_message_cost();
            if (ac.waver_busy_refinement_config.verbose || ac.waver_merge_refinement_config.verbose) {
                std::cout << "[vcycle] message cost " << cost_before << " -> " << cost_after << std::endl;
            }
            if (cost_after >= cost_before) { load_partition(saved); }
            timings.vcycle_ms += get_milli_seconds(p, get_time_point());
        }

        void load_partition(const std::vector<partition_t> &part) {
            structures_fresh = false;
            for (vertex_t u = 0; u < graphs.back().n; ++u) { p_manager.partition[u] = part[u]; }
            p_manager.recalculate_weights(graphs.back());
        }

        f64 current_message_cost() {
            return SpectralTracker::message_cost(graphs.back(), p_manager.partition, d_oracle, ac.alpha_end(),
                                                 ac.beta_end(), ac.alpha_link, ac.beta_link);
        }

        // ------------------------------------------------------------------------------------------------
        // partition helpers
        // ------------------------------------------------------------------------------------------------

        std::vector<partition_t> snapshot() const {
            std::vector<partition_t> part(graphs.back().n);
            for (vertex_t u = 0; u < part.size(); ++u) { part[u] = p_manager[u]; }
            return part;
        }

        // The boundary matches the partition.
        // p_manager.size() is not kept up to date by the swap refinements, so the vertices per PE are not checked.
        bool state_is_consistent() const {
            return assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, ac.k, 1, false);
        }
    };
}

#endif //HEIPROMAP_WAVERMULTILEVEL_SOLVER_H
