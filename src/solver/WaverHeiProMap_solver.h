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
#include <iomanip>
#include <sstream>
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
#include "../utility/hungarian.h"
#include "../utility/qap.h"
#include "../utility/spectral_embedding.h"
#include "../utility/spectral_tracker.h"
#include "../utility/utils.h"

#include "../coarsening/global_path_algorithm.h"
#include "../coarsening/size_constrained_lp.h"
#include "../coarsening/heavy_edge_matching.h"
#include "../rebalance/rebalancer.h"
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
        WaverHeiProMapConfiguration ac;
        RandomEngine random_engine;

        std::chrono::high_resolution_clock::time_point sp;

        std::vector<graph_t> graphs;

        p_manager_t p_manager;
        bv_manager_t bv_manager;
        LargeQuotientGraph q_graph;
        block_conn_t block_conn;
        GridDistanceOracle d_oracle;
        SpectralTracker spectral;

        // matching
        std::vector<Mapping> mappings;
        GlobalPathAlgorithmMatcher gpa_matcher;
        SizeConstrainedLP size_constrained_lp;
        HeavyEdgeMatching heavy_edge_matching;

        Rebalancer rebalancer;

        // refinement
        WaverLabelPropagationRefinement<false> lp_refine;

        f64 io_ms = 0.0;
        f64 msg_refine_ms = 0.0;
        f64 misc_ms = 0.0;
        f64 coarsening_ms = 0.0;
        f64 contraction_ms = 0.0;
        f64 initial_partitioning_ms = 0.0;
        f64 uncontraction_ms = 0.0;
        f64 rebalance_ms = 0.0;
        f64 refinement_ms = 0.0;
        f64 lp_refine_ms = 0.0;

        void activate_all_blocks() {
            for (partition_t id = 0; id < ac.k; ++id) {
                p_manager.set_hierarchy_level(id, 0);
            }
        }

        void update_lmax_array() {
            const CSRGraph &topo = d_oracle.current_topology();
            std::vector<weight_t> lmax_vec(ac.k);
            if (topo.uniform_v_weights) {
                weight_t lmax = std::ceil((1.0 + ac.imbalance) * ((f64) graphs[0].g_weight / (f64) ac.k));
                for (partition_t i = 0; i < ac.k; ++i) {
                    lmax_vec[i] = lmax;
                }
            } else {
                f64 topo_total_w = topo.g_weight > 0 ? (f64) topo.g_weight : (f64) ac.k;
                for (partition_t i = 0; i < ac.k; ++i) {
                    double fair_share = (double) topo.v_weight(i) / topo_total_w;
                    lmax_vec[i] = std::ceil((1.0 + ac.imbalance) * graphs[0].g_weight * fair_share);
                }
            }
            weight_t global_lmin = static_cast<weight_t>(ac.min_imbalance * ((f64) graphs[0].g_weight / (f64) ac.k));
            if (graphs[0].n >= ac.k && global_lmin < 1) global_lmin = 1;
            weight_t lmax_min = lmax_vec[0], lmax_max = lmax_vec[0];
            weight_t lmin_min = std::min(global_lmin, lmax_vec[0] / 2);
            weight_t lmin_max = lmin_min;
            for (partition_t i = 0; i < ac.k; ++i) {
                const weight_t lmin_i = std::min(global_lmin, lmax_vec[i] / 2);
                p_manager.set_lmax(i, lmax_vec[i]);
                p_manager.set_lmin(i, lmin_i);
                lmax_min = std::min(lmax_min, lmax_vec[i]);
                lmax_max = std::max(lmax_max, lmax_vec[i]);
                lmin_min = std::min(lmin_min, lmin_i);
                lmin_max = std::max(lmin_max, lmin_i);
            }
        }

        void dump_heatmap(const std::string &tag) {
            const size_t W = static_cast<size_t>(d_oracle.current_W());
            const size_t H = (ac.k + W - 1) / W;
            const size_t scale = std::max(size_t(1), 512 / std::max(W, H));
            const size_t pw = W * scale, ph = H * scale;

            weight_t max_bw = 1;
            for (partition_t b = 0; b < ac.k; ++b) {
                max_bw = std::max(max_bw, p_manager.bweights[b]);
            }

            auto heat_color = [](double t, uint8_t &r, uint8_t &g, uint8_t &b) {
                if (t <= 0.0) {
                    r = 0;
                    g = 0;
                    b = 0;
                    return;
                }
                t = std::min(t, 1.0);
                if (t < 0.25) {
                    r = 0;
                    g = uint8_t(255 * t * 4);
                    b = 255;
                } else if (t < 0.50) {
                    r = 0;
                    g = 255;
                    b = uint8_t(255 * (2 - t * 4));
                } else if (t < 0.75) {
                    r = uint8_t(255 * (t * 4 - 2));
                    g = 255;
                    b = 0;
                } else {
                    r = 255;
                    g = uint8_t(255 * (4 - t * 4));
                    b = 0;
                }
            };

            std::string fname = "/tmp/heatmap_" + tag + ".ppm";
            std::ofstream out(fname, std::ios::binary);
            out << "P6\n" << pw << " " << ph << "\n255\n";
            for (size_t py = 0; py < ph; ++py) {
                for (size_t px = 0; px < pw; ++px) {
                    partition_t b = static_cast<partition_t>((py / scale) * W + (px / scale));
                    double t = (b < ac.k) ? static_cast<double>(p_manager.bweights[b]) / max_bw : 0.0;
                    uint8_t r, g, bl;
                    heat_color(t, r, g, bl);
                    out.put(r);
                    out.put(g);
                    out.put(bl);
                }
            }
        }

        void diffuse_into_empties() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "diffuse_into_empties");
            const size_t W = static_cast<size_t>(d_oracle.current_W());
            const partition_t k = ac.k;
            const auto &g = graphs.back();

            auto get_neighbors = [&](partition_t b, partition_t *nbrs) -> size_t {
                size_t cnt = 0;
                const size_t c = b % W;
                if (b >= static_cast<partition_t>(W)) nbrs[cnt++] = b - static_cast<partition_t>(W);
                if (b + static_cast<partition_t>(W) < k) nbrs[cnt++] = b + static_cast<partition_t>(W);
                if (c > 0) nbrs[cnt++] = b - 1;
                if (c + 1 < W && b + 1 < k) nbrs[cnt++] = b + 1;
                return cnt;
            };

            // Target weight: fair share per PE
            const weight_t avg_bw = g.g_weight / k;

            // Build per-block vertex lists
            std::vector<std::vector<vertex_t> > blk_v(k);
            for (vertex_t u = 0; u < g.n; ++u) {
                const partition_t b = p_manager.partition[u];
                if (b != NO_ID && b < k) { blk_v[b].push_back(u); }
            }

            // Phase 1: aggressively fill empties.
            // Pull from ANY non-empty neighbor, allow donor to go as low
            // as 1 (just don't empty it).  Repeat until no empties touch
            // a non-empty neighbor.
            size_t round = 0;
            bool progress = true;
            while (progress) {
                progress = false;
                size_t filled_this_round = 0;
                size_t still_empty_count = 0;
                size_t skipped_cap = 0, skipped_floor = 0, skipped_no_donor = 0;
                for (partition_t e = 0; e < k; ++e) {
                    if (p_manager.bweights[e] > 0) { continue; }
                    still_empty_count++;

                    partition_t nbrs[4];
                    const size_t n_nbrs = get_neighbors(e, nbrs);

                    partition_t donor = NO_ID;
                    weight_t donor_bw = 0;
                    for (size_t i = 0; i < n_nbrs; ++i) {
                        if (p_manager.bweights[nbrs[i]] > donor_bw) {
                            donor_bw = p_manager.bweights[nbrs[i]];
                            donor = nbrs[i];
                        }
                    }
                    const weight_t donor_lmin = p_manager.get_lmin(donor == NO_ID ? 0 : donor);
                    if (donor == NO_ID || donor_bw <= donor_lmin) {
                        skipped_no_donor++;
                        continue;
                    }

                    const weight_t e_cap = p_manager.get_lmax(e);
                    const weight_t donor_floor = std::max(donor_lmin, donor_bw / 2);
                    bool moved_any = false;
                    for (const vertex_t u: blk_v[donor]) {
                        if (p_manager.partition[u] != donor) { continue; }
                        const weight_t uw = g.v_weight(u);
                        if (p_manager.bweights[e] + uw > e_cap) {
                            skipped_cap++;
                            continue;
                        }
                        if (p_manager.bweights[donor] - uw < donor_floor) {
                            skipped_floor++;
                            continue;
                        }

                        p_manager.bweights[donor] -= uw;
                        p_manager.partition[u] = e;
                        p_manager.bweights[e] += uw;
                        blk_v[e].push_back(u);
                        progress = true;
                        moved_any = true;
                    }
                    if (moved_any) { filled_this_round++; }
                }
                round++;
            }

            // Phase 2: equalize — pull from above-average neighbors into
            // below-average blocks.
            progress = true;
            while (progress) {
                progress = false;
                for (partition_t e = 0; e < k; ++e) {
                    if (p_manager.bweights[e] >= avg_bw) { continue; }

                    const weight_t e_cap = p_manager.get_lmax(e);

                    partition_t nbrs[4];
                    const size_t n_nbrs = get_neighbors(e, nbrs);

                    partition_t donor = NO_ID;
                    weight_t donor_bw = avg_bw;
                    for (size_t i = 0; i < n_nbrs; ++i) {
                        if (p_manager.bweights[nbrs[i]] > donor_bw) {
                            donor_bw = p_manager.bweights[nbrs[i]];
                            donor = nbrs[i];
                        }
                    }
                    if (donor == NO_ID) { continue; }

                    for (const vertex_t u: blk_v[donor]) {
                        if (p_manager.partition[u] != donor) { continue; }
                        const weight_t uw = g.v_weight(u);
                        if (p_manager.bweights[e] + uw > e_cap) { continue; }
                        if (p_manager.bweights[donor] - uw < avg_bw) { continue; }

                        p_manager.bweights[donor] -= uw;
                        p_manager.partition[u] = e;
                        p_manager.bweights[e] += uw;
                        blk_v[e].push_back(u);
                        progress = true;

                        if (p_manager.bweights[e] >= avg_bw) { break; }
                    }
                }
            }
        }

        void reinit_structures_for_new_k() {
            p_manager.initialize(graphs[0].n, ac.k, graphs[0].g_weight);
            activate_all_blocks();
            update_lmax_array();
            bv_manager.initialize(graphs[0].n, ac.k);

            rebalancer.initialize(graphs[0].n, graphs[0].m, ac.k, random_engine.get_u64());
        }

        void initialize_common(CSRGraph topology) {
            random_engine = RandomEngine(ac.seed);

            d_oracle.initialize(std::move(topology), ac.threads);
            ac.k = d_oracle.current_k();

            p_manager.initialize(graphs[0].n, ac.k, graphs[0].g_weight);
            activate_all_blocks();
            update_lmax_array();
            bv_manager.initialize(graphs[0].n, ac.k);
            HEAVYASSERT(assert_state_pre_partitioning(graphs[0], p_manager, ac.k));

            gpa_matcher.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine, ac.global_path_algorithm_config);
            size_constrained_lp.initialize(graphs[0].n, graphs[0].m, ac.k, random_engine.get_u64(), ac.size_constrained_lp_config);

            rebalancer.initialize(graphs[0].n, graphs[0].m, ac.k, random_engine.get_u64());

            if (ac.waver_label_propagation_config.enabled) {
                lp_refine.initialize(graphs[0].n, graphs[0].m, ac.k, ac.threads, random_engine.get_u64(), ac.waver_label_propagation_config);
            }
        }

    public:
        explicit WaverHeiProMapSolver(const WaverHeiProMapConfiguration &t_ac) {
            sp = std::chrono::high_resolution_clock::now();

            ac = t_ac;

            auto sp_io = get_time_point();
            graphs.emplace_back(ac.graph_in);
            io_ms += get_milli_seconds(sp_io, get_time_point());

            if (ac.topology_in.empty()) {
                std::cerr << "Topology required for WaverMap!" << std::endl;
                exit(1);
            }

            initialize_common(CSRGraph(ac.topology_in));
        }

        explicit WaverHeiProMapSolver(graph_t &&g, graph_t &&topology, const WaverHeiProMapConfiguration &t_ac) {
            sp = std::chrono::high_resolution_clock::now();

            ac = t_ac;
            graphs.emplace_back(std::move(g));

            initialize_common(std::move(topology));
        }

        std::vector<vertex_t> solve() {
            internal_solve();

            std::vector<partition_t> p(graphs[0].n);
            for (vertex_t u = 0; u < graphs[0].n; ++u) { p[u] = p_manager[u]; }

            const auto ep = std::chrono::high_resolution_clock::now();
            f64 duration = get_seconds(sp, ep);

            q_graph.compute_from_scratch(graphs[0], p_manager);
            block_conn.compute_from_scratch(graphs[0], p_manager);
            HEAVYASSERT(assert_state_after_partitioning(graphs[0], p_manager, bv_manager, q_graph, block_conn, ac.k));
            MessageCostBreakdown cost = get_message_cost(q_graph, ac.k, d_oracle,
                                                         ac.alpha_end(), ac.beta_end(),
                                                         ac.alpha_link, ac.beta_link);
            f64 max_link_load = get_max_link_load(q_graph, ac.k, d_oracle.current_W(),
                                                  ac.alpha_link, ac.beta_link);
            weight_t lmax = std::ceil((1.0 + ac.imbalance) * ((f64) graphs[0].g_weight / (f64) ac.k));

            size_t n_empty_partitions = 0;
            size_t n_overloaded_partitions = 0;
            weight_t sum_too_much = 0;
            for (partition_t id = 0; id < ac.k; ++id) {
                n_empty_partitions += p_manager.get_bweight(id) == 0;
                n_overloaded_partitions += p_manager.get_bweight(id) > lmax;
                sum_too_much += std::max((weight_t) 0, p_manager.get_bweight(id) - lmax);
            }

            std::cout << "Graph                 : " << ac.graph_in << std::endl;
            std::cout << "Total time (s)        : " << duration << std::endl;
            std::cout << "#Nodes                : " << graphs[0].n << std::endl;
            std::cout << "#Edges                : " << graphs[0].m << std::endl;
            std::cout << "k                     : " << ac.k << std::endl;
            std::cout << "Lmax                  : " << lmax << std::endl;
            std::cout << "Min imbalance         : " << ac.min_imbalance << std::endl;
            std::cout << "Threads               : " << ac.threads << std::endl;
            std::cout << "alpha_send            : " << ac.alpha_send << std::endl;
            std::cout << "alpha_recv            : " << ac.alpha_recv << std::endl;
            std::cout << "alpha_link            : " << ac.alpha_link << std::endl;
            std::cout << "beta_send             : " << ac.beta_send << std::endl;
            std::cout << "beta_recv             : " << ac.beta_recv << std::endl;
            std::cout << "beta_link             : " << ac.beta_link << std::endl;
            std::cout << "-----------------------" << std::endl;
            cost.print(std::cout);
            std::cout << std::left << std::setw(36) << "Max link load"
                    << ": " << std::right << std::setw(16) << static_cast<u64>(max_link_load) << "\n";
            std::cout << "max block w           : " << max(p_manager.get_bweights()) << std::endl;
            std::cout << "#empty partitions     : " << n_empty_partitions << std::endl;
            std::cout << "#oload partitions     : " << n_overloaded_partitions << std::endl;
            std::cout << "Sum oload weights     : " << sum_too_much << std::endl;
            std::cout << "------- Time -------" << std::endl;
            std::cout << "IO                    : " << io_ms << std::endl;
            std::cout << "Coarsening            : " << coarsening_ms << std::endl;
            std::cout << "Contraction           : " << contraction_ms << std::endl;
            std::cout << "Init. Part.           : " << initial_partitioning_ms << std::endl;
            std::cout << "Uncontraction         : " << uncontraction_ms << std::endl;
            std::cout << "Rebalance             : " << rebalance_ms << std::endl;
            std::cout << "Refinement            : " << refinement_ms << std::endl;
            std::cout << "  Label Prop.         : " << lp_refine_ms << std::endl;
            std::cout << "  Msg Refinement      : " << msg_refine_ms << std::endl;
            std::cout << "Misc                  : " << misc_ms << std::endl;
            std::cout << "ALL                   : " << io_ms + misc_ms + coarsening_ms + contraction_ms + initial_partitioning_ms + uncontraction_ms + rebalance_ms + refinement_ms << std::endl;

            return p;
        }


        // Compact balance summary: "over=2 under=0 empty=0 [3717..4095]/4096".
        std::string balance_summary() {
            partition_t n_empty = 0, n_over = 0, n_under = 0;
            weight_t max_bw = 0, min_bw = std::numeric_limits<weight_t>::max(), lmax_max = 0;
            for (partition_t id = 0; id < ac.k; ++id) {
                if (!p_manager.is_active(id)) continue;
                const weight_t bw = p_manager.get_bweight(id);
                const weight_t lmin = p_manager.get_lmin(id);
                lmax_max = std::max(lmax_max, p_manager.get_lmax(id));
                if (bw == 0) { ++n_empty; } else if (bw > p_manager.get_lmax(id)) { ++n_over; } else if (lmin > 0 && bw < lmin) { ++n_under; }
                max_bw = std::max(max_bw, bw);
                if (bw > 0) min_bw = std::min(min_bw, bw);
            }
            std::ostringstream os;
            os << "over=" << n_over << " under=" << n_under << " empty=" << n_empty
               << " bw=[" << min_bw << ".." << max_bw << "]/" << lmax_max;
            return os.str();
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

        // One line per stage. If the tracked spectral placement beats the multilevel partition, adopt it
        // (rebalanced and refined). The spectral tracker itself is never overwritten.
        void print_level_stats(const std::string &tag, u64 level) {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "print_level_stats");
            const graph_t &g = graphs.back();
            std::ostringstream os;
            os << std::left << std::setw(9) << tag << std::right
               << " n=" << std::setw(6) << g.n << " k=" << std::setw(6) << ac.k;

            if (ac.spectral_solution && spectral.size() == g.n) {
                std::vector<partition_t> sp_part;
                spectral.place(g, d_oracle.current_topology(), d_oracle.current_W(), sp_part);
                const f64 ml = SpectralTracker::message_cost(g, p_manager.partition, d_oracle,
                                                             ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
                const f64 sp = SpectralTracker::message_cost(g, sp_part, d_oracle,
                                                             ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
                os << " | ml=" << std::fixed << std::setprecision(2) << ml / 1e6 << "M sp=" << sp / 1e6 << "M";
                if (sp < ml) {
                    std::vector<partition_t> ml_part(g.n);
                    for (vertex_t u = 0; u < g.n; ++u) ml_part[u] = p_manager[u];
                    load_partition(sp_part);
                    recompute_datastructures();
                    rebalance(level);
                    diffuse_into_empties();
                    recompute_datastructures();
                    refinement();
                    const f64 refined = current_message_cost();
                    if (refined < ml) {
                        os << " -> sp refined=" << refined / 1e6 << "M";
                    } else {
                        load_partition(ml_part);
                        recompute_datastructures();
                        os << " -> ml (sp refined=" << refined / 1e6 << "M)";
                    }
                } else {
                    os << " -> ml";
                }
            }

            u64 n_pairs = 0;
            for (partition_t id = 0; id < ac.k; ++id) { n_pairs += q_graph.neighbors(id).size(); }
            os << " | |Eq|=" << n_pairs << " " << balance_summary();
            std::cout << os.str() << std::endl;
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
            update_lmax_array();
            return true;
        }

        void optimize_block_assignment() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "optimize_block_assignment");
            const auto &g = graphs.back();
            const partition_t k = ac.k;

            // Build flow matrix: flow[b1][b2] = sum of (alpha_link + beta_link * w) over edges between b1 and b2
            std::vector<f64> flow(k * k, 0.0);
            for (vertex_t u = 0; u < g.n; ++u) {
                const partition_t bu = p_manager.partition[u];
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const partition_t bv = p_manager.partition[g.edges_v[j]];
                    if (bu != bv) {
                        const f64 w = static_cast<f64>(g.edge_w(j));
                        flow[bu * k + bv] += ac.alpha_link + ac.beta_link * w;
                    }
                }
            }

            // Build distance matrix for all cell pairs
            std::vector<f64> dist(k * k, 0.0);
            for (partition_t i = 0; i < k; ++i) {
                for (partition_t j = i + 1; j < k; ++j) {
                    f64 d = static_cast<f64>(d_oracle.get(i, j));
                    dist[i * k + j] = d;
                    dist[j * k + i] = d;
                }
            }

            // Evaluate QAP cost for a given permutation (block b -> cell perm[b])
            auto qap_cost = [&](const std::vector<partition_t> &perm) -> f64 {
                f64 cost = 0.0;
                for (partition_t b1 = 0; b1 < k; ++b1) {
                    for (partition_t b2 = b1 + 1; b2 < k; ++b2) {
                        f64 f = flow[b1 * k + b2] + flow[b2 * k + b1];
                        if (f > 0.0) {
                            cost += f * dist[perm[b1] * k + perm[b2]];
                        }
                    }
                }
                return cost;
            };

            // Evaluate swap delta: cost change if we swap perm[a] and perm[b]
            auto swap_delta = [&](const std::vector<partition_t> &perm, partition_t a, partition_t b) -> f64 {
                f64 delta = 0.0;
                const partition_t ca = perm[a], cb = perm[b];
                for (partition_t x = 0; x < k; ++x) {
                    if (x == a || x == b) continue;
                    const partition_t cx = perm[x];
                    f64 fa = flow[a * k + x] + flow[x * k + a];
                    f64 fb = flow[b * k + x] + flow[x * k + b];
                    if (fa > 0.0) delta += fa * (dist[cb * k + cx] - dist[ca * k + cx]);
                    if (fb > 0.0) delta += fb * (dist[ca * k + cx] - dist[cb * k + cx]);
                }
                f64 fab = flow[a * k + b] + flow[b * k + a];
                if (fab > 0.0) delta += fab * (dist[cb * k + ca] - dist[ca * k + cb]);
                return delta;
            };

            // Local search: repeatedly find best improving swap
            auto local_search = [&](std::vector<partition_t> &perm) {
                bool improved = true;
                while (improved) {
                    improved = false;
                    partition_t best_a = 0, best_b = 0;
                    f64 best_delta = -1e-9;
                    for (partition_t a = 0; a < k; ++a) {
                        for (partition_t b = a + 1; b < k; ++b) {
                            f64 d = swap_delta(perm, a, b);
                            if (d < best_delta) {
                                best_delta = d;
                                best_a = a;
                                best_b = b;
                            }
                        }
                    }
                    if (best_delta < -1e-9) {
                        std::swap(perm[best_a], perm[best_b]);
                        improved = true;
                    }
                }
            };

            // Current assignment: block b is at cell b (identity permutation after greedy)
            std::vector<partition_t> best_perm(k);
            std::iota(best_perm.begin(), best_perm.end(), 0);
            local_search(best_perm);
            f64 best_cost = qap_cost(best_perm);

            // Try spectral assignment as starting permutation
            {
                // Build a dense CSRGraph from the flow matrix for spectral embedding
                CSRGraph flow_graph;
                std::vector<std::vector<std::pair<vertex_t, weight_t>>> adj(k);
                for (partition_t i = 0; i < k; ++i) {
                    for (partition_t j = 0; j < k; ++j) {
                        if (i == j) continue;
                        f64 f = flow[i * k + j];
                        if (f > 0.0) {
                            adj[i].emplace_back(j, static_cast<weight_t>(std::max(1.0, f)));
                        }
                    }
                }
                size_t total_edges = 0;
                for (partition_t i = 0; i < k; ++i) total_edges += adj[i].size();
                flow_graph.n = k;
                flow_graph.m = total_edges;
                flow_graph.neighborhoods.initialize(k + 1);
                flow_graph.edges_v.initialize(total_edges);
                flow_graph.edges_w.initialize(total_edges);
                flow_graph.v_weights.initialize(k, 1);
                flow_graph.g_weight = k;
                flow_graph.neighborhoods[0] = 0;
                for (partition_t i = 0; i < k; ++i) {
                    flow_graph.neighborhoods[i + 1] = flow_graph.neighborhoods[i] + adj[i].size();
                    for (size_t e = 0; e < adj[i].size(); ++e) {
                        size_t idx = flow_graph.neighborhoods[i] + e;
                        flow_graph.edges_v[idx] = adj[i][e].first;
                        flow_graph.edges_w[idx] = adj[i][e].second;
                    }
                }

                auto spectral_embedding = SpectralEmbedding::compute_laplacian_eigenvectors(flow_graph, 2, ac.seed);

                // Normalize embedding to [0,1]
                const size_t W = static_cast<size_t>(d_oracle.current_W());
                const size_t H = (k + W - 1) / W;
                for (size_t d = 0; d < 2; ++d) {
                    f64 mn = spectral_embedding[d], mx = spectral_embedding[d];
                    for (partition_t i = 0; i < k; ++i) {
                        mn = std::min(mn, spectral_embedding[i * 2 + d]);
                        mx = std::max(mx, spectral_embedding[i * 2 + d]);
                    }
                    f64 range = mx - mn;
                    if (range < 1e-12) range = 1.0;
                    for (partition_t i = 0; i < k; ++i) {
                        spectral_embedding[i * 2 + d] = (spectral_embedding[i * 2 + d] - mn) / range;
                    }
                }

                // Grid cell coordinates normalized to [0,1]
                std::vector<f64> grid_x(k), grid_y(k);
                for (partition_t c = 0; c < k; ++c) {
                    grid_x[c] = static_cast<f64>(c % W) / static_cast<f64>(std::max(W - 1, static_cast<size_t>(1)));
                    grid_y[c] = static_cast<f64>(c / W) / static_cast<f64>(std::max(H - 1, static_cast<size_t>(1)));
                }

                // Try 8 orientation variants
                for (int variant = 0; variant < 8; ++variant) {
                    bool flip_x = (variant & 1) != 0;
                    bool flip_y = (variant & 2) != 0;
                    bool swap_axes = (variant & 4) != 0;

                    std::vector<std::vector<weight_t>> cost_mat(k, std::vector<weight_t>(k));
                    for (partition_t b = 0; b < k; ++b) {
                        f64 bx = spectral_embedding[b * 2 + 0];
                        f64 by = spectral_embedding[b * 2 + 1];
                        if (swap_axes) std::swap(bx, by);
                        if (flip_x) bx = 1.0 - bx;
                        if (flip_y) by = 1.0 - by;
                        for (partition_t c = 0; c < k; ++c) {
                            f64 dx = bx - grid_x[c];
                            f64 dy = by - grid_y[c];
                            cost_mat[b][c] = static_cast<weight_t>((dx * dx + dy * dy) * 1e6);
                        }
                    }

                    auto assignment = solve_hungarian(cost_mat);
                    std::vector<partition_t> perm(k);
                    for (partition_t b = 0; b < k; ++b) perm[b] = static_cast<partition_t>(assignment[b]);
                    local_search(perm);
                    f64 cost = qap_cost(perm);
                    if (cost < best_cost) {
                        best_cost = cost;
                        best_perm = perm;
                    }
                }
            }

            // Random restarts
            const size_t n_restarts = std::max(static_cast<size_t>(1000), static_cast<size_t>(k * k));
            std::mt19937_64 rng(ac.seed + 12345);
            for (size_t r = 0; r < n_restarts; ++r) {
                std::vector<partition_t> perm(k);
                std::iota(perm.begin(), perm.end(), 0);
                std::shuffle(perm.begin(), perm.end(), rng);
                local_search(perm);
                f64 cost = qap_cost(perm);
                if (cost < best_cost) {
                    best_cost = cost;
                    best_perm = perm;
                }
            }

            // Check if permutation is non-trivial
            bool changed = false;
            for (partition_t b = 0; b < k; ++b) {
                if (best_perm[b] != b) { changed = true; break; }
            }
            if (!changed) return;

            // Apply permutation: remap all vertex assignments
            // Build inverse: inv_perm[cell] = block that should go there
            // We need: vertex u was in block b, now block b goes to cell perm[b]
            // But cells ARE blocks in our representation, so we need to reassign
            // vertex u from block b to block perm[b]
            std::vector<weight_t> new_bweights(k, 0);
            for (vertex_t u = 0; u < g.n; ++u) {
                partition_t old_block = p_manager.partition[u];
                partition_t new_block = best_perm[old_block];
                p_manager.partition[u] = new_block;
                new_bweights[new_block] += g.v_weight(u);
            }
            for (partition_t b = 0; b < k; ++b) {
                p_manager.bweights[b] = new_bweights[b];
            }

        }

    private:
        void internal_solve() {
            u64 level = 0;

            // Pre-coarsen the topology oracle as far as possible.
            while (d_oracle.can_coarsen()) {
                d_oracle.coarsen();
                ac.k = d_oracle.current_k();
                reinit_structures_for_new_k();
            }

            update_lmax_array();

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
            optimize_block_assignment();
            recompute_datastructures();
            rebalance(level);
            refinement();
            print_level_stats("init", level);

            // Graph uncoarsening loop.
            while (!mappings.empty()) {
                level -= 1;

                // Try to uncoarsen topology levels before graph uncoarsening.
                while (d_oracle.depth() > 0) {
                    if (!try_uncoarsen_topology()) break;

                    rebalance(level);
                    diffuse_into_empties();
                    recompute_datastructures();
                    refinement();
                }

                if (ac.spectral_solution) spectral.prolong(mappings.back(), graphs[graphs.size() - 2]);
                uncoarsening();

                rebalance(level);
                diffuse_into_empties();
                recompute_datastructures();
                refinement();
                print_level_stats("L" + std::to_string(level), level);
            }

            // Uncoarsen remaining topology levels after graph is fully uncoarsened.
            while (d_oracle.depth() > 0) {
                if (!try_uncoarsen_topology()) {
                    rebalance(0);
                    diffuse_into_empties();
                    recompute_datastructures();
                    refinement();
                    print_level_stats("topo_blocked", 0);
                    break;
                }

                rebalance(0);
                diffuse_into_empties();
                recompute_datastructures();
                refinement();
                print_level_stats("topo", 0);
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

            initial_partitioning_ms += get_milli_seconds(p, get_time_point());

            q_graph.compute_from_scratch(graphs.back(), p_manager);
            block_conn.compute_from_scratch(graphs.back(), p_manager);
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

            coarsening_ms += get_milli_seconds(p, get_time_point());
        }

        void contraction() {
            auto p = get_time_point();
            graphs.emplace_back();
            graphs.back().contract(graphs[graphs.size() - 2], mappings.back(), ac.threads);
            p_manager.contract(mappings.back());

            contraction_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_pre_partitioning(graphs.back(), p_manager, ac.k));
        }

        void uncoarsen_topology() {
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "uncoarsen_topology");
            auto p = get_time_point();

            // Build coarse->fine map for partition expansion before popping oracle level
            partition_t coarse_k = d_oracle.current_k();
            partition_t fine_k = d_oracle.fine_k();
            const Mapping &topo_map = d_oracle.top_mapping();

            std::vector<std::vector<partition_t> > coarse_to_fine(coarse_k);
            for (partition_t f = 0; f < fine_k; ++f) {
                coarse_to_fine[topo_map.get(f)].push_back(f);
            }

            d_oracle.uncoarsen();
            ac.k = d_oracle.current_k();

            // Expand partition: assign each graph vertex from coarse block to least-loaded fine block
            std::vector<partition_t> old_partition(graphs.back().n);
            for (vertex_t u = 0; u < graphs.back().n; ++u) {
                old_partition[u] = p_manager[u];
            }

            reinit_structures_for_new_k();
            p_manager.reset_weights();

            // mark all unassigned so cost computation can skip them
            for (vertex_t u = 0; u < graphs.back().n; ++u) {
                p_manager.partition[u] = NO_ID;
            }

            // blk_verts[b] = vertices assigned to fine block b; used by second pass
            std::vector<std::vector<vertex_t> > blk_verts(fine_k);

            for (vertex_t u = 0; u < graphs.back().n; ++u) {
                const auto &fine_blocks = coarse_to_fine[old_partition[u]];
                const weight_t u_weight = graphs.back().v_weight(u);

                // collect distinct assigned neighbor fine blocks (for alpha_end)
                partition_t neighbor_blocks[256];
                size_t n_neighbor_blocks = 0;
                for (size_t j = graphs.back().neighborhoods[u]; j < graphs.back().neighborhoods[u + 1]; ++j) {
                    const partition_t v_block = p_manager.partition[graphs.back().edges_v[j]];
                    if (v_block == NO_ID) { continue; }
                    bool found = false;
                    for (size_t x = 0; x < n_neighbor_blocks; ++x) {
                        if (neighbor_blocks[x] == v_block) {
                            found = true;
                            break;
                        }
                    }
                    if (!found && n_neighbor_blocks < 256) { neighbor_blocks[n_neighbor_blocks++] = v_block; }
                }

                partition_t best = NO_ID;
                f64 min_cost = std::numeric_limits<f64>::max();

                for (const partition_t block: fine_blocks) {
                    if (p_manager.bweights[block] + u_weight > p_manager.get_lmax(block)) { continue; }

                    f64 cost = 0.0;
                    for (size_t x = 0; x < n_neighbor_blocks; ++x) {
                        if (neighbor_blocks[x] != block) { cost += ac.alpha_end(); }
                    }
                    for (size_t j = graphs.back().neighborhoods[u]; j < graphs.back().neighborhoods[u + 1]; ++j) {
                        const vertex_t v = graphs.back().edges_v[j];
                        const partition_t v_block = p_manager.partition[v];
                        if (v_block == NO_ID) { continue; }
                        const f64 w = static_cast<f64>(graphs.back().edge_w(j));
                        const f64 dist = static_cast<f64>(d_oracle.get(block, v_block));
                        cost += ac.beta_end() * w * static_cast<f64>(block != v_block)
                                + ac.alpha_link * dist
                                + ac.beta_link * w * dist;
                    }

                    // Penalize blocks above lmin to spread load across siblings
                    const weight_t block_lmin = p_manager.get_lmin(block);
                    if (block_lmin > 0 && p_manager.bweights[block] >= block_lmin) {
                        f64 over_ratio = static_cast<f64>(p_manager.bweights[block]) / static_cast<f64>(block_lmin);
                        cost += ac.alpha_end() * over_ratio;
                    }

                    if (best == NO_ID || cost < min_cost) {
                        min_cost = cost;
                        best = block;
                    }
                }

                // fallback: all fine blocks overloaded — pick least loaded
                if (best == NO_ID) {
                    best = fine_blocks[0];
                    weight_t min_w = p_manager.bweights[fine_blocks[0]];
                    for (size_t i = 1; i < fine_blocks.size(); ++i) {
                        if (p_manager.bweights[fine_blocks[i]] < min_w) {
                            min_w = p_manager.bweights[fine_blocks[i]];
                            best = fine_blocks[i];
                        }
                    }
                }

                p_manager.set(u, u_weight, best);
                blk_verts[best].push_back(u);
            }

            // Second pass: guarantee every fine block gets at least one vertex.
            // For each empty fine block e, steal the cheapest-to-move vertex from
            // its most-loaded sibling (same coarse parent).  "Cheapest" = minimum
            // increase in hop cost when switching from donor block s to e.
            const graph_t &g = graphs.back();
            for (partition_t c = 0; c < coarse_k; ++c) {
                const auto &siblings = coarse_to_fine[c];
                if (siblings.size() <= 1) { continue; }

                for (const partition_t e: siblings) {
                    if (p_manager.bweights[e] > 0) { continue; }

                    vertex_t best_u = NO_ID;
                    partition_t best_donor = NO_ID;
                    f64 best_delta = std::numeric_limits<f64>::max();

                    for (const partition_t s: siblings) {
                        if (s == e) { continue; }
                        for (const vertex_t u: blk_verts[s]) {
                            if (p_manager.partition[u] != s) { continue; } // already moved
                            if (g.v_weight(u) > p_manager.get_lmax(e)) { continue; }

                            // hop-cost delta for moving u from s to e
                            f64 delta = 0.0;
                            for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                                const partition_t vb = p_manager.partition[g.edges_v[j]];
                                if (vb == NO_ID || vb == s || vb == e) { continue; }
                                const f64 w = static_cast<f64>(g.edge_w(j));
                                delta += (ac.alpha_link + ac.beta_link * w)
                                        * (static_cast<f64>(d_oracle.get(e, vb))
                                           - static_cast<f64>(d_oracle.get(s, vb)));
                            }

                            if (delta < best_delta) {
                                best_delta = delta;
                                best_u = u;
                                best_donor = s;
                            }
                        }
                    }

                    if (best_u == NO_ID) { continue; } // coarse block has zero vertices — nothing to steal

                    // Perform the move (bv_manager/q_graph not built yet, use raw arrays)
                    const weight_t uw = g.v_weight(best_u);
                    p_manager.bweights[best_donor] -= uw;
                    p_manager.partition[best_u] = e;
                    p_manager.bweights[e] += uw;
                    blk_verts[e].push_back(best_u);
                }
            }

            // Third pass: avalanche steal along density gradient.
            // For each empty block, find the adjacent PE with highest weight,
            // pick a sender group from it (all vertices sending to the same
            // target), and move the whole group.  If that empties the donor,
            // the "hole" cascades along the gradient until it reaches a block
            // that can absorb the loss.  Each move is one grid hop.
            {
                const size_t W = static_cast<size_t>(d_oracle.current_W());
                auto get_neighbors = [&](partition_t b, partition_t *nbrs) -> size_t {
                    size_t cnt = 0;
                    const size_t c = b % W;
                    if (b >= static_cast<partition_t>(W)) nbrs[cnt++] = b - static_cast<partition_t>(W);
                    if (b + static_cast<partition_t>(W) < fine_k) nbrs[cnt++] = b + static_cast<partition_t>(W);
                    if (c > 0) nbrs[cnt++] = b - 1;
                    if (c + 1 < W && b + 1 < fine_k) nbrs[cnt++] = b + 1;
                    return cnt;
                };

                // scratch for building sender groups
                std::vector<size_t> target_idx(fine_k, SIZE_MAX);
                std::vector<partition_t> target_dirty;

                std::vector<partition_t> still_empty;
                for (partition_t e = 0; e < fine_k; ++e) {
                    if (p_manager.bweights[e] == 0) {
                        still_empty.push_back(e);
                    }
                }

                for (const partition_t e: still_empty) {
                    if (p_manager.bweights[e] > 0) { continue; }

                    partition_t hole = e;
                    constexpr size_t MAX_CHAIN = 64;
                    for (size_t step = 0; step < MAX_CHAIN; ++step) {
                        partition_t nbrs[4];
                        const size_t n_nbrs = get_neighbors(hole, nbrs);

                        // Pick neighbor with highest block weight
                        partition_t donor = NO_ID;
                        weight_t donor_bw = 0;
                        for (size_t i = 0; i < n_nbrs; ++i) {
                            if (p_manager.bweights[nbrs[i]] > donor_bw) {
                                donor_bw = p_manager.bweights[nbrs[i]];
                                donor = nbrs[i];
                            }
                        }
                        if (donor == NO_ID) { break; }

                        // Build sender groups in donor
                        std::vector<std::vector<vertex_t> > sg_verts;
                        std::vector<partition_t> sg_target;
                        std::vector<weight_t> sg_weight;

                        for (const vertex_t u: blk_verts[donor]) {
                            if (p_manager.partition[u] != donor) { continue; }
                            const weight_t uw = g.v_weight(u);
                            partition_t ext_tgts[256];
                            size_t n_ext = 0;
                            for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                                const partition_t vb = p_manager.partition[g.edges_v[j]];
                                if (vb == NO_ID || vb == donor) { continue; }
                                bool dup = false;
                                for (size_t x = 0; x < n_ext; ++x) {
                                    if (ext_tgts[x] == vb) {
                                        dup = true;
                                        break;
                                    }
                                }
                                if (!dup && n_ext < 256) { ext_tgts[n_ext++] = vb; }
                            }
                            for (size_t x = 0; x < n_ext; ++x) {
                                const partition_t D = ext_tgts[x];
                                if (target_idx[D] == SIZE_MAX) {
                                    target_idx[D] = sg_verts.size();
                                    sg_verts.emplace_back();
                                    sg_target.push_back(D);
                                    sg_weight.push_back(0);
                                    target_dirty.push_back(D);
                                }
                                const size_t idx = target_idx[D];
                                sg_verts[idx].push_back(u);
                                sg_weight[idx] += uw;
                            }
                        }
                        for (const partition_t D: target_dirty) { target_idx[D] = SIZE_MAX; }
                        target_dirty.clear();

                        // Pick lightest group that fits in hole
                        const weight_t hole_cap = p_manager.get_lmax(hole);
                        size_t pick = SIZE_MAX;
                        weight_t pick_w = std::numeric_limits<weight_t>::max();
                        for (size_t i = 0; i < sg_verts.size(); ++i) {
                            if (sg_weight[i] > hole_cap) { continue; }
                            if (sg_weight[i] < pick_w) {
                                pick_w = sg_weight[i];
                                pick = i;
                            }
                        }
                        if (pick == SIZE_MAX) { break; }

                        // Move the group
                        for (const vertex_t u: sg_verts[pick]) {
                            if (p_manager.partition[u] != donor) { continue; }
                            const weight_t uw = g.v_weight(u);
                            p_manager.bweights[donor] -= uw;
                            p_manager.partition[u] = hole;
                            p_manager.bweights[hole] += uw;
                            blk_verts[hole].push_back(u);
                        }

                        if (p_manager.bweights[donor] > 0) { break; }
                        hole = donor;
                    }
                }
            }

            // Fourth pass: evacuate small islands of non-empty PEs.
            // BFS on the PE grid to find connected components of non-empty
            // blocks.  Small components get their vertices moved to the
            // nearest block in the largest component.
            {
                const size_t W = static_cast<size_t>(d_oracle.current_W());
                auto get_neighbors = [&](partition_t b, partition_t *nbrs) -> size_t {
                    size_t cnt = 0;
                    const size_t c = b % W;
                    if (b >= static_cast<partition_t>(W)) nbrs[cnt++] = b - static_cast<partition_t>(W);
                    if (b + static_cast<partition_t>(W) < fine_k) nbrs[cnt++] = b + static_cast<partition_t>(W);
                    if (c > 0) nbrs[cnt++] = b - 1;
                    if (c + 1 < W && b + 1 < fine_k) nbrs[cnt++] = b + 1;
                    return cnt;
                };

                // BFS to label connected components of non-empty blocks
                std::vector<partition_t> comp_id(fine_k, NO_ID);
                std::vector<std::vector<partition_t> > components;
                std::vector<partition_t> bfs_queue;

                for (partition_t b = 0; b < fine_k; ++b) {
                    if (p_manager.bweights[b] == 0 || comp_id[b] != NO_ID) { continue; }
                    const partition_t cid = static_cast<partition_t>(components.size());
                    components.emplace_back();
                    bfs_queue.clear();
                    bfs_queue.push_back(b);
                    comp_id[b] = cid;
                    for (size_t head = 0; head < bfs_queue.size(); ++head) {
                        const partition_t cur = bfs_queue[head];
                        components[cid].push_back(cur);
                        partition_t nbrs[4];
                        const size_t n_nbrs = get_neighbors(cur, nbrs);
                        for (size_t i = 0; i < n_nbrs; ++i) {
                            if (p_manager.bweights[nbrs[i]] > 0 && comp_id[nbrs[i]] == NO_ID) {
                                comp_id[nbrs[i]] = cid;
                                bfs_queue.push_back(nbrs[i]);
                            }
                        }
                    }
                }

                // Find the largest component
                size_t largest_idx = 0;
                for (size_t i = 1; i < components.size(); ++i) {
                    if (components[i].size() > components[largest_idx].size()) {
                        largest_idx = i;
                    }
                }

                // Evacuate all small components into the largest
                for (size_t ci = 0; ci < components.size(); ++ci) {
                    if (ci == largest_idx) { continue; }
                    for (const partition_t b: components[ci]) {
                        if (p_manager.bweights[b] == 0) { continue; }
                        // Find nearest block in largest component with capacity
                        const weight_t b_weight = p_manager.bweights[b];
                        partition_t target = NO_ID;
                        weight_t best_dist = std::numeric_limits<weight_t>::max();
                        for (const partition_t s: components[largest_idx]) {
                            if (p_manager.bweights[s] + b_weight > p_manager.get_lmax(s)) { continue; }
                            weight_t d = d_oracle.get(b, s);
                            if (d < best_dist) {
                                best_dist = d;
                                target = s;
                            }
                        }
                        if (target == NO_ID) { continue; }
                        for (const vertex_t u: blk_verts[b]) {
                            if (p_manager.partition[u] != b) { continue; }
                            const weight_t uw = g.v_weight(u);
                            p_manager.bweights[b] -= uw;
                            p_manager.partition[u] = target;
                            p_manager.bweights[target] += uw;
                            blk_verts[target].push_back(u);
                        }
                    }
                }
            }

            // Fifth pass: brute-force sibling permutation optimization.
            // For each coarse parent with 2-4 children, try all permutations
            // of which vertex group maps to which child PE; keep the best.
            {
                auto perm_cost = [&](const std::vector<vertex_t> &verts, partition_t block) -> f64 {
                    f64 cost = 0.0;
                    for (const vertex_t u : verts) {
                        for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                            const partition_t nb = p_manager.partition[g.edges_v[j]];
                            if (nb == NO_ID) continue;
                            const f64 w = static_cast<f64>(g.edge_w(j));
                            const f64 dist = static_cast<f64>(d_oracle.get(block, nb));
                            if (block != nb) {
                                cost += ac.alpha_end() + ac.beta_end() * w;
                            }
                            cost += ac.alpha_link * dist + ac.beta_link * w * dist;
                        }
                    }
                    return cost;
                };

                for (partition_t c = 0; c < coarse_k; ++c) {
                    const auto &siblings = coarse_to_fine[c];
                    const size_t ns = siblings.size();
                    if (ns <= 1 || ns > 4) continue;

                    // Collect vertex sets for each sibling
                    std::vector<std::vector<vertex_t>> sib_verts(ns);
                    std::vector<weight_t> sib_weights(ns);
                    for (size_t i = 0; i < ns; ++i) {
                        for (const vertex_t u : blk_verts[siblings[i]]) {
                            if (p_manager.partition[u] == siblings[i]) {
                                sib_verts[i].push_back(u);
                            }
                        }
                        sib_weights[i] = p_manager.bweights[siblings[i]];
                    }

                    // Generate all permutations of indices [0..ns-1]
                    std::vector<size_t> perm(ns);
                    std::iota(perm.begin(), perm.end(), 0);

                    std::vector<size_t> best_perm = perm;
                    f64 best_cost = std::numeric_limits<f64>::max();

                    // Temporarily assign identity permutation to compute baseline
                    // (partition already has identity mapping)
                    do {
                        // Check capacity: vertex group i goes to PE siblings[perm[i]]
                        bool feasible = true;
                        for (size_t i = 0; i < ns; ++i) {
                            if (sib_weights[i] > p_manager.get_lmax(siblings[perm[i]])) {
                                feasible = false;
                                break;
                            }
                        }
                        if (!feasible) continue;

                        // Apply permutation temporarily
                        for (size_t i = 0; i < ns; ++i) {
                            for (const vertex_t u : sib_verts[i]) {
                                p_manager.partition[u] = siblings[perm[i]];
                            }
                        }

                        f64 cost = 0.0;
                        for (size_t i = 0; i < ns; ++i) {
                            cost += perm_cost(sib_verts[i], siblings[perm[i]]);
                        }

                        if (cost < best_cost) {
                            best_cost = cost;
                            best_perm = perm;
                        }
                    } while (std::next_permutation(perm.begin(), perm.end()));

                    // Apply best permutation
                    for (size_t i = 0; i < ns; ++i) {
                        for (const vertex_t u : sib_verts[i]) {
                            p_manager.partition[u] = siblings[best_perm[i]];
                        }
                        p_manager.bweights[siblings[best_perm[i]]] = sib_weights[i];
                    }

                    // Rebuild blk_verts for affected siblings
                    for (size_t i = 0; i < ns; ++i) {
                        blk_verts[siblings[i]].clear();
                    }
                    for (size_t i = 0; i < ns; ++i) {
                        for (const vertex_t u : sib_verts[i]) {
                            blk_verts[siblings[best_perm[i]]].push_back(u);
                        }
                    }
                }
            }

            // Sixth pass: diffuse vertices from dense borders into empty regions
            diffuse_into_empties();

            uncontraction_ms += get_milli_seconds(p, get_time_point());
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

            uncontraction_ms += get_milli_seconds(p, get_time_point());

            q_graph.compute_from_scratch(graphs.back(), p_manager);
            block_conn.compute_from_scratch(graphs.back(), p_manager);
            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
        }

        void rebalance(const u64 level) {
            auto p = get_time_point();

            if (level == 0) {
                rebalancer.rebalance_last_layer<true, true>(graphs.back(), p_manager, bv_manager, q_graph, d_oracle, block_conn);
            } else {
                rebalancer.rebalance<true, true>(graphs.back(), p_manager, bv_manager, q_graph, d_oracle, block_conn);
            }

            rebalance_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
        }

        void refinement() {
            auto p = get_time_point();

            if (ac.waver_label_propagation_config.enabled) {
                auto sp_local = get_time_point();
                lp_refine.refine(graphs.back(), d_oracle, bv_manager, p_manager, q_graph, block_conn);
                lp_refine_ms += get_milli_seconds(sp_local, get_time_point());
                {
                    partition_t ne = 0;
                    for (partition_t id = 0; id < ac.k; ++id)
                        if (p_manager.get_bweight(id) == 0 && p_manager.is_active(id)) ++ne;
                    if (ne > 0) std::cout << "    [after LP] empty=" << ne << std::endl;
                }
            }

            if (ac.waver_message_refinement_config.enabled) {
                auto sp_local = get_time_point();
                WaverMessageRefinement::refine(graphs.back(), p_manager, bv_manager, q_graph, d_oracle, block_conn,
                                               ac.waver_message_refinement_config,
                                               ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
                msg_refine_ms += get_milli_seconds(sp_local, get_time_point());
                {
                    partition_t ne = 0;
                    for (partition_t id = 0; id < ac.k; ++id)
                        if (p_manager.get_bweight(id) == 0 && p_manager.is_active(id)) ++ne;
                    if (ne > 0) std::cout << "    [after MsgRef] empty=" << ne << std::endl;
                }
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

            refinement_ms += get_milli_seconds(p, get_time_point());
            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), p_manager, bv_manager, q_graph, block_conn, ac.k));
        }
    };
}

#endif //HEIPROMAP_WAVER_MAP_SOLVER_H
