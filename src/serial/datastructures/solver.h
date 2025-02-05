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
#include "boundary_vertex_manager.h"
#include "graph_csr.h"
#include "partition_manager.h"
#include "quotient_graph.h"
#include "statistic_collector.h"
#include "../../definitions.h"
#include "../../macros.h"
#include "../coarsening/heavy_edge_matcher.h"
#include "../partitioning/global_multisection.h"
#include "../partitioning/kaffpa_partitioner.h"
#include "../refinement/label_propagation_refinement_Faraj20.h"
#include "../refinement/quotient_graph_refinement_Faraj20.h"
#include "../utility/assert_state.h"
#include "../utility/qap.h"
#include "../utility/utils.h"
#include "../coarsening/global_path_algorithm.h"

namespace HeiProMap {
    /**
     * Solver for serial Process Mapping.
     */
    class Solver {
        std::vector<GraphCSR> graphs;
        ActiveVertexManager   av_manager;
        PartitionManager      p_manager;
        BoundaryVertexManager bv_manager;
        QuotientGraph         q_graph;

        // distance
        std::vector<partition_t> hierarchy;
        std::vector<weight_t>    distance;
        partition_t              k;
        DistanceOracle           d_oracle;

        // balance
        f64      m_imbalance = 0.0;
        weight_t lmax        = 0;

        // multilevel
        vertex_t threshold;

        // matching
        std::vector<std::vector<EdgeUV>> matches;
        HeavyEdgeMatcher                 he_matcher;
        GlobalPathAlgorithmMatcher       gpa_matcher;

        // refinement
        LabelPropagationRefinementFaraj20 lp_refine_faraj20;
        QuotientGraphRefinementFaraj20    qg_refine_faraj20;

        // statistics
        StatisticCollector stat_collect;

    public:
        Solver(const std::string &t_graph_in,
               const std::vector<partition_t> &t_hierarchy,
               const std::vector<weight_t> &t_distance,
               const f64 t_imbalance) {
            const auto sp_graph_io = std::chrono::high_resolution_clock::now();
            graphs.emplace_back(t_graph_in);
            const auto ep_graph_io = std::chrono::high_resolution_clock::now();

            const auto sp_io = std::chrono::high_resolution_clock::now();

            hierarchy = t_hierarchy;
            distance  = t_distance;
            k         = prod<partition_t>(hierarchy);

            // manager
            av_manager.initialize(graphs.back().get_n());
            p_manager.initialize(graphs.back().get_n(), k);
            bv_manager.initialize(graphs.back().get_n(), k);
            q_graph.initialize(k);
            HEAVYASSERT(assert_state_pre_partitioning(graphs.back(), av_manager));

            // distance
            d_oracle.initialize(hierarchy, distance);

            // balance
            m_imbalance = t_imbalance;
            lmax        = std::ceil((1.0 + t_imbalance) * ((f64) graphs.back().get_weight() / (f64) k));

            // multilevel
            threshold = 0; // graphs.back().get_n() / 500;

            // matching
            // ge_matcher.initialize(graphs.back().get_n());
            he_matcher.initialize(graphs.back().get_n(), lmax);
            gpa_matcher.initialize(graphs.back().get_n(), lmax);

            // refinement
            lp_refine_faraj20.initialize(graphs.back().get_n(), hierarchy, distance, lmax);
            qg_refine_faraj20.initialize(graphs.back().get_n(), hierarchy, distance, lmax);

            const auto ep_io = std::chrono::high_resolution_clock::now();
            stat_collect.set_io(get_seconds(sp_graph_io, ep_graph_io), get_seconds(sp_io, ep_io));
        }

        std::vector<vertex_t> solve() {
            internal_solve();

#if STATISTICCOLLECTOR
            weight_t              qap = get_qap(graphs.back(), av_manager, p_manager, d_oracle);
            std::vector<weight_t> pweights;
            for (partition_t      id  = 0; id < k; ++id) { pweights.push_back(p_manager.get_bweight(id)); }
            stat_collect.set_final(qap, pweights, lmax);
#endif
            stat_collect.finalize();

            std::cout << stat_collect.to_JSON() << std::endl;

            std::vector<partition_t> p(graphs.back().get_n());
            for (vertex_t            u = 0; u < graphs.back().get_n(); ++u) { p[u] = p_manager[u]; }

            return p;
        }

    private:
        void internal_solve() {
            s32 level = 0;

            while (av_manager.get_n_active() > k * 64) {
                matching(level);
                coarsening(level);

                weight_t      max_w = 0;
                for (vertex_t u: av_manager) {
                    max_w = std::max(max_w, graphs.back().get_weight(u));
                }

                std::cout << level << " " << av_manager.get_n_active() << " " << graphs.back().get_m() << " heavy vertex " << max_w << std::endl;

                level += 1;
            }

            partition();

            while (level > 0) {
                level -= 1;
                uncoarsening(level);
                refinement(level);

                std::cout << level << " " << av_manager.get_n_active() << " " << graphs.back().get_m() << std::endl;
            }
        }

        void partition() {
            const auto sp_partition = std::chrono::high_resolution_clock::now();

            // KaffpaPartitioner partitioner;
            GlobalMultisectionPartitioner partitioner;
            partitioner.partition(graphs.back(), av_manager, p_manager, hierarchy, distance, m_imbalance);

            // initialize boundary vertices and quotient graph
            for (const vertex_t u: av_manager) {
                for (size_t i = 0; i < graphs.back().size(u); ++i) {
                    const vertex_t    v    = graphs.back().neighbor(u, i);
                    const weight_t    w    = graphs.back().get_weight(u, i);
                    const partition_t u_id = p_manager[u];
                    const partition_t v_id = p_manager[v];

                    if (u_id != v_id) {
                        bv_manager.add(u, u_id); // boundary vertex
                        q_graph.add_edge(u_id, v_id, w); // quotient graph
                    }
                }
            }

            const auto ep_partition = std::chrono::high_resolution_clock::now();
            stat_collect.set_partition_time(get_seconds(sp_partition, ep_partition));

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), av_manager, p_manager, bv_manager, k));
#if STATISTICCOLLECTOR
            stat_collect.set_partition_stats(get_qap(graphs.back(), av_manager, p_manager, d_oracle), p_manager.get_bweights(), lmax);
#endif
        }

        void matching(const s32 level) {
            const auto sp_match = std::chrono::high_resolution_clock::now();

            matches.emplace_back();
            matches.back().reserve(av_manager.get_n_active() / 2);

            // ge_matcher.match(graphs.back(), av_manager, matches.back());
            // gpa_matcher.match(graphs.back(), av_manager, matches.back());
            he_matcher.match(graphs.back(), av_manager, matches.back());

            const auto ep_match = std::chrono::high_resolution_clock::now();
            stat_collect.set_matching_time(get_seconds(sp_match, ep_match), level);

#if STATISTICCOLLECTOR
            stat_collect.set_matching_stats(level, matches.back().size());
#endif
        }

        void coarsening(const s32 level) {
            const auto sp_coarse = std::chrono::high_resolution_clock::now();

            graphs.emplace_back(graphs.back(), matches.back()); // coarse the graph
            av_manager.contract(matches.back());

            const auto ep_coarse = std::chrono::high_resolution_clock::now();
            stat_collect.set_coarsening_time(get_seconds(sp_coarse, ep_coarse), level);

            HEAVYASSERT(assert_state_pre_partitioning(graphs.back(), av_manager));
#if STATISTICCOLLECTOR
            stat_collect.set_coarsening_stats(av_manager.get_n_active(), level);
#endif
        }

        void uncoarsening(const s32 level) {
            const auto sp_uncoarse = std::chrono::high_resolution_clock::now();

            p_manager.uncontract(matches.back());
            av_manager.uncontract(matches.back());
            bv_manager.uncontract(matches.back(), graphs[graphs.size() - 2], graphs[graphs.size() - 1], av_manager, p_manager);
            graphs.pop_back(); // this is doing uncontraction

            matches.pop_back(); // throw away the matching, not needed anymore

            const auto ep_uncoarse = std::chrono::high_resolution_clock::now();
            stat_collect.set_uncoarsening_time(get_seconds(sp_uncoarse, ep_uncoarse), level);

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), av_manager, p_manager, bv_manager, k));
#if STATISTICCOLLECTOR
            stat_collect.set_uncoarsening_stats(level, av_manager.get_n_active());
#endif
        }

        void refinement(const s32 level) {
            const auto sp_refinement = std::chrono::high_resolution_clock::now();

            // lp_refine_faraj20.refine(graphs.back(), av_manager, bv_manager, p_manager, d_oracle, q_graph);
            // qg_refine_faraj20.refine(graphs.back(), av_manager, bv_manager, p_manager, d_oracle, q_graph);

            const auto ep_refinement = std::chrono::high_resolution_clock::now();
            stat_collect.set_refinement_time(get_seconds(sp_refinement, ep_refinement), level);

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), av_manager, p_manager, bv_manager, k));
#if STATISTICCOLLECTOR
            stat_collect.set_refinement_stats(level, get_qap(graphs.back(), av_manager, p_manager, d_oracle));
#endif
        }
    };
}

#endif //HEIPROMAP_SOLVER_H
