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

#ifndef HEIPROMAP_PARALLEL_SOLVER_H
#define HEIPROMAP_PARALLEL_SOLVER_H

#include "parallel_boundary_vertex_manager.h"
#include "parallel_partition_manager.h"
#include "parallel_active_vertex_manager.h"
#include "parallel_graph.h"
#include "parallel_distance_oracle.h"
#include "parallel_quotient_graph.h"
#include "../utility/parallel_assert_states.h"

namespace HeiProMap {
    /**
     * Solver for parallel Process Mapping.
     */
    class ParallelSolver {
    private:
        ParallelAlgorithmConfiguration ac;

        std::vector<ParallelGraph> graphs;

        ParallelActiveVertexManager   av_manager;
        ParallelPartitionManager      p_manager;
        ParallelBoundaryVertexManager bv_manager;
        ParallelQuotientGraph         q_graph;
        ParallelDistanceOracle        d_oracle;

        // balance
        weight_t lmax = 0;

        // threads
        u64 threads = 1;

        // matching
        std::vector<EdgeUV *>    matches;
        std::vector<size_t>      matches_size;
        ParallelHeavyEdgeMatcher parallel_he_matcher;

        // refinement

        // statistics

    public:
        explicit ParallelSolver(ParallelAlgorithmConfiguration &t_ac) {
            ac = t_ac;

            // threads
            threads = ac.threads;

            const auto sp_graph_io = std::chrono::high_resolution_clock::now();
            graphs.emplace_back(ac.graph_in, threads);
            const auto ep_graph_io = std::chrono::high_resolution_clock::now();

            const auto sp_io = std::chrono::high_resolution_clock::now();

            // balance
            lmax = ceil((1.0 + ac.imbalance) * ((f64) graphs[0].get_weight() / (f64) ac.k));

            // manager
            av_manager.initialize(graphs[0].get_n());
            p_manager.initialize(graphs[0].get_n(), ac.k, lmax);
            bv_manager.initialize(graphs[0].get_n(), ac.k);
            q_graph.initialize(ac.k);
            HEAVYASSERT(assert_state_pre_partitioning(graphs[0], av_manager, threads));

            // distance
            d_oracle.initialize(ac.hierarchy, ac.distance, threads);

            // matching
            parallel_he_matcher.initialize(graphs[0].get_n(), graphs[0].get_m(), ac.k, lmax, ac.seed);

            // refinement

            const auto ep_io = std::chrono::high_resolution_clock::now();
        }

        std::vector<vertex_t> solve() {
            internal_solve();

            std::vector<partition_t> p(graphs.back().get_n());
            for (vertex_t            u = 0; u < graphs.back().get_n(); ++u) { p[u] = p_manager[u]; }

            return p;
        }

    private:
        void internal_solve() {
            s32 level = 0;

            std::cout << "level " << level << " size " << av_manager.get_n_active() << std::endl;

            while (av_manager.get_n_active() > ac.k * 64) {
                matching(level);
                coarsening(level);
                level += 1;
                std::cout << "level " << level << " size " << av_manager.get_n_active() << std::endl;
            }

            std::cout << "partition" << std::endl;
            partition();

            while (level > 0) {
                level -= 1;
                uncoarsening(level);
                refinement(level);
                std::cout << "level " << level << " size " << av_manager.get_n_active() << std::endl;
            }
        }

        void partition() {
            const auto sp_partition = std::chrono::high_resolution_clock::now();

            if (ac.parallel_partitioning_algorithm_id == PARALLEL_PARTITIONING_ALG_KAFFPA) {
                KaffpaPartitioner partitioner;
                partitioner.partition(ac.parallel_kaffpa_partitioner_config, graphs.back(), av_manager, p_manager, ac.hierarchy, ac.distance, ac.imbalance, ac.seed);
            } else {
                std::cout << "Parallel Partitioning algorithm " << parallel_partitioning_algorithm_to_string(ac.parallel_partitioning_algorithm_id) << " with id " << ac.parallel_partitioning_algorithm_id << " not known!" << std::endl;
                exit(EXIT_FAILURE);
            }

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
            std::cout << "Partitioning : " << get_seconds(sp_partition, ep_partition) << std::endl;

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), av_manager, p_manager, bv_manager, q_graph, ac.k, threads));
        }

        void matching(s32 level) {
            const auto sp_match = std::chrono::high_resolution_clock::now();

            size_t t_n_64       = round_up_64(av_manager.get_n_active() / 2);
            EdgeUV *matches_arr = (EdgeUV *) aligned_alloc(64, t_n_64 * sizeof(EdgeUV));
            matches.push_back(matches_arr);
            matches_size.push_back(0);

            if (ac.parallel_coarsening_algorithm_id == PARALLEL_COARSENING_ALG_HEAVY_MATCHING) {
                parallel_he_matcher.match(level, ac.parallel_heavy_edge_matcher_config, graphs.back(), av_manager, matches.back(), matches_size.back());
            } else {
                std::cout << "Coarsening algorithm " << coarsening_algorithm_to_string(ac.parallel_coarsening_algorithm_id) << " with id " << ac.parallel_coarsening_algorithm_id << " not known!" << std::endl;
                exit(EXIT_FAILURE);
            }

            const auto ep_match = std::chrono::high_resolution_clock::now();
            std::cout << "Matching : " << get_seconds(sp_match, ep_match) << std::endl;
        }

        void coarsening(s32 level) {
            const auto sp_coarse = std::chrono::high_resolution_clock::now();

            graphs.emplace_back(graphs.back(), matches.back(), matches_size.back(), threads); // coarse the graph
            av_manager.contract(matches.back(), matches_size.back());

            const auto ep_coarse = std::chrono::high_resolution_clock::now();
            std::cout << "Coarsening : " << get_seconds(sp_coarse, ep_coarse) << std::endl;

            HEAVYASSERT(assert_state_pre_partitioning(graphs.back(), av_manager, threads));
        }

        void uncoarsening(s32 level) {
            const auto sp_uncoarse = std::chrono::high_resolution_clock::now();

            p_manager.uncontract(matches.back(), matches_size.back());
            av_manager.uncontract(matches.back(), matches_size.back());
            bv_manager.uncontract(matches.back(), matches_size.back(), graphs[graphs.size() - 2], graphs[graphs.size() - 1], av_manager, p_manager);
            graphs.pop_back(); // this is doing uncontraction

            free(matches.back());
            matches.pop_back(); // throw away the matching, not needed anymore
            matches_size.pop_back();

            const auto ep_uncoarse = std::chrono::high_resolution_clock::now();
            std::cout << "Uncoarsening : " << get_seconds(sp_uncoarse, ep_uncoarse) << std::endl;

            HEAVYASSERT(assert_state_after_partitioning(graphs.back(), av_manager, p_manager, bv_manager, q_graph, ac.k, threads));
        }

        void refinement(s32 level) {

        }
    };
}
#endif //HEIPROMAP_PARALLEL_SOLVER_H
