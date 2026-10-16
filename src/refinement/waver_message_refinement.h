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

#ifndef HEIPROMAP_WAVER_MESSAGE_REFINEMENT_H
#define HEIPROMAP_WAVER_MESSAGE_REFINEMENT_H

#include <limits>
#include <vector>

#include "../definitions.h"
#include "../datastructures/block_conn.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/large_quotient_graph.h"
#include "../datastructures/partition_manager.h"
#include "../utility/profiler.h"

namespace HeiProMap {

    // Configuration for WaverMessageRefinement.
    struct WaverMessageRefinementConfiguration {
        std::string name;
        bool enabled  = false;
        u64  max_iters = 5000;

        explicit WaverMessageRefinementConfiguration(std::string n) : name(std::move(n)) {}
    };

    // Message-group refinement.
    //
    // For each quotient edge (A, B), the "sender group" S_AB is the set of all
    // vertices in block A that have at least one application edge to block B.
    // S_AB collectively generates the single coalesced message A→B.
    //
    // We try to relocate S_AB as a unit to a candidate block C, which:
    //   - eliminates the A→B message pair (saving alpha_end)
    //   - either combines with an existing C→B message (no new pair cost)
    //     or creates a shorter-distance C→B message (saving hop cost)
    //
    // Candidates: block B itself, and all current neighbors of B in the quotient graph.
    //
    // Cost delta for moving group G from A to C (positive = improvement):
    //   alpha_end terms: track which block pairs appear/disappear via q_graph weights
    //   beta_end/alpha_link/beta_link terms: per-edge, skipping intra-group edges
    //
    class WaverMessageRefinement {
        using graph_t      = CSRGraph;
        using p_manager_t  = PartitionManager;
        using bv_manager_t = BoundaryVertexManager;
        using block_conn_t = BlockConn;
        using q_graph_t    = LargeQuotientGraph;

        static void apply_move(const graph_t      &g,
                               bv_manager_t       &bv_manager,
                               p_manager_t        &p_manager,
                               const vertex_t      u,
                               const weight_t      u_weight,
                               const partition_t   from_id,
                               const partition_t   to_id) {
            bv_manager.move(g, p_manager, u, from_id, to_id);
            p_manager.move_serial(u, u_weight, from_id, to_id);
        }

        // Compute the cost delta (positive = improvement) for moving all vertices in
        // `group` from block A to block C.
        //
        // `weight_to_D` : sum of edge weights from group to each external block D (D != A,
        //                 intra-group edges excluded).
        // `weight_to_A` : sum of edge weights from group to non-group vertices in A.
        template<typename DistanceOracleT>
        static f64 compute_group_delta(const graph_t                         &g,
                                       const partition_t                       A,
                                       const partition_t                       C,
                                       const std::vector<vertex_t>            &group,
                                       const std::vector<bool>                &in_group,
                                       const std::vector<std::pair<partition_t, weight_t>> &ext_neighbors,
                                       const weight_t                          weight_to_A,
                                       const p_manager_t                      &p_manager,
                                       const q_graph_t                        &q_graph,
                                       const DistanceOracleT                  &d_oracle,
                                       const block_conn_t                     &block_conn,
                                       const f64 alpha_end,
                                       const f64 beta_end,
                                       const f64 alpha_link,
                                       const f64 beta_link) {
            f64 delta = 0.0;

            // ---------------------------------------------------------------
            // alpha_end: track block-pair appearances / disappearances
            // ---------------------------------------------------------------
            // For each external block D that group connects to:
            //   A→D may disappear (if group carried all A→D weight)
            //   C→D may appear   (if C had no prior connection to D)
            for (const auto &[D, w_group_to_D] : ext_neighbors) {
                // A→D: disappears if group was the only A→D sender
                weight_t w_A_D = q_graph.get_weight(A, D);
                if (w_A_D > 0 && w_A_D == w_group_to_D) {
                    delta += alpha_end; // A→D pair removed
                }
                // C→D: appears if C had no prior connection to D (and D != C)
                if (D != C) {
                    weight_t w_C_D = q_graph.get_weight(C, D);
                    if (w_C_D == 0) {
                        delta -= alpha_end; // new C→D pair
                    }
                }
            }

            // Group edges to non-group A-vertices become C→A (intra → inter)
            if (weight_to_A > 0) {
                weight_t w_C_A = q_graph.get_weight(C, A);
                if (w_C_A == 0) {
                    delta -= alpha_end; // new C→A pair
                }
            }

            // ---------------------------------------------------------------
            // beta_end, alpha_link, beta_link: per edge, skip intra-group
            // ---------------------------------------------------------------
            for (const vertex_t u : group) {
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const vertex_t    v = g.edges_v[j];
                    if (in_group[v]) { continue; } // both endpoints move — no change

                    const partition_t D = p_manager[v];
                    const f64         w = static_cast<f64>(g.edge_w(j));

                    if (D == A) {
                        // intra (A→A) becomes inter (C→A)
                        delta -= beta_end  * w;
                        delta -= (alpha_link + beta_link * w) * static_cast<f64>(d_oracle.get(C, A));
                    } else if (D == C) {
                        // inter (A→C) becomes intra (C→C)
                        delta += beta_end  * w;
                        delta += (alpha_link + beta_link * w) * static_cast<f64>(d_oracle.get(A, C));
                    } else {
                        // A→D becomes C→D: only distance changes, both cut
                        const f64 d_old = static_cast<f64>(d_oracle.get(A, D));
                        const f64 d_new = static_cast<f64>(d_oracle.get(C, D));
                        delta += (alpha_link + beta_link * w) * (d_old - d_new);
                    }
                }
            }

            return delta;
        }

    public:
        template<typename DistanceOracleT>
        static void refine(graph_t                                   &g,
                           p_manager_t                               &p_manager,
                           bv_manager_t                              &bv_manager,
                           q_graph_t                                 &q_graph,
                           const DistanceOracleT                     &d_oracle,
                           block_conn_t                              &block_conn,
                           const WaverMessageRefinementConfiguration &config,
                           const f64 alpha_end,
                           const f64 beta_end,
                           const f64 alpha_link,
                           const f64 beta_link) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverMessageRefinement", "refine");

            const partition_t k = p_manager.get_k();

            std::vector<bool>   in_group(g.n, false);
            std::vector<vertex_t> group;

            // per-run accumulator for external neighbor weights
            std::vector<std::pair<partition_t, weight_t>> ext_neighbors;
            std::vector<weight_t> tmp_weight(k, 0); // scratch: D → total group→D weight
            std::vector<partition_t> dirty; // blocks touched in tmp_weight

            bool any_move = true;
            for (u64 iter = 0; iter < config.max_iters && any_move; ++iter) {
                HEIPROMAP_PROFILE_SCOPE("refinement", "WaverMessageRefinement", "iteration");
                any_move = false;

                // Snapshot quotient edges to avoid iterator invalidation during moves
                std::vector<std::pair<partition_t, partition_t>> edges;
                for (partition_t A = 0; A < k; ++A) {
                    for (const auto &he : q_graph.neighbors(A)) {
                        if (A < he.target) { // process each undirected pair once
                            edges.emplace_back(A, he.target);
                        }
                    }
                }

                for (const auto &[A, B] : edges) {
                    // Skip if pair no longer exists (prior move may have dissolved it)
                    if (q_graph.get_weight(A, B) == 0) { continue; }

                    // -----------------------------------------------------------
                    // Collect S_AB: boundary vertices in A connecting to B
                    // -----------------------------------------------------------
                    group.clear();
                    for (size_t i = 0; i < bv_manager.size(A); ++i) {
                        const vertex_t u = bv_manager.get(A, i);
                        for (size_t j = block_conn.start(u); j < block_conn.end(u); ++j) {
                            if (block_conn.get_id(j) == B) {
                                group.push_back(u);
                                in_group[u] = true;
                                break;
                            }
                        }
                    }
                    if (group.empty()) {
                        for (vertex_t u : group) { in_group[u] = false; }
                        continue;
                    }

                    // Also consider S_BA (vertices in B connecting to A) as a separate group
                    // — handled in the next iteration when (B,A) is processed as (A,B) with A<B swapped.
                    // Here we only move S_AB (group in A).

                    // -----------------------------------------------------------
                    // Build ext_neighbors and weight_to_A for this group
                    // -----------------------------------------------------------
                    weight_t group_weight  = 0;
                    weight_t weight_to_A   = 0;

                    for (const vertex_t u : group) {
                        group_weight += g.v_weight(u);
                        for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                            const vertex_t    v = g.edges_v[j];
                            if (in_group[v]) { continue; }
                            const partition_t D = p_manager[v];
                            const weight_t    w = g.edge_w(j);
                            if (D == A) {
                                weight_to_A += w;
                            } else {
                                if (tmp_weight[D] == 0) { dirty.push_back(D); }
                                tmp_weight[D] += w;
                            }
                        }
                    }
                    ext_neighbors.clear();
                    for (const partition_t D : dirty) {
                        ext_neighbors.emplace_back(D, tmp_weight[D]);
                        tmp_weight[D] = 0;
                    }
                    dirty.clear();

                    // -----------------------------------------------------------
                    // Try candidates: B itself + all current neighbors of B +
                    // empty blocks that are closer to B than A is
                    // -----------------------------------------------------------
                    partition_t best_C     = NO_ID;
                    f64         best_delta = 0.0;

                    auto try_candidate = [&](const partition_t C) {
                        if (C == A) { return; }
                        if (p_manager.bweights[C] + group_weight > p_manager.get_lmax(C)) { return; }
                        if (p_manager.bweights[A] - group_weight < p_manager.get_lmin(A)) { return; }
                        f64 d = compute_group_delta(g, A, C, group, in_group,
                                                    ext_neighbors, weight_to_A,
                                                    p_manager, q_graph, d_oracle, block_conn,
                                                    alpha_end, beta_end, alpha_link, beta_link);
                        if (d > best_delta) { best_delta = d; best_C = C; }
                    };

                    try_candidate(B); // move senders INTO the receiver block
                    for (const auto &he : q_graph.neighbors(B)) {
                        try_candidate(he.target);
                    }


                    // -----------------------------------------------------------
                    // Apply best move
                    // -----------------------------------------------------------
                    if (best_C != NO_ID) {
                        for (const vertex_t u : group) {
                            q_graph.move(g, p_manager, u, A, best_C);
                            block_conn.move(g, u, A, best_C);
                            apply_move(g, bv_manager, p_manager, u, g.v_weight(u), A, best_C);
                        }
                        any_move = true;
                    }

                    for (vertex_t u : group) { in_group[u] = false; }
                }
            }
        }
    };

} // namespace HeiProMap

#endif // HEIPROMAP_WAVER_MESSAGE_REFINEMENT_H
