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

// Message-group refinement (Waver). For every quotient edge (A, B), the sender group S_AB is the set of vertices of
// block A with an edge to block B; it makes up the single coalesced message A -> B. The refiner tries to move S_AB as a unit
// to B or to a quotient neighbor of B, with room, which can remove the A -> B message or shorten it, and applies the move with
// the largest drop of the message cost (alpha_end per pair, beta_end per unit of volume, alpha_link / beta_link per hop).
// Both orientations of every edge are visited, and only pairs whose inputs changed are evaluated again. It keeps its own
// quotient rows, so no quotient graph has to be maintained outside. Unit of work: whole sender groups; runs on every level.


#ifndef HEIPROMAP_WAVER_MESSAGE_REFINEMENT_H
#define HEIPROMAP_WAVER_MESSAGE_REFINEMENT_H

#include <algorithm>
#include <chrono>
#include <iostream>
#include <limits>
#include <vector>

#include "../definitions.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../utility/profiler.h"

namespace HeiProMap {
    // Configuration for WaverMessageRefinement.
    struct WaverMessageRefinementConfiguration {
        std::string name;
        bool enabled = false;
        u64 max_iters = 5000;
        bool verbose = false; // print the statistics of every iteration

        explicit WaverMessageRefinementConfiguration(std::string n) : name(std::move(n)) {
        }
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
    //   alpha_end terms: track which block pairs appear/disappear via the quotient weights
    //   beta_end/alpha_link/beta_link terms: per-edge, skipping intra-group edges
    //
    class WaverMergeRefinement;

    class WaverMessageRefinement {
        friend class WaverMergeRefinement; // uses the quotient rows, the move and its cost delta

        // A move must gain more than the floating point noise of the delta, or zero-gain moves cycle.
        static constexpr f64 kMinGain = 1e-6;

        using graph_t = CSRGraph;
        using p_manager_t = PartitionManager;
        using bv_manager_t = BoundaryVertexManager;

        // The quotient graph, kept as one row per block: the blocks it exchanges messages with (sorted by id) and the
        // edge weight to each. A row is built on first use from the boundary vertices of its block and dropped when a
        // move changes it, so no quotient graph has to be maintained outside this refinement.
        class QuotientRows {
        public:
            struct Entry {
                partition_t target;
                weight_t weight;
            };

            void initialize(const graph_t &t_g, const p_manager_t &t_pm, const bv_manager_t &t_bv, const partition_t k) {
                g = &t_g;
                pm = &t_pm;
                bv = &t_bv;
                rows.assign(k, {});
                valid.assign(k, 0);
                acc.assign(k, 0);
            }

            const std::vector<Entry> &neighbors(const partition_t x) {
                if (!valid[x]) { build(x); }
                return rows[x];
            }

            // Edge weight between the blocks x and y (x != y), 0 if they exchange no message.
            weight_t weight(const partition_t x, const partition_t y) {
                const std::vector<Entry> &row = neighbors(x);
                const auto it = std::lower_bound(row.begin(), row.end(), y,
                                                 [](const Entry &e, const partition_t t) { return e.target < t; });
                return it != row.end() && it->target == y ? it->weight : 0;
            }

            // Applies a change of the weight between the blocks x and y to the row of x. A row that has not been built
            // yet is left alone: it is built from the partition when it is first used.
            void add(const partition_t x, const partition_t y, const weight_t delta) {
                if (!valid[x]) { return; }
                std::vector<Entry> &row = rows[x];
                auto it = std::lower_bound(row.begin(), row.end(), y,
                                           [](const Entry &e, const partition_t t) { return e.target < t; });
                if (it != row.end() && it->target == y) {
                    it->weight += delta;
                    if (it->weight == 0) { row.erase(it); }
                } else {
                    row.insert(it, Entry{y, delta});
                }
            }

            // A group moved from the block a to the block c, after the move. ext holds the weight (d, w) from the group to
            // every block d != a, to_a the weight from the group to the vertices left in a.
            void apply_group_move(const partition_t a, const partition_t c,
                                  const std::vector<std::pair<partition_t, weight_t> > &ext, const weight_t to_a) {
                for (const auto &[d, w]: ext) {
                    add(a, d, -w);
                    if (d == c) {
                        add(c, a, -w); // the edges between the group and c are inside c now
                    } else {
                        add(c, d, w);
                        add(d, a, -w);
                        add(d, c, w);
                    }
                }
                if (to_a > 0) {
                    // the edges between the group and the rest of a are cut now
                    add(a, c, to_a);
                    add(c, a, to_a);
                }
            }

        private:
            const graph_t *g = nullptr;
            const p_manager_t *pm = nullptr;
            const bv_manager_t *bv = nullptr;
            std::vector<std::vector<Entry> > rows;
            std::vector<char> valid;
            std::vector<weight_t> acc; // scratch: weight from the block being built to each block
            std::vector<partition_t> dirty; // blocks with a non-zero entry in acc

            void build(const partition_t x) {
                for (size_t i = 0; i < bv->size(x); ++i) {
                    const vertex_t u = bv->get(x, i);
                    for (size_t j = g->neighborhoods[u]; j < g->neighborhoods[u + 1]; ++j) {
                        const partition_t d = (*pm)[g->edges_v[j]];
                        if (d == x) { continue; }
                        if (acc[d] == 0) { dirty.push_back(d); }
                        acc[d] += g->edge_w(j);
                    }
                }
                std::sort(dirty.begin(), dirty.end());
                std::vector<Entry> &row = rows[x];
                row.clear();
                for (const partition_t d: dirty) {
                    row.push_back({d, acc[d]});
                    acc[d] = 0;
                }
                dirty.clear();
                valid[x] = 1;
            }
        };

        static void apply_move(const graph_t &g,
                               bv_manager_t &bv_manager,
                               p_manager_t &p_manager,
                               const vertex_t u,
                               const weight_t u_weight,
                               const partition_t from_id,
                               const partition_t to_id) {
            bv_manager.move(g, p_manager, u, from_id, to_id);
            p_manager.move_serial(u, u_weight, from_id, to_id);
        }


        // Compute the cost delta (positive = improvement) for moving all vertices of a group from block A to block C.
        //
        // `ext_neighbors` : (D, w) for every block D != A that the group has edges to: their weight sum w (edges
        //                   inside the group are not counted: both endpoints move).
        // `weight_to_A`   : sum of edge weights from the group to vertices of A outside the group.
        template<typename DistanceOracleT>
        static f64 compute_group_delta(const partition_t A,
                                       const partition_t C,
                                       const std::vector<std::pair<partition_t, weight_t> > &ext_neighbors,
                                       const weight_t weight_to_A,
                                       QuotientRows &quotient,
                                       const DistanceOracleT &d_oracle,
                                       const f64 alpha_end,
                                       const f64 beta_end,
                                       const f64 alpha_link,
                                       const f64 beta_link) {
            f64 delta = 0.0;
            bool pair_A_C_removed = false;

            // ---------------------------------------------------------------
            // alpha_end, alpha_link: per block pair (one message), so they are charged when a pair appears or
            // disappears. A pair at distance h costs alpha_end + alpha_link * h.
            // ---------------------------------------------------------------
            // For each external block D that group connects to:
            //   A→D may disappear (if group carried all A→D weight)
            //   C→D may appear   (if C had no prior connection to D)
            for (const auto &[D, w_group_to_D]: ext_neighbors) {
                // A→D: disappears if group was the only A→D sender
                weight_t w_A_D = quotient.weight(A, D);
                if (w_A_D > 0 && w_A_D == w_group_to_D) {
                    delta += alpha_end + alpha_link * static_cast<f64>(d_oracle.get(A, D)); // A→D pair removed
                    if (D == C) { pair_A_C_removed = true; }
                }
                // C→D: appears if C had no prior connection to D (and D != C)
                if (D != C) {
                    weight_t w_C_D = quotient.weight(C, D);
                    if (w_C_D == 0) {
                        delta -= alpha_end + alpha_link * static_cast<f64>(d_oracle.get(C, D)); // new C→D pair
                    }
                }
            }

            // Group edges to non-group A-vertices become C→A (intra → inter)
            if (weight_to_A > 0) {
                weight_t w_C_A = quotient.weight(C, A);
                if (w_C_A == 0 || pair_A_C_removed) {
                    // new C→A pair, or the A→C pair removed above is recreated by these edges
                    delta -= alpha_end + alpha_link * static_cast<f64>(d_oracle.get(C, A));
                }
            }

            // ---------------------------------------------------------------
            // beta_end, beta_link: per unit of volume, summed per block instead of per edge
            // ---------------------------------------------------------------
            for (const auto &[D, w_group_to_D]: ext_neighbors) {
                const f64 w = static_cast<f64>(w_group_to_D);
                if (D == C) {
                    // inter (A→C) becomes intra (C→C)
                    delta += beta_end * w;
                    delta += beta_link * w * static_cast<f64>(d_oracle.get(A, C));
                } else {
                    // A→D becomes C→D: only distance changes, both cut
                    delta += beta_link * w * (static_cast<f64>(d_oracle.get(A, D)) - static_cast<f64>(d_oracle.get(C, D)));
                }
            }
            if (weight_to_A > 0) {
                // intra (A→A) becomes inter (C→A)
                const f64 w = static_cast<f64>(weight_to_A);
                delta -= beta_end * w;
                delta -= beta_link * w * static_cast<f64>(d_oracle.get(C, A));
            }

            return delta;
        }

    public:
        template<typename DistanceOracleT>
        static void refine(graph_t &g,
                           p_manager_t &p_manager,
                           bv_manager_t &bv_manager,
                           const DistanceOracleT &d_oracle,
                           const WaverMessageRefinementConfiguration &config,
                           const f64 alpha_end,
                           const f64 beta_end,
                           const f64 alpha_link,
                           const f64 beta_link) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverMessageRefinement", "refine");

            const partition_t k = p_manager.get_k();

            QuotientRows quotient;
            std::vector<bool> in_group;
            std::vector<std::vector<vertex_t> > groups; // groups[i]: S_AB of the i-th quotient neighbor B of A
            std::vector<partition_t> targets; // quotient neighbors of A
            std::vector<int32_t> slot_of; // block -> index in targets, -1 if not a neighbor
            std::vector<char> group_valid; // groups[i] is built and up to date
            std::vector<u64> seen_stamp; // block was already seen for the current vertex
            u64 stamp = 0;
            weight_t min_vw = 1; // no group is lighter than this

            // Only pairs whose inputs changed are evaluated again. The evaluation of (A, B) reads the blocks A, B and
            // the quotient neighbors of B (the candidates), and the quotient weights between them and the blocks
            // adjacent to the group. A move from A0 to C0 changes the rows of A0, C0 and of the blocks adjacent to the
            // group, so those are marked as touched. A pair is evaluated again in the next iteration if A is touched,
            // or B or a neighbor of B is touched.
            std::vector<u64> last_touched; // 1 + the iteration in which the block was last touched, 0 = never
            std::vector<char> touched_prev; // touched in the previous iteration
            std::vector<char> near_touched; // the block or one of its quotient neighbors is in touched_prev

            // per-run accumulator for external neighbor weights
            std::vector<std::pair<partition_t, weight_t> > ext_neighbors;
            std::vector<weight_t> tmp_weight; // scratch: D → total group→D weight
            std::vector<char> tmp_seen; // scratch: D is in dirty
            std::vector<partition_t> dirty; // blocks touched in tmp_weight

            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverMessageRefinement", "setup");
            quotient.initialize(g, p_manager, bv_manager, k);
            in_group.assign(g.n, false);
            tmp_weight.assign(k, 0);
            tmp_seen.assign(k, 0);
            slot_of.assign(k, -1);
            seen_stamp.assign(k, 0);
            last_touched.assign(k, 0);
            min_vw = 1;
            if (!g.uniform_v_weights && g.n > 0) {
                min_vw = g.v_weight(0);
                for (vertex_t u = 1; u < g.n; ++u) { min_vw = std::min(min_vw, g.v_weight(u)); }
            }
            touched_prev.assign(k, 0);
            near_touched.assign(k, 0);


            bool any_move = true;
            for (u64 iter = 0; iter < config.max_iters && any_move; ++iter) {
                HEIPROMAP_PROFILE_SCOPE("refinement", "WaverMessageRefinement", "iteration");
                any_move = false;

                // A group can only move into a block with room left. Without any, nothing can be accepted.
                bool has_slack = false;
                for (partition_t C = 0; C < k && !has_slack; ++C) {
                    has_slack = p_manager.bweights[C] < p_manager.get_lmax(C);
                }
                if (!has_slack) { break; }

                // The first iteration evaluates every pair.
                if (iter > 0) {
                    std::fill(near_touched.begin(), near_touched.end(), 0);
                    for (partition_t X = 0; X < k; ++X) {
                        touched_prev[X] = last_touched[X] >= iter; // touched in iteration iter - 1 (stored as iter)
                        if (!touched_prev[X]) { continue; }
                        near_touched[X] = 1;
                        for (const auto &he: quotient.neighbors(X)) { near_touched[he.target] = 1; }
                    }
                }


                for (partition_t A = 0; A < k; ++A) {
                    // Snapshot the neighbors: moves change the quotient graph. Both orientations of a quotient edge
                    // are visited, so the sender group of either block can move.
                    targets.clear();
                    for (const auto &he: quotient.neighbors(A)) {
                        if (iter == 0 || touched_prev[A] || near_touched[he.target]) { targets.push_back(he.target); }
                    }
                    if (targets.empty()) { continue; }

                    group_valid.assign(targets.size(), 0);
                    if (groups.size() < targets.size()) { groups.resize(targets.size()); }

                    // A necessary condition for a move of S_AB (whose weight is at least min_vw): A can lose that
                    // much, and B or a neighbor of B can take it. Most pairs fail this, so their groups are not built.
                    auto may_move = [&](const partition_t B) {
                        if (p_manager.bweights[A] - min_vw < p_manager.get_lmin(A)) { return false; }
                        if (p_manager.bweights[B] + min_vw <= p_manager.get_lmax(B)) { return true; }
                        for (const auto &he: quotient.neighbors(B)) {
                            if (he.target != A && p_manager.bweights[he.target] + min_vw <= p_manager.get_lmax(he.target)) {
                                return true;
                            }
                        }
                        return false;
                    };

                    for (size_t idx = 0; idx < targets.size(); ++idx) {
                        const partition_t B = targets[idx];
                        // Skip if pair no longer exists (prior move may have dissolved it)
                        if (quotient.weight(A, B) == 0) { continue; }
                        if (!may_move(B)) { continue; }

                        // -------------------------------------------------------
                        // Collect S_AB: boundary vertices in A connecting to B. One pass over the boundary of A
                        // builds the groups of this and all later neighbors that may move; a move invalidates them.
                        // -------------------------------------------------------
                        if (!group_valid[idx]) {
                            for (size_t i = idx; i < targets.size(); ++i) {
                                group_valid[i] = i == idx || may_move(targets[i]);
                                if (group_valid[i]) {
                                    groups[i].clear();
                                    slot_of[targets[i]] = static_cast<int32_t>(i);
                                }
                            }
                            for (size_t i = 0; i < bv_manager.size(A); ++i) {
                                const vertex_t u = bv_manager.get(A, i);
                                ++stamp;
                                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                                    const partition_t D = p_manager[g.edges_v[j]];
                                    if (seen_stamp[D] == stamp) { continue; }
                                    seen_stamp[D] = stamp;
                                    const int32_t slot = slot_of[D];
                                    if (slot >= 0) { groups[slot].push_back(u); }
                                }
                            }
                            for (size_t i = idx; i < targets.size(); ++i) { slot_of[targets[i]] = -1; }
                        }
                        const std::vector<vertex_t> &group = groups[idx];
                        if (group.empty()) {
                            continue;
                        }
                        for (const vertex_t u: group) { in_group[u] = true; }

                        // Only S_AB (the group in A) moves here; S_BA is tried when B is the outer block.

                        // -----------------------------------------------------------
                        // Build ext_neighbors and weight_to_A for this group
                        // -----------------------------------------------------------
                        weight_t group_weight = 0;
                        weight_t weight_to_A = 0;

                        for (const vertex_t u: group) { group_weight += g.v_weight(u); }

                        // With the exact weight: is there any candidate that can take the group? Else skip the scan.
                        bool any_candidate = p_manager.bweights[A] - group_weight >= p_manager.get_lmin(A)
                                             && p_manager.bweights[B] + group_weight <= p_manager.get_lmax(B);
                        if (!any_candidate && p_manager.bweights[A] - group_weight >= p_manager.get_lmin(A)) {
                            for (const auto &he: quotient.neighbors(B)) {
                                if (he.target != A
                                    && p_manager.bweights[he.target] + group_weight <= p_manager.get_lmax(he.target)) {
                                    any_candidate = true;
                                    break;
                                }
                            }
                        }
                        if (!any_candidate) {
                            for (const vertex_t u: group) { in_group[u] = false; }
                            continue;
                        }

                        // The weight from the group to every other block, and to the vertices of A outside the group.
                        for (const vertex_t u: group) {
                            for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                                const vertex_t v = g.edges_v[j];
                                const partition_t D = p_manager[v];
                                const weight_t w = g.edge_w(j);
                                if (D == A) {
                                    if (!in_group[v]) { weight_to_A += w; }
                                    continue;
                                }
                                if (!tmp_seen[D]) {
                                    tmp_seen[D] = 1;
                                    dirty.push_back(D);
                                }
                                tmp_weight[D] += w;
                            }
                        }
                        ext_neighbors.clear();
                        for (const partition_t D: dirty) {
                            ext_neighbors.emplace_back(D, tmp_weight[D]);
                            tmp_weight[D] = 0;
                            tmp_seen[D] = 0;
                        }
                        dirty.clear();

                        // -----------------------------------------------------------
                        // Try candidates: B itself + all current neighbors of B +
                        // empty blocks that are closer to B than A is
                        // -----------------------------------------------------------
                        partition_t best_C = NO_ID;
                        f64 best_delta = 0.0;

                        auto try_candidate = [&](const partition_t C) {
                            if (C == A) { return; }
                            if (p_manager.bweights[C] + group_weight > p_manager.get_lmax(C)) { return; }
                            if (p_manager.bweights[A] - group_weight < p_manager.get_lmin(A)) { return; }
                            f64 d = compute_group_delta(A, C, ext_neighbors, weight_to_A, quotient, d_oracle,
                                                        alpha_end, beta_end, alpha_link, beta_link);
                            if (d > best_delta + kMinGain) {
                                best_delta = d;
                                best_C = C;
                            }
                        };

                        try_candidate(B); // move senders INTO the receiver block
                        for (const auto &he: quotient.neighbors(B)) {
                            try_candidate(he.target);
                        }

                        // -----------------------------------------------------------
                        // Apply best move
                        // -----------------------------------------------------------
                        if (best_C != NO_ID) {
                            for (const vertex_t u: group) {
                                apply_move(g, bv_manager, p_manager, u, g.v_weight(u), A, best_C);
                            }
                            quotient.apply_group_move(A, best_C, ext_neighbors, weight_to_A);
                            last_touched[A] = last_touched[best_C] = iter + 1;
                            for (const auto &[D, w_group_to_D]: ext_neighbors) { last_touched[D] = iter + 1; }
                            any_move = true;
                            std::fill(group_valid.begin(), group_valid.end(), 0);
                        }

                        for (vertex_t u: group) { in_group[u] = false; }
                    }
                }
            }
        }
    };
} // namespace HeiProMap

#endif //HEIPROMAP_WAVER_MESSAGE_REFINEMENT_H
