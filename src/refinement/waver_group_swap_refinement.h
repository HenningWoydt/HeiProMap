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

// Group and vertex swap refinement (Waver). The swap refinement on the group level: two PEs exchange one sender group each,
// instead of the whole contents. On the vertex level the two PEs exchange one single vertex each: a vertex of A with an edge
// to C for a vertex of C with an edge to A (the best few of each are tried). For a PE A and a PE C within the radius on the grid, a sender group of A (the vertices of A with an
// edge to a block E, where C is E or a neighbor of E) moves to C while a sender group of C moves to A, if both PEs stay
// within their weight limits (a PE outside them already may not get worse). The pairs are ranked by the cost change of the
// two moves alone; the best few are applied one after the other with their exact cost change and undone again, and the
// best that pays is kept. Rounds repeat until no exchange pays. Options: --max-group-swaps and, with single vertices, --max-vertex-swaps.

#ifndef HEIPROMAP_WAVER_GROUP_SWAP_REFINEMENT_H
#define HEIPROMAP_WAVER_GROUP_SWAP_REFINEMENT_H

#include <algorithm>
#include <vector>

#include "../definitions.h"
#include "../utility/profiler.h"
#include "waver_group_moves.h"
#include "waver_swap_refinement.h"

namespace HeiProMap {
    class WaverGroupSwapRefinement {
        using Base = WaverGroupMoves;
        using graph_t = Base::graph_t;
        using p_manager_t = Base::p_manager_t;
        using bv_manager_t = Base::bv_manager_t;
        using Candidate = Base::Candidate;
        using Move = Base::Move;

        // the exchanges whose exact cost change is computed per pair of PEs
        static constexpr size_t kExact = 6;

        template<typename DistanceOracleT>
        class GroupSwapEngine : public Base::Engine<DistanceOracleT> {
            using E = Base::Engine<DistanceOracleT>;
            using E::g; using E::pm; using E::oracle; using E::quotient; using E::locked;
            using E::collect_group; using E::measure; using E::apply; using E::delta; using E::undo; using E::targets_of;

            // A sender group of a block: the vertices with an edge to `toward`.
            struct Group {
                partition_t toward;
                std::vector<vertex_t> vertices;
                weight_t weight;
            };

            // One exchange: group i of block a for group j of block c.
            struct Exchange {
                f64 approx;
                partition_t a, c;
                size_t i, j;
            };

            const WaverSwapRefinementConfiguration &cfg;
            std::vector<std::vector<Group> > groups;
            std::vector<char> fresh;

        public:
            GroupSwapEngine(graph_t &t_g, p_manager_t &t_pm, bv_manager_t &t_bv, const DistanceOracleT &t_oracle,
                            const WaverSwapRefinementConfiguration &t_cfg,
                            const f64 t_alpha_end, const f64 t_beta_end, const f64 t_alpha_link, const f64 t_beta_link)
                : E(t_g, t_pm, t_bv, t_oracle, t_alpha_end, t_beta_end, t_alpha_link, t_beta_link), cfg(t_cfg),
                  groups(this->k), fresh(this->k, 0) {}

            // Exchanges until a round finds none or max_swaps are made. Returns the number of exchanges.
            u64 run(const bool vertices, const u64 max_swaps) {
                const partition_t k = this->k;
                u64 total = 0;
                for (u64 round = 0; round < cfg.max_rounds && total < max_swaps; ++round) {
                    u64 swaps = 0;
                    for (partition_t a = 0; a < k && total + swaps < max_swaps; ++a) {
                        swaps += (vertices ? improve_vertices(a) : improve(a)) ? 1 : 0;
                    }
                    total += swaps;
                    if (swaps == 0) { break; }
                }
                return total;
            }

        private:
            // The sender groups of a block towards each of its quotient neighbors, kept until a swap changes them.
            const std::vector<Group> &groups_of(const partition_t b) {
                if (!fresh[b]) {
                    groups[b].clear();
                    std::vector<partition_t> toward;
                    for (const auto &he: quotient.neighbors(b)) { toward.push_back(he.target); }
                    for (const partition_t e: toward) {
                        Group group;
                        group.toward = e;
                        collect_group(b, e, group.vertices);
                        if (group.vertices.empty()) { continue; }
                        group.weight = 0;
                        for (const vertex_t u: group.vertices) { group.weight += g.v_weight(u); }
                        groups[b].push_back(std::move(group));
                    }
                    fresh[b] = 1;
                }
                return groups[b];
            }

            // A block may end in its limits, widened to the weight it has now (it may not get worse).
            bool in_band(const partition_t b, const weight_t end) const {
                return end >= std::min(pm.get_lmin(b), pm.bweights[b]) && end <= std::max(pm.get_lmax(b), pm.bweights[b]);
            }

            // The groups of a toward e can go to the blocks e and the neighbors of e; the group of c that comes back is one
            // that may go to a: its target block is a, or a neighbor of the block it has edges to.
            bool may_go_to(const partition_t c, const partition_t toward, const partition_t a) {
                return toward == a || quotient.weight(toward, a) > 0;
            }

            // Applies the exchange and returns its exact cost change (positive = better); the moves are undone if keep is false.
            f64 exchange(const partition_t a, const partition_t c, const Group &first, const Group &second, const bool keep) {
                std::vector<Move> done;
                f64 total = 0.0;
                Candidate one;
                one.group = first.vertices;
                measure(a, one);
                total += delta(a, c, one);
                apply(one, a, c);
                done.push_back(Move{one.group, a, c});
                Candidate two;
                two.group = second.vertices;
                measure(c, two);
                total += delta(c, a, two);
                apply(two, c, a);
                done.push_back(Move{two.group, c, a});
                if (!keep) {
                    for (size_t i = done.size(); i-- > 0;) { undo(done[i]); }
                }
                return total;
            }

            // The best exchange that a can take part in; returns whether one was made.
            bool improve(const partition_t a) {
                std::vector<Exchange> candidates;
                const std::vector<Group> mine = groups_of(a);
                std::vector<partition_t> targets;
                for (size_t i = 0; i < mine.size(); ++i) {
                    Candidate one;
                    one.group = mine[i].vertices;
                    measure(a, one);
                    targets_of(a, mine[i].toward, targets);
                    for (const partition_t c: targets) {
                        if (static_cast<u64>(oracle.get(a, c)) > cfg.radius) { continue; }
                        const std::vector<Group> theirs = groups_of(c);
                        for (size_t j = 0; j < theirs.size(); ++j) {
                            if (!may_go_to(c, theirs[j].toward, a)) { continue; }
                            if (!in_band(a, pm.bweights[a] - mine[i].weight + theirs[j].weight)) { continue; }
                            if (!in_band(c, pm.bweights[c] - theirs[j].weight + mine[i].weight)) { continue; }
                            Candidate two;
                            two.group = theirs[j].vertices;
                            measure(c, two);
                            candidates.push_back({delta(a, c, one) + delta(c, a, two), a, c, i, j});
                        }
                    }
                }
                if (candidates.empty()) { return false; }
                const size_t exact = std::min(kExact, candidates.size());
                std::partial_sort(candidates.begin(), candidates.begin() + exact, candidates.end(),
                                  [](const Exchange &x, const Exchange &y) { return x.approx > y.approx; });
                f64 best = cfg.min_gain;
                size_t best_index = candidates.size();
                for (size_t t = 0; t < exact; ++t) {
                    const Exchange &x = candidates[t];
                    const f64 total = exchange(x.a, x.c, mine[x.i], groups_of(x.c)[x.j], false);
                    if (total > best) {
                        best = total;
                        best_index = t;
                    }
                }
                if (best_index == candidates.size()) { return false; }
                const Exchange x = candidates[best_index];
                const Group second = groups_of(x.c)[x.j];
                exchange(x.a, x.c, mine[x.i], second, true);
                // the groups of both blocks and of the blocks around them have changed
                for (const partition_t b: {x.a, x.c}) {
                    fresh[b] = 0;
                    for (const auto &he: quotient.neighbors(b)) { fresh[he.target] = 0; }
                }
                return true;
            }

            // A vertex of a and a vertex of a neighboring block c change places. For every neighbor c within the radius, the
            // best few vertices of a towards c and of c towards a (by the cost change of their move alone) are paired and the
            // pairs are applied one after the other with their exact cost change.
            bool improve_vertices(const partition_t a) {
                static constexpr size_t kTop = 4;
                struct Ranked {
                    f64 gain;
                    vertex_t vertex;
                };
                const std::vector<Group> mine = groups_of(a);
                f64 best = cfg.min_gain;
                partition_t best_c = NO_ID;
                vertex_t best_u = 0, best_v = 0;
                for (const Group &toward_c: mine) {
                    const partition_t c = toward_c.toward;
                    if (static_cast<u64>(oracle.get(a, c)) > cfg.radius) { continue; }
                    const std::vector<Group> theirs = groups_of(c);
                    const Group *toward_a = nullptr;
                    for (const Group &x: theirs) {
                        if (x.toward == a) { toward_a = &x; }
                    }
                    if (toward_a == nullptr) { continue; }
                    std::vector<Ranked> from_a, from_c;
                    for (const vertex_t u: toward_c.vertices) {
                        Candidate one;
                        one.group = {u};
                        measure(a, one);
                        from_a.push_back({delta(a, c, one), u});
                    }
                    for (const vertex_t v: toward_a->vertices) {
                        Candidate two;
                        two.group = {v};
                        measure(c, two);
                        from_c.push_back({delta(c, a, two), v});
                    }
                    const auto better = [](const Ranked &x, const Ranked &y) { return x.gain > y.gain; };
                    std::partial_sort(from_a.begin(), from_a.begin() + std::min(kTop, from_a.size()), from_a.end(), better);
                    std::partial_sort(from_c.begin(), from_c.begin() + std::min(kTop, from_c.size()), from_c.end(), better);
                    for (size_t i = 0; i < std::min(kTop, from_a.size()); ++i) {
                        for (size_t j = 0; j < std::min(kTop, from_c.size()); ++j) {
                            const vertex_t u = from_a[i].vertex, v = from_c[j].vertex;
                            const weight_t wu = g.v_weight(u), wv = g.v_weight(v);
                            if (!in_band(a, pm.bweights[a] - wu + wv) || !in_band(c, pm.bweights[c] - wv + wu)) { continue; }
                            const f64 total = exchange(a, c, Group{c, {u}, wu}, Group{a, {v}, wv}, false);
                            if (total > best) {
                                best = total;
                                best_c = c;
                                best_u = u;
                                best_v = v;
                            }
                        }
                    }
                }
                if (best_c == NO_ID) { return false; }
                exchange(a, best_c, Group{best_c, {best_u}, g.v_weight(best_u)}, Group{a, {best_v}, g.v_weight(best_v)}, true);
                for (const partition_t b: {a, best_c}) {
                    fresh[b] = 0;
                    for (const auto &he: quotient.neighbors(b)) { fresh[he.target] = 0; }
                }
                return true;
            }
        };

    public:
        // Exchanges of a sender group of each of two PEs, at most config.max_group_swaps. Returns their number.
        template<typename DistanceOracleT>
        static u64 refine_groups(graph_t &g, p_manager_t &p_manager, bv_manager_t &bv_manager,
                                 const DistanceOracleT &d_oracle, const WaverSwapRefinementConfiguration &config,
                                 const f64 alpha_end, const f64 beta_end, const f64 alpha_link, const f64 beta_link) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverGroupSwapRefinement", "refine_groups");
            return GroupSwapEngine<DistanceOracleT>(g, p_manager, bv_manager, d_oracle, config,
                                                    alpha_end, beta_end, alpha_link, beta_link)
                    .run(false, config.max_group_swaps);
        }

        // Exchanges of a vertex of each of two PEs, at most config.max_vertex_swaps. Returns their number.
        template<typename DistanceOracleT>
        static u64 refine_vertices(graph_t &g, p_manager_t &p_manager, bv_manager_t &bv_manager,
                                   const DistanceOracleT &d_oracle, const WaverSwapRefinementConfiguration &config,
                                   const f64 alpha_end, const f64 beta_end, const f64 alpha_link, const f64 beta_link) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverGroupSwapRefinement", "refine_vertices");
            return GroupSwapEngine<DistanceOracleT>(g, p_manager, bv_manager, d_oracle, config,
                                                    alpha_end, beta_end, alpha_link, beta_link)
                    .run(true, config.max_vertex_swaps);
        }
    };
} // namespace HeiProMap

#endif //HEIPROMAP_WAVER_GROUP_SWAP_REFINEMENT_H
