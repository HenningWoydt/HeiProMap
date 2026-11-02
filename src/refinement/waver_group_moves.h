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

// Group moves (Waver). The engine shared by the refinements that move sender groups: the sender group S_AB is the
// set of vertices of block A with an edge to block B, and it makes up the single coalesced message A -> B. The engine
// keeps the quotient rows of the partition (WaverMessageRefinement::QuotientRows) and offers what such a refinement
// needs: collecting and measuring a group, the cost change of moving it to one or many blocks, applying a move,
// undoing it, a chain of moves that is kept or undone as a whole (finish), and the ejection of a group out of an
// overloaded block (eject). The chain, scatter and split refinements, and the busy, cycle and crowd refinements, derive
// their engines from it.

#ifndef HEIPROMAP_WAVER_GROUP_MOVES_H
#define HEIPROMAP_WAVER_GROUP_MOVES_H

#include <algorithm>
#include <utility>
#include <vector>

#include "../definitions.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "waver_message_refinement.h"

namespace HeiProMap {
    class WaverGroupMoves {
    public:
        using graph_t = CSRGraph;
        using p_manager_t = PartitionManager;
        using bv_manager_t = BoundaryVertexManager;
        using Quotient = WaverMessageRefinement::QuotientRows;
        using Ext = std::vector<std::pair<partition_t, weight_t> >;

        struct Move {
            std::vector<vertex_t> group;
            partition_t from = NO_ID;
            partition_t to = NO_ID;
        };

        // A group of one block together with its edge weights to the other blocks.
        struct Candidate {
            std::vector<vertex_t> group;
            Ext ext;
            weight_t weight_to_own = 0;
            weight_t weight = 0;
        };

        template<typename DistanceOracleT>
        class Engine {
        public:
            Engine(graph_t &t_g, p_manager_t &t_pm, bv_manager_t &t_bv, const DistanceOracleT &t_oracle,
                   const f64 t_alpha_end, const f64 t_beta_end, const f64 t_alpha_link, const f64 t_beta_link)
                : g(t_g), pm(t_pm), bv(t_bv), oracle(t_oracle),
                  alpha_end(t_alpha_end), beta_end(t_beta_end), alpha_link(t_alpha_link), beta_link(t_beta_link) {
                k = pm.get_k();
                quotient.initialize(g, pm, bv, k);
                in_group.assign(g.n, 0);
                locked.assign(g.n, 0);
                tmp_weight.assign(k, 0);
                tmp_seen.assign(k, 0);
                gw_weight.assign(k, 0);
                gw_stamp.assign(k, 0);
                cache_weight.assign(k, 0);
                slot_of.assign(k, -1);
            }

        protected:
            graph_t &g;
            p_manager_t &pm;
            bv_manager_t &bv;
            const DistanceOracleT &oracle;
            f64 alpha_end, beta_end, alpha_link, beta_link;
            partition_t k = 0;

            Quotient quotient;
            std::vector<char> in_group; // scratch: vertex is in the group being measured
            std::vector<char> locked;   // vertex was moved by the chain being built
            std::vector<weight_t> tmp_weight;
            std::vector<char> tmp_seen;
            std::vector<partition_t> dirty;
            std::vector<weight_t> gw_weight;   // scratch of group_weights(): block -> weight of the group
            std::vector<u32> gw_stamp;         // scratch of group_weights(): last vertex counted per block
            u32 gw_clock = 0;
            // The sender-group weights of one block, valid until a chain is kept (see weight_for()).
            std::vector<weight_t> cache_weight;
            std::vector<partition_t> cache_blocks;
            partition_t cache_block = NO_ID;
            u64 cache_epoch = 0, epoch = 1;
            Ext ext_sorted;                // scratch of deltas()
            std::vector<char> vanishes;    // scratch of deltas()
            std::vector<partition_t> ext_x, ext_y; // scratch of deltas()
            std::vector<f64> target_deltas; // scratch of deltas()
            std::vector<int32_t> slot_of; // scratch: block -> index in the targets of a split, -1 if none

            // The blocks with the most messages first: they decide the busiest endpoint.
            std::vector<partition_t> blocks_by_degree() {
                std::vector<partition_t> order(k);
                std::vector<size_t> degree(k);
                for (partition_t a = 0; a < k; ++a) {
                    order[a] = a;
                    degree[a] = quotient.neighbors(a).size();
                }
                std::stable_sort(order.begin(), order.end(),
                                 [&](const partition_t x, const partition_t y) { return degree[x] > degree[y]; });
                return order;
            }

            // The vertices of block x that have an edge to block e and have not been moved by the current chain.
            void collect_group(const partition_t x, const partition_t e, std::vector<vertex_t> &out) const {
                out.clear();
                for (size_t i = 0; i < bv.size(x); ++i) {
                    const vertex_t u = bv.get(x, i);
                    if (locked[u]) { continue; }
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        if (pm[g.edges_v[j]] == e) {
                            out.push_back(u);
                            break;
                        }
                    }
                }
            }

            // The weight of the sender group of x towards every neighbor block e (the vertices of x with an edge to e, not
            // moved by the current chain), for all e at once: one pass over the boundary of x. Returns the blocks e.
            void group_weights(const partition_t x, std::vector<partition_t> &blocks) {
                blocks.clear();
                for (size_t i = 0; i < bv.size(x); ++i) {
                    const vertex_t u = bv.get(x, i);
                    if (locked[u]) { continue; }
                    ++gw_clock;
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const partition_t e = pm[g.edges_v[j]];
                        if (e == x || gw_stamp[e] == gw_clock) { continue; }
                        gw_stamp[e] = gw_clock;
                        if (gw_weight[e] == 0) { blocks.push_back(e); }
                        gw_weight[e] += g.v_weight(u);
                    }
                }
            }

            // The weight of the sender group of x towards e, from a table of all groups of x that is kept until a chain is
            // committed (the moves of a chain that is undone leave the partition as it was). Saves collecting a group
            // for the pairs that cannot be moved anyway.
            weight_t weight_for(const partition_t x, const partition_t e) {
                if (cache_block != x || cache_epoch != epoch) {
                    for (const partition_t b: cache_blocks) { cache_weight[b] = 0; }
                    group_weights(x, cache_blocks);
                    for (const partition_t b: cache_blocks) {
                        cache_weight[b] = gw_weight[b];
                        gw_weight[b] = 0;
                    }
                    cache_block = x;
                    cache_epoch = epoch;
                }
                return cache_weight[e];
            }

            // The weight of the group, its edge weight to every other block, and to the vertices of x outside it.
            void measure(const partition_t x, Candidate &c) {
                c.weight = 0;
                c.weight_to_own = 0;
                for (const vertex_t u: c.group) {
                    in_group[u] = 1;
                    c.weight += g.v_weight(u);
                }
                for (const vertex_t u: c.group) {
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const vertex_t v = g.edges_v[j];
                        const partition_t d = pm[v];
                        const weight_t w = g.edge_w(j);
                        if (d == x) {
                            if (!in_group[v]) { c.weight_to_own += w; }
                            continue;
                        }
                        if (!tmp_seen[d]) {
                            tmp_seen[d] = 1;
                            dirty.push_back(d);
                        }
                        tmp_weight[d] += w;
                    }
                }
                c.ext.clear();
                for (const partition_t d: dirty) {
                    c.ext.emplace_back(d, tmp_weight[d]);
                    tmp_weight[d] = 0;
                    tmp_seen[d] = 0;
                }
                dirty.clear();
                for (const vertex_t u: c.group) { in_group[u] = 0; }
            }

            f64 delta(const partition_t from, const partition_t to, const Candidate &c) {
                return WaverMessageRefinement::compute_group_delta(from, to, c.ext, c.weight_to_own, quotient, oracle,
                                                                   alpha_end, beta_end, alpha_link, beta_link);
            }

            // The cost change of moving a group from `from` to every block in `targets` (the same values as delta()).
            // What does not depend on the target is computed once, and the quotient row of a target is joined with the
            // external blocks of the group in one pass instead of one search per block.
            void deltas(const partition_t from, const std::vector<partition_t> &targets, const Candidate &c,
                        std::vector<f64> &out) {
                ext_sorted = c.ext;
                std::sort(ext_sorted.begin(), ext_sorted.end());
                const size_t ne = ext_sorted.size();
                // grid coordinates, so that a distance is two subtractions instead of four divisions
                const partition_t W = oracle.current_W();
                const auto dist = [W](const partition_t xp, const partition_t yp, const partition_t xq, const partition_t yq) {
                    return static_cast<f64>((xp > xq ? xp - xq : xq - xp) + (yp > yq ? yp - yq : yq - yp));
                };
                const partition_t from_x = from % W, from_y = from / W;
                ext_x.resize(ne);
                ext_y.resize(ne);
                f64 removed = 0.0, beta_from = 0.0;
                vanishes.assign(ne, 0);
                for (size_t i = 0; i < ne; ++i) {
                    const auto &[d, w] = ext_sorted[i];
                    ext_x[i] = d % W;
                    ext_y[i] = d / W;
                    const f64 dist_a = dist(from_x, from_y, ext_x[i], ext_y[i]);
                    const weight_t w_a_d = quotient.weight(from, d);
                    if (w_a_d > 0 && w_a_d == w) {
                        vanishes[i] = 1;
                        removed += alpha_end + alpha_link * dist_a;
                    }
                    beta_from += beta_link * static_cast<f64>(w) * dist_a;
                }
                out.resize(targets.size());
                for (size_t t = 0; t < targets.size(); ++t) {
                    const partition_t to = targets[t];
                    const partition_t to_x = to % W, to_y = to / W;
                    const auto &row = quotient.neighbors(to);
                    f64 delta = removed + beta_from;
                    bool pair_removed = false;
                    size_t j = 0;
                    for (size_t i = 0; i < ne; ++i) {
                        const auto &[d, w] = ext_sorted[i];
                        if (d == to) {
                            pair_removed = vanishes[i] != 0;
                            delta += beta_end * static_cast<f64>(w);
                            continue;
                        }
                        while (j < row.size() && row[j].target < d) { ++j; }
                        const f64 dist_c = dist(to_x, to_y, ext_x[i], ext_y[i]);
                        if (j >= row.size() || row[j].target != d) { delta -= alpha_end + alpha_link * dist_c; }
                        delta -= beta_link * static_cast<f64>(w) * dist_c;
                    }
                    if (c.weight_to_own > 0) {
                        const f64 w = static_cast<f64>(c.weight_to_own);
                        const f64 dist_ca = dist(to_x, to_y, from_x, from_y);
                        if (quotient.weight(to, from) == 0 || pair_removed) { delta -= alpha_end + alpha_link * dist_ca; }
                        delta -= beta_end * w + beta_link * w * dist_ca;
                    }
                    out[t] = delta;
                }
            }

            void apply(const Candidate &c, const partition_t from, const partition_t to) {
                for (const vertex_t u: c.group) {
                    WaverMessageRefinement::apply_move(g, bv, pm, u, g.v_weight(u), from, to);
                }
                quotient.apply_group_move(from, to, c.ext, c.weight_to_own);
            }

            // The part of a group of x that goes to the block to: the vertices that have more edge weight to `to` than to
            // x first, until at least `need` weight is taken, and never more than `room`. False if that is not possible.
            bool partial_group(const partition_t x, const partition_t to, const std::vector<vertex_t> &group,
                               const weight_t need, const weight_t room, Candidate &out) {
                std::vector<std::pair<f64, vertex_t> > scored;
                scored.reserve(group.size());
                for (const vertex_t u: group) {
                    f64 score = 0.0;
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const partition_t d = pm[g.edges_v[j]];
                        if (d == to) {
                            score += static_cast<f64>(g.edge_w(j));
                        } else if (d == x) {
                            score -= static_cast<f64>(g.edge_w(j));
                        }
                    }
                    scored.emplace_back(-score, u);
                }
                std::sort(scored.begin(), scored.end());
                out.group.clear();
                weight_t taken = 0;
                for (const auto &[score, u]: scored) {
                    const weight_t uw = g.v_weight(u);
                    if (taken + uw > room) { continue; }
                    out.group.push_back(u);
                    taken += uw;
                    if (taken >= need) { break; }
                }
                if (taken < need) { return false; }
                measure(x, out);
                return true;
            }

            // Moves a group back to the block it came from.
            void undo(Move &m) {
                Candidate c;
                c.group = std::move(m.group);
                measure(m.to, c);
                apply(c, m.to, m.from);
            }

            // The blocks a group of x with edges to e can go to: e and the quotient neighbors of e.
            void targets_of(const partition_t x, const partition_t e, std::vector<partition_t> &out) {
                out.clear();
                out.push_back(e);
                for (const auto &he: quotient.neighbors(e)) {
                    if (he.target != x) { out.push_back(he.target); }
                }
            }

            // Keeps the moves of a chain if it is valid and pays off, otherwise undoes them all.
            int finish(std::vector<Move> &chain, const bool resolved, const f64 total, const f64 min_gain = 1e-6) {
                for (const Move &m: chain) {
                    for (const vertex_t u: m.group) { locked[u] = 0; }
                }
                if (resolved && total > min_gain) { return 1; }
                for (size_t i = chain.size(); i-- > 0;) { undo(chain[i]); }
                return -1;
            }

            // The edge weight of every vertex of a group to every block in targets (row i of the result is vertex i).
            void compute_affinity(const std::vector<vertex_t> &group, const std::vector<partition_t> &targets,
                                  std::vector<f64> &affinity) {
                const size_t nt = targets.size();
                for (size_t i = 0; i < nt; ++i) { slot_of[targets[i]] = static_cast<int32_t>(i); }
                affinity.assign(group.size() * nt, 0.0);
                for (size_t i = 0; i < group.size(); ++i) {
                    const vertex_t u = group[i];
                    f64 *row = &affinity[i * nt];
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const int32_t s = slot_of[pm[g.edges_v[j]]];
                        if (s >= 0) { row[s] += static_cast<f64>(g.edge_w(j)); }
                    }
                }
                for (size_t i = 0; i < nt; ++i) { slot_of[targets[i]] = -1; }
            }

            // Moves a sender group out of the overloaded block cur, so that cur is within its limit again. The group
            // and its target are chosen by the cost (groups that fit the target are preferred); on the last step of a
            // chain only a move into a block with room is useful. Updates cur to the block that received the group.
            bool eject(const partition_t from, const bool last, std::vector<Move> &chain, f64 &total, partition_t &cur,
                       const bool partial = false) {
                const weight_t over = pm.bweights[from] - pm.get_lmax(from);
                std::vector<partition_t> neighbors;
                for (const auto &he: quotient.neighbors(from)) { neighbors.push_back(he.target); }

                Candidate cand, best;
                partition_t best_to = NO_ID;
                bool best_fits = false;
                f64 best_delta = 0.0;
                std::vector<partition_t> targets;
                std::vector<std::vector<vertex_t> > groups(neighbors.size());
                // The weight of every group first: only groups that are heavy enough, and (on the last step) have a
                // block with room to go to, are collected and priced.
                std::vector<partition_t> weighed;
                if (!partial) {
                    group_weights(from, weighed);
                    neighbors = weighed;
                    groups.resize(neighbors.size());
                }
                for (size_t idx = 0; idx < neighbors.size(); ++idx) {
                    const partition_t e = neighbors[idx];
                    if (!partial) {
                        const weight_t group_weight = gw_weight[e];
                        gw_weight[e] = 0;
                        if (group_weight < over || pm.bweights[from] - group_weight < pm.get_lmin(from)) { continue; }
                        if (last) {
                            bool room = pm.bweights[e] + group_weight <= pm.get_lmax(e);
                            for (const auto &hn: quotient.neighbors(e)) {
                                if (room) { break; }
                                room = hn.target != from && pm.bweights[hn.target] + group_weight <= pm.get_lmax(hn.target);
                            }
                            if (!room) { continue; }
                        }
                    }
                    collect_group(from, e, cand.group);
                    if (cand.group.empty()) { continue; }
                    if (partial) { groups[idx] = cand.group; }
                    measure(from, cand);
                    if (cand.weight < over) { continue; }
                    if (pm.bweights[from] - cand.weight < pm.get_lmin(from)) { continue; }
                    targets_of(from, e, targets);
                    if (last) {
                        // only a move into a block with room is useful on the last step
                        size_t kept = 0;
                        for (const partition_t to: targets) {
                            if (pm.bweights[to] + cand.weight <= pm.get_lmax(to)) { targets[kept++] = to; }
                        }
                        targets.resize(kept);
                        if (targets.empty()) { continue; }
                    }
                    deltas(from, targets, cand, target_deltas);
                    for (size_t t = 0; t < targets.size(); ++t) {
                        const partition_t to = targets[t];
                        const bool fits = pm.bweights[to] + cand.weight <= pm.get_lmax(to);
                        const f64 d = target_deltas[t];
                        const bool better = best_to == NO_ID || (fits && !best_fits) || (fits == best_fits && d > best_delta);
                        if (better) {
                            best = cand;
                            best_to = to;
                            best_fits = fits;
                            best_delta = d;
                        }
                    }
                }
                // No whole group fits a block with room: give up just the weight that is over, to a block with room.
                if (partial && !best_fits) {
                    for (size_t idx = 0; idx < neighbors.size(); ++idx) {
                        if (groups[idx].empty()) { continue; }
                        targets_of(from, neighbors[idx], targets);
                        for (const partition_t to: targets) {
                            const weight_t room = pm.get_lmax(to) - pm.bweights[to];
                            if (room < over) { continue; }
                            if (!partial_group(from, to, groups[idx], over, room, cand)) { continue; }
                            if (pm.bweights[from] - cand.weight < pm.get_lmin(from)) { continue; }
                            const f64 d = delta(from, to, cand);
                            if (!best_fits || d > best_delta) {
                                best = cand;
                                best_to = to;
                                best_fits = true;
                                best_delta = d;
                            }
                        }
                    }
                }
                if (best_to == NO_ID) { return false; }

                apply(best, from, best_to);
                for (const vertex_t u: best.group) { locked[u] = 1; }
                total += best_delta;
                chain.push_back(Move{best.group, from, best_to});
                cur = best_to;
                return true;
            }
        };
    };
} // namespace HeiProMap

#endif //HEIPROMAP_WAVER_GROUP_MOVES_H
