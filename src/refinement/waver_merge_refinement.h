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

// Merge refinement (Waver). Combines messages by moving a sender group into a block that has no room for it. The
// message refinement only moves a group into a block with room left, and the blocks are nearly always full. Here the move is
// made anyway and the block may exceed its upper limit for a moment; it then passes load on by moving one of its own sender
// groups to another block, and so on, along a chain. The chain is kept if every block is within its limits again and the
// message cost went down in total, otherwise all its moves are undone. Variants: scatter (spread the group over the receiver
// and its neighbors within their room), partial groups, and splitting the group among several blocks. This is also the
// engine (quotient rows, group measuring, delta, undo) the busy, cycle and crowd refinements are built on.


#ifndef HEIPROMAP_WAVER_MERGE_REFINEMENT_H
#define HEIPROMAP_WAVER_MERGE_REFINEMENT_H

#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../definitions.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../utility/profiler.h"
#include "waver_message_refinement.h"

namespace HeiProMap {
    // Hyperparameters of the merge refinement.
    class WaverMergeRefinementConfiguration {
    public:
        explicit WaverMergeRefinementConfiguration(std::string t_name) : name(std::move(t_name)) {}

        std::string name;
        bool enabled = false;
        u64 max_depth = 3;       // moves in one chain
        u64 max_rounds = 3;      // rounds over all quotient edges; stops earlier once a round commits no chain
        f64 min_gain = 1e-6;     // a chain must improve the cost by more than this
        bool partial = false;    // an overloaded block may also give up only part of a sender group
        u64 split_parts = 0;     // a sender group may also be split among up to this many blocks (below 2: never)
        bool scatter = false;    // no chains: the sender group is spread over B and all its neighbors, within their room
        bool verbose = false;    // print the statistics of every round
    };

    // Merge refinement: combines messages by moving a sender group into a block that has no room for it.
    //
    // The message refinement moves the sender group S_AB (the vertices of A with an edge to B) to B or to a neighbor C
    // of B when that saves messages, but only into a block with room left, and with the blocks as full as they are
    // there is hardly ever one. Here the move is made anyway, and C may exceed its upper weight limit for a moment. C
    // then has to give up weight: one of its own sender groups S_CE moves to E or to a neighbor of E, which can merge
    // messages again and may overload that block in turn. The chain goes on until a block has room for what it
    // receives. It is kept if every block is within its limits again and the cost went down in total (the single
    // moves may cost something), otherwise all its moves are undone.
    //
    // With `scatter`, no chain is built and no block is ever overloaded. The sender group S_AB is spread over B and all
    // the quotient neighbors of B at once: every vertex goes to the block it has the most edge weight to among those
    // that still have room for it, so each part merges into a message that block already sends. The parts are priced
    // one after the other, each against the quotient left by the previous ones, and the move is kept if the total
    // cost went down. It is dropped when a vertex finds no block with room.
    //
    // With `split_parts`, a chain can also start by splitting the sender group S_AB among up to that many blocks (B and
    // the neighbors of B): every vertex goes to the block it has the most edge weight to, so that each part merges into
    // a message that block already sends. The parts are priced one after the other, each against the quotient left by
    // the previous ones, and every block that ends up over its limit passes load on like in a plain chain.
    //
    // A group that is too heavy for any block with room cannot be passed on, and the chain ends. With `partial`, an
    // overloaded block then gives up just enough weight instead: the vertices of the group that are tied most closely
    // to the receiving block, rather than to the block they leave.
    class WaverBusyRefinement;
    class WaverCycleRefinement;
    class WaverCrowdRefinement;

    class WaverMergeRefinement {
        friend class WaverBusyRefinement; // derive their engines from the merge engine
        friend class WaverCycleRefinement;
        friend class WaverCrowdRefinement;
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
                   const WaverMergeRefinementConfiguration &t_cfg,
                   const f64 t_alpha_end, const f64 t_beta_end, const f64 t_alpha_link, const f64 t_beta_link)
                : g(t_g), pm(t_pm), bv(t_bv), oracle(t_oracle), cfg(t_cfg),
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

            void run() {
                std::vector<partition_t> order(k);
                std::vector<partition_t> targets;
                for (u64 round = 0; round < cfg.max_rounds; ++round) {
                    HEIPROMAP_PROFILE_SCOPE("refinement", "WaverMergeRefinement", "round");
                    // The blocks with the most messages first: they decide the busiest endpoint.
                    for (partition_t a = 0; a < k; ++a) { order[a] = a; }
                    std::vector<size_t> degree(k);
                    for (partition_t a = 0; a < k; ++a) { degree[a] = quotient.neighbors(a).size(); }
                    std::stable_sort(order.begin(), order.end(),
                                     [&](const partition_t x, const partition_t y) { return degree[x] > degree[y]; });

                    u64 attempts = 0, committed = 0, rolled_back = 0;
                    split_tried = split_kept = 0;
                    f64 gain = 0.0;
                    for (const partition_t a: order) {
                        targets.clear();
                        for (const auto &he: quotient.neighbors(a)) { targets.push_back(he.target); }
                        for (const partition_t b: targets) {
                            if (quotient.weight(a, b) == 0) { continue; }
                            f64 chain_gain = 0.0;
                            if (cfg.scatter) {
                                const int result = try_scatter(a, b, chain_gain);
                                if (result == 0) { continue; }
                                ++attempts;
                                if (result > 0) {
                                    ++committed;
                                    ++epoch;
                                    gain += chain_gain;
                                } else {
                                    ++rolled_back;
                                }
                                continue;
                            }
                            int result = try_chain(a, b, chain_gain);
                            if (result <= 0 && cfg.split_parts >= 2) {
                                f64 split_gain = 0.0;
                                const int split = try_split(a, b, split_gain);
                                if (split > 0) {
                                    result = 1;
                                    chain_gain = split_gain;
                                } else if (split < 0 && result == 0) {
                                    result = -1;
                                }
                            }
                            if (result == 0) { continue; }
                            ++attempts;
                            if (result > 0) {
                                ++committed;
                                ++epoch;
                                gain += chain_gain;
                            } else {
                                ++rolled_back;
                            }
                        }
                    }
                    if (cfg.verbose) {
                        std::cout << "[merge] round " << round << ": " << attempts << " chains tried, " << committed
                                  << " kept, " << rolled_back << " undone, gain " << gain << "; splits: " << split_tried
                                  << " tried, " << split_kept << " kept" << std::endl;
                    }
                    if (committed == 0) { break; }
                }
            }

        protected:
            graph_t &g;
            p_manager_t &pm;
            bv_manager_t &bv;
            const DistanceOracleT &oracle;
            const WaverMergeRefinementConfiguration &cfg;
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
            u64 split_tried = 0, split_kept = 0;

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

            // False if a chain that starts by moving a group of the given weight from a to c cannot end within the
            // limits, which is found without making any move. If c has room for the group the chain is done at once.
            // Otherwise c has to give up weight `over`: with one move left that is only possible into a block with at
            // least that much room (a itself has room for the group it just lost), found among the blocks that a group
            // of c can go to. With more moves left it is not decided here.
            bool can_resolve(const partition_t a, const partition_t c, const weight_t group_weight) {
                const weight_t over = pm.bweights[c] + group_weight - pm.get_lmax(c);
                if (over <= 0) { return true; }
                if (cfg.max_depth < 2) { return false; }
                if (cfg.max_depth > 2) { return true; }
                const auto room = [&](const partition_t t) {
                    const weight_t w = t == a ? pm.bweights[a] - group_weight : pm.bweights[t];
                    return pm.get_lmax(t) - w;
                };
                for (const auto &he: quotient.neighbors(c)) {
                    const partition_t e = he.target;
                    if (e != c && room(e) >= over) { return true; }
                    for (const auto &hn: quotient.neighbors(e)) {
                        if (hn.target != c && room(hn.target) >= over) { return true; }
                    }
                }
                return false;
            }

            // One chain that starts with the sender group of a towards b. Returns 0 if there is no move to start
            // with, +1 if a chain was kept (its gain is in total), -1 if one was tried and undone.
            int try_chain(const partition_t a, const partition_t b, f64 &total) {
                const weight_t group_weight = weight_for(a, b);
                if (group_weight == 0 || pm.bweights[a] - group_weight < pm.get_lmin(a)) { return 0; }
                Candidate first;
                collect_group(a, b, first.group);
                if (first.group.empty()) { return 0; }
                measure(a, first);

                std::vector<partition_t> targets;
                targets_of(a, b, targets);
                partition_t best_c = NO_ID;
                f64 best_delta = cfg.min_gain;
                deltas(a, targets, first, target_deltas);
                for (size_t t = 0; t < targets.size(); ++t) {
                    if (target_deltas[t] > best_delta) {
                        best_delta = target_deltas[t];
                        best_c = targets[t];
                    }
                }
                if (best_c == NO_ID) { return 0; }
                if (!can_resolve(a, best_c, first.weight)) { return 0; }

                std::vector<Move> chain;
                total = best_delta;
                apply(first, a, best_c);
                for (const vertex_t u: first.group) { locked[u] = 1; }
                chain.push_back(Move{first.group, a, best_c});

                partition_t cur = best_c;
                bool resolved = pm.bweights[cur] <= pm.get_lmax(cur);
                for (u64 depth = 1; !resolved && depth < cfg.max_depth; ++depth) {
                    const bool last = depth + 1 == cfg.max_depth;
                    if (!eject(cur, last, chain, total, cur)) { break; }
                    resolved = pm.bweights[cur] <= pm.get_lmax(cur);
                }

                return finish(chain, resolved, total);
            }

            // Keeps the moves of a chain if it is valid and pays off, otherwise undoes them all.
            int finish(std::vector<Move> &chain, const bool resolved, const f64 total) {
                for (const Move &m: chain) {
                    for (const vertex_t u: m.group) { locked[u] = 0; }
                }
                if (resolved && total > cfg.min_gain) { return 1; }
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

            // The sender group of a towards b, spread over b and all its neighbors within the room they have. Returns
            // like try_chain, without any overload: 0 if there is nothing to try, +1 if the move was kept, -1 if it
            // was made and undone.
            int try_scatter(const partition_t a, const partition_t b, f64 &total) {
                const weight_t group_weight = weight_for(a, b);
                if (group_weight == 0 || pm.bweights[a] - group_weight < pm.get_lmin(a)) { return 0; }
                std::vector<partition_t> targets;
                targets_of(a, b, targets);
                const size_t nt = targets.size();
                std::vector<weight_t> room(nt);
                weight_t total_room = 0;
                for (size_t t = 0; t < nt; ++t) {
                    room[t] = std::max<weight_t>(0, pm.get_lmax(targets[t]) - pm.bweights[targets[t]]);
                    total_room += room[t];
                }
                if (total_room < group_weight) { return 0; }
                Candidate whole;
                collect_group(a, b, whole.group);
                if (whole.group.empty()) { return 0; }
                measure(a, whole);
                const size_t nv = whole.group.size();

                std::vector<f64> affinity;
                compute_affinity(whole.group, targets, affinity);

                // The vertices with the strongest tie to some block choose first.
                std::vector<f64> best_affinity(nv, 0.0);
                for (size_t i = 0; i < nv; ++i) {
                    for (size_t t = 0; t < nt; ++t) { best_affinity[i] = std::max(best_affinity[i], affinity[i * nt + t]); }
                }
                std::vector<size_t> order(nv);
                for (size_t i = 0; i < nv; ++i) { order[i] = i; }
                std::stable_sort(order.begin(), order.end(),
                                 [&](const size_t x, const size_t y) { return best_affinity[x] > best_affinity[y]; });

                std::vector<std::vector<vertex_t> > parts(nt);
                for (const size_t i: order) {
                    const weight_t vw = g.v_weight(whole.group[i]);
                    size_t best = nt;
                    for (size_t t = 0; t < nt; ++t) {
                        if (room[t] < vw || affinity[i * nt + t] <= 0.0) { continue; }
                        if (best == nt || affinity[i * nt + t] > affinity[i * nt + best]) { best = t; }
                    }
                    if (best == nt) { return 0; } // this vertex has no block with room that it is tied to
                    room[best] -= vw;
                    parts[best].push_back(whole.group[i]);
                }

                ++split_tried;
                std::vector<Move> chain;
                total = 0.0;
                for (size_t t = 0; t < nt; ++t) {
                    if (parts[t].empty()) { continue; }
                    Candidate part;
                    part.group = std::move(parts[t]);
                    measure(a, part);
                    total += delta(a, targets[t], part);
                    apply(part, a, targets[t]);
                    for (const vertex_t u: part.group) { locked[u] = 1; }
                    chain.push_back(Move{part.group, a, targets[t]});
                }
                const int result = finish(chain, true, total);
                if (result > 0) { ++split_kept; }
                return result;
            }

            // A chain that starts by splitting the sender group of a towards b among several blocks. Returns like
            // try_chain.
            int try_split(const partition_t a, const partition_t b, f64 &total) {
                Candidate whole;
                collect_group(a, b, whole.group);
                if (whole.group.size() < 2) { return 0; }
                measure(a, whole);
                if (pm.bweights[a] - whole.weight < pm.get_lmin(a)) { return 0; }
                std::vector<partition_t> targets;
                targets_of(a, b, targets);
                if (targets.size() < 2) { return 0; }

                // The edge weight of every vertex to every candidate block, and its best block.
                const size_t nt = targets.size(), nv = whole.group.size();
                for (size_t i = 0; i < nt; ++i) { slot_of[targets[i]] = static_cast<int32_t>(i); }
                std::vector<f64> affinity(nv * nt, 0.0);
                std::vector<size_t> choice(nv, 0);
                std::vector<f64> part_affinity(nt, 0.0);
                for (size_t i = 0; i < nv; ++i) {
                    const vertex_t u = whole.group[i];
                    f64 *row = &affinity[i * nt];
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const int32_t s = slot_of[pm[g.edges_v[j]]];
                        if (s >= 0) { row[s] += static_cast<f64>(g.edge_w(j)); }
                    }
                    for (size_t t = 1; t < nt; ++t) {
                        if (row[t] > row[choice[i]]) { choice[i] = t; }
                    }
                    part_affinity[choice[i]] += row[choice[i]];
                }
                for (size_t i = 0; i < nt; ++i) { slot_of[targets[i]] = -1; }

                // Keep the parts with the most affinity, at most split_parts of them.
                std::vector<size_t> by_affinity;
                for (size_t t = 0; t < nt; ++t) {
                    if (part_affinity[t] > 0.0) { by_affinity.push_back(t); }
                }
                if (by_affinity.size() < 2) { return 0; }
                std::sort(by_affinity.begin(), by_affinity.end(),
                          [&](const size_t x, const size_t y) { return part_affinity[x] > part_affinity[y]; });
                if (by_affinity.size() > cfg.split_parts) { by_affinity.resize(cfg.split_parts); }
                std::vector<std::vector<vertex_t> > parts(nt);
                for (size_t i = 0; i < nv; ++i) {
                    size_t t = choice[i];
                    if (std::find(by_affinity.begin(), by_affinity.end(), t) == by_affinity.end()) {
                        t = by_affinity[0];
                        for (const size_t x: by_affinity) {
                            if (affinity[i * nt + x] > affinity[i * nt + t]) { t = x; }
                        }
                    }
                    parts[t].push_back(whole.group[i]);
                }
                if (std::count_if(parts.begin(), parts.end(), [](const auto &p) { return !p.empty(); }) < 2) { return 0; }

                ++split_tried;
                std::vector<Move> chain;
                std::vector<partition_t> pending;
                total = 0.0;
                for (const size_t t: by_affinity) {
                    if (parts[t].empty()) { continue; }
                    Candidate part;
                    part.group = std::move(parts[t]);
                    measure(a, part);
                    total += delta(a, targets[t], part);
                    apply(part, a, targets[t]);
                    for (const vertex_t u: part.group) { locked[u] = 1; }
                    chain.push_back(Move{part.group, a, targets[t]});
                    if (pm.bweights[targets[t]] > pm.get_lmax(targets[t])) { pending.push_back(targets[t]); }
                }

                // Every overloaded receiver passes load on, one move each, as in a plain chain.
                bool resolved = true;
                u64 ejections = 0;
                while (!pending.empty()) {
                    const partition_t x = pending.back();
                    pending.pop_back();
                    if (pm.bweights[x] <= pm.get_lmax(x)) { continue; }
                    if (ejections >= cfg.max_depth) {
                        resolved = false;
                        break;
                    }
                    partition_t received = NO_ID;
                    if (!eject(x, ejections + 1 == cfg.max_depth, chain, total, received)) {
                        resolved = false;
                        break;
                    }
                    ++ejections;
                    if (pm.bweights[received] > pm.get_lmax(received)) { pending.push_back(received); }
                }
                const int result = finish(chain, resolved, total);
                if (result > 0) { ++split_kept; }
                return result;
            }

            // Moves a sender group out of the overloaded block cur, so that cur is within its limit again. The group
            // and its target are chosen by the cost (groups that fit the target are preferred); on the last step of a
            // chain only a move into a block with room is useful. Updates cur to the block that received the group.
            bool eject(const partition_t from, const bool last, std::vector<Move> &chain, f64 &total, partition_t &cur) {
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
                if (!cfg.partial) {
                    group_weights(from, weighed);
                    neighbors = weighed;
                    groups.resize(neighbors.size());
                }
                for (size_t idx = 0; idx < neighbors.size(); ++idx) {
                    const partition_t e = neighbors[idx];
                    if (!cfg.partial) {
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
                    if (cfg.partial) { groups[idx] = cand.group; }
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
                if (cfg.partial && !best_fits) {
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

    public:
        template<typename DistanceOracleT>
        static void refine(graph_t &g,
                           p_manager_t &p_manager,
                           bv_manager_t &bv_manager,
                           const DistanceOracleT &d_oracle,
                           const WaverMergeRefinementConfiguration &config,
                           const f64 alpha_end,
                           const f64 beta_end,
                           const f64 alpha_link,
                           const f64 beta_link) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverMergeRefinement", "refine");
            Engine<DistanceOracleT>(g, p_manager, bv_manager, d_oracle, config,
                                    alpha_end, beta_end, alpha_link, beta_link).run();
        }
    };
} // namespace HeiProMap

#endif //HEIPROMAP_WAVER_MERGE_REFINEMENT_H
