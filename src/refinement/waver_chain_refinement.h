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

// Chain refinement (Waver). Combines messages by moving a sender group into a block that has no room for it. The
// message refinement only moves a group into a block with room left, and the blocks are nearly always full. Here the move is
// made anyway and the block may exceed its upper limit for a moment; it then passes load on by moving one of its own sender
// groups to another block, and so on, along a chain. The chain is kept if every block is within its limits again and the
// message cost went down in total, otherwise all its moves are undone. With `partial`, an overloaded block that has no
// whole group to give away gives up just enough weight. Options: --chain-refinement, --chain-depth, --chain-rounds,
// --chain-round-gain, --chain-partial.

#ifndef HEIPROMAP_WAVER_CHAIN_REFINEMENT_H
#define HEIPROMAP_WAVER_CHAIN_REFINEMENT_H

#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../definitions.h"
#include "../utility/profiler.h"
#include "waver_group_moves.h"

namespace HeiProMap {
    // Hyperparameters of the chain refinement.
    class WaverChainRefinementConfiguration {
    public:
        explicit WaverChainRefinementConfiguration(std::string t_name) : name(std::move(t_name)) {}

        std::string name;
        bool enabled = false;
        u64 max_depth = 3;        // moves in one chain
        u64 max_rounds = 3;       // rounds over all quotient edges; stops earlier once a round commits no chain
        f64 min_round_gain = 0.0; // stop when a round gains less than this share of the gain of the first round (0: off)
        f64 min_gain = 1e-6;      // a chain must improve the cost by more than this
        bool partial = false;     // an overloaded block may also give up only part of a sender group
        bool verbose = false;     // print the statistics of every round
    };

    // Chain refinement: combines messages by moving a sender group into a block that has no room for it.
    //
    // The message refinement moves the sender group S_AB (the vertices of A with an edge to B) to B or to a neighbor C
    // of B when that saves messages, but only into a block with room left, and with the blocks as full as they are
    // there is hardly ever one. Here the move is made anyway, and C may exceed its upper weight limit for a moment. C
    // then has to give up weight: one of its own sender groups S_CE moves to E or to a neighbor of E, which can merge
    // messages again and may overload that block in turn. The chain goes on until a block has room for what it
    // receives. It is kept if every block is within its limits again and the cost went down in total (the single
    // moves may cost something), otherwise all its moves are undone.
    //
    // A group that is too heavy for any block with room cannot be passed on, and the chain ends. With `partial`, an
    // overloaded block then gives up just enough weight instead: the vertices of the group that are tied most closely
    // to the receiving block, rather than to the block they leave.
    class WaverChainRefinement {
        using Base = WaverGroupMoves;
        using graph_t = Base::graph_t;
        using p_manager_t = Base::p_manager_t;
        using bv_manager_t = Base::bv_manager_t;
        using Candidate = Base::Candidate;
        using Move = Base::Move;

        template<typename DistanceOracleT>
        class ChainEngine : public Base::Engine<DistanceOracleT> {
            using E = Base::Engine<DistanceOracleT>;
            using E::pm; using E::blocks_by_degree; using E::quotient; using E::k; using E::epoch; using E::locked; using E::target_deltas;
            using E::weight_for; using E::collect_group; using E::measure; using E::targets_of; using E::deltas;
            using E::apply; using E::eject; using E::finish; using E::alpha_end; using E::beta_end;
            const WaverChainRefinementConfiguration &cfg;

        public:
            ChainEngine(graph_t &t_g, p_manager_t &t_pm, bv_manager_t &t_bv, const DistanceOracleT &t_oracle,
                        const WaverChainRefinementConfiguration &t_cfg,
                        const f64 t_alpha_end, const f64 t_beta_end, const f64 t_alpha_link, const f64 t_beta_link)
                : E(t_g, t_pm, t_bv, t_oracle, t_alpha_end, t_beta_end, t_alpha_link, t_beta_link), cfg(t_cfg) {}

            void run() {
                f64 first_gain = 0.0;
                std::vector<partition_t> targets;
                for (u64 round = 0; round < cfg.max_rounds; ++round) {
                    HEIPROMAP_PROFILE_SCOPE("refinement", "WaverChainRefinement", "round");
                    const std::vector<partition_t> order = blocks_by_degree();
                    u64 attempts = 0, committed = 0, rolled_back = 0;
                    f64 gain = 0.0;
                    for (const partition_t a: order) {
                        targets.clear();
                        for (const auto &he: quotient.neighbors(a)) { targets.push_back(he.target); }
                        for (const partition_t b: targets) {
                            if (quotient.weight(a, b) == 0) { continue; }
                            f64 chain_gain = 0.0;
                            const int result = try_chain(a, b, chain_gain);
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
                        std::cout << "[chain] round " << round << ": " << attempts << " chains tried, " << committed
                                  << " kept, " << rolled_back << " undone, gain " << gain << std::endl;
                    }
                    if (committed == 0) { break; }
                    if (round == 0) { first_gain = gain; }
                    if (cfg.min_round_gain > 0.0 && gain < cfg.min_round_gain * first_gain) { break; }
                }
            }

        private:
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
                    if (!eject(cur, last, chain, total, cur, cfg.partial)) { break; }
                    resolved = pm.bweights[cur] <= pm.get_lmax(cur);
                }

                return finish(chain, resolved, total, cfg.min_gain);
            }
        };

    public:
        template<typename DistanceOracleT>
        static void refine(graph_t &g,
                           p_manager_t &p_manager,
                           bv_manager_t &bv_manager,
                           const DistanceOracleT &d_oracle,
                           const WaverChainRefinementConfiguration &config,
                           const f64 alpha_end,
                           const f64 beta_end,
                           const f64 alpha_link,
                           const f64 beta_link) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverChainRefinement", "refine");
            ChainEngine<DistanceOracleT>(g, p_manager, bv_manager, d_oracle, config,
                                         alpha_end, beta_end, alpha_link, beta_link).run();
        }
    };
} // namespace HeiProMap

#endif //HEIPROMAP_WAVER_CHAIN_REFINEMENT_H
