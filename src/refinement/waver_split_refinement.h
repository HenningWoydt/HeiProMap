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

// Split refinement (Waver). Splits a sender group among several blocks. The sender group S_AB (the vertices of block A
// with an edge to block B) is divided among up to `parts` of the blocks B and the neighbors of B: every vertex goes to
// the block it has the most edge weight to, so that each part merges into a message that block already sends. The
// parts are priced one after the other, each against the quotient left by the previous ones, and every block that ends
// up over its weight limit passes load on like in a chain (see WaverChainRefinement). The move is kept if every block
// is within its limits again and the cost went down in total, otherwise all its moves are undone. Options:
// --group-split-parts, --group-split-depth, --group-split-rounds, --group-split-partial.

#ifndef HEIPROMAP_WAVER_SPLIT_REFINEMENT_H
#define HEIPROMAP_WAVER_SPLIT_REFINEMENT_H

#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../definitions.h"
#include "../utility/profiler.h"
#include "waver_group_moves.h"

namespace HeiProMap {
    // Hyperparameters of the split refinement.
    class WaverSplitRefinementConfiguration {
    public:
        explicit WaverSplitRefinementConfiguration(std::string t_name) : name(std::move(t_name)) {}

        std::string name;
        bool enabled = false;
        u64 parts = 3;            // a sender group is split among up to this many blocks (below 2: never)
        u64 max_depth = 3;        // moves that pass the load of overloaded blocks on
        u64 max_rounds = 3;       // rounds over all quotient edges; stops earlier once a round commits no move
        f64 min_round_gain = 0.0; // stop when a round gains less than this share of the gain of the first round (0: off)
        f64 min_gain = 1e-6;      // a split must improve the cost by more than this
        bool partial = false;     // an overloaded block may also give up only part of a sender group
        bool verbose = false;     // print the statistics of every round
    };

    class WaverSplitRefinement {
        using Base = WaverGroupMoves;
        using graph_t = Base::graph_t;
        using p_manager_t = Base::p_manager_t;
        using bv_manager_t = Base::bv_manager_t;
        using Candidate = Base::Candidate;
        using Move = Base::Move;

        template<typename DistanceOracleT>
        class SplitEngine : public Base::Engine<DistanceOracleT> {
            using E = Base::Engine<DistanceOracleT>;
            using E::g; using E::pm; using E::blocks_by_degree; using E::quotient; using E::k; using E::epoch; using E::locked; using E::slot_of;
            using E::collect_group; using E::measure; using E::targets_of; using E::delta; using E::apply;
            using E::eject; using E::finish;
            const WaverSplitRefinementConfiguration &cfg;
            u64 split_tried = 0, split_kept = 0;

        public:
            SplitEngine(graph_t &t_g, p_manager_t &t_pm, bv_manager_t &t_bv, const DistanceOracleT &t_oracle,
                        const WaverSplitRefinementConfiguration &t_cfg,
                        const f64 t_alpha_end, const f64 t_beta_end, const f64 t_alpha_link, const f64 t_beta_link)
                : E(t_g, t_pm, t_bv, t_oracle, t_alpha_end, t_beta_end, t_alpha_link, t_beta_link), cfg(t_cfg) {}

            void run() {
                f64 first_gain = 0.0;
                std::vector<partition_t> targets;
                for (u64 round = 0; round < cfg.max_rounds; ++round) {
                    HEIPROMAP_PROFILE_SCOPE("refinement", "WaverSplitRefinement", "round");
                    const std::vector<partition_t> order = blocks_by_degree();
                    u64 attempts = 0, committed = 0, rolled_back = 0;
                    split_tried = split_kept = 0;
                    f64 gain = 0.0;
                    for (const partition_t a: order) {
                        targets.clear();
                        for (const auto &he: quotient.neighbors(a)) { targets.push_back(he.target); }
                        for (const partition_t b: targets) {
                            if (quotient.weight(a, b) == 0) { continue; }
                            f64 split_gain = 0.0;
                            const int result = try_split(a, b, split_gain);
                            if (result == 0) { continue; }
                            ++attempts;
                            if (result > 0) {
                                ++committed;
                                ++epoch;
                                gain += split_gain;
                            } else {
                                ++rolled_back;
                            }
                        }
                    }
                    if (cfg.verbose) {
                        std::cout << "[split] round " << round << ": " << attempts << " tried, " << committed
                                  << " kept, " << rolled_back << " undone, gain " << gain << std::endl;
                    }
                    if (committed == 0) { break; }
                    if (round == 0) { first_gain = gain; }
                    if (cfg.min_round_gain > 0.0 && gain < cfg.min_round_gain * first_gain) { break; }
                }
            }

        private:
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

                // Keep the parts with the most affinity, at most `parts` of them.
                std::vector<size_t> by_affinity;
                for (size_t t = 0; t < nt; ++t) {
                    if (part_affinity[t] > 0.0) { by_affinity.push_back(t); }
                }
                if (by_affinity.size() < 2) { return 0; }
                std::sort(by_affinity.begin(), by_affinity.end(),
                          [&](const size_t x, const size_t y) { return part_affinity[x] > part_affinity[y]; });
                if (by_affinity.size() > cfg.parts) { by_affinity.resize(cfg.parts); }
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
                    if (!eject(x, ejections + 1 == cfg.max_depth, chain, total, received, cfg.partial)) {
                        resolved = false;
                        break;
                    }
                    ++ejections;
                    if (pm.bweights[received] > pm.get_lmax(received)) { pending.push_back(received); }
                }
                const int result = finish(chain, resolved, total, cfg.min_gain);
                if (result > 0) { ++split_kept; }
                return result;
            }
        };

    public:
        template<typename DistanceOracleT>
        static void refine(graph_t &g,
                           p_manager_t &p_manager,
                           bv_manager_t &bv_manager,
                           const DistanceOracleT &d_oracle,
                           const WaverSplitRefinementConfiguration &config,
                           const f64 alpha_end,
                           const f64 beta_end,
                           const f64 alpha_link,
                           const f64 beta_link) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverSplitRefinement", "refine");
            SplitEngine<DistanceOracleT>(g, p_manager, bv_manager, d_oracle, config,
                                         alpha_end, beta_end, alpha_link, beta_link).run();
        }
    };
} // namespace HeiProMap

#endif //HEIPROMAP_WAVER_SPLIT_REFINEMENT_H
