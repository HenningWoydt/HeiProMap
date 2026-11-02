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

// Scatter refinement (Waver). Spreads a sender group over the receiver and its neighbors. The sender group S_AB (the
// vertices of block A with an edge to block B) is not moved to one block but divided among B and all the quotient
// neighbors of B that still have room: every vertex goes to the block it has the most edge weight to, so each part
// merges into a message that block already sends. No block is ever overloaded. The parts are priced one after the
// other, each against the quotient left by the previous ones, and the move is kept if the total cost went down. It is
// dropped when a vertex finds no block with room. Options: --scatter-refinement, --scatter-rounds.

#ifndef HEIPROMAP_WAVER_SCATTER_REFINEMENT_H
#define HEIPROMAP_WAVER_SCATTER_REFINEMENT_H

#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../definitions.h"
#include "../utility/profiler.h"
#include "waver_group_moves.h"

namespace HeiProMap {
    // Hyperparameters of the scatter refinement.
    class WaverScatterRefinementConfiguration {
    public:
        explicit WaverScatterRefinementConfiguration(std::string t_name) : name(std::move(t_name)) {}

        std::string name;
        bool enabled = false;
        u64 max_rounds = 3;       // rounds over all quotient edges; stops earlier once a round commits no move
        f64 min_round_gain = 0.0; // stop when a round gains less than this share of the gain of the first round (0: off)
        f64 min_gain = 1e-6;      // a scatter must improve the cost by more than this
        bool verbose = false;     // print the statistics of every round
    };

    class WaverScatterRefinement {
        using Base = WaverGroupMoves;
        using graph_t = Base::graph_t;
        using p_manager_t = Base::p_manager_t;
        using bv_manager_t = Base::bv_manager_t;
        using Candidate = Base::Candidate;
        using Move = Base::Move;

        template<typename DistanceOracleT>
        class ScatterEngine : public Base::Engine<DistanceOracleT> {
            using E = Base::Engine<DistanceOracleT>;
            using E::g; using E::pm; using E::blocks_by_degree; using E::quotient; using E::k; using E::epoch; using E::locked;
            using E::weight_for; using E::collect_group; using E::measure; using E::targets_of;
            using E::compute_affinity; using E::delta; using E::apply; using E::finish;
            const WaverScatterRefinementConfiguration &cfg;

        public:
            ScatterEngine(graph_t &t_g, p_manager_t &t_pm, bv_manager_t &t_bv, const DistanceOracleT &t_oracle,
                          const WaverScatterRefinementConfiguration &t_cfg,
                          const f64 t_alpha_end, const f64 t_beta_end, const f64 t_alpha_link, const f64 t_beta_link)
                : E(t_g, t_pm, t_bv, t_oracle, t_alpha_end, t_beta_end, t_alpha_link, t_beta_link), cfg(t_cfg) {}

            void run() {
                f64 first_gain = 0.0;
                std::vector<partition_t> targets;
                for (u64 round = 0; round < cfg.max_rounds; ++round) {
                    HEIPROMAP_PROFILE_SCOPE("refinement", "WaverScatterRefinement", "round");
                    const std::vector<partition_t> order = blocks_by_degree();
                    u64 attempts = 0, committed = 0, rolled_back = 0;
                    f64 gain = 0.0;
                    for (const partition_t a: order) {
                        targets.clear();
                        for (const auto &he: quotient.neighbors(a)) { targets.push_back(he.target); }
                        for (const partition_t b: targets) {
                            if (quotient.weight(a, b) == 0) { continue; }
                            f64 scatter_gain = 0.0;
                            const int result = try_scatter(a, b, scatter_gain);
                            if (result == 0) { continue; }
                            ++attempts;
                            if (result > 0) {
                                ++committed;
                                ++epoch;
                                gain += scatter_gain;
                            } else {
                                ++rolled_back;
                            }
                        }
                    }
                    if (cfg.verbose) {
                        std::cout << "[scatter] round " << round << ": " << attempts << " tried, " << committed
                                  << " kept, " << rolled_back << " undone, gain " << gain << std::endl;
                    }
                    if (committed == 0) { break; }
                    if (round == 0) { first_gain = gain; }
                    if (cfg.min_round_gain > 0.0 && gain < cfg.min_round_gain * first_gain) { break; }
                }
            }

        private:
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
                return finish(chain, true, total, cfg.min_gain);
            }
        };

    public:
        template<typename DistanceOracleT>
        static void refine(graph_t &g,
                           p_manager_t &p_manager,
                           bv_manager_t &bv_manager,
                           const DistanceOracleT &d_oracle,
                           const WaverScatterRefinementConfiguration &config,
                           const f64 alpha_end,
                           const f64 beta_end,
                           const f64 alpha_link,
                           const f64 beta_link) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverScatterRefinement", "refine");
            ScatterEngine<DistanceOracleT>(g, p_manager, bv_manager, d_oracle, config,
                                           alpha_end, beta_end, alpha_link, beta_link).run();
        }
    };
} // namespace HeiProMap

#endif //HEIPROMAP_WAVER_SCATTER_REFINEMENT_H
