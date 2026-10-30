#ifndef HEIPROMAP_WAVER_CYCLE_REFINEMENT_H
#define HEIPROMAP_WAVER_CYCLE_REFINEMENT_H

#include <algorithm>
#include <iostream>
#include <vector>

#include "waver_merge_refinement.h"

namespace HeiProMap {
    // Hyperparameters of the cycle refinement.
    class WaverCycleRefinementConfiguration {
    public:
        explicit WaverCycleRefinementConfiguration(std::string t_name) : name(std::move(t_name)) {}

        std::string name;
        bool enabled = false;
        u64 max_length = 4;   // blocks on one cycle
        u64 max_rounds = 3;   // rounds over all blocks; stops earlier once a round commits no cycle
        u64 max_tries = 8;    // first moves tried per block and round
        f64 min_gain = 1e-6;  // a cycle must improve the cost by more than this
        bool verbose = false;
    };

    // Cycle refinement: moves sender groups around a cycle of blocks.
    //
    // With the blocks as full as they are, no single move has room to go to. Along a cycle every block gives a group
    // away and receives one from the previous block, so the weights stay (nearly) the same: block A1 gives a sender
    // group to A2, A2 gives one to A3, and so on until the last block gives one back to A1. The cycle is built
    // greedily from A1: each step moves the group of the current block that pays off most (a single step may cost
    // something) to a block that is not yet on the cycle, or back to A1, which closes it. Every block that gives a group
    // away has to stay within its limits after giving, and A1 within its limits after receiving. The cycle is kept if it
    // closes and the total cost went down, otherwise all its moves are undone.
    class WaverCycleRefinement {
        using graph_t = CSRGraph;
        using p_manager_t = PartitionManager;
        using bv_manager_t = BoundaryVertexManager;
        using Base = WaverMergeRefinement;

        template<typename DistanceOracleT>
        class CycleEngine : public Base::Engine<DistanceOracleT> {
            using E = Base::Engine<DistanceOracleT>;
            using E::g; using E::pm; using E::quotient; using E::locked;
            using E::collect_group; using E::measure; using E::apply; using E::delta; using E::targets_of;
            using E::finish;

        public:
            CycleEngine(graph_t &t_g, p_manager_t &t_pm, bv_manager_t &t_bv, const DistanceOracleT &t_oracle,
                        const WaverCycleRefinementConfiguration &t_cycle,
                        const f64 t_alpha_end, const f64 t_beta_end, const f64 t_alpha_link, const f64 t_beta_link)
                : E(t_g, t_pm, t_bv, t_oracle, unused_config(), t_alpha_end, t_beta_end, t_alpha_link, t_beta_link),
                  cycle(t_cycle) {}

            void run() {
                const partition_t k = pm.get_k();
                std::vector<partition_t> order(k);
                for (u64 round = 0; round < cycle.max_rounds; ++round) {
                    for (partition_t a = 0; a < k; ++a) { order[a] = a; }
                    std::vector<size_t> degree(k);
                    for (partition_t a = 0; a < k; ++a) { degree[a] = quotient.neighbors(a).size(); }
                    std::stable_sort(order.begin(), order.end(),
                                     [&](const partition_t x, const partition_t y) { return degree[x] > degree[y]; });
                    u64 tried = 0, kept = 0, length_sum = 0;
                    f64 gain = 0.0;
                    for (const partition_t a: order) {
                        std::vector<partition_t> firsts;
                        for (const auto &he: quotient.neighbors(a)) { firsts.push_back(he.target); }
                        if (firsts.size() > cycle.max_tries) { firsts.resize(cycle.max_tries); }
                        for (const partition_t e: firsts) {
                            if (quotient.weight(a, e) == 0) { continue; }
                            f64 total = 0.0;
                            u64 length = 0;
                            const int result = try_cycle(a, e, total, length);
                            if (result == 0) { continue; }
                            ++tried;
                            if (result > 0) {
                                ++kept;
                                gain += total;
                                length_sum += length;
                                break;
                            }
                        }
                    }
                    if (cycle.verbose) {
                        std::cout << "[cycle] round " << round << ": " << tried << " cycles tried, " << kept
                                  << " kept (mean length " << (kept ? static_cast<f64>(length_sum) / kept : 0.0)
                                  << "), gain " << gain << std::endl;
                    }
                    if (kept == 0) { break; }
                }
            }

        private:
            const WaverCycleRefinementConfiguration &cycle;

            static const WaverMergeRefinementConfiguration &unused_config() {
                static const WaverMergeRefinementConfiguration c("unused");
                return c;
            }

            // The best move of a group of `from` with edges to e: its target, delta and group. A move to `home`
            // (closing the cycle) needs room in home; a move to another block is only allowed if that block is not yet
            // on the cycle. The block `from` must stay within its limits.
            bool best_move(const partition_t from, const partition_t e, const partition_t home,
                           const std::vector<partition_t> &on_cycle, const bool allow_open, const bool allow_close,
                           Base::Candidate &out, partition_t &to, f64 &best_delta) {
                Base::Candidate group;
                collect_group(from, e, group.group);
                if (group.group.empty()) { return false; }
                measure(from, group);
                if (pm.bweights[from] - group.weight < pm.get_lmin(from)) { return false; }
                std::vector<partition_t> targets;
                targets_of(from, e, targets);
                bool found = false;
                for (const partition_t c: targets) {
                    if (c == home) {
                        if (!allow_close || pm.bweights[home] + group.weight > pm.get_lmax(home)) { continue; }
                    } else {
                        if (!allow_open || std::find(on_cycle.begin(), on_cycle.end(), c) != on_cycle.end()) { continue; }
                    }
                    const f64 d = delta(from, c, group);
                    if (!found || d > best_delta) {
                        found = true;
                        best_delta = d;
                        to = c;
                    }
                }
                if (found) { out = std::move(group); }
                return found;
            }

            // A cycle from a that starts with the sender group of a towards e. Returns 0 if there is no move to start
            // with, +1 if a cycle was kept, -1 if one was tried and undone.
            int try_cycle(const partition_t a, const partition_t e, f64 &total, u64 &length) {
                std::vector<Base::Move> path;
                std::vector<partition_t> on_cycle{a};
                Base::Candidate group;
                partition_t to = NO_ID;
                f64 d = 0.0;
                if (!best_move(a, e, a, on_cycle, true, false, group, to, d)) { return 0; }
                total = d;
                apply(group, a, to);
                for (const vertex_t u: group.group) { locked[u] = 1; }
                path.push_back(Base::Move{group.group, a, to});
                on_cycle.push_back(to);

                bool closed = false;
                partition_t cur = to;
                for (u64 step = 1; step < cycle.max_length && !closed; ++step) {
                    const bool last = step + 1 == cycle.max_length;
                    // the next move: the best over all groups of cur, closing the cycle if that is possible
                    bool found = false;
                    Base::Candidate best;
                    partition_t best_to = NO_ID;
                    f64 best_delta = 0.0;
                    bool best_closes = false;
                    std::vector<partition_t> nbrs;
                    for (const auto &he: quotient.neighbors(cur)) { nbrs.push_back(he.target); }
                    for (const partition_t n: nbrs) {
                        Base::Candidate cand;
                        partition_t cand_to = NO_ID;
                        f64 cand_delta = 0.0;
                        if (!best_move(cur, n, a, on_cycle, !last, true, cand, cand_to, cand_delta)) { continue; }
                        const bool closes = cand_to == a;
                        // closing moves win over open ones
                        if (!found || (closes && !best_closes) || (closes == best_closes && cand_delta > best_delta)) {
                            found = true;
                            best = std::move(cand);
                            best_to = cand_to;
                            best_delta = cand_delta;
                            best_closes = closes;
                        }
                    }
                    if (!found) { break; }
                    apply(best, cur, best_to);
                    for (const vertex_t u: best.group) { locked[u] = 1; }
                    total += best_delta;
                    path.push_back(Base::Move{best.group, cur, best_to});
                    if (best_closes) {
                        closed = true;
                    } else {
                        on_cycle.push_back(best_to);
                        cur = best_to;
                    }
                }
                length = path.size();
                const int result = finish(path, closed && total > cycle.min_gain, total);
                return result;
            }
        };

    public:
        template<typename DistanceOracleT>
        static void refine(graph_t &g, p_manager_t &p_manager, bv_manager_t &bv_manager,
                           const DistanceOracleT &d_oracle, const WaverCycleRefinementConfiguration &config,
                           const f64 alpha_end, const f64 beta_end, const f64 alpha_link, const f64 beta_link) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverCycleRefinement", "refine");
            CycleEngine<DistanceOracleT>(g, p_manager, bv_manager, d_oracle, config,
                                         alpha_end, beta_end, alpha_link, beta_link).run();
        }
    };
} // namespace HeiProMap

#endif //HEIPROMAP_WAVER_CYCLE_REFINEMENT_H
