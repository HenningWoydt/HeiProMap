#ifndef HEIPROMAP_WAVER_BUSY_REFINEMENT_H
#define HEIPROMAP_WAVER_BUSY_REFINEMENT_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "waver_merge_refinement.h"
#include "waver_region_load.h"

namespace HeiProMap {
    // Hyperparameters of the busy-PE refinement.
    class WaverBusyRefinementConfiguration {
    public:
        explicit WaverBusyRefinementConfiguration(std::string t_name) : name(std::move(t_name)) {}

        std::string name;
        bool enabled = false;
        u64 max_rounds = 5;     // rounds; stops earlier once a round commits no move
        f64 top_fraction = 0.1; // the share of blocks, busiest first, that is worked on in a round
        f64 power = 4.0;        // the loads of the blocks touched by a move are compared as sum(load^power)
        f64 min_gain = 1e-6;    // relative improvement of that sum a move needs
        f64 alpha = 0.0;        // load per message; 0: the alpha_end of the cost model
        f64 beta = 0.0;         // load per unit of volume; 0: the beta_end of the cost model
        u64 max_tries = 8;      // sender groups tried per busy block and round
        bool exact = false;     // the load of a block is the engine cost model (WaverRegionLoad with `costs`), not only its messages
        WaverRegionLoad::Weights costs; // the engine cost model of the exact load
        bool chain = false;     // a group may also go to a block without room if that block gives a group back (exact load only)
        bool wide = false;      // a group may also go to any block it has edges to, not only to the receiver and its neighbors
        bool verbose = false;
    };

    // Busy-PE refinement: the makespan follows the busiest engine, not the total number of messages.
    //
    // The load of a block is what its engine handles: for each message it exchanges, alpha_end plus beta_end times the
    // volume. In every round the busiest blocks are visited. For such a block a, the sender groups S_ae that are
    // cheapest to get rid of are tried, one by one, to every block that can take it (e and the neighbors of e, with room
    // left). A move changes the load of a, of the receiver and of every block the group has edges to, and it is kept
    // when the p-norm of those loads goes down (a high power makes the busiest of them count most), so it may
    // add messages elsewhere as long as the busy block ends up relieved.
    class WaverBusyRefinement {
        using graph_t = CSRGraph;
        using p_manager_t = PartitionManager;
        using bv_manager_t = BoundaryVertexManager;
        using Base = WaverMergeRefinement;

        template<typename DistanceOracleT>
        class BusyEngine : public Base::Engine<DistanceOracleT> {
            using E = Base::Engine<DistanceOracleT>;
            using E::g; using E::pm; using E::bv; using E::cfg; using E::k; using E::quotient;
            using E::alpha_end; using E::beta_end;
            using E::collect_group; using E::measure; using E::apply; using E::undo; using E::targets_of;
            using E::locked; using E::eject;

        public:
            BusyEngine(graph_t &t_g, p_manager_t &t_pm, bv_manager_t &t_bv, const DistanceOracleT &t_oracle,
                       const WaverBusyRefinementConfiguration &t_busy,
                       const f64 t_alpha_end, const f64 t_beta_end, const f64 t_alpha_link, const f64 t_beta_link)
                : E(t_g, t_pm, t_bv, t_oracle, unused_config(), t_alpha_end, t_beta_end, t_alpha_link, t_beta_link),
                  busy(t_busy) {
                if (busy.exact) { tracker.initialize(t_g, t_pm, this->k, busy.costs); }
            }

            void run() {
                for (u64 round = 0; round < busy.max_rounds; ++round) {
                    std::vector<f64> loads(k);
                    for (partition_t x = 0; x < k; ++x) { loads[x] = load(x); }
                    std::vector<partition_t> order(k);
                    for (partition_t x = 0; x < k; ++x) { order[x] = x; }
                    const size_t top = std::max<size_t>(1, static_cast<size_t>(busy.top_fraction * static_cast<f64>(k)));
                    std::partial_sort(order.begin(), order.begin() + std::min(top, order.size()), order.end(),
                                      [&](const partition_t x, const partition_t y) { return loads[x] > loads[y]; });
                    order.resize(std::min(top, order.size()));
                    const f64 threshold = loads[order.back()];
                    const f64 max_before = loads[order.front()];

                    u64 tried = 0, kept = 0;
                    for (const partition_t a: order) {
                        if (load(a) < threshold) { continue; }
                        const auto result = relieve(a);
                        tried += result.first;
                        kept += result.second;
                    }
                    if (busy.verbose) {
                        f64 max_after = 0.0;
                        for (partition_t x = 0; x < k; ++x) { max_after = std::max(max_after, load(x)); }
                        std::cout << "[busy] round " << round << ": max load " << max_before << " -> " << max_after
                                  << ", " << tried << " moves tried, " << kept << " kept" << std::endl;
                    }
                    if (kept == 0) { break; }
                }
            }

        private:
            const WaverBusyRefinementConfiguration &busy;
            WaverRegionLoad tracker;
            std::vector<char> in_set_mark;

            static const WaverMergeRefinementConfiguration &unused_config() {
                static const WaverMergeRefinementConfiguration c("unused");
                return c;
            }

            f64 load(const partition_t x) {
                if (busy.exact) { return tracker.load(x, quotient.neighbors(x).size()); }
                f64 l = 0.0;
                const f64 a = busy.alpha > 0.0 ? busy.alpha : alpha_end;
                const f64 b = busy.beta > 0.0 ? busy.beta : beta_end;
                for (const auto &he: quotient.neighbors(x)) { l += a + b * static_cast<f64>(he.weight); }
                return l;
            }

            f64 cost(const std::vector<partition_t> &blocks) {
                f64 c = 0.0;
                for (const partition_t x: blocks) { c += std::pow(load(x), busy.power); }
                return c;
            }

            // Moves the group of a to the full block c, which then gives one of its own groups back (to a block with
            // room, often a itself), and keeps both moves if the loads of the blocks they touch go down (sum of the
            // loads to the power). Otherwise everything is restored. Returns whether the moves were kept.
            bool try_exchange(const partition_t a, typename Base::Candidate &group, const partition_t c) {
                std::vector<typename Base::Move> moves;
                apply(group, a, c);
                for (const vertex_t u: group.group) { locked[u] = 1; }
                moves.push_back(typename Base::Move{group.group, a, c});
                bool ok = true;
                if (pm.bweights[c] > pm.get_lmax(c)) {
                    f64 total = 0.0;
                    partition_t received = NO_ID;
                    ok = eject(c, true, moves, total, received) && pm.bweights[received] <= pm.get_lmax(received);
                }
                for (const auto &m: moves) { for (const vertex_t u: m.group) { locked[u] = 0; } }
                const auto restore = [&]() {
                    for (size_t i = moves.size(); i-- > 0;) {
                        typename Base::Move copy = moves[i];
                        undo(copy);
                    }
                };
                if (!ok) {
                    restore();
                    return false;
                }
                // the blocks the moves touch: where the groups were and are, and the blocks of their neighbors
                std::vector<partition_t> touched;
                for (const auto &m: moves) {
                    touched.push_back(m.from);
                    touched.push_back(m.to);
                    for (const vertex_t u: m.group) {
                        for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) { touched.push_back(pm[g.edges_v[j]]); }
                    }
                }
                tracker.clear_dirty();
                for (const auto &m: moves) { for (const vertex_t u: m.group) { tracker.move(u, m.from, m.to); } }
                for (const partition_t x: tracker.dirty_blocks()) { touched.push_back(x); }
                std::sort(touched.begin(), touched.end());
                touched.erase(std::unique(touched.begin(), touched.end()), touched.end());
                const f64 after = cost(touched);
                for (size_t i = moves.size(); i-- > 0;) { for (const vertex_t u: moves[i].group) { tracker.move(u, moves[i].to, moves[i].from); } }
                restore();
                const f64 before = cost(touched);
                if (!(before - after > busy.min_gain * before)) { return false; }
                // keep: make the moves again
                for (const auto &m: moves) {
                    typename Base::Candidate again;
                    again.group = m.group;
                    measure(m.from, again);
                    apply(again, m.from, m.to);
                    for (const vertex_t u: m.group) { tracker.move(u, m.from, m.to); }
                }
                return true;
            }

            // Tries to lower the load of a. Returns (moves tried, moves kept).
            std::pair<u64, u64> relieve(const partition_t a) {
                u64 tried = 0, kept = 0;
                // The cheapest messages of a first: giving up a small sender group removes a message outright.
                std::vector<std::pair<weight_t, partition_t> > cand;
                for (const auto &he: quotient.neighbors(a)) { cand.emplace_back(he.weight, he.target); }
                std::sort(cand.begin(), cand.end());
                if (cand.size() > busy.max_tries) { cand.resize(busy.max_tries); }

                std::vector<partition_t> targets, affected, full;
                for (const auto &[w, e]: cand) {
                    if (quotient.weight(a, e) == 0) { continue; }
                    typename Base::Candidate group;
                    collect_group(a, e, group.group);
                    if (group.group.empty()) { continue; }
                    measure(a, group);
                    if (pm.bweights[a] - group.weight < pm.get_lmin(a)) { continue; }
                    targets_of(a, e, targets);
                    if (busy.wide) {
                        // the blocks the group already touches: moving it there can end the message to e outright
                        for (const auto &[d, dw]: group.ext) {
                            if (std::find(targets.begin(), targets.end(), d) == targets.end()) { targets.push_back(d); }
                        }
                    }

                    partition_t best = NO_ID;
                    f64 best_gain = 0.0;
                    full.clear();
                    for (const partition_t c: targets) {
                        if (pm.bweights[c] + group.weight > pm.get_lmax(c)) {
                            if (busy.chain && busy.exact && c != a) { full.push_back(c); }
                            continue;
                        }
                        f64 before, after;
                        if (busy.exact) {
                            // the blocks whose counters change are only known once the move is made
                            tracker.clear_dirty();
                            apply(group, a, c);
                            for (const vertex_t u: group.group) { tracker.move(u, a, c); }
                            affected.assign({a, c});
                            for (const auto &[d, dw]: group.ext) { affected.push_back(d); }
                            for (const partition_t x: tracker.dirty_blocks()) { affected.push_back(x); }
                            std::sort(affected.begin(), affected.end());
                            affected.erase(std::unique(affected.begin(), affected.end()), affected.end());
                            after = cost(affected);
                            for (const vertex_t u: group.group) { tracker.move(u, c, a); }
                            typename Base::Move m{group.group, a, c};
                            undo(m);
                            before = cost(affected);
                        } else {
                            affected.assign({a, c});
                            for (const auto &[d, dw]: group.ext) { affected.push_back(d); }
                            std::sort(affected.begin(), affected.end());
                            affected.erase(std::unique(affected.begin(), affected.end()), affected.end());
                            before = cost(affected);
                            apply(group, a, c);
                            after = cost(affected);
                            typename Base::Move m{group.group, a, c};
                            undo(m);
                        }
                        ++tried;
                        const f64 gain = before - after;
                        if (gain > best_gain && gain > busy.min_gain * before) {
                            best_gain = gain;
                            best = c;
                        }
                    }
                    if (best == NO_ID) {
                        // no block with room: a block that is full may take the group if it gives one back
                        for (size_t i = 0; i < full.size() && i < 3; ++i) {
                            ++tried;
                            if (try_exchange(a, group, full[i])) {
                                ++kept;
                                break;
                            }
                        }
                        continue;
                    }
                    // undo() consumed nothing of `group` but re-measure: the quotient rows changed back and forth
                    measure(a, group);
                    apply(group, a, best);
                    if (busy.exact) { for (const vertex_t u: group.group) { tracker.move(u, a, best); } }
                    ++kept;

                }
                return {tried, kept};
            }
        };

    public:
        template<typename DistanceOracleT>
        static void refine(graph_t &g, p_manager_t &p_manager, bv_manager_t &bv_manager,
                           const DistanceOracleT &d_oracle, const WaverBusyRefinementConfiguration &config,
                           const f64 alpha_end, const f64 beta_end, const f64 alpha_link, const f64 beta_link) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverBusyRefinement", "refine");
            BusyEngine<DistanceOracleT>(g, p_manager, bv_manager, d_oracle, config,
                                        alpha_end, beta_end, alpha_link, beta_link).run();
        }
    };
} // namespace HeiProMap

#endif //HEIPROMAP_WAVER_BUSY_REFINEMENT_H
