#ifndef HEIPROMAP_WAVER_SMOOTH_REFINEMENT_H
#define HEIPROMAP_WAVER_SMOOTH_REFINEMENT_H

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../definitions.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../utility/profiler.h"
#include "waver_region_load.h"

namespace HeiProMap {
    // Hyperparameters of the smooth refinement.
    class WaverSmoothRefinementConfiguration {
    public:
        explicit WaverSmoothRefinementConfiguration(std::string t_name) : name(std::move(t_name)) {}

        std::string name;
        bool enabled = false;
        u64 max_rounds = 5;       // sweeps over the boundary vertices; stops earlier once a sweep moves nothing
        f64 hop_weight = 1.0;     // cost of one hop of one cut edge, next to the engine load weights
        f64 min_gain = 1e-6;      // a move must lower the cost by more than this
        WaverRegionLoad::Weights costs; // the engine cost model of the load terms
        bool verbose = false;
    };

    // Smooth refinement: single boundary vertices move to another cell where they have neighbors, if that lowers
    //   (neighbor cells) * message weight + (send regions, receive regions, values) * their weights + hops of cut edges.
    // The first part is what the group moves of the other refinements cannot reach: a cell that touches one block through
    // a single vertex loses that neighbor by giving the vertex away, which no cut-edge objective rewards enough.
    class WaverSmoothRefinement {
        using graph_t = CSRGraph;
        using p_manager_t = PartitionManager;

    public:
        template<typename DistanceOracleT>
        static void refine(graph_t &g, p_manager_t &p_manager, const DistanceOracleT &d_oracle,
                           const WaverSmoothRefinementConfiguration &config) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverSmoothRefinement", "refine");
            const partition_t k = p_manager.get_k();
            const int64_t W = static_cast<int64_t>(d_oracle.current_W());
            const WaverRegionLoad::Weights &weights = config.costs;

            WaverRegionLoad tracker;
            tracker.initialize(g, p_manager, k, weights);

            auto key = [k](const partition_t x, const partition_t y) -> u64 {
                const u64 lo = std::min(x, y), hi = std::max(x, y);
                return lo * static_cast<u64>(k) + hi;
            };
            auto dist = [&](const partition_t p, const partition_t q) -> f64 {
                const int64_t dx = static_cast<int64_t>(p) % W - static_cast<int64_t>(q) % W;
                const int64_t dy = static_cast<int64_t>(p) / W - static_cast<int64_t>(q) / W;
                return static_cast<f64>((dx < 0 ? -dx : dx) + (dy < 0 ? -dy : dy));
            };

            // cut[(x, y)]: the number of edges between the blocks x and y
            std::unordered_map<u64, int64_t> cut;
            for (vertex_t u = 0; u < g.n; ++u) {
                const partition_t x = p_manager.partition[u];
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const partition_t y = p_manager.partition[g.edges_v[j]];
                    if (y > x) { ++cut[key(x, y)]; }
                }
            }
            auto cut_of = [&](const partition_t x, const partition_t y) -> int64_t {
                const auto it = cut.find(key(x, y));
                return it == cut.end() ? 0 : it->second;
            };

            std::vector<std::pair<partition_t, int64_t> > around; // blocks of the neighbors of u (not a) and their counts
            std::vector<partition_t> touched;
            u64 total_moves = 0;
            for (u64 round = 0; round < config.max_rounds; ++round) {
                u64 moves = 0;
                for (vertex_t u = 0; u < g.n; ++u) {
                    const partition_t a = p_manager.partition[u];
                    const weight_t w = g.v_weight(u);
                    if (p_manager.n_vertices[a] <= 1 || p_manager.bweights[a] - w < p_manager.lmin[a]) { continue; }

                    around.clear();
                    bool boundary = false;
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const partition_t c = p_manager.partition[g.edges_v[j]];
                        if (c == a) { continue; }
                        boundary = true;
                        bool found = false;
                        for (auto &e: around) {
                            if (e.first == c) {
                                ++e.second;
                                found = true;
                                break;
                            }
                        }
                        if (!found) { around.emplace_back(c, 1); }
                    }
                    if (!boundary) { continue; }

                    f64 best_delta = -config.min_gain;
                    partition_t best = NO_ID;
                    for (const auto &[b, to_b]: around) {
                        if (p_manager.bweights[b] + w > p_manager.lmax[b]) { continue; }
                        // neighbor cells and hops of the cut edges
                        f64 delta = 0.0;
                        int64_t pairs = 0;
                        // edges to a: become cut with b; edges to b: stop being cut; edges to c: (a,c) -> (b,c)
                        int64_t to_a = 0;
                        for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                            if (p_manager.partition[g.edges_v[j]] == a) { ++to_a; }
                        }
                        // the pair (a,b)
                        {
                            const int64_t before = cut_of(a, b), after = before - to_b + to_a;
                            pairs += (after > 0) - (before > 0);
                            delta += config.hop_weight * dist(a, b) * static_cast<f64>(to_a - to_b);
                        }
                        for (const auto &[c, n_c]: around) {
                            if (c == b) { continue; }
                            const int64_t before_a = cut_of(a, c), after_a = before_a - n_c;
                            const int64_t before_b = cut_of(b, c), after_b = before_b + n_c;
                            pairs += (after_a > 0) - (before_a > 0) + (after_b > 0) - (before_b > 0);
                            delta += config.hop_weight * static_cast<f64>(n_c) * (dist(b, c) - dist(a, c));
                        }
                        delta += 2.0 * weights.message * static_cast<f64>(pairs);

                        // regions and values, from the tracker
                        tracker.clear_dirty();
                        tracker.move(u, a, b);
                        touched = tracker.dirty_blocks();
                        f64 after = 0.0;
                        for (const partition_t x: touched) { after += tracker.load(x, 0); }
                        tracker.move(u, b, a);
                        f64 before = 0.0;
                        for (const partition_t x: touched) { before += tracker.load(x, 0); }
                        delta += after - before;

                        if (delta < best_delta) {
                            best_delta = delta;
                            best = b;
                        }
                    }
                    if (best == NO_ID) { continue; }

                    // apply
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const partition_t c = p_manager.partition[g.edges_v[j]];
                        if (c != a) { --cut[key(a, c)]; }
                    }
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const partition_t c = p_manager.partition[g.edges_v[j]];
                        if (c != best) { ++cut[key(best, c)]; }
                    }
                    tracker.move(u, a, best);
                    p_manager.move_serial(u, w, a, best);
                    ++moves;
                }
                total_moves += moves;
                if (config.verbose) {
                    std::cout << "[smooth] round " << round << ": moved " << moves << " vertices" << std::endl;
                }
                if (moves == 0) { break; }
            }
            if (config.verbose) { std::cout << "[smooth] total moves " << total_moves << std::endl; }
        }
    };
} // namespace HeiProMap

#endif //HEIPROMAP_WAVER_SMOOTH_REFINEMENT_H
