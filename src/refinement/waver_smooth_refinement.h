// Smooth refinement (Waver). Single boundary vertices move to another block where they have neighbors, if that lowers
// (neighbor cells) * message weight + (send regions, receive regions, values) * their weights + hops of cut edges. This is
// what the group moves cannot reach: a block that touches another through a single vertex loses that neighbor by giving the
// vertex away. Uses WaverRegionLoad for the region and value terms. Option: --smooth-refinement.

#ifndef HEIPROMAP_WAVER_SMOOTH_REFINEMENT_H
#define HEIPROMAP_WAVER_SMOOTH_REFINEMENT_H

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
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
    class WaverSmoothRefinement {
        using graph_t = CSRGraph;
        using p_manager_t = PartitionManager;

        // Fast open-addressing hash table for cut edges between blocks (u64 key -> int64_t count)
        struct FastCutTable {
            struct Entry {
                u64 key = 0;
                int64_t count = 0;
            };
            std::vector<Entry> table;
            u64 mask = 0;

            void init(size_t expected_elements) {
                size_t cap = 16;
                while (cap < expected_elements * 2) { cap <<= 1; }
                mask = cap - 1;
                table.assign(cap, {0, 0});
            }

            static inline u64 hash(u64 x) {
                x ^= x >> 30;
                x *= 0xbf58476d1ce4e5b9ULL;
                x ^= x >> 27;
                x *= 0x94d049bb133111ebULL;
                x ^= x >> 31;
                return x;
            }

            inline int64_t get(u64 k) const {
                u64 idx = hash(k) & mask;
                while (true) {
                    if (table[idx].key == k) { return table[idx].count; }
                    if (table[idx].key == 0) { return 0; }
                    idx = (idx + 1) & mask;
                }
            }

            inline void add(u64 k, int64_t delta) {
                u64 idx = hash(k) & mask;
                while (true) {
                    if (table[idx].key == k) {
                        table[idx].count += delta;
                        return;
                    }
                    if (table[idx].key == 0) {
                        table[idx].key = k;
                        table[idx].count = delta;
                        return;
                    }
                    idx = (idx + 1) & mask;
                }
            }
        };

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
                return lo * static_cast<u64>(k) + hi + 1; // +1 to ensure non-zero key
            };
            auto dist = [&](const partition_t p, const partition_t q) -> f64 {
                const int64_t dx = static_cast<int64_t>(p) % W - static_cast<int64_t>(q) % W;
                const int64_t dy = static_cast<int64_t>(p) / W - static_cast<int64_t>(q) / W;
                return static_cast<f64>((dx < 0 ? -dx : dx) + (dy < 0 ? -dy : dy));
            };

            // 1. Initial scan: count cut edges and build active boundary list
            FastCutTable cut;
            size_t estimated_cut_pairs = 0;
            std::vector<vertex_t> boundary_vertices;
            boundary_vertices.reserve(g.n / 4);

            for (vertex_t u = 0; u < g.n; ++u) {
                const partition_t x = p_manager.partition[u];
                bool is_boundary = false;
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const partition_t y = p_manager.partition[g.edges_v[j]];
                    if (y != x) {
                        is_boundary = true;
                        if (y > x) { ++estimated_cut_pairs; }
                    }
                }
                if (is_boundary) {
                    boundary_vertices.push_back(u);
                }
            }

            cut.init(estimated_cut_pairs + 1024);
            for (const vertex_t u: boundary_vertices) {
                const partition_t x = p_manager.partition[u];
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const partition_t y = p_manager.partition[g.edges_v[j]];
                    if (y > x) { cut.add(key(x, y), 1); }
                }
            }

            std::vector<std::pair<partition_t, int64_t> > around; // blocks of the neighbors of u (not a) and their counts
            std::vector<partition_t> touched;
            std::vector<vertex_t> next_boundary;
            std::vector<char> is_in_boundary(g.n, 0);
            for (const vertex_t u: boundary_vertices) { is_in_boundary[u] = 1; }

            u64 total_moves = 0;
            for (u64 round = 0; round < config.max_rounds; ++round) {
                u64 moves = 0;
                next_boundary.clear();

                for (const vertex_t u: boundary_vertices) {
                    const partition_t a = p_manager.partition[u];
                    const weight_t w = g.v_weight(u);
                    if (p_manager.n_vertices[a] <= 1 || p_manager.bweights[a] - w < p_manager.lmin[a]) {
                        next_boundary.push_back(u);
                        continue;
                    }

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
                    if (!boundary) {
                        is_in_boundary[u] = 0;
                        continue; // No longer boundary
                    }

                    f64 best_delta = -config.min_gain;
                    partition_t best = NO_ID;
                    for (const auto &[b, to_b]: around) {
                        if (p_manager.bweights[b] + w > p_manager.lmax[b]) { continue; }
                        // neighbor cells and hops of the cut edges
                        f64 delta = 0.0;
                        int64_t pairs = 0;
                        int64_t to_a = 0;
                        for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                            if (p_manager.partition[g.edges_v[j]] == a) { ++to_a; }
                        }
                        // the pair (a,b)
                        {
                            const int64_t before = cut.get(key(a, b)), after = before - to_b + to_a;
                            pairs += (after > 0) - (before > 0);
                            delta += config.hop_weight * dist(a, b) * static_cast<f64>(to_a - to_b);
                        }
                        for (const auto &[c, n_c]: around) {
                            if (c == b) { continue; }
                            const int64_t before_a = cut.get(key(a, c)), after_a = before_a - n_c;
                            const int64_t before_b = cut.get(key(b, c)), after_b = before_b + n_c;
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
                    if (best == NO_ID) {
                        next_boundary.push_back(u);
                        continue;
                    }

                    // apply move
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const partition_t c = p_manager.partition[g.edges_v[j]];
                        if (c != a) { cut.add(key(a, c), -1); }
                    }
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const partition_t c = p_manager.partition[g.edges_v[j]];
                        if (c != best) { cut.add(key(best, c), 1); }
                    }
                    tracker.move(u, a, best);
                    p_manager.move_serial(u, w, a, best);
                    ++moves;

                    next_boundary.push_back(u);
                    // Add neighbors of u to boundary if they weren't already
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const vertex_t v = g.edges_v[j];
                        if (!is_in_boundary[v]) {
                            is_in_boundary[v] = 1;
                            next_boundary.push_back(v);
                        }
                    }
                }
                total_moves += moves;
                boundary_vertices = std::move(next_boundary);
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
