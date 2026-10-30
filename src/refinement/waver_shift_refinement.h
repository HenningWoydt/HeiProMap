#ifndef HEIPROMAP_WAVER_SHIFT_REFINEMENT_H
#define HEIPROMAP_WAVER_SHIFT_REFINEMENT_H

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

namespace HeiProMap {
    // Hyperparameters of the shift refinement.
    class WaverShiftRefinementConfiguration {
    public:
        explicit WaverShiftRefinementConfiguration(std::string t_name) : name(std::move(t_name)) {}

        std::string name;
        bool enabled = false;
        u64 max_rounds = 5;        // rounds over the longest messages; stops earlier once a round shifts nothing
        u64 max_messages = 5000;   // only this many longest messages are considered per round
        f64 min_gain = 1e-6;       // a shift must improve the cost by more than this
        bool verbose = false;
    };

    // Shift refinement: places the two blocks of a long message next to each other.
    //
    // One block of the message is moved to a grid cell next to the other. The block that was in that cell and every
    // block on a path from there to the cell the moved block left shift by one cell along the path, so that the cell
    // the moved block left is filled again. The result is one long move and many moves of a single hop. Both ends of
    // the message, the four cells around the other end and both L-shaped paths (horizontal first or vertical first) are
    // tried, and the shift with the largest drop of the hop cost is made. Only the placement of the blocks changes: the
    // messages, their volume and the contents of the blocks stay as they are.
    class WaverShiftRefinement {
        using graph_t = CSRGraph;
        using p_manager_t = PartitionManager;

        struct Row {
            partition_t target;
            f64 weight;
        };

    public:
        template<typename DistanceOracleT>
        static void refine(graph_t &g, p_manager_t &p_manager, const DistanceOracleT &d_oracle,
                           const WaverShiftRefinementConfiguration &config,
                           const f64 alpha_link, const f64 beta_link) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverShiftRefinement", "refine");
            const partition_t k = p_manager.get_k();
            const int64_t W = static_cast<int64_t>(d_oracle.current_W());

            // the quotient rows: block -> (block, hop cost per unit of distance)
            std::vector<std::vector<vertex_t> > blk_verts(k);
            for (vertex_t u = 0; u < g.n; ++u) { blk_verts[p_manager.partition[u]].push_back(u); }
            std::vector<std::vector<Row> > rows(k);
            {
                std::vector<f64> acc(k, 0.0);
                std::vector<partition_t> dirty;
                for (partition_t x = 0; x < k; ++x) {
                    for (const vertex_t u: blk_verts[x]) {
                        for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                            const partition_t y = p_manager.partition[g.edges_v[j]];
                            if (y == x) { continue; }
                            if (acc[y] == 0.0) { dirty.push_back(y); }
                            acc[y] += static_cast<f64>(g.edge_w(j));
                        }
                    }
                    for (const partition_t y: dirty) {
                        rows[x].push_back({y, alpha_link + beta_link * acc[y]});
                        acc[y] = 0.0;
                    }
                    dirty.clear();
                }
            }

            // pos[block]: the cell of a block; occ[cell]: the block in it (or NO_ID)
            std::vector<partition_t> pos(k), occ(k);
            for (partition_t x = 0; x < k; ++x) { pos[x] = x; occ[x] = x; }
            auto dist = [&](const partition_t p, const partition_t q) -> f64 {
                const int64_t dx = static_cast<int64_t>(p) % W - static_cast<int64_t>(q) % W;
                const int64_t dy = static_cast<int64_t>(p) / W - static_cast<int64_t>(q) / W;
                return static_cast<f64>((dx < 0 ? -dx : dx) + (dy < 0 ? -dy : dy));
            };

            std::vector<char> in_set(k, 0);
            // The hop cost of the blocks in `set` at their current cells, with the others fixed.
            auto set_cost = [&](const std::vector<partition_t> &set) -> f64 {
                f64 c = 0.0;
                for (const partition_t x: set) { in_set[x] = 1; }
                for (const partition_t x: set) {
                    for (const Row &r: rows[x]) {
                        const f64 d = r.weight * dist(pos[x], pos[r.target]);
                        c += in_set[r.target] ? 0.5 * d : d;
                    }
                }
                for (const partition_t x: set) { in_set[x] = 0; }
                return c;
            };

            std::vector<partition_t> path, set, moved_to;
            // The cells of an L-shaped path from t to s (inclusive); false if it leaves the grid.
            auto make_path = [&](const partition_t t, const partition_t s, const bool horizontal_first) -> bool {
                path.clear();
                int64_t x = static_cast<int64_t>(t) % W, y = static_cast<int64_t>(t) / W;
                const int64_t sx = static_cast<int64_t>(s) % W, sy = static_cast<int64_t>(s) / W;
                auto push = [&]() {
                    const int64_t c = y * W + x;
                    if (c < 0 || c >= static_cast<int64_t>(k)) { return false; }
                    path.push_back(static_cast<partition_t>(c));
                    return true;
                };
                if (!push()) { return false; }
                for (int phase = 0; phase < 2; ++phase) {
                    const bool horiz = (phase == 0) == horizontal_first;
                    if (horiz) {
                        while (x != sx) { x += sx > x ? 1 : -1; if (!push()) { return false; } }
                    } else {
                        while (y != sy) { y += sy > y ? 1 : -1; if (!push()) { return false; } }
                    }
                }
                return true;
            };

            size_t total = 0;
            for (u64 round = 0; round < config.max_rounds; ++round) {
                // the longest messages (by hop cost), once per pair
                std::vector<std::pair<f64, std::pair<partition_t, partition_t> > > msgs;
                for (partition_t x = 0; x < k; ++x) {
                    for (const Row &r: rows[x]) {
                        if (x < r.target && dist(pos[x], pos[r.target]) > 1.0) {
                            msgs.push_back({r.weight * dist(pos[x], pos[r.target]), {x, r.target}});
                        }
                    }
                }
                std::sort(msgs.begin(), msgs.end(), [](const auto &a, const auto &b) { return a.first > b.first; });
                if (msgs.size() > config.max_messages) { msgs.resize(config.max_messages); }

                size_t shifts = 0;
                for (const auto &entry: msgs) {
                    for (int side = 0; side < 2; ++side) {
                        const partition_t a = side == 0 ? entry.second.first : entry.second.second;
                        const partition_t b = side == 0 ? entry.second.second : entry.second.first;
                        if (dist(pos[a], pos[b]) <= 1.0) { break; }
                        // move b next to a
                        const int64_t ax = static_cast<int64_t>(pos[a]) % W, ay = static_cast<int64_t>(pos[a]) / W;
                        const int64_t dx[4] = {1, -1, 0, 0}, dy[4] = {0, 0, 1, -1};
                        f64 best_gain = config.min_gain;
                        std::vector<partition_t> best_path;
                        for (int dir = 0; dir < 4; ++dir) {
                            const int64_t tx = ax + dx[dir], ty = ay + dy[dir];
                            if (tx < 0 || tx >= W || ty < 0) { continue; }
                            const int64_t tc = ty * W + tx;
                            if (tc >= static_cast<int64_t>(k)) { continue; }
                            const partition_t t = static_cast<partition_t>(tc);
                            if (t == pos[b] || t == pos[a]) { continue; }
                            for (int hv = 0; hv < 2; ++hv) {
                                if (!make_path(t, pos[b], hv == 0)) { continue; }
                                // path[0] = t, ..., path.back() = pos[b]; b goes to t, the block at path[i] to
                                // path[i + 1]
                                set.clear();
                                bool fits = true;
                                for (size_t i = 0; i + 1 < path.size(); ++i) {
                                    const partition_t blk = occ[path[i]];
                                    if (blk == NO_ID || blk == a) { fits = false; break; } // keep a where it is
                                    set.push_back(blk);
                                }
                                if (!fits) { continue; }
                                set.push_back(b);
                                // the weight limits of the cells the blocks move to
                                for (size_t i = 0; i + 1 < path.size() && fits; ++i) {
                                    const partition_t blk = occ[path[i]];
                                    fits = p_manager.bweights[blk] <= p_manager.get_lmax(path[i + 1]) &&
                                           p_manager.bweights[blk] >= p_manager.get_lmin(path[i + 1]);
                                }
                                fits = fits && p_manager.bweights[b] <= p_manager.get_lmax(t) &&
                                       p_manager.bweights[b] >= p_manager.get_lmin(t);
                                if (!fits) { continue; }

                                const f64 before = set_cost(set);
                                std::vector<partition_t> old_pos(set.size());
                                for (size_t i = 0; i < set.size(); ++i) { old_pos[i] = pos[set[i]]; }
                                for (size_t i = 0; i + 1 < set.size(); ++i) { pos[set[i]] = path[i + 1]; }
                                pos[b] = t;
                                const f64 after = set_cost(set);
                                for (size_t i = 0; i < set.size(); ++i) { pos[set[i]] = old_pos[i]; }
                                const f64 gain = before - after;
                                if (gain > best_gain) {
                                    best_gain = gain;
                                    best_path = path;
                                }
                            }
                        }
                        if (best_path.empty()) { continue; }
                        // make the shift
                        std::vector<partition_t> blocks;
                        for (size_t i = 0; i + 1 < best_path.size(); ++i) { blocks.push_back(occ[best_path[i]]); }
                        const partition_t b_cell_old = pos[b];
                        for (size_t i = 0; i < blocks.size(); ++i) {
                            pos[blocks[i]] = best_path[i + 1];
                            occ[best_path[i + 1]] = blocks[i];
                        }
                        pos[b] = best_path[0];
                        occ[best_path[0]] = b;
                        (void) b_cell_old;
                        ++shifts;
                        break;
                    }
                }
                total += shifts;
                if (config.verbose) {
                    std::cout << "[shift] round " << round << ": " << msgs.size() << " long messages, " << shifts
                              << " shifts" << std::endl;
                }
                if (shifts == 0) { break; }
            }
            if (total == 0) { return; }

            // relabel: the vertices and the weights of every block follow it to its new cell
            std::vector<weight_t> new_weights(k, 0);
            for (partition_t x = 0; x < k; ++x) { new_weights[pos[x]] = p_manager.bweights[x]; }
            for (vertex_t u = 0; u < g.n; ++u) { p_manager.partition[u] = pos[p_manager.partition[u]]; }
            for (partition_t x = 0; x < k; ++x) { p_manager.bweights[x] = new_weights[x]; }
        }
    };
} // namespace HeiProMap

#endif //HEIPROMAP_WAVER_SHIFT_REFINEMENT_H
