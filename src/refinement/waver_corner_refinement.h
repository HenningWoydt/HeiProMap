// Corner refinement (Waver). Run once on the finished mapping of the input graph. Single vertices move, or swap with
// an adjacent vertex of a full PE, between neighboring PEs when that lowers sum(engine load^power) over the PEs (the
// busiest count most). The load of a PE is the engine cost model: a cost per message (neighbor PE), per cut edge and per
// internal edge. The number of neighbors enters smoothly, as sum over neighbors of 1 - exp(-w / tau) for a cut weight w,
// so every vertex moved out of a small contact lowers the objective a little and a corner contact of a few edges, which
// costs a whole message on two engines, shrinks until it is gone. Unit of work: single vertices (all weights 1, PE limit
// lmax). Options: --corner-rounds, --corner-power, --corner-tau, --corner-top, --corner-swaps.

#ifndef HEIPROMAP_WAVER_CORNER_REFINEMENT_H
#define HEIPROMAP_WAVER_CORNER_REFINEMENT_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "../definitions.h"
#include "../datastructures/csr_graph.h"

namespace HeiProMap {

    struct WaverCornerRefinementConfiguration {
        bool enabled = false;
        u64 max_rounds = 10;
        f64 power = 6.0;        // the engine loads are compared as sum(load^power): the busiest PEs count most
        f64 tau = 2.0;          // a contact of w cut edges costs a message 1 - exp(-w / tau): small contacts shrink first
        f64 message = 655.0;    // engine cost of a message (send + receive)
        f64 cut_edge = 20.4;    // engine cost of a cut edge (a value sent and received)
        f64 internal_edge = 9.0; // engine cost of an edge inside a PE (two local matrix entries)
        f64 min_gain = 1e-9;    // relative improvement of the objective a move needs
        f64 top_fraction = 1.0; // the share of PEs, busiest first, whose vertices are tried in a round
        const std::vector<f64> *weights = nullptr; // per PE: weight of its load in the objective (null: 1 for all), e.g. how
                                                  // often the PE is on the critical chain of a schedule
        u64 swap_tries = 4;     // neighbors in the target block a vertex may be swapped with if the block is full
        bool verbose = false;
    };

    // Corner refinement: single vertices move between neighboring PEs so that the busiest engines get less work. The
    // work of a PE is the engine cost model: a cost per message (neighbor PE), per cut edge and per internal edge. The
    // number of neighbors enters smoothly, as sum over neighbors of 1 - exp(-w / tau) with w the cut weight to that
    // neighbor, so that every vertex moved out of a small contact lowers the objective a little; a contact of a few
    // edges at a corner, which costs a whole message, shrinks until it is gone. The objective is the sum of the loads
    // to a power. A move into a PE without room is a swap with a vertex of that PE that is adjacent to it. Works on the
    // finest graph (one vertex per row), where all vertex weights are 1 and the PE limits are known.
    class WaverCornerRefinement {
    public:
        // part: the PE of every vertex (changed in place). lmax: weight limit of a PE.
        static u64 refine(const CSRGraph &g, std::vector<partition_t> &part, const u32 K, const weight_t lmax,
                          const WaverCornerRefinementConfiguration &cfg) {
            State s(g, part, K, lmax, cfg);
            u64 total_moves = 0;
            for (u64 round = 0; round < cfg.max_rounds; ++round) {
                const f64 before = s.objective();
                s.rebuild_members();
                u64 moves = 0;
                // the busiest PEs first
                std::vector<u32> order(K);
                for (u32 p = 0; p < K; ++p) { order[p] = p; }
                std::sort(order.begin(), order.end(), [&](u32 a, u32 b) { return s.load[a] > s.load[b]; });
                const u32 limit = std::max<u32>(1, static_cast<u32>(std::ceil(cfg.top_fraction * static_cast<f64>(K))));
                for (u32 i = 0; i < limit && i < K; ++i) { moves += s.sweep_pe(order[i]); }
                total_moves += moves;
                if (cfg.verbose) {
                    std::cout << "[corner] round " << round << ": " << moves << " moves, objective " << before << " -> "
                              << s.objective() << ", max load " << s.max_load() << std::endl;
                }
                if (moves == 0) { break; }
            }
            return total_moves;
        }

    private:
        struct State {
            const CSRGraph &g;
            std::vector<partition_t> &part;
            const u32 K;
            const weight_t lmax;
            const WaverCornerRefinementConfiguration &cfg;
            std::vector<std::vector<vertex_t> > members; // vertices of a PE (may hold stale entries)
            std::unordered_map<u64, s64> pair_weight;    // cut weight between two PEs (key a * K + b, a < b)
            std::vector<f64> msgsum, cutw, internal, load;
            std::vector<weight_t> bw;
            f64 sum_power = 0.0;
            std::vector<partition_t> blocks;             // scratch
            std::vector<s64> counts;                     // scratch: edges of the moved vertex to each block

            State(const CSRGraph &t_g, std::vector<partition_t> &t_part, const u32 t_K, const weight_t t_lmax,
                  const WaverCornerRefinementConfiguration &t_cfg)
                : g(t_g), part(t_part), K(t_K), lmax(t_lmax), cfg(t_cfg) {
                members.resize(K);
                msgsum.assign(K, 0.0);
                cutw.assign(K, 0.0);
                internal.assign(K, 0.0);
                load.assign(K, 0.0);
                bw.assign(K, 0);
                counts.assign(K, 0);
                for (vertex_t u = 0; u < g.n; ++u) {
                    const partition_t a = part[u];
                    members[a].push_back(u);
                    bw[a] += g.v_weight(u);
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const vertex_t v = g.edges_v[j];
                        const partition_t b = part[v];
                        const f64 w = static_cast<f64>(g.edge_w(j));
                        if (b == a) {
                            internal[a] += 0.5 * w;
                        } else {
                            cutw[a] += w;
                            if (a < b) { pair_weight[key(a, b)] += static_cast<s64>(g.edge_w(j)); }
                        }
                    }
                }
                for (const auto &[k, w]: pair_weight) {
                    const u32 a = static_cast<u32>(k / K), b = static_cast<u32>(k % K);
                    const f64 m = smooth(w);
                    msgsum[a] += m;
                    msgsum[b] += m;
                }
                for (u32 p = 0; p < K; ++p) {
                    load[p] = compute_load(p);
                    sum_power += omega(p) * std::pow(load[p], cfg.power);
                }
            }

            void rebuild_members() {
                for (auto &m: members) { m.clear(); }
                for (vertex_t u = 0; u < g.n; ++u) { members[part[u]].push_back(u); }
            }

            f64 omega(const u32 p) const { return cfg.weights ? (*cfg.weights)[p] : 1.0; }

            u64 key(const u32 a, const u32 b) const { return static_cast<u64>(std::min(a, b)) * K + std::max(a, b); }
            f64 smooth(const s64 w) const { return w <= 0 ? 0.0 : 1.0 - std::exp(-static_cast<f64>(w) / cfg.tau); }
            f64 compute_load(const u32 p) const {
                return cfg.message * msgsum[p] + cfg.cut_edge * cutw[p] + cfg.internal_edge * internal[p];
            }
            f64 objective() const { return sum_power; }
            f64 max_load() const { return *std::max_element(load.begin(), load.end()); }

            void change_pair(const u32 a, const u32 b, const s64 delta, std::vector<u32> &touched) {
                s64 &w = pair_weight[key(a, b)];
                const f64 before = smooth(w);
                w += delta;
                const f64 diff = smooth(w) - before;
                msgsum[a] += diff;
                msgsum[b] += diff;
                if (w == 0) { pair_weight.erase(key(a, b)); }
                touched.push_back(a);
                touched.push_back(b);
            }

            // Moves u to the PE c, updating the pair weights, the loads and the objective.
            void move(const vertex_t u, const partition_t c, std::vector<u32> &touched) {
                const partition_t a = part[u];
                // edges of u to each block
                blocks.clear();
                f64 deg = 0.0;
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const partition_t d = part[g.edges_v[j]];
                    if (counts[d] == 0) { blocks.push_back(d); }
                    counts[d] += static_cast<s64>(g.edge_w(j));
                    deg += static_cast<f64>(g.edge_w(j));
                }
                const f64 to_a = static_cast<f64>(counts[a]), to_c = static_cast<f64>(counts[c]);
                touched.clear();
                for (const partition_t d: blocks) {
                    if (d != a) { change_pair(a, d, -counts[d], touched); }
                    if (d != c) { change_pair(c, d, counts[d], touched); }
                }
                cutw[a] += -(deg - to_a) + to_a;
                cutw[c] += (deg - to_c) - to_c;
                internal[a] -= to_a;
                internal[c] += to_c;
                touched.push_back(a);
                touched.push_back(c);
                for (const partition_t d: blocks) { counts[d] = 0; }
                bw[a] -= g.v_weight(u);
                bw[c] += g.v_weight(u);
                part[u] = c;
                std::sort(touched.begin(), touched.end());
                touched.erase(std::unique(touched.begin(), touched.end()), touched.end());
                for (const u32 p: touched) {
                    sum_power -= omega(p) * std::pow(load[p], cfg.power);
                    load[p] = compute_load(p);
                    sum_power += omega(p) * std::pow(load[p], cfg.power);
                }
            }

            // Tries the boundary vertices of the PE a: a move to a neighboring PE with room, or a swap with an adjacent
            // vertex of a full one. Kept if the objective goes down.
            u64 sweep_pe(const u32 a) {
                u64 moves = 0;
                std::vector<u32> touched;
                std::vector<vertex_t> mine = members[a];
                std::vector<partition_t> targets;
                for (const vertex_t u: mine) {
                    if (part[u] != a) { continue; }
                    targets.clear();
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const partition_t d = part[g.edges_v[j]];
                        if (d != a && std::find(targets.begin(), targets.end(), d) == targets.end()) { targets.push_back(d); }
                    }
                    for (const partition_t c: targets) {
                        const f64 before = sum_power;
                        const weight_t uw = g.v_weight(u);
                        if (bw[c] + uw <= lmax && bw[a] - uw > 0) {
                            move(u, c, touched);
                            if (sum_power < before * (1.0 - cfg.min_gain)) {
                                ++moves;
                                break;
                            }
                            move(u, a, touched);
                            continue;
                        }
                        // a full block: swap u with an adjacent vertex v of c
                        u64 tried = 0;
                        bool done = false;
                        for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1] && tried < cfg.swap_tries; ++j) {
                            const vertex_t v = g.edges_v[j];
                            if (part[v] != c) { continue; }
                            ++tried;
                            move(u, c, touched);
                            move(v, a, touched);
                            if (sum_power < before * (1.0 - cfg.min_gain)) {
                                ++moves;
                                done = true;
                                break;
                            }
                            move(v, c, touched);
                            move(u, a, touched);
                        }
                        if (done) { break; }
                    }
                }
                return moves;
            }
        };
    };

}

#endif //HEIPROMAP_WAVER_CORNER_REFINEMENT_H
