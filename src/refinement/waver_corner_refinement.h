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
        f64 round_gain = 1e-3;  // stop when a round lowers the objective by less than this share
        f64 power = 6.0;        // the engine loads are compared as sum(load^power): the busiest PEs count most
        f64 tau = 2.0;          // a contact of w cut edges costs a message 1 - exp(-w / tau): small contacts shrink first
        f64 message = 655.0;    // engine cost of a message (send + receive)
        f64 cut_edge = 20.4;    // engine cost of a cut edge (a value sent and received)
        f64 internal_edge = 9.0; // engine cost of an edge inside a PE (two local matrix entries)
        f64 min_gain = 1e-9;    // relative improvement of the objective a move needs
        f64 top_fraction = 1.0; // the share of PEs, busiest first, whose vertices are tried in a round
        const std::vector<f64> *weights = nullptr; // per PE: weight of its load in the objective (null: 1 for all), e.g. how
                                                  // often the PE is on the critical chain of a schedule
        bool pair_objective = false; // the objective is sum over message pairs (a, b) of (load_a + load_b)^power: a message
                                     // needs both engines, so two busy PEs next to each other are worse than a busy one
                                     // next to a light one
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
                if (moves == 0 || s.objective() > before * (1.0 - cfg.round_gain)) { break; }
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
            std::vector<std::vector<u32> > adj;          // the PEs each PE has a message with
            std::vector<std::vector<s64> > adjw;         // the cut weight to each of them (same order as adj)
            std::vector<u64> pair_scratch;
            std::vector<weight_t> bw;
            f64 sum_power = 0.0;
            std::vector<f64> smooth_table;               // smooth(w) for small w
            int int_power = 0;                           // cfg.power when it is a small integer, else 0
            std::vector<partition_t> blocks;             // scratch
            std::vector<s64> counts;                     // scratch: edges of the moved vertex to each block
            std::vector<s64> wa, wc;                     // scratch: cut weight of the pairs of a and of c, by the other PE
            std::vector<f64> dmsg;                       // scratch: change of the message sum of a PE

            State(const CSRGraph &t_g, std::vector<partition_t> &t_part, const u32 t_K, const weight_t t_lmax,
                  const WaverCornerRefinementConfiguration &t_cfg)
                : g(t_g), part(t_part), K(t_K), lmax(t_lmax), cfg(t_cfg) {
                members.resize(K);
                if (cfg.power >= 1.0 && cfg.power <= 16.0 && cfg.power == std::floor(cfg.power)) { int_power = static_cast<int>(cfg.power); }
                smooth_table.resize(4096);
                for (size_t w = 1; w < smooth_table.size(); ++w) { smooth_table[w] = 1.0 - std::exp(-static_cast<f64>(w) / cfg.tau); }
                msgsum.assign(K, 0.0);
                cutw.assign(K, 0.0);
                internal.assign(K, 0.0);
                load.assign(K, 0.0);
                bw.assign(K, 0);
                counts.assign(K, 0);
                wa.assign(K, 0);
                wc.assign(K, 0);
                dmsg.assign(K, 0.0);
                adj.assign(K, {});
                adjw.assign(K, {});
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
                    adj[a].push_back(b);
                    adj[b].push_back(a);
                    adjw[a].push_back(w);
                    adjw[b].push_back(w);
                }
                for (u32 p = 0; p < K; ++p) {
                    load[p] = compute_load(p);
                    if (!cfg.pair_objective) { sum_power += omega(p) * powl_(load[p]); }
                }
                if (cfg.pair_objective) {
                    for (const auto &[k, w]: pair_weight) {
                        sum_power += powl_(load[static_cast<u32>(k / K)] + load[static_cast<u32>(k % K)]);
                    }
                }
                pair_weight.clear();
            }

            void rebuild_members() {
                for (auto &m: members) { m.clear(); }
                for (vertex_t u = 0; u < g.n; ++u) { members[part[u]].push_back(u); }
            }

            f64 omega(const u32 p) const { return cfg.weights ? (*cfg.weights)[p] : 1.0; }

            u64 key(const u32 a, const u32 b) const { return static_cast<u64>(std::min(a, b)) * K + std::max(a, b); }
            f64 smooth(const s64 w) const {
                if (w <= 0) { return 0.0; }
                if (w < static_cast<s64>(smooth_table.size())) { return smooth_table[static_cast<size_t>(w)]; }
                return 1.0 - std::exp(-static_cast<f64>(w) / cfg.tau);
            }
            // x^power; an integer power (the default 6) is a few multiplications instead of a pow call
            f64 powl_(const f64 x) const {
                if (int_power > 0) {
                    f64 r = x;
                    for (int i = 1; i < int_power; ++i) { r *= x; }
                    return r;
                }
                return std::pow(x, cfg.power);
            }
            f64 compute_load(const u32 p) const {
                return cfg.message * msgsum[p] + cfg.cut_edge * cutw[p] + cfg.internal_edge * internal[p];
            }
            f64 objective() const { return sum_power; }
            f64 max_load() const { return *std::max_element(load.begin(), load.end()); }

            // sum of (load_a + load_b)^power over the pairs with a message that have an end in `pes` (each pair once)
            f64 pair_sum(const std::vector<u32> &pes) {
                pair_scratch.clear();
                for (const u32 p: pes) { for (const u32 q: adj[p]) { pair_scratch.push_back(key(p, q)); } }
                std::sort(pair_scratch.begin(), pair_scratch.end());
                pair_scratch.erase(std::unique(pair_scratch.begin(), pair_scratch.end()), pair_scratch.end());
                f64 total = 0.0;
                for (const u64 k: pair_scratch) { total += powl_(load[static_cast<u32>(k / K)] + load[static_cast<u32>(k % K)]); }
                return total;
            }

            void change_pair(const u32 a, const u32 b, const s64 delta, std::vector<u32> &touched) {
                const auto ia = std::find(adj[a].begin(), adj[a].end(), b) - adj[a].begin();
                const s64 old_w = ia < static_cast<std::ptrdiff_t>(adj[a].size()) ? adjw[a][ia] : 0;
                const s64 new_w = old_w + delta;
                const f64 diff = smooth(new_w) - smooth(old_w);
                if (old_w <= 0 && new_w > 0) {
                    adj[a].push_back(b);
                    adjw[a].push_back(new_w);
                    adj[b].push_back(a);
                    adjw[b].push_back(new_w);
                } else {
                    const auto ib = std::find(adj[b].begin(), adj[b].end(), a) - adj[b].begin();
                    if (new_w <= 0) {
                        adj[a][ia] = adj[a].back(); adj[a].pop_back();
                        adjw[a][ia] = adjw[a].back(); adjw[a].pop_back();
                        adj[b][ib] = adj[b].back(); adj[b].pop_back();
                        adjw[b][ib] = adjw[b].back(); adjw[b].pop_back();
                    } else {
                        adjw[a][ia] = new_w;
                        adjw[b][ib] = new_w;
                    }
                }
                msgsum[a] += diff;
                msgsum[b] += diff;
                touched.push_back(a);
                touched.push_back(b);
            }

            // The objective after moving u to the PE c, without changing the state (not for the pair objective).
            f64 eval_move(const vertex_t u, const partition_t c) {
                const partition_t a = part[u];
                blocks.clear();
                f64 deg = 0.0;
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const partition_t d = part[g.edges_v[j]];
                    if (counts[d] == 0) { blocks.push_back(d); }
                    counts[d] += static_cast<s64>(g.edge_w(j));
                    deg += static_cast<f64>(g.edge_w(j));
                }
                const f64 to_a = static_cast<f64>(counts[a]), to_c = static_cast<f64>(counts[c]);
                for (size_t i = 0; i < adj[a].size(); ++i) { wa[adj[a][i]] = adjw[a][i]; }
                for (size_t i = 0; i < adj[c].size(); ++i) { wc[adj[c][i]] = adjw[c][i]; }
                // the pairs of a and c change by the edges of u to the other PEs, the pair (a, c) by both directions
                f64 dmsg_a = 0.0, dmsg_c = 0.0;
                for (const partition_t d: blocks) {
                    if (d == a || d == c) { continue; }
                    const f64 da = smooth(wa[d] - counts[d]) - smooth(wa[d]);
                    const f64 dc = smooth(wc[d] + counts[d]) - smooth(wc[d]);
                    dmsg_a += da;
                    dmsg_c += dc;
                    dmsg[d] = da + dc;
                }
                const s64 net_ac = counts[a] - counts[c];
                if (net_ac != 0) {
                    const f64 dac = smooth(wa[c] + net_ac) - smooth(wa[c]);
                    dmsg_a += dac;
                    dmsg_c += dac;
                }
                f64 result = sum_power;
                for (const partition_t d: blocks) {
                    if (d == a || d == c) { continue; }
                    const f64 nl = cfg.message * (msgsum[d] + dmsg[d]) + cfg.cut_edge * cutw[d] + cfg.internal_edge * internal[d];
                    result += omega(d) * (powl_(nl) - powl_(load[d]));
                    dmsg[d] = 0.0;
                }
                const f64 la = cfg.message * (msgsum[a] + dmsg_a) + cfg.cut_edge * (cutw[a] - (deg - to_a) + to_a)
                               + cfg.internal_edge * (internal[a] - to_a);
                const f64 lc = cfg.message * (msgsum[c] + dmsg_c) + cfg.cut_edge * (cutw[c] + (deg - to_c) - to_c)
                               + cfg.internal_edge * (internal[c] + to_c);
                result += omega(a) * (powl_(la) - powl_(load[a])) + omega(c) * (powl_(lc) - powl_(load[c]));
                for (size_t i = 0; i < adj[a].size(); ++i) { wa[adj[a][i]] = 0; }
                for (size_t i = 0; i < adj[c].size(); ++i) { wc[adj[c][i]] = 0; }
                for (const partition_t d: blocks) { counts[d] = 0; }
                return result;
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
                std::vector<u32> around;
                f64 old_pairs = 0.0;
                if (cfg.pair_objective) {
                    // every pair with an end in a PE whose load changes: those of a, c and the blocks u has edges to
                    around.assign(blocks.begin(), blocks.end());
                    around.push_back(a);
                    around.push_back(c);
                    std::sort(around.begin(), around.end());
                    around.erase(std::unique(around.begin(), around.end()), around.end());
                    old_pairs = pair_sum(around);
                }
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
                if (cfg.pair_objective) {
                    for (const u32 p: touched) { load[p] = compute_load(p); }
                    sum_power += pair_sum(around) - old_pairs;
                    return;
                }
                for (const u32 p: touched) {
                    sum_power -= omega(p) * powl_(load[p]);
                    load[p] = compute_load(p);
                    sum_power += omega(p) * powl_(load[p]);
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
                            if (!cfg.pair_objective) {
                                if (eval_move(u, c) < before * (1.0 - cfg.min_gain)) {
                                    move(u, c, touched);
                                    ++moves;
                                    break;
                                }
                                continue;
                            }
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
                        if (!cfg.pair_objective) {
                            // u goes to c once, every candidate v is only evaluated and the first improving one moves
                            move(u, c, touched);
                            for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1] && tried < cfg.swap_tries; ++j) {
                                const vertex_t v = g.edges_v[j];
                                if (part[v] != c) { continue; }
                                ++tried;
                                if (eval_move(v, a) < before * (1.0 - cfg.min_gain)) {
                                    move(v, a, touched);
                                    ++moves;
                                    done = true;
                                    break;
                                }
                            }
                            if (done) { break; }
                            move(u, a, touched);
                            continue;
                        }
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
