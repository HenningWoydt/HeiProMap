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

#ifndef HEIPROMAP_WAVER_ENGINE_REFINEMENT_H
#define HEIPROMAP_WAVER_ENGINE_REFINEMENT_H

#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../definitions.h"
#include "../datastructures/block_conn.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/large_quotient_graph.h"
#include "../datastructures/partition_manager.h"
#include "../utility/engine_load.h"
#include "../utility/profiler.h"

namespace HeiProMap {

    struct WaverEngineRefinementConfiguration {
        std::string name;
        bool enabled = false;
        u64 max_passes = 30;  // passes over the busiest PEs
        f64 ratio = 0.9;      // PEs with at least this fraction of the maximum engine load are worked on
        bool verbose = false;

        // Pair-volume mode: minimizes the total engine work plus pair_penalty * (|S_pq| - pair_threshold)^2 for
        // every directed pair (p, q) whose message carries more than pair_threshold values. The last messages of the
        // critical chains of the schedules are large ones, exchanged in sequence between two neighboring PEs.
        f64 pair_penalty = 0.0; // 0 = off
        u64 pair_threshold = 60;
        u64 pair_passes = 10;

        explicit WaverEngineRefinementConfiguration(std::string n) : name(std::move(n)) {}
    };

    // Lowers the load of the busiest compute engine (see engine_load.h), which bounds the makespan from below.
    //
    // Single-vertex moves out of the PEs with a load close to the maximum. A move of u from A to C changes the
    // rows and local nonzeros of A and C, and the sets of rows that are sent: u itself, and every neighbor v of u
    // that starts or stops having a neighbor in A or C. All of those PEs are re-evaluated; the move is made if it
    // lowers the largest load among them. Then the sorted vector of all loads decreases lexicographically, so the
    // search ends.
    class WaverEngineRefinement {
        using graph_t = CSRGraph;
        using p_manager_t = PartitionManager;
        using bv_manager_t = BoundaryVertexManager;
        using block_conn_t = BlockConn;
        using q_graph_t = LargeQuotientGraph;

        struct Change {
            partition_t p, q; // the set of rows of p with a neighbor in q
            int delta;        // +1: a row joins, -1: a row leaves
        };

        // The state the engine loads are computed from.
        struct State {
            EngineCosts c;
            std::vector<u64> rows, local_nnz;
            std::vector<u32> receives;
            std::vector<u64> sent_values, received_values;
            std::vector<std::vector<std::pair<partition_t, u32>>> out; // out[p]: (q, |S_pq|)

            f64 load(const partition_t p) const {
                return c.startup_fixed + c.row * static_cast<f64>(rows[p]) + c.local_nnz * static_cast<f64>(local_nnz[p])
                       + c.send_msg * static_cast<f64>(out[p].size()) + c.send_value * static_cast<f64>(sent_values[p])
                       + c.recv_msg * static_cast<f64>(receives[p]) + c.recv_value * static_cast<f64>(received_values[p]);
            }

            void add(const partition_t p, const partition_t q, const int delta) {
                auto &o = out[p];
                size_t i = 0;
                while (i < o.size() && o[i].first != q) { ++i; }
                if (i == o.size()) {
                    o.emplace_back(q, 0u); // delta > 0 here
                    ++receives[q];
                }
                o[i].second = static_cast<u32>(static_cast<int>(o[i].second) + delta);
                sent_values[p] = static_cast<u64>(static_cast<long long>(sent_values[p]) + delta);
                received_values[q] = static_cast<u64>(static_cast<long long>(received_values[q]) + delta);
                if (o[i].second == 0) {
                    o[i] = o.back();
                    o.pop_back();
                    --receives[q];
                }
            }
        };

        static weight_t conn_to(const block_conn_t &block_conn, const vertex_t v, const partition_t b) {
            for (size_t j = block_conn.start(v); j < block_conn.end(v); ++j) {
                if (block_conn.get_id(j) == b) { return block_conn.get_w(j); }
            }
            return 0;
        }

        // The effect of moving u from A to C on the sets of sent rows (changes) and on rows / local nonzeros.
        static void collect_changes(const graph_t &g, const p_manager_t &p_manager, const block_conn_t &block_conn,
                                    const vertex_t u, const partition_t A, const partition_t C,
                                    std::vector<Change> &changes, u64 &n_in_A, u64 &n_in_C) {
            changes.clear();
            n_in_A = n_in_C = 0;
            for (size_t j = block_conn.start(u); j < block_conn.end(u); ++j) {
                const partition_t b = block_conn.get_id(j);
                if (b != A) { changes.push_back({A, b, -1}); }
                if (b != C) { changes.push_back({C, b, +1}); }
            }
            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t v = g.edges_v[i];
                const weight_t w = g.edge_w(i);
                const partition_t D = p_manager[v];
                n_in_A += D == A;
                n_in_C += D == C;
                if (D != A && conn_to(block_conn, v, A) - w == 0) { changes.push_back({D, A, -1}); }
                if (D != C && conn_to(block_conn, v, C) == 0) { changes.push_back({D, C, +1}); }
            }
        }

        static void apply_changes(State &s, const std::vector<Change> &changes, const partition_t A, const partition_t C,
                                  const u64 vw, const u64 n_in_A, const u64 n_in_C, const int sign) {
            for (const Change &ch : changes) { s.add(ch.p, ch.q, sign * ch.delta); }
            const auto d = static_cast<long long>(sign);
            s.rows[A] = static_cast<u64>(static_cast<long long>(s.rows[A]) - d * static_cast<long long>(vw));
            s.rows[C] = static_cast<u64>(static_cast<long long>(s.rows[C]) + d * static_cast<long long>(vw));
            s.local_nnz[A] = static_cast<u64>(static_cast<long long>(s.local_nnz[A]) - d * static_cast<long long>(1 + 2 * n_in_A));
            s.local_nnz[C] = static_cast<u64>(static_cast<long long>(s.local_nnz[C]) + d * static_cast<long long>(1 + 2 * n_in_C));
        }

        static State build_state(const graph_t &g, const p_manager_t &p_manager, const block_conn_t &block_conn,
                                 const partition_t k) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverEngineRefinement", "setup");
            State s;
            s.rows.assign(k, 0);
            s.local_nnz.assign(k, 0);
            s.receives.assign(k, 0);
            s.sent_values.assign(k, 0);
            s.received_values.assign(k, 0);
            s.out.assign(k, {});
            for (vertex_t u = 0; u < g.n; ++u) {
                const partition_t a = p_manager[u];
                s.rows[a] += static_cast<u64>(g.v_weight(u));
                ++s.local_nnz[a];
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    s.local_nnz[a] += p_manager[g.edges_v[i]] == a;
                }
                for (size_t j = block_conn.start(u); j < block_conn.end(u); ++j) {
                    const partition_t b = block_conn.get_id(j);
                    if (b != a) { s.add(a, b, +1); }
                }
            }
            return s;
        }

        static u32 pair_count(const State &s, const partition_t p, const partition_t q) {
            for (const auto &e : s.out[p]) {
                if (e.first == q) { return e.second; }
            }
            return 0;
        }

    public:
        static void refine(const graph_t &g,
                           p_manager_t &p_manager,
                           bv_manager_t &bv_manager,
                           q_graph_t &q_graph,
                           block_conn_t &block_conn,
                           const WaverEngineRefinementConfiguration &config) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverEngineRefinement", "refine");

            const partition_t k = p_manager.get_k();
            State s = build_state(g, p_manager, block_conn, k);

            std::vector<f64> loads(k);
            std::vector<partition_t> order(k);
            std::vector<Change> changes, best_changes;
            std::vector<vertex_t> candidates;
            std::vector<partition_t> touched;
            std::vector<u64> stamp(k, 0);
            u64 stamp_now = 0;
            u64 total_moves = 0;
            f64 first_max = 0.0, last_max = 0.0;

            for (u64 pass = 0; pass < config.max_passes; ++pass) {
                f64 max_load = 0.0;
                for (partition_t p = 0; p < k; ++p) {
                    loads[p] = s.load(p);
                    max_load = std::max(max_load, loads[p]);
                }
                if (pass == 0) { first_max = max_load; }
                last_max = max_load;

                order.clear();
                for (partition_t p = 0; p < k; ++p) {
                    if (loads[p] >= config.ratio * max_load) { order.push_back(p); }
                }
                std::sort(order.begin(), order.end(), [&](const partition_t a, const partition_t b) {
                    return loads[a] != loads[b] ? loads[a] > loads[b] : a < b;
                });

                u64 moves = 0;
                for (const partition_t A : order) {
                    // The vertices of A that have a neighbor outside; moves change the boundary, so copy it.
                    candidates.clear();
                    for (size_t i = 0; i < bv_manager.size(A); ++i) { candidates.push_back(bv_manager.get(A, i)); }

                    for (const vertex_t u : candidates) {
                        if (p_manager[u] != A) { continue; }
                        const f64 load_A = s.load(A);
                        if (load_A < config.ratio * max_load) { break; }
                        const weight_t vw = g.v_weight(u);

                        partition_t best_C = NO_ID;
                        f64 best_new = load_A - 1e-6; // the largest load of the touched PEs must end up below this
                        u64 best_nA = 0, best_nC = 0;

                        for (size_t j = block_conn.start(u); j < block_conn.end(u); ++j) {
                            const partition_t C = block_conn.get_id(j);
                            if (C == A) { continue; }
                            if (p_manager.bweights[C] + vw > p_manager.get_lmax(C)) { continue; }
                            if (p_manager.bweights[A] - vw < p_manager.get_lmin(A)) { continue; }

                            u64 nA = 0, nC = 0;
                            collect_changes(g, p_manager, block_conn, u, A, C, changes, nA, nC);

                            // Largest load among the touched PEs before the move ...
                            ++stamp_now;
                            touched.clear();
                            auto touch = [&](const partition_t p) {
                                if (stamp[p] != stamp_now) { stamp[p] = stamp_now; touched.push_back(p); }
                            };
                            touch(A);
                            touch(C);
                            for (const Change &ch : changes) { touch(ch.p); touch(ch.q); }
                            f64 old_max = 0.0;
                            for (const partition_t p : touched) { old_max = std::max(old_max, s.load(p)); }

                            // ... and after it.
                            apply_changes(s, changes, A, C, static_cast<u64>(vw), nA, nC, +1);
                            f64 new_max = 0.0;
                            for (const partition_t p : touched) { new_max = std::max(new_max, s.load(p)); }
                            apply_changes(s, changes, A, C, static_cast<u64>(vw), nA, nC, -1);

                            if (new_max < old_max - 1e-6 && new_max < best_new) {
                                best_new = new_max;
                                best_C = C;
                                best_nA = nA;
                                best_nC = nC;
                                best_changes = changes;
                            }
                        }

                        if (best_C == NO_ID) { continue; }
                        apply_changes(s, best_changes, A, best_C, static_cast<u64>(vw), best_nA, best_nC, +1);
                        q_graph.move(g, p_manager, u, A, best_C);
                        block_conn.move(g, u, A, best_C);
                        bv_manager.move(g, p_manager, u, A, best_C);
                        p_manager.move_serial(u, vw, A, best_C);
                        ++moves;
                    }
                }
                total_moves += moves;
                if (config.verbose) {
                    std::cout << "[WaverEngineRefinement] pass " << pass << ": max load " << max_load << ", PEs worked on "
                              << order.size() << ", moves " << moves << std::endl;
                }
                if (moves == 0) { break; }
            }

            if (config.verbose) {
                f64 final_max = 0.0;
                for (partition_t p = 0; p < k; ++p) { final_max = std::max(final_max, s.load(p)); }
                std::cout << "[WaverEngineRefinement] max load " << first_max << " -> " << final_max << " in "
                          << total_moves << " moves" << std::endl;
            }
            (void) last_max;
        }

        // Pair-volume mode (see the configuration). Moves a vertex of a block with an oversized outgoing message
        // if that lowers (sum of the engine loads of the touched PEs) + (penalty of the touched pairs).
        static void refine_pairs(const graph_t &g,
                                 p_manager_t &p_manager,
                                 bv_manager_t &bv_manager,
                                 q_graph_t &q_graph,
                                 block_conn_t &block_conn,
                                 const WaverEngineRefinementConfiguration &config) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverEngineRefinement", "refine_pairs");

            const partition_t k = p_manager.get_k();
            State s = build_state(g, p_manager, block_conn, k);
            const f64 lambda = config.pair_penalty;
            const auto threshold = static_cast<f64>(config.pair_threshold);
            auto penalty = [&](const u32 count) {
                const f64 over = static_cast<f64>(count) - threshold;
                return over > 0.0 ? lambda * over * over : 0.0;
            };

            std::vector<Change> changes, best_changes;
            std::vector<vertex_t> candidates;
            std::vector<partition_t> touched, big_blocks;
            std::vector<std::pair<partition_t, partition_t>> pairs;
            std::vector<u64> stamp(k, 0);
            u64 stamp_now = 0;
            u64 total_moves = 0;
            f64 first_penalty = -1.0, last_penalty = 0.0;

            // Sum over the unique pairs of the changes of the penalty of their current count.
            auto pairs_penalty = [&]() {
                f64 sum = 0.0;
                for (const auto &pq : pairs) { sum += penalty(pair_count(s, pq.first, pq.second)); }
                return sum;
            };

            for (u64 pass = 0; pass < config.pair_passes; ++pass) {
                big_blocks.clear();
                f64 total_penalty = 0.0;
                u64 big_pairs = 0;
                for (partition_t p = 0; p < k; ++p) {
                    bool big = false;
                    for (const auto &e : s.out[p]) {
                        if (static_cast<f64>(e.second) > threshold) { big = true; ++big_pairs; total_penalty += penalty(e.second); }
                    }
                    if (big) { big_blocks.push_back(p); }
                }
                if (first_penalty < 0.0) { first_penalty = total_penalty; }
                last_penalty = total_penalty;
                if (big_blocks.empty()) { break; }

                u64 moves = 0;
                for (const partition_t A : big_blocks) {
                    candidates.clear();
                    for (size_t i = 0; i < bv_manager.size(A); ++i) { candidates.push_back(bv_manager.get(A, i)); }

                    for (const vertex_t u : candidates) {
                        if (p_manager[u] != A) { continue; }
                        bool in_big_pair = false;
                        for (size_t j = block_conn.start(u); j < block_conn.end(u) && !in_big_pair; ++j) {
                            const partition_t q = block_conn.get_id(j);
                            in_big_pair = q != A && static_cast<f64>(pair_count(s, A, q)) > threshold;
                        }
                        if (!in_big_pair) { continue; }
                        const weight_t vw = g.v_weight(u);

                        partition_t best_C = NO_ID;
                        f64 best_delta = -1e-6;
                        u64 best_nA = 0, best_nC = 0;
                        for (size_t j = block_conn.start(u); j < block_conn.end(u); ++j) {
                            const partition_t C = block_conn.get_id(j);
                            if (C == A) { continue; }
                            if (p_manager.bweights[C] + vw > p_manager.get_lmax(C)) { continue; }
                            if (p_manager.bweights[A] - vw < p_manager.get_lmin(A)) { continue; }

                            u64 nA = 0, nC = 0;
                            collect_changes(g, p_manager, block_conn, u, A, C, changes, nA, nC);

                            ++stamp_now;
                            touched.clear();
                            auto touch = [&](const partition_t p) {
                                if (stamp[p] != stamp_now) { stamp[p] = stamp_now; touched.push_back(p); }
                            };
                            touch(A);
                            touch(C);
                            pairs.clear();
                            for (const Change &ch : changes) {
                                touch(ch.p);
                                touch(ch.q);
                                pairs.emplace_back(ch.p, ch.q);
                            }
                            std::sort(pairs.begin(), pairs.end());
                            pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());

                            f64 before = pairs_penalty();
                            for (const partition_t p : touched) { before += s.load(p); }
                            apply_changes(s, changes, A, C, static_cast<u64>(vw), nA, nC, +1);
                            f64 after = pairs_penalty();
                            for (const partition_t p : touched) { after += s.load(p); }
                            apply_changes(s, changes, A, C, static_cast<u64>(vw), nA, nC, -1);

                            if (after - before < best_delta) {
                                best_delta = after - before;
                                best_C = C;
                                best_nA = nA;
                                best_nC = nC;
                                best_changes = changes;
                            }
                        }

                        if (best_C == NO_ID) { continue; }
                        apply_changes(s, best_changes, A, best_C, static_cast<u64>(vw), best_nA, best_nC, +1);
                        q_graph.move(g, p_manager, u, A, best_C);
                        block_conn.move(g, u, A, best_C);
                        bv_manager.move(g, p_manager, u, A, best_C);
                        p_manager.move_serial(u, vw, A, best_C);
                        ++moves;
                    }
                }
                total_moves += moves;
                if (config.verbose) {
                    std::cout << "[WaverEngineRefinement] pair pass " << pass << ": " << big_pairs << " pairs above "
                              << config.pair_threshold << " values, penalty " << total_penalty << ", moves " << moves
                              << std::endl;
                }
                if (moves == 0) { break; }
            }
            if (config.verbose) {
                std::cout << "[WaverEngineRefinement] pair penalty " << first_penalty << " -> " << last_penalty
                          << " (before the last pass) in " << total_moves << " moves" << std::endl;
            }
        }
    };
}

#endif //HEIPROMAP_WAVER_ENGINE_REFINEMENT_H
