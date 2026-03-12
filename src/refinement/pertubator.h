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

#ifndef HEIPROMAP_PERTURBATOR_H
#define HEIPROMAP_PERTURBATOR_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "../definitions.h"
#include "../definitions_1.h"
#include "../definitions_2.h"
#include "../definitions_3.h"
#include "../utility/qap.h"
#include "../utility/random_engine.h"

namespace HeiProMap {
    namespace perturbator_detail {
        struct PertMove {
            vertex_t u;
            partition_t from;
            partition_t to;
            weight_t weight;
            s64 gain; // positive = improvement, negative = deterioration
        };

        inline std::vector<weight_t> compute_block_weights(const graph_t &g,
                                                           bv_manager_t &bv_manager,
                                                           p_manager_t &p_manager) {
            std::vector<weight_t> block_weights(p_manager.k, 0);

            for (partition_t b = 0; b < p_manager.k; ++b) {
                forall_bv_id_iu(bv_manager, b, i, u)
                    {
                        block_weights[b] += g.v_weights[u];
                    }
                endfor
            }

            return block_weights;
        }

        inline void apply_move(graph_t &g,
                               bv_manager_t &bv_manager,
                               p_manager_t &p_manager,
                               q_graph_t &q_graph,
                               block_conn_t &block_conn,
                               std::vector<weight_t> &block_weights,
                               const PertMove &m) {
            bv_manager.move(g, p_manager, m.u, m.from, m.to);
            q_graph.move(g, p_manager, m.u, m.from, m.to);
            block_conn.move(g, m.u, m.from, m.to);
            p_manager.move(m.u, m.weight, m.from, m.to);

            block_weights[m.from] -= m.weight;
            block_weights[m.to] += m.weight;
        }

        inline void undo_move(graph_t &g,
                              bv_manager_t &bv_manager,
                              p_manager_t &p_manager,
                              q_graph_t &q_graph,
                              block_conn_t &block_conn,
                              std::vector<weight_t> &block_weights,
                              const PertMove &m) {
            bv_manager.move(g, p_manager, m.u, m.to, m.from);
            q_graph.move(g, p_manager, m.u, m.to, m.from);
            block_conn.move(g, m.u, m.to, m.from);
            p_manager.move(m.u, m.weight, m.to, m.from);

            block_weights[m.to] -= m.weight;
            block_weights[m.from] += m.weight;
        }

        inline void collect_candidates_from_block(graph_t &g,
                                                  d_oracle_t &d_oracle,
                                                  bv_manager_t &bv_manager,
                                                  p_manager_t &p_manager,
                                                  partition_t from_block,
                                                  const std::vector<weight_t> &block_weights,
                                                  weight_t lmax,
                                                  const std::vector<u8> &used_vertices,
                                                  std::vector<PertMove> &out) {
            out.clear();

            std::vector<u8> seen_block(p_manager.k, false);

            forall_bv_id_iu(bv_manager, from_block, i, u)
                {
                    if (used_vertices[u]) {
                        continue;
                    }

                    const weight_t u_weight = g.v_weights[u];

                    std::fill(seen_block.begin(), seen_block.end(), false);

                    forall_guiv(g, u, j, v)
                        {
                            const partition_t to_block = p_manager[v];

                            if (to_block == from_block) {
                                continue;
                            }

                            if (seen_block[to_block]) {
                                continue;
                            }
                            seen_block[to_block] = true;

                            // Strict balance preservation:
                            // only consider moves that keep the target block feasible immediately.
                            if (block_weights[to_block] + u_weight > lmax) {
                                continue;
                            }

                            const s64 gain = get_u_qap_delta(g, u, from_block, to_block, p_manager, d_oracle);

                            out.push_back({u, from_block, to_block, u_weight, gain});
                        }
                    endfor
                }
            endfor
        }

        inline void shuffle_candidates(RandomEngine &random_engine,
                                       std::vector<PertMove> &moves) {
            if (moves.empty()) {
                return;
            }

            for (u64 i = moves.size() - 1; i > 0; --i) {
                const u64 j = random_engine.get_u64() % (i + 1);
                std::swap(moves[i], moves[j]);
            }
        }

        inline bool better_chain(s64 candidate_gain,
                                 const std::vector<PertMove> &candidate_chain,
                                 s64 best_gain,
                                 const std::vector<PertMove> &best_chain) {
            if (candidate_chain.empty()) {
                return false;
            }

            if (best_chain.empty()) {
                return true;
            }

            if (candidate_gain > best_gain) {
                return true;
            }

            if (candidate_gain < best_gain) {
                return false;
            }

            // Tie-breaker: prefer longer chains as stronger perturbations.
            return candidate_chain.size() > best_chain.size();
        }

        inline bool find_best_chain_from_start_block(graph_t &g,
                                                     d_oracle_t &d_oracle,
                                                     bv_manager_t &bv_manager,
                                                     p_manager_t &p_manager,
                                                     q_graph_t &q_graph,
                                                     block_conn_t &block_conn,
                                                     partition_t start_block,
                                                     weight_t lmax,
                                                     s64 max_total_loss,
                                                     u64 max_depth,
                                                     u64 max_branching,
                                                     RandomEngine &random_engine,
                                                     std::vector<PertMove> &best_chain_out) {
            std::vector<weight_t> block_weights = compute_block_weights(g, bv_manager, p_manager);

            // Assumption: vertex ids are in [0, g.n).
            // This is consistent with your graph structure style.
            std::vector<u8> used_vertices(g.n, false);

            std::vector<PertMove> current_chain;
            std::vector<PertMove> best_chain;
            std::vector<PertMove> candidates;

            s64 current_gain = 0;
            s64 best_gain = std::numeric_limits<s64>::lowest();

            auto dfs = [&](auto &&self, partition_t current_block, u64 depth) -> void {
                if (!current_chain.empty()) {
                    if (better_chain(current_gain, current_chain, best_gain, best_chain)) {
                        best_gain = current_gain;
                        best_chain = current_chain;
                    }
                }

                if (depth >= max_depth) {
                    return;
                }

                collect_candidates_from_block(g,
                                              d_oracle,
                                              bv_manager,
                                              p_manager,
                                              current_block,
                                              block_weights,
                                              lmax,
                                              used_vertices,
                                              candidates);

                if (candidates.empty()) {
                    return;
                }

                // Randomize first, then stable-sort by gain descending.
                // This gives randomized tie-breaking among equal-gain moves.
                shuffle_candidates(random_engine, candidates);

                std::stable_sort(candidates.begin(), candidates.end(),
                                 [](const PertMove &a, const PertMove &b) {
                                     return a.gain > b.gain;
                                 });

                const u64 limit = std::min<u64>(max_branching, candidates.size());

                for (u64 idx = 0; idx < limit; ++idx) {
                    const PertMove m = candidates[idx];

                    // Prevent chains whose accumulated loss becomes too large.
                    // gain >= 0 is good, gain < 0 is bad.
                    if (current_gain + m.gain < -max_total_loss) {
                        continue;
                    }

                    used_vertices[m.u] = true;
                    apply_move(g, bv_manager, p_manager, q_graph, block_conn, block_weights, m);
                    current_chain.push_back(m);
                    current_gain += m.gain;

                    self(self, m.to, depth + 1);

                    current_gain -= m.gain;
                    current_chain.pop_back();
                    undo_move(g, bv_manager, p_manager, q_graph, block_conn, block_weights, m);
                    used_vertices[m.u] = false;
                }
            };

            dfs(dfs, start_block, 0);

            if (best_chain.empty()) {
                return false;
            }

            best_chain_out = std::move(best_chain);
            return true;
        }
    } // namespace perturbator_detail

    inline void perturbate(graph_t &g,
                           d_oracle_t &d_oracle,
                           bv_manager_t &bv_manager,
                           p_manager_t &p_manager,
                           q_graph_t &q_graph,
                           block_conn_t &block_conn,
                           f64 imbalance) {
        ScopedTimer _t("refinement", "Pertubator", "perturbate");
        using namespace perturbator_detail;

        // Tunable parameters
        const u64 max_iterations = 8;
        const u64 max_depth = 4; // maximum chain length
        const u64 max_branching = 3; // number of candidate moves explored per depth
        RandomEngine random_engine(0);

        const weight_t lmax =
                static_cast<weight_t>(std::ceil((1.0 + imbalance) *
                                                (static_cast<f64>(g.g_weight) /
                                                 static_cast<f64>(p_manager.k))));

        // Allow a small objective deterioration over the whole chain.
        // This is intentionally conservative.
        //
        // Interpretation:
        //   total_chain_gain >= -max_total_loss
        //
        // You will likely want to tune this.
        const s64 max_total_loss =
                std::max<s64>(1, static_cast<s64>(std::ceil(0.01 * static_cast<f64>(g.g_weight))));

        // Random start offset to avoid always scanning blocks in the same order.
        partition_t start_offset = static_cast<partition_t>(random_engine.get_u64() % p_manager.k);

        for (u64 iter = 0; iter < max_iterations; ++iter) {
            std::vector<PertMove> best_chain_this_iteration;
            s64 best_gain_this_iteration = std::numeric_limits<s64>::lowest();

            // Try all blocks once in rotated order and keep the best chain.
            for (partition_t step = 0; step < p_manager.k; ++step) {
                const partition_t start_block =
                        static_cast<partition_t>((start_offset + step) % p_manager.k);

                std::vector<PertMove> candidate_chain;
                const bool found = find_best_chain_from_start_block(g,
                                                                    d_oracle,
                                                                    bv_manager,
                                                                    p_manager,
                                                                    q_graph,
                                                                    block_conn,
                                                                    start_block,
                                                                    lmax,
                                                                    max_total_loss,
                                                                    max_depth,
                                                                    max_branching,
                                                                    random_engine,
                                                                    candidate_chain);

                if (!found) {
                    continue;
                }

                s64 candidate_gain = 0;
                for (const PertMove &m: candidate_chain) {
                    candidate_gain += m.gain;
                }

                if (better_chain(candidate_gain,
                                 candidate_chain,
                                 best_gain_this_iteration,
                                 best_chain_this_iteration)) {
                    best_gain_this_iteration = candidate_gain;
                    best_chain_this_iteration = std::move(candidate_chain);
                }
            }

            if (best_chain_this_iteration.empty()) {
                break;
            }

            // Commit the best chain found in this iteration.
            for (const PertMove &m: best_chain_this_iteration) {
                bv_manager.move(g, p_manager, m.u, m.from, m.to);
                q_graph.move(g, p_manager, m.u, m.from, m.to);
                block_conn.move(g, m.u, m.from, m.to);
                p_manager.move(m.u, m.weight, m.from, m.to);
            }

            // Change block scan start for next iteration.
            start_offset = static_cast<partition_t>(random_engine.get_u64() % p_manager.k);
        }
    }
} // namespace HeiProMap

#endif // HEIPROMAP_PERTURBATOR_H
