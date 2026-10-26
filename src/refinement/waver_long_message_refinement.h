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

#ifndef HEIPROMAP_WAVER_LONG_MESSAGE_REFINEMENT_H
#define HEIPROMAP_WAVER_LONG_MESSAGE_REFINEMENT_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "../definitions.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../utility/profiler.h"

namespace HeiProMap {

    // Hyperparameters of the long-message refinement: the longest communication edges are shortened by swapping
    // an endpoint's block into a rectangle around the other endpoint.
    class WaverLongMessageRefinementConfiguration {
    public:
        explicit WaverLongMessageRefinementConfiguration(std::string t_name) : name(std::move(t_name)) {}

        std::string name;
        bool enabled = false;
        u64 padding = 1;                // half-side of the searched rectangle around the other endpoint
        u64 distance_percentile = 90;   // edges at least this long (percentile of all edge distances) count as long
        u64 max_rounds = 3000;          // rounds over the long edges; stops earlier once a round swaps nothing
        u64 max_edges = 5000;           // only the this many costliest long edges are considered
        f64 min_gain = 1e-6;            // a swap must improve the cost by more than this
    };

    class WaverLongMessageRefinement {
        using graph_t     = CSRGraph;
        using p_manager_t = PartitionManager;

    public:
        template<typename DistanceOracleT>
        // Returns the number of block swaps performed.
        static size_t refine(graph_t              &g,
                            p_manager_t          &p_manager,
                            const DistanceOracleT &d_oracle,
                            const WaverLongMessageRefinementConfiguration &config,
                            f64                   alpha_end,
                            f64                   beta_end,
                            f64                   alpha_link,
                            f64                   beta_link) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverLongMessageRefinement", "refine");
            const size_t W = static_cast<size_t>(d_oracle.current_W());
            const size_t H = (p_manager.get_k() + W - 1) / W;
            const partition_t k = p_manager.get_k();

            std::vector<std::vector<vertex_t>> blk_verts(k);
            for (vertex_t u = 0; u < g.n; ++u) {
                blk_verts[p_manager.partition[u]].push_back(u);
            }

            auto block_cost = [&](const std::vector<vertex_t> &verts, partition_t block) -> f64 {
                f64 cost = 0.0;
                for (const vertex_t u : verts) {
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const partition_t nb = p_manager.partition[g.edges_v[j]];
                        const f64 w = static_cast<f64>(g.edge_w(j));
                        const f64 dist = static_cast<f64>(d_oracle.get(block, nb));
                        if (block != nb) {
                            cost += alpha_end + beta_end * w;
                        }
                        cost += alpha_link * dist + beta_link * w * dist;
                    }
                }
                return cost;
            };

            // Collect long edges with accumulated weighted distance
            struct LongEdge {
                partition_t p, q;
                f64 weighted_dist;
            };

            // P90 distance as threshold
            std::vector<weight_t> all_dists;
            for (vertex_t u = 0; u < g.n; ++u) {
                const partition_t ub = p_manager.partition[u];
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const partition_t vb = p_manager.partition[g.edges_v[j]];
                    if (ub < vb) {
                        all_dists.push_back(d_oracle.get(ub, vb));
                    }
                }
            }
            if (all_dists.empty()) return 0;
            std::sort(all_dists.begin(), all_dists.end());
            const weight_t dist_threshold = all_dists[all_dists.size() * config.distance_percentile / 100];
            if (dist_threshold <= 1) return 0;

            std::vector<LongEdge> long_edges;
            std::unordered_map<u64, size_t> pair_to_idx;
            for (vertex_t u = 0; u < g.n; ++u) {
                const partition_t ub = p_manager.partition[u];
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const partition_t vb = p_manager.partition[g.edges_v[j]];
                    if (ub >= vb) continue;
                    const weight_t dist = d_oracle.get(ub, vb);
                    if (dist < dist_threshold) continue;
                    const u64 key = static_cast<u64>(ub) * k + vb;
                    const f64 w = static_cast<f64>(g.edge_w(j));
                    const f64 wd = static_cast<f64>(dist) * (alpha_link + beta_link * w);
                    auto it = pair_to_idx.find(key);
                    if (it == pair_to_idx.end()) {
                        pair_to_idx[key] = long_edges.size();
                        long_edges.push_back({ub, vb, wd});
                    } else {
                        long_edges[it->second].weighted_dist += wd;
                    }
                }
            }

            std::sort(long_edges.begin(), long_edges.end(),
                      [](const LongEdge &a, const LongEdge &b) {
                          return a.weighted_dist > b.weighted_dist;
                      });

            const size_t max_rounds = config.max_rounds;
            const size_t max_edges = std::min(long_edges.size(), static_cast<size_t>(config.max_edges));
            size_t total_swaps = 0;

            for (size_t round = 0; round < max_rounds; ++round) {
                size_t swaps = 0;

                for (size_t ei = 0; ei < max_edges; ++ei) {
                    const auto &le = long_edges[ei];
                    const partition_t pa = le.p;
                    const partition_t pb = le.q;
                    if (blk_verts[pa].empty() || blk_verts[pb].empty()) continue;
                    if (d_oracle.get(pa, pb) < dist_threshold) continue;

                    const size_t ax = pa % W, ay = pa / W;
                    const size_t bx = pb % W, by = pb / W;

                    // Try swapping each endpoint toward the other:
                    // search a rectangle of side 2*padding+1 centered on the OTHER endpoint
                    for (int side = 0; side < 2; ++side) {
                        const partition_t src = (side == 0) ? pa : pb;
                        const size_t cx = (side == 0) ? bx : ax;
                        const size_t cy = (side == 0) ? by : ay;

                        const size_t x_lo = (cx >= config.padding) ? cx - config.padding : 0;
                        const size_t y_lo = (cy >= config.padding) ? cy - config.padding : 0;
                        const size_t x_hi = std::min(cx + config.padding, W - 1);
                        const size_t y_hi = std::min(cy + config.padding, H - 1);

                        const f64 cost_src = block_cost(blk_verts[src], src);
                        partition_t best_q = NO_ID;
                        f64 best_gain = config.min_gain;

                        for (size_t y = y_lo; y <= y_hi; ++y) {
                            for (size_t x = x_lo; x <= x_hi; ++x) {
                                const partition_t q = static_cast<partition_t>(y * W + x);
                                if (q >= k || q == pa || q == pb) continue;
                                if (blk_verts[q].empty()) continue;
                                if (p_manager.bweights[q] > p_manager.get_lmax(src) ||
                                    p_manager.bweights[src] > p_manager.get_lmax(q)) continue;

                                const f64 cost_q = block_cost(blk_verts[q], q);
                                const f64 cost_before = cost_src + cost_q;

                                for (const vertex_t u : blk_verts[src]) p_manager.partition[u] = q;
                                for (const vertex_t v : blk_verts[q]) p_manager.partition[v] = src;

                                const f64 cost_after = block_cost(blk_verts[src], q) +
                                                       block_cost(blk_verts[q], src);

                                for (const vertex_t u : blk_verts[src]) p_manager.partition[u] = src;
                                for (const vertex_t v : blk_verts[q]) p_manager.partition[v] = q;

                                const f64 gain = cost_before - cost_after;
                                if (gain > best_gain) {
                                    best_gain = gain;
                                    best_q = q;
                                }
                            }
                        }

                        if (best_q != NO_ID) {
                            for (const vertex_t u : blk_verts[src]) p_manager.partition[u] = best_q;
                            for (const vertex_t v : blk_verts[best_q]) p_manager.partition[v] = src;
                            const weight_t ws = p_manager.bweights[src];
                            const weight_t wq = p_manager.bweights[best_q];
                            p_manager.bweights[src] = wq;
                            p_manager.bweights[best_q] = ws;
                            std::swap(blk_verts[src], blk_verts[best_q]);
                            swaps++;
                            break;
                        }
                    }
                }

                total_swaps += swaps;
                if (swaps == 0) break;
            }
            return total_swaps;
        }
    };

} // namespace HeiProMap

#endif //HEIPROMAP_WAVER_LONG_MESSAGE_REFINEMENT_H
