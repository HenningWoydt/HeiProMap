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

// Long-message refinement (Waver). Shortens the longest messages by swapping whole blocks (the contents of two
// PEs) on the PE grid: for the cut edges whose endpoints are farthest apart (a percentile of the distances), each endpoint's
// block is tried at every position in a rectangle around the other endpoint, and the swap is made if the message cost drops.
// Only blocks change places, so the contents and the edges between them never change; the cost of a position needs one
// distance per neighbor content. Changes the PE labels (p_manager only), so it must run after the refiners that need
// the boundary and quotient structures. Option: --long-message-padding.


#ifndef HEIPROMAP_WAVER_LONG_MESSAGE_REFINEMENT_H
#define HEIPROMAP_WAVER_LONG_MESSAGE_REFINEMENT_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
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

            // Only whole blocks change places, so the blocks (the "contents", named by their block at the start) and
            // the edges between them never change: for every content the neighbor contents with the number of edges
            // and their weight sum. A block's cost then needs one distance per neighbor instead of one per edge.
            struct Nb {
                partition_t c;
                f64 cnt;   // edges to the neighbor content
                f64 wsum;  // their weight sum
            };
            std::vector<std::vector<Nb>> nbs(k);
            std::vector<u64> n_in(k, 0);
            for (vertex_t u = 0; u < g.n; ++u) {
                const partition_t cu = p_manager.partition[u];
                ++n_in[cu];
                auto &lst = nbs[cu];
                size_t last = 0; // consecutive edges of a vertex often go to the same block
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const partition_t cv = p_manager.partition[g.edges_v[j]];
                    if (cv == cu) continue;
                    const f64 w = static_cast<f64>(g.edge_w(j));
                    size_t i = last;
                    if (i >= lst.size() || lst[i].c != cv) {
                        for (i = 0; i < lst.size() && lst[i].c != cv; ++i) {}
                        if (i == lst.size()) lst.push_back({cv, 0.0, 0.0});
                    }
                    lst[i].cnt += 1.0;
                    lst[i].wsum += w;
                    last = i;
                }
            }

            // Edge distance threshold: the distance_percentile of the distances of all (undirected) cut edges, from
            // a histogram. A cut edge between the blocks c < c' is counted once, from the side of c.
            std::vector<u64> hist;
            u64 n_cut = 0;
            for (partition_t c = 0; c < k; ++c) {
                for (const Nb &nb : nbs[c]) {
                    if (nb.c < c) continue;
                    const auto d = static_cast<size_t>(d_oracle.get(c, nb.c));
                    if (d >= hist.size()) hist.resize(d + 1, 0);
                    hist[d] += static_cast<u64>(nb.cnt);
                    n_cut += static_cast<u64>(nb.cnt);
                }
            }
            if (n_cut == 0) return 0;
            const u64 rank = n_cut * config.distance_percentile / 100; // index into the sorted distances
            weight_t dist_threshold = 0;
            for (u64 d = 0, seen = 0; d < hist.size(); ++d) {
                seen += hist[d];
                if (seen > rank) { dist_threshold = static_cast<weight_t>(d); break; }
            }
            if (dist_threshold <= 1) return 0;

            // Long edges: pairs of blocks at least that far apart, with their accumulated weighted distance.
            struct LongEdge {
                partition_t p, q;
                f64 weighted_dist;
            };
            std::vector<LongEdge> long_edges;
            for (partition_t c = 0; c < k; ++c) {
                for (const Nb &nb : nbs[c]) {
                    if (nb.c < c) continue;
                    const weight_t dist = d_oracle.get(c, nb.c);
                    if (dist < dist_threshold) continue;
                    long_edges.push_back({c, nb.c, static_cast<f64>(dist) * (alpha_link * nb.cnt + beta_link * nb.wsum)});
                }
            }
            std::sort(long_edges.begin(), long_edges.end(),
                      [](const LongEdge &a, const LongEdge &b) {
                          if (a.weighted_dist != b.weighted_dist) return a.weighted_dist > b.weighted_dist;
                          return a.p != b.p ? a.p < b.p : a.q < b.q;
                      });

            // The placement: content at[p] sits at block (position) p, content c at pos[c]; bw[p] is the weight there.
            std::vector<partition_t> at(k), pos(k);
            std::vector<weight_t> bw(k);
            for (partition_t p = 0; p < k; ++p) {
                at[p] = p;
                pos[p] = p;
                bw[p] = p_manager.bweights[p];
            }

            // Cost of content c placed at position P, while the content `moved` is placed at `moved_pos` instead of
            // pos[moved] (pass NO_ID for no change).
            auto cost_of = [&](const partition_t c, const partition_t P, const partition_t moved,
                               const partition_t moved_pos) -> f64 {
                f64 cost = 0.0;
                for (const Nb &nb : nbs[c]) {
                    const partition_t Pn = nb.c == moved ? moved_pos : pos[nb.c];
                    const f64 dist = static_cast<f64>(d_oracle.get(P, Pn));
                    cost += nb.cnt * (alpha_end + alpha_link * dist) + nb.wsum * (beta_end + beta_link * dist);
                }
                return cost;
            };

            const size_t max_rounds = config.max_rounds;
            const size_t max_edges = std::min(long_edges.size(), static_cast<size_t>(config.max_edges));
            size_t total_swaps = 0;

            for (size_t round = 0; round < max_rounds; ++round) {
                size_t swaps = 0;

                for (size_t ei = 0; ei < max_edges; ++ei) {
                    const auto &le = long_edges[ei];
                    const partition_t pa = le.p;
                    const partition_t pb = le.q;
                    if (n_in[at[pa]] == 0 || n_in[at[pb]] == 0) continue;
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

                        const partition_t c_src = at[src];
                        const f64 cost_src = cost_of(c_src, src, NO_ID, NO_ID);
                        partition_t best_q = NO_ID;
                        f64 best_gain = config.min_gain;

                        for (size_t y = y_lo; y <= y_hi; ++y) {
                            for (size_t x = x_lo; x <= x_hi; ++x) {
                                const partition_t q = static_cast<partition_t>(y * W + x);
                                if (q >= k || q == pa || q == pb) continue;
                                const partition_t c_q = at[q];
                                if (n_in[c_q] == 0) continue;
                                if (bw[q] > p_manager.get_lmax(src) || bw[src] > p_manager.get_lmax(q)) continue;

                                const f64 cost_before = cost_src + cost_of(c_q, q, NO_ID, NO_ID);
                                // the two contents trade places
                                const f64 cost_after = cost_of(c_src, q, c_q, src) + cost_of(c_q, src, c_src, q);

                                const f64 gain = cost_before - cost_after;
                                if (gain > best_gain) {
                                    best_gain = gain;
                                    best_q = q;
                                }
                            }
                        }

                        if (best_q != NO_ID) {
                            const partition_t c_q = at[best_q];
                            at[src] = c_q;
                            at[best_q] = c_src;
                            pos[c_src] = best_q;
                            pos[c_q] = src;
                            std::swap(bw[src], bw[best_q]);
                            swaps++;
                            break;
                        }
                    }
                }

                total_swaps += swaps;
                if (swaps == 0) break;
            }

            if (total_swaps > 0) {
                for (vertex_t u = 0; u < g.n; ++u) { p_manager.partition[u] = pos[p_manager.partition[u]]; }
                for (partition_t p = 0; p < k; ++p) { p_manager.bweights[p] = bw[p]; }
            }
            return total_swaps;
        }
    };

} // namespace HeiProMap

#endif //HEIPROMAP_WAVER_LONG_MESSAGE_REFINEMENT_H
