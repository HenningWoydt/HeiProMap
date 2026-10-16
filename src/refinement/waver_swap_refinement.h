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

#ifndef HEIPROMAP_WAVER_SWAP_REFINEMENT_H
#define HEIPROMAP_WAVER_SWAP_REFINEMENT_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "../definitions.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../utility/profiler.h"

namespace HeiProMap {

    class WaverSwapRefinement {
        using graph_t     = CSRGraph;
        using p_manager_t = PartitionManager;

        static void add_radius_nbrs(size_t cx, size_t cy, size_t r,
                                    partition_t p, size_t W, size_t H, partition_t k,
                                    std::vector<partition_t> &nbrs,
                                    std::vector<bool> &seen) {
            for (size_t dy = 0; dy <= r; ++dy) {
                const size_t dx_max = r - dy;
                for (size_t dx = 0; dx <= dx_max; ++dx) {
                    if (dx == 0 && dy == 0) continue;
                    size_t coords[4][2];
                    size_t n_coords = 0;
                    if (cx + dx < W && cy + dy < H) coords[n_coords][0] = cx + dx, coords[n_coords][1] = cy + dy, n_coords++;
                    if (dx > 0 && cx >= dx && cy + dy < H) coords[n_coords][0] = cx - dx, coords[n_coords][1] = cy + dy, n_coords++;
                    if (dy > 0 && cx + dx < W && cy >= dy) coords[n_coords][0] = cx + dx, coords[n_coords][1] = cy - dy, n_coords++;
                    if (dx > 0 && dy > 0 && cx >= dx && cy >= dy) coords[n_coords][0] = cx - dx, coords[n_coords][1] = cy - dy, n_coords++;
                    if (dx == 0 && cy + dy < H) { coords[0][0] = cx; coords[0][1] = cy + dy; n_coords = 1; }
                    if (dx == 0 && dy > 0 && cy >= dy) { coords[n_coords][0] = cx; coords[n_coords][1] = cy - dy; n_coords++; }
                    for (size_t c = 0; c < n_coords; ++c) {
                        partition_t q = static_cast<partition_t>(coords[c][1] * W + coords[c][0]);
                        if (q < k && q != p && !seen[q]) {
                            seen[q] = true;
                            nbrs.push_back(q);
                        }
                    }
                }
            }
        }

    public:
        template<typename DistanceOracleT>
        static void swap_refine(graph_t              &g,
                                p_manager_t          &p_manager,
                                const DistanceOracleT &d_oracle,
                                size_t                neighbor_radius,
                                bool                  extended,
                                f64                   alpha_end,
                                f64                   beta_end,
                                f64                   alpha_link,
                                f64                   beta_link) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverSwapRefinement", "swap_refine");
            const size_t W = static_cast<size_t>(d_oracle.current_W());
            const size_t H = (p_manager.get_k() + W - 1) / W;
            const partition_t k = p_manager.get_k();

            std::vector<std::vector<vertex_t>> blk_verts(k);
            for (vertex_t u = 0; u < g.n; ++u) {
                blk_verts[p_manager.partition[u]].push_back(u);
            }

            auto get_candidates = [&](partition_t p, std::vector<partition_t> &nbrs,
                                      std::vector<bool> &seen) {
                nbrs.clear();
                add_radius_nbrs(p % W, p / W, neighbor_radius, p, W, H, k, nbrs, seen);

                if (extended) {
                    for (const vertex_t u : blk_verts[p]) {
                        for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                            const partition_t nb = p_manager.partition[g.edges_v[j]];
                            if (nb != p) {
                                add_radius_nbrs(nb % W, nb / W, neighbor_radius, p, W, H, k, nbrs, seen);
                            }
                        }
                    }
                }

                for (const partition_t q : nbrs) seen[q] = false;
            };

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

            const size_t max_rounds = 100;
            size_t total_swaps = 0;
            bool progress = true;
            size_t round = 0;
            std::vector<partition_t> nbrs;
            std::vector<bool> seen(k, false);
            while (progress && round < max_rounds) {
                progress = false;
                size_t swaps = 0;
                for (partition_t p = 0; p < k; ++p) {
                    if (blk_verts[p].empty()) continue;

                    get_candidates(p, nbrs, seen);

                    partition_t best_q = NO_ID;
                    f64 best_gain = 1e-6;
                    const f64 cost_p = block_cost(blk_verts[p], p);

                    for (const partition_t q : nbrs) {
                        if (blk_verts[q].empty()) continue;

                        const weight_t wp = p_manager.bweights[p];
                        const weight_t wq = p_manager.bweights[q];
                        if (wq > p_manager.get_lmax(p) || wp > p_manager.get_lmax(q)) continue;

                        const f64 cost_q = block_cost(blk_verts[q], q);
                        const f64 cost_before = cost_p + cost_q;

                        for (const vertex_t u : blk_verts[p]) p_manager.partition[u] = q;
                        for (const vertex_t v : blk_verts[q]) p_manager.partition[v] = p;

                        const f64 cost_after = block_cost(blk_verts[p], q) + block_cost(blk_verts[q], p);

                        for (const vertex_t u : blk_verts[p]) p_manager.partition[u] = p;
                        for (const vertex_t v : blk_verts[q]) p_manager.partition[v] = q;

                        const f64 gain = cost_before - cost_after;
                        if (gain > best_gain) {
                            best_gain = gain;
                            best_q = q;
                        }
                    }

                    if (best_q != NO_ID) {
                        for (const vertex_t u : blk_verts[p]) p_manager.partition[u] = best_q;
                        for (const vertex_t v : blk_verts[best_q]) p_manager.partition[v] = p;
                        const weight_t wp = p_manager.bweights[p];
                        const weight_t wq = p_manager.bweights[best_q];
                        p_manager.bweights[p] = wq;
                        p_manager.bweights[best_q] = wp;
                        std::swap(blk_verts[p], blk_verts[best_q]);
                        swaps++;
                        progress = true;
                    }
                }
                total_swaps += swaps;
                round++;
            }
        }

        template<typename DistanceOracleT>
        static void vertex_move_refine(graph_t              &g,
                                       p_manager_t          &p_manager,
                                       const DistanceOracleT &d_oracle,
                                       size_t                neighbor_radius,
                                       bool                  extended,
                                       f64                   alpha_end,
                                       f64                   beta_end,
                                       f64                   alpha_link,
                                       f64                   beta_link) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "WaverSwapRefinement", "vertex_move_refine");
            const size_t W = static_cast<size_t>(d_oracle.current_W());
            const size_t H = (p_manager.get_k() + W - 1) / W;
            const partition_t k = p_manager.get_k();

            std::vector<std::vector<vertex_t>> blk_verts(k);
            for (vertex_t u = 0; u < g.n; ++u) {
                blk_verts[p_manager.partition[u]].push_back(u);
            }

            auto get_candidates = [&](partition_t p, vertex_t u,
                                      std::vector<partition_t> &nbrs,
                                      std::vector<bool> &seen) {
                nbrs.clear();
                add_radius_nbrs(p % W, p / W, neighbor_radius, p, W, H, k, nbrs, seen);

                if (extended) {
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const partition_t nb = p_manager.partition[g.edges_v[j]];
                        if (nb != p) {
                            add_radius_nbrs(nb % W, nb / W, neighbor_radius, p, W, H, k, nbrs, seen);
                        }
                    }
                }

                for (const partition_t q : nbrs) seen[q] = false;
            };

            auto vertex_cost = [&](vertex_t u, partition_t block) -> f64 {
                f64 cost = 0.0;
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const partition_t nb = p_manager.partition[g.edges_v[j]];
                    const f64 w = static_cast<f64>(g.edge_w(j));
                    const f64 dist = static_cast<f64>(d_oracle.get(block, nb));
                    if (block != nb) {
                        cost += alpha_end + beta_end * w;
                    }
                    cost += alpha_link * dist + beta_link * w * dist;
                }
                return cost;
            };

            const size_t max_rounds = 100;
            size_t total_moves = 0;
            bool progress = true;
            size_t round = 0;
            std::vector<partition_t> nbrs;
            std::vector<bool> seen(k, false);
            while (progress && round < max_rounds) {
                progress = false;
                size_t moves = 0;
                for (partition_t p = 0; p < k; ++p) {
                    if (blk_verts[p].size() <= 1) continue;

                    for (size_t vi = 0; vi < blk_verts[p].size(); ++vi) {
                        const vertex_t u = blk_verts[p][vi];
                        const weight_t u_w = g.v_weight(u);

                        if (p_manager.bweights[p] - u_w < p_manager.get_lmin(p)) continue;

                        get_candidates(p, u, nbrs, seen);

                        const f64 cost_here = vertex_cost(u, p);
                        partition_t best_q = NO_ID;
                        f64 best_gain = 1e-6;

                        for (const partition_t q : nbrs) {
                            if (p_manager.bweights[q] + u_w > p_manager.get_lmax(q)) continue;

                            const f64 cost_there = vertex_cost(u, q);
                            const f64 gain = cost_here - cost_there;
                            if (gain > best_gain) {
                                best_gain = gain;
                                best_q = q;
                            }
                        }

                        if (best_q != NO_ID) {
                            p_manager.partition[u] = best_q;
                            p_manager.bweights[p] -= u_w;
                            p_manager.bweights[best_q] += u_w;
                            blk_verts[best_q].push_back(u);
                            blk_verts[p][vi] = blk_verts[p].back();
                            blk_verts[p].pop_back();
                            --vi;
                            moves++;
                            progress = true;
                        }
                    }
                }
                total_moves += moves;
                round++;
            }
        }
    };

} // namespace HeiProMap

#endif //HEIPROMAP_WAVER_SWAP_REFINEMENT_H
