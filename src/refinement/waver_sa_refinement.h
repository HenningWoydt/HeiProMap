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

#ifndef HEIPROMAP_WAVER_SA_REFINEMENT_H
#define HEIPROMAP_WAVER_SA_REFINEMENT_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "../definitions.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"

namespace HeiProMap {

    class WaverSARefinement {
        using graph_t     = CSRGraph;
        using p_manager_t = PartitionManager;

    public:
        template<typename DistanceOracleT>
        static void refine(graph_t              &g,
                            p_manager_t          &p_manager,
                            const DistanceOracleT &d_oracle,
                            size_t                max_iterations,
                            size_t                max_radius,
                            u64                   seed,
                            f64                   alpha_end,
                            f64                   beta_end,
                            f64                   alpha_link,
                            f64                   beta_link) {
            const size_t W = static_cast<size_t>(d_oracle.current_W());
            const partition_t k = p_manager.get_k();
            const size_t H = (k + W - 1) / W;

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

            auto pick_neighbor = [&](partition_t p, std::mt19937_64 &rng_local,
                                     size_t radius) -> partition_t {
                const size_t px = p % W, py = p / W;
                std::uniform_int_distribution<int> offset_dist(
                    -static_cast<int>(radius), static_cast<int>(radius));
                for (size_t attempt = 0; attempt < 10; ++attempt) {
                    const int dx = offset_dist(rng_local);
                    const int dy = offset_dist(rng_local);
                    if (dx == 0 && dy == 0) continue;
                    const int qx = static_cast<int>(px) + dx;
                    const int qy = static_cast<int>(py) + dy;
                    if (qx < 0 || qx >= static_cast<int>(W) ||
                        qy < 0 || qy >= static_cast<int>(H)) continue;
                    const partition_t q = static_cast<partition_t>(qy * W + qx);
                    if (q < k && !blk_verts[q].empty()) return q;
                }
                return NO_ID;
            };

            // Build cost-weighted block sampling: blocks with higher cost are
            // picked more often so SA focuses on the worst-mapped blocks
            std::vector<f64> blk_costs(k, 0.0);
            for (partition_t b = 0; b < k; ++b) {
                blk_costs[b] = block_cost(blk_verts[b], b);
            }
            std::vector<f64> cumulative(k);
            cumulative[0] = blk_costs[0];
            for (partition_t b = 1; b < k; ++b) {
                cumulative[b] = cumulative[b - 1] + blk_costs[b];
            }

            std::mt19937_64 rng(seed);
            std::uniform_real_distribution<f64> unit_dist(0.0, 1.0);

            auto sample_block = [&]() -> partition_t {
                f64 r = unit_dist(rng) * cumulative[k - 1];
                auto it = std::lower_bound(cumulative.begin(), cumulative.end(), r);
                return static_cast<partition_t>(it - cumulative.begin());
            };

            // Calibrate: sample worsening deltas
            std::vector<f64> sample_deltas;
            for (size_t i = 0; i < 2000; ++i) {
                const partition_t p = sample_block();
                if (blk_verts[p].empty()) continue;
                const partition_t q = pick_neighbor(p, rng, max_radius);
                if (q == NO_ID) continue;

                const f64 cost_before = block_cost(blk_verts[p], p) + block_cost(blk_verts[q], q);
                for (const vertex_t u : blk_verts[p]) p_manager.partition[u] = q;
                for (const vertex_t v : blk_verts[q]) p_manager.partition[v] = p;
                const f64 cost_after = block_cost(blk_verts[p], q) + block_cost(blk_verts[q], p);
                for (const vertex_t u : blk_verts[p]) p_manager.partition[u] = p;
                for (const vertex_t v : blk_verts[q]) p_manager.partition[v] = q;

                const f64 delta = cost_after - cost_before;
                if (delta > 0) sample_deltas.push_back(delta);
            }

            f64 T;
            if (sample_deltas.size() < 10) {
                T = 1.0;
            } else {
                std::sort(sample_deltas.begin(), sample_deltas.end());
                // Start from the very smallest worsening delta (P1),
                // accept it with ~20% probability
                f64 p1_delta = sample_deltas[std::max(static_cast<size_t>(1),
                                             sample_deltas.size() / 100)];
                // exp(-p1/T) = 0.2 => T = p1 / ln(5) = p1 / 1.609
                T = p1_delta / 1.609;
            }

            // Cool based on actual evaluated swaps, not loop iterations
            const f64 T_min = T * 1e-6;
            const f64 cooling = std::pow(T_min / T, 1.0 / static_cast<f64>(max_iterations));

            // Save best solution
            std::vector<partition_t> best_partition(g.n);
            std::vector<weight_t> best_bweights(k);
            f64 best_total_cost = 0.0;
            for (partition_t b = 0; b < k; ++b) {
                best_total_cost += block_cost(blk_verts[b], b);
            }
            f64 current_total_cost = best_total_cost;
            for (vertex_t u = 0; u < g.n; ++u) best_partition[u] = p_manager.partition[u];
            for (partition_t b = 0; b < k; ++b) best_bweights[b] = p_manager.bweights[b];

            size_t accepted = 0, improved = 0, evaluated = 0;

            for (size_t iter = 0; iter < max_iterations; ++iter) {
                const partition_t p = sample_block();
                if (blk_verts[p].empty()) continue;
                const partition_t q = pick_neighbor(p, rng, max_radius);
                if (q == NO_ID) continue;

                if (p_manager.bweights[q] > p_manager.get_lmax(p) ||
                    p_manager.bweights[p] > p_manager.get_lmax(q)) continue;

                const f64 cost_p = block_cost(blk_verts[p], p);
                const f64 cost_q = block_cost(blk_verts[q], q);
                const f64 cost_before = cost_p + cost_q;

                for (const vertex_t u : blk_verts[p]) p_manager.partition[u] = q;
                for (const vertex_t v : blk_verts[q]) p_manager.partition[v] = p;
                const f64 cost_after = block_cost(blk_verts[p], q) + block_cost(blk_verts[q], p);

                const f64 delta = cost_after - cost_before;
                evaluated++;

                bool accept;
                if (delta <= 0) {
                    accept = true;
                    improved++;
                } else {
                    accept = unit_dist(rng) < std::exp(-delta / T);
                }

                if (accept) {
                    const weight_t wp = p_manager.bweights[p];
                    const weight_t wq = p_manager.bweights[q];
                    p_manager.bweights[p] = wq;
                    p_manager.bweights[q] = wp;
                    std::swap(blk_verts[p], blk_verts[q]);
                    current_total_cost += delta;
                    accepted++;

                    if (current_total_cost < best_total_cost) {
                        best_total_cost = current_total_cost;
                        for (vertex_t u = 0; u < g.n; ++u) best_partition[u] = p_manager.partition[u];
                        for (partition_t b = 0; b < k; ++b) best_bweights[b] = p_manager.bweights[b];
                    }
                } else {
                    for (const vertex_t u : blk_verts[p]) p_manager.partition[u] = p;
                    for (const vertex_t v : blk_verts[q]) p_manager.partition[v] = q;
                }

                T *= cooling;
            }

            std::cout << "    [SA] k=" << k << " eval=" << evaluated
                      << " accepted=" << accepted << " improved=" << improved
                      << " T_init=" << (T / std::pow(cooling, evaluated))
                      << " T_final=" << T << std::endl;

            for (vertex_t u = 0; u < g.n; ++u) p_manager.partition[u] = best_partition[u];
            for (partition_t b = 0; b < k; ++b) p_manager.bweights[b] = best_bweights[b];
        }
    };

} // namespace HeiProMap

#endif //HEIPROMAP_WAVER_SA_REFINEMENT_H
