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

#ifndef HEIPROMAP_WAVER_GREEDY_PARTITIONER_H
#define HEIPROMAP_WAVER_GREEDY_PARTITIONER_H

#include <algorithm>
#include <limits>
#include <numeric>
#include <vector>

#include "../definitions.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../utility/random_engine.h"
#include "../utility/aligned_array.h"

namespace HeiProMap {

    // Greedy initial partitioner using the full four-term cost model:
    //
    //   Cost(u -> block) = alpha_end  * Δ|E_π|
    //                    + beta_end   * ΔΣS
    //                    + alpha_link * ΔΣh
    //                    + beta_link  * ΔΣh·S
    //
    // For each already-assigned neighbor v of u with block v_block:
    //   Δ|E_π|  += (block != v_block)                           [new cut edge]
    //   ΔΣS     += w(u,v) * (block != v_block)                  [cut volume]
    //   ΔΣh     += dist(block, v_block)                         [hop count]
    //   ΔΣh·S   += w(u,v) * dist(block, v_block)               [hop volume]
    //
    // We minimize total cost over all candidate blocks.
    template<typename DistanceOracleT>
    inline void waver_greedy_partition(const graph_t &g,
                                       const DistanceOracleT &d_oracle,
                                       const f64 imbalance,
                                       const u64 seed,
                                       PartitionManager &p_manager,
                                       AlignedArray<vertex_t> &vertices,
                                       const f64 alpha_end,
                                       const f64 beta_end,
                                       const f64 alpha_link,
                                       const f64 beta_link) {
        const vertex_t    n                = g.n;
        const partition_t k                = d_oracle.get_k();

        p_manager.reset_weights();
        for (vertex_t u = 0; u < n; ++u) {
            p_manager.partition[u] = NO_ID;
        }

        for (vertex_t u = 0; u < n; ++u) {
            vertices[u] = u;
        }

        if (seed != 0) {
            RandomEngine re(seed);
            std::shuffle(vertices.get_ptr(), vertices.get_ptr() + n, re.generator);
        }

        for (vertex_t i = 0; i < n; ++i) {
            const vertex_t u        = vertices[i];
            const weight_t u_weight = g.v_weight(u);

            partition_t best_block = NO_ID;
            f64         min_cost   = std::numeric_limits<f64>::max();

            // collect distinct assigned neighbor blocks once per u (reused across block candidates)
            partition_t neighbor_blocks[256];
            size_t      n_neighbor_blocks = 0;
            for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                const partition_t v_block = p_manager.partition[g.edges_v[j]];
                if (v_block == NO_ID) { continue; }
                bool found = false;
                for (size_t x = 0; x < n_neighbor_blocks; ++x) {
                    if (neighbor_blocks[x] == v_block) { found = true; break; }
                }
                if (!found && n_neighbor_blocks < 256) {
                    neighbor_blocks[n_neighbor_blocks++] = v_block;
                }
            }

            for (partition_t block = 0; block < k; ++block) {
                if (p_manager.bweights[block] + u_weight > p_manager.get_lmax(block)) { continue; }

                // alpha_end: once per distinct neighbor block (not per edge)
                f64 cost = 0.0;
                for (size_t x = 0; x < n_neighbor_blocks; ++x) {
                    if (neighbor_blocks[x] != block) {
                        cost += alpha_end;
                    }
                }

                // beta_end, alpha_link, beta_link: per edge
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const vertex_t    v       = g.edges_v[j];
                    const partition_t v_block = p_manager.partition[v];
                    if (v_block == NO_ID) { continue; }

                    const f64 w    = static_cast<f64>(g.edge_w(j));
                    const f64 dist = static_cast<f64>(d_oracle.get(block, v_block));

                    cost += beta_end   * w * static_cast<f64>(block != v_block)
                          + alpha_link * dist
                          + beta_link  * w * dist;
                }

                // Penalize blocks above lmin to steer vertices toward under-filled blocks
                const weight_t block_lmin = p_manager.get_lmin(block);
                if (block_lmin > 0 && p_manager.bweights[block] >= block_lmin) {
                    f64 over_ratio = static_cast<f64>(p_manager.bweights[block]) / static_cast<f64>(block_lmin);
                    cost += alpha_end * over_ratio;
                }

                if (cost < min_cost ||
                    (cost == min_cost && p_manager.bweights[block] < p_manager.bweights[best_block])) {
                    min_cost   = cost;
                    best_block = block;
                }
            }

            if (best_block == NO_ID) {
                best_block = 0;
                weight_t min_w = p_manager.bweights[0];
                for (partition_t block = 1; block < k; ++block) {
                    if (p_manager.bweights[block] < min_w) {
                        min_w      = p_manager.bweights[block];
                        best_block = block;
                    }
                }
            }

            p_manager.set(u, u_weight, best_block);
        }
    }

}

#endif //HEIPROMAP_WAVER_GREEDY_PARTITIONER_H
