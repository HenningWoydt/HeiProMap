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

#ifndef HEIPROMAP_EMPTY_BLOCK_DIFFUSION_H
#define HEIPROMAP_EMPTY_BLOCK_DIFFUSION_H

#include <algorithm>
#include <vector>

#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../definitions.h"
#include "../utility/profiler.h"

namespace HeiProMap {
    // Balancing step for grid mappings: fills empty PEs by moving vertices from their grid neighbors, then
    // equalizes below-average PEs from above-average neighbors. Respects lmax and keeps donors above lmin.
    class EmptyBlockDiffusion {
    public:
        static void diffuse(const CSRGraph &g, PartitionManager &p_manager, const size_t W) {
            HEIPROMAP_PROFILE_SCOPE("rebalance", "EmptyBlockDiffusion", "diffuse");
            const partition_t k = p_manager.get_k();

            auto get_neighbors = [&](partition_t b, partition_t *nbrs) -> size_t {
                size_t cnt = 0;
                const size_t c = b % W;
                if (b >= static_cast<partition_t>(W)) nbrs[cnt++] = b - static_cast<partition_t>(W);
                if (b + static_cast<partition_t>(W) < k) nbrs[cnt++] = b + static_cast<partition_t>(W);
                if (c > 0) nbrs[cnt++] = b - 1;
                if (c + 1 < W && b + 1 < k) nbrs[cnt++] = b + 1;
                return cnt;
            };

            // Target weight: fair share per PE
            const weight_t avg_bw = g.g_weight / k;

            // Build per-block vertex lists
            std::vector<std::vector<vertex_t> > blk_v(k);
            for (vertex_t u = 0; u < g.n; ++u) {
                const partition_t b = p_manager.partition[u];
                if (b != NO_ID && b < k) { blk_v[b].push_back(u); }
            }

            // Phase 1: aggressively fill empties.
            // Pull from ANY non-empty neighbor, allow donor to go as low
            // as 1 (just don't empty it).  Repeat until no empties touch
            // a non-empty neighbor.
            bool progress = true;
            while (progress) {
                progress = false;
                for (partition_t e = 0; e < k; ++e) {
                    if (p_manager.bweights[e] > 0) { continue; }

                    partition_t nbrs[4];
                    const size_t n_nbrs = get_neighbors(e, nbrs);

                    partition_t donor = NO_ID;
                    weight_t donor_bw = 0;
                    for (size_t i = 0; i < n_nbrs; ++i) {
                        if (p_manager.bweights[nbrs[i]] > donor_bw) {
                            donor_bw = p_manager.bweights[nbrs[i]];
                            donor = nbrs[i];
                        }
                    }
                    const weight_t donor_lmin = p_manager.get_lmin(donor == NO_ID ? 0 : donor);
                    if (donor == NO_ID || donor_bw <= donor_lmin) {
                        continue;
                    }

                    const weight_t e_cap = p_manager.get_lmax(e);
                    const weight_t donor_floor = std::max(donor_lmin, donor_bw / 2);
                    for (const vertex_t u: blk_v[donor]) {
                        if (p_manager.partition[u] != donor) { continue; }
                        const weight_t uw = g.v_weight(u);
                        if (p_manager.bweights[e] + uw > e_cap) { continue; }
                        if (p_manager.bweights[donor] - uw < donor_floor) { continue; }

                        p_manager.bweights[donor] -= uw;
                        p_manager.partition[u] = e;
                        p_manager.bweights[e] += uw;
                        blk_v[e].push_back(u);
                        progress = true;
                    }
                }
            }

            // Phase 2: equalize — pull from above-average neighbors into
            // below-average blocks.
            progress = true;
            while (progress) {
                progress = false;
                for (partition_t e = 0; e < k; ++e) {
                    if (p_manager.bweights[e] >= avg_bw) { continue; }

                    const weight_t e_cap = p_manager.get_lmax(e);

                    partition_t nbrs[4];
                    const size_t n_nbrs = get_neighbors(e, nbrs);

                    partition_t donor = NO_ID;
                    weight_t donor_bw = avg_bw;
                    for (size_t i = 0; i < n_nbrs; ++i) {
                        if (p_manager.bweights[nbrs[i]] > donor_bw) {
                            donor_bw = p_manager.bweights[nbrs[i]];
                            donor = nbrs[i];
                        }
                    }
                    if (donor == NO_ID) { continue; }

                    for (const vertex_t u: blk_v[donor]) {
                        if (p_manager.partition[u] != donor) { continue; }
                        const weight_t uw = g.v_weight(u);
                        if (p_manager.bweights[e] + uw > e_cap) { continue; }
                        if (p_manager.bweights[donor] - uw < avg_bw) { continue; }

                        p_manager.bweights[donor] -= uw;
                        p_manager.partition[u] = e;
                        p_manager.bweights[e] += uw;
                        blk_v[e].push_back(u);
                        progress = true;

                        if (p_manager.bweights[e] >= avg_bw) { break; }
                    }
                }
            }
        }
    };
}

#endif //HEIPROMAP_EMPTY_BLOCK_DIFFUSION_H
