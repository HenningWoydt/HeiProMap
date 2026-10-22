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

#ifndef HEIPROMAP_PARTITION_ALIGNER_H
#define HEIPROMAP_PARTITION_ALIGNER_H

#include <vector>
#include <map>
#include <algorithm>

#include "../definitions.h"
#include "../datastructures/dyn_graph.h"
#include "../utility/hungarian.h"
#include "../utility/aligned_array.h"

namespace HeiProMap {
    inline std::vector<int> match_blocks_general(const DynGraph &g,
                                                 const AlignedArray<partition_t> &previous_partition,
                                                 const std::vector<partition_t> &new_partition,
                                                 const std::vector<u64> &hierarchy,
                                                 int level, partition_t p_old, partition_t p_new,
                                                 const std::vector<partition_t> &D,
                                                 const std::vector<char> &is_new) {
        int count = (int) hierarchy[level];
        std::vector<std::vector<weight_t> > cost(count, std::vector<weight_t>(count, 0));
        int num_levels = (int) hierarchy.size();

        for (vertex_t v = 0; v < g.n; ++v) {
            if (is_new[v]) continue;
            partition_t id_old = previous_partition[v];
            partition_t id_new = new_partition[v];

            if (level < num_levels - 1) {
                if (id_old / D[level + 1] != p_old || id_new / D[level + 1] != p_new) continue;
            }

            cost[(id_old / D[level]) % count][(id_new / D[level]) % count] += g.v_weights[v];
        }

        weight_t max_w = 0;
        for (int i = 0; i < count; ++i)
            for (int j = 0; j < count; ++j) max_w = std::max(max_w, cost[i][j]);

        std::vector<std::vector<weight_t> > neg_cost(count, std::vector<weight_t>(count));
        for (int i = 0; i < count; ++i)
            for (int j = 0; j < count; ++j) neg_cost[i][j] = max_w - cost[i][j];

        return solve_hungarian(neg_cost);
    }

    inline void run_alignment_recursive(const DynGraph &g,
                                        const AlignedArray<partition_t> &previous_partition,
                                        std::vector<partition_t> &new_partition,
                                        const std::vector<u64> &hierarchy,
                                        int level, partition_t p_old, partition_t p_new,
                                        const std::vector<partition_t> &D,
                                        std::vector<std::map<partition_t, std::vector<int> > > &inv_m,
                                        const std::vector<char> &is_new) {
        if (level < 0) return;

        std::vector<int> m = match_blocks_general(g, previous_partition, new_partition, hierarchy, level, p_old, p_new, D, is_new);

        std::vector<int> inv(m.size());
        for (size_t i = 0; i < m.size(); ++i) inv[m[i]] = (int) i;
        inv_m[level][p_new] = inv;

        int count = (int) hierarchy[level];
        for (int i = 0; i < count; ++i) {
            run_alignment_recursive(g, previous_partition, new_partition, hierarchy, level - 1, p_old * count + i, p_new * count + m[i], D, inv_m, is_new);
        }
    }

    /**
     * Hierarchically aligns new_partition to match previous_partition minimizing migration cost.
     */
    inline void align_partitions(const DynGraph &g,
                                 const AlignedArray<partition_t> &previous_partition,
                                 std::vector<partition_t> &new_partition,
                                 const std::vector<u64> &hierarchy) {
        if (g.n == 0) return;

        bool all_zeros = true;
        for (vertex_t v = 0; v < g.n; ++v) {
            if (previous_partition[v] != 0) {
                all_zeros = false;
                break;
            }
        }
        if (all_zeros) return;

        int num_levels = (int) hierarchy.size();
        std::vector<partition_t> D(num_levels + 1, 1);
        for (int i = 0; i < num_levels; ++i) D[i + 1] = D[i] * hierarchy[i];

        // New vertices have no previous block, so they must not influence the matching.
        std::vector<char> is_new = g.new_vertex_mask();

        std::vector<std::map<partition_t, std::vector<int> > > inv_m(num_levels);
        run_alignment_recursive(g, previous_partition, new_partition, hierarchy, num_levels - 1, 0, 0, D, inv_m, is_new);

        for (vertex_t v = 0; v < g.n; ++v) {
            partition_t old_new_id = new_partition[v];
            partition_t aligned_id = 0;
            for (int l = num_levels - 1; l >= 0; --l) {
                partition_t prefix_new = old_new_id / D[l + 1];
                int comp_new = (int) ((old_new_id / D[l]) % hierarchy[l]);
                int comp_old = inv_m[l][prefix_new][comp_new];
                aligned_id = aligned_id * hierarchy[l] + comp_old;
            }
            new_partition[v] = aligned_id;
        }
    }
}

#endif // HEIPROMAP_PARTITION_ALIGNER_H
