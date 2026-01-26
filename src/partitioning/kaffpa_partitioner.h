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

#ifndef HEIPROMAP_KAFFPA_PARTITIONER_H
#define HEIPROMAP_KAFFPA_PARTITIONER_H

#include "../definitions.h"
#include "../utility/assert_state.h"
#include "../../extern/local/kahip/include/kaHIP_interface.h"

namespace HeiProMap {
    enum KaffpaPartitionMode {
        KAFFPA_PARTITION_STRONG,
        KAFFPA_PARTITION_ECO,
        KAFFPA_PARTITION_FAST
    };

    inline void kaffpa_partition(graph_t &g,
                                 partition_t k,
                                 f64 imb,
                                 KaffpaPartitionMode kaffpa_mode,
                                 u64 seed,
                                 AlignedArray<partition_t> &partition,
                                 u64 kappa) {
        ASSERT(assert_graph(g));

        int n = (int) g.n;
        int m = (int) g.m;
        int *vwgt = (int *) malloc(n * sizeof(int));
        int *xadj = (int *) malloc((n + 1) * sizeof(int));
        int *adjcwgt = (int *) malloc(m * sizeof(int));
        int *adjncy = (int *) malloc(m * sizeof(int));
        int nparts = (int) k;
        double imbalance = std::max(imb, 0.0);
        bool suppress_output = true;
        int mode = FAST;
        int edge_cut;
        int best_edge_cut = std::numeric_limits<int>::max();
        int *part = (int *) malloc(n * sizeof(int));

        if (kaffpa_mode == KAFFPA_PARTITION_STRONG) { mode = STRONG; }
        if (kaffpa_mode == KAFFPA_PARTITION_ECO) { mode = ECO; }
        if (kaffpa_mode == KAFFPA_PARTITION_FAST) { mode = FAST; }

        // set vertex weights
        for (int i = 0; i < n; i++) { vwgt[i] = (int) g.v_weights[i]; }

        // set xadj
        for (int i = 0; i < n + 1; i++) { xadj[i] = (int) g.neighborhoods[i]; }

        // set adjncy
        for (int i = 0; i < m; i++) { adjncy[i] = (int) g.edges_v[i]; }

        // set adjcwgt
        for (int i = 0; i < m; i++) { adjcwgt[i] = (int) g.edges_w[i]; }

        for (u64 j = 0; j < kappa; j++) {
            kaffpa(&n, vwgt, xadj, adjcwgt, adjncy, &nparts, &imbalance, suppress_output, (int) seed + j, mode, &edge_cut, part);

            if (edge_cut < best_edge_cut) {
                best_edge_cut = edge_cut;
                for (int i = 0; i < n; i++) { partition[i] = part[i]; }
            }
        }

        free(vwgt);
        free(xadj);
        free(adjcwgt);
        free(adjncy);
        free(part);
    }
}

#endif //HEIPROMAP_KAFFPA_PARTITIONER_H
