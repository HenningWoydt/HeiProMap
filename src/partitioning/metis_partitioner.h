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

#ifndef HEIPROMAP_METIS_PARTITIONER_H
#define HEIPROMAP_METIS_PARTITIONER_H

#include "../definitions.h"
#include "../utility/assert_state.h"
#include "../../extern/local/metis/include/metis.h"

namespace HeiProMap {
    enum MetisPartitionMode {
        METIS_PARTITION_RECURSIVE, // METIS_PartGraphRecursive (ptype=RB)
        METIS_PARTITION_KWAY       // METIS_PartGraphKway (ptype=KWAY)
    };

    inline void metis_partition(graph_t &g,
                                partition_t k,
                                f64 imb,
                                MetisPartitionMode metis_mode,
                                u64 seed,
                                AlignedArray<partition_t> &partition) {
        ASSERT(assert_graph(g));

        // METIS uses idx_t/real_t types (configured by IDXTYPEWIDTH/REALTYPEWIDTH in metis.h)
        idx_t nvtxs = (idx_t) g.n;
        idx_t ncon = (idx_t) 1;
        idx_t nparts = (idx_t) k;

        // Your graph CSR:
        // - neighborhoods: size n+1
        // - edges_v/edges_w: size m
        // METIS expects xadj size n+1, adjncy size xadj[n].
        idx_t *xadj = (idx_t *) std::malloc((size_t) (nvtxs + 1) * sizeof(idx_t));
        idx_t *adjncy = (idx_t *) std::malloc((size_t) g.m * sizeof(idx_t));

        // Optional weights (we provide them like KaHIP wrapper does).
        idx_t *vwgt = (idx_t *) std::malloc((size_t) nvtxs * sizeof(idx_t));
        idx_t *adjwgt = (idx_t *) std::malloc((size_t) g.m * sizeof(idx_t));

        idx_t *part = (idx_t *) std::malloc((size_t) nvtxs * sizeof(idx_t));

        ASSERT(xadj && adjncy && vwgt && adjwgt && part);

        for (idx_t i = 0; i < nvtxs; ++i) {
            vwgt[i] = (idx_t) g.v_weights[(u64) i];
        }
        for (idx_t i = 0; i < nvtxs + 1; ++i) {
            xadj[i] = (idx_t) g.neighborhoods[(u64) i];
        }
        for (idx_t i = 0; i < (idx_t) g.m; ++i) {
            adjncy[i] = (idx_t) g.edges_v[(u64) i];
            adjwgt[i] = (idx_t) g.edges_w[(u64) i];
        }

        // METIS options
        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options); // sets sane defaults :contentReference[oaicite:1]{index=1}

        // Choose scheme
        if (metis_mode == METIS_PARTITION_RECURSIVE) {
            options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;
        } else {
            options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
        }

        // Use edge-cut objective (default usually, but be explicit)
        options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;

        // Seed
        options[METIS_OPTION_SEED] = (idx_t) seed;

        // Imbalance tolerance:
        // METIS uses ubvec as a multiplicative factor (> 1.0).
        // We map your imb (epsilon) to ubvec = 1.0 + imb.
        real_t ubvec[1];
        ubvec[0] = (real_t) (1.0 + imb); // must be > 1.0 :contentReference[oaicite:2]{index=2}
        ubvec[0] = std::max(ubvec[0], (real_t) 1.0);

        // Let METIS assume equal target weights per part by passing nullptr for tpwgts.
        idx_t edgecut = 0;

        if (metis_mode == METIS_PARTITION_RECURSIVE) {
            METIS_PartGraphRecursive(
                &nvtxs, &ncon, xadj, adjncy,
                vwgt,
                /*vsize=*/nullptr,
                adjwgt,
                &nparts,
                /*tpwgts=*/nullptr,
                ubvec,
                options,
                &edgecut,
                part
            );
        } else {
            METIS_PartGraphKway(
                &nvtxs, &ncon, xadj, adjncy,
                vwgt,
                /*vsize=*/nullptr,
                adjwgt,
                &nparts,
                /*tpwgts=*/nullptr,
                ubvec,
                options,
                &edgecut,
                part
            );
        }

        for (idx_t i = 0; i < nvtxs; ++i) {
            partition[(u64) i] = (partition_t) part[i];
        }

        std::free(xadj);
        std::free(adjncy);
        std::free(vwgt);
        std::free(adjwgt);
        std::free(part);
    }
}

#endif //HEIPROMAP_METIS_PARTITIONER_H
