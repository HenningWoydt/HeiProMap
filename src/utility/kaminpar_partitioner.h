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

#ifndef HEIPROMAP_KAMINPAR_PARTITIONER_H
#define HEIPROMAP_KAMINPAR_PARTITIONER_H

#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

#include <kaminpar.h>

#include "../definitions.h"
#include "../datastructures/csr_graph.h"
#include "mapping.h"
#include "profiler.h"
#include "utils.h"

namespace HeiProMap {

    // Thin wrapper around the KaMinPar shared-memory partitioner (extern/KaMinPar).
    // Partitions g into k blocks with the given imbalance and returns part[u] in [0, k).
    inline std::vector<partition_t> kaminpar_partition(const CSRGraph &g,
                                                       const partition_t k,
                                                       const f64 epsilon,
                                                       const u64 threads,
                                                       const u64 seed,
                                                       const std::string &preset = "default") {
        HEIPROMAP_PROFILE_SCOPE("partition", "KaMinPar", "partition");
        using NodeID = kaminpar::shm::NodeID;
        using EdgeID = kaminpar::shm::EdgeID;
        using NodeWeight = kaminpar::shm::NodeWeight;
        using EdgeWeight = kaminpar::shm::EdgeWeight;
        using BlockID = kaminpar::shm::BlockID;

        static_assert(sizeof(EdgeID) == sizeof(u64), "KaMinPar must be built with 64-bit edge IDs");

        if (g.n == 0) { return {}; }
        if (k <= 1) { return std::vector<partition_t>(g.n, 0); }

        const size_t m = g.neighborhoods[g.n];
        std::vector<EdgeID> xadj(g.n + 1);
        for (vertex_t u = 0; u <= g.n; ++u) { xadj[u] = static_cast<EdgeID>(g.neighborhoods[u]); }
        std::vector<NodeID> adjncy(m);
        std::vector<EdgeWeight> adjwgt(m);
        for (size_t e = 0; e < m; ++e) {
            adjncy[e] = static_cast<NodeID>(g.edges_v[e]);
            adjwgt[e] = g.uniform_e_weights ? EdgeWeight{1} : static_cast<EdgeWeight>(g.edges_w[e]);
        }
        std::vector<NodeWeight> vwgt;
        if (!g.uniform_v_weights) {
            vwgt.resize(g.n);
            for (vertex_t u = 0; u < g.n; ++u) { vwgt[u] = static_cast<NodeWeight>(g.v_weights[u]); }
        }

        kaminpar::KaMinPar shm(static_cast<int>(std::max<u64>(threads, 1)),
                               kaminpar::shm::create_context_by_preset_name(preset));
        shm.set_output_level(kaminpar::OutputLevel::QUIET);
        kaminpar::KaMinPar::reseed(static_cast<int>(seed));
        shm.copy_graph(xadj, adjncy, vwgt, adjwgt);

        std::vector<BlockID> part(g.n);
        shm.compute_partition(static_cast<BlockID>(k), static_cast<double>(epsilon), part);

        std::vector<partition_t> out(g.n);
        for (vertex_t u = 0; u < g.n; ++u) { out[u] = static_cast<partition_t>(part[u]); }
        return out;
    }

    // Result of partition_into_quotient: part[u] is the block of original vertex u, quotient has one
    // unit-weight vertex per block and the summed edge weights between blocks (so quotient.n == k).
    struct QuotientResult {
        CSRGraph quotient;
        std::vector<partition_t> part;
    };

    // Partitions g into k blocks with KaMinPar and contracts the blocks into the quotient graph.
    // kaminpar_threads == 0 uses all hardware threads.
    inline QuotientResult partition_into_quotient(const CSRGraph &g,
                                                  const partition_t k,
                                                  const f64 epsilon,
                                                  const u64 kaminpar_threads,
                                                  const u64 contraction_threads,
                                                  const u64 seed,
                                                  const std::string &preset) {
        const u64 threads = kaminpar_threads > 0 ? kaminpar_threads
                                                 : std::max<u64>(1, std::thread::hardware_concurrency());
        const auto sp = get_time_point();

        QuotientResult result;
        result.part = kaminpar_partition(g, k, epsilon, threads, seed, preset);

        Mapping mapping;
        mapping.initialize(g.n);
        mapping.set_coarse_n(k);
        for (vertex_t u = 0; u < g.n; ++u) { mapping.set(u, result.part[u]); }

        result.quotient.contract(g, mapping, contraction_threads);
        result.quotient.uniform_v_weights = true;
        result.quotient.v_weights.free_memory();
        result.quotient.g_weight = k;

        std::cout << "[kaminpar] n=" << g.n << " -> quotient n=" << result.quotient.n << " m=" << result.quotient.m
                  << " (" << get_seconds(sp, get_time_point()) << " s, " << threads << " threads)" << std::endl;
        return result;
    }
}

#endif //HEIPROMAP_KAMINPAR_PARTITIONER_H
