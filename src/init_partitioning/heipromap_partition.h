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

#ifndef HEIPROMAP_DYN_PARTITIONER_H
#define HEIPROMAP_DYN_PARTITIONER_H

#include <omp.h>
#include "../datastructures/dyn_graph.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../utility/random_engine.h"
#include "../solver/HeiProMap_solver.h"

namespace HeiProMap {
    inline void heipromap_partition(const DynGraph &g,
                                    const std::vector<partition_t> &hierarchy,
                                    const std::vector<weight_t> &distance,
                                    f64 imbalance,
                                    u64 seed,
                                    u64 n_threads,
                                    const std::string &mode_str,
                                    std::vector<partition_t> &partition) {
        const vertex_t invalid_v = std::numeric_limits<vertex_t>::max();
        std::vector<vertex_t> sparse_to_dense(g.n, invalid_v);
        std::vector<vertex_t> dense_to_sparse;
        dense_to_sparse.reserve(g.n);

        for (vertex_t u = 0; u < g.n; ++u) {
            if (g.vertex_exists(u)) {
                sparse_to_dense[u] = (vertex_t) dense_to_sparse.size();
                dense_to_sparse.push_back(u);
            }
        }

        const vertex_t n_dense = (vertex_t) dense_to_sparse.size();
        if (n_dense == 0) {
            partition.assign(g.n, 0);
            return;
        }

        // 2. Build dense CSRGraph
        CSRGraph csr_g(n_dense, g.m, g.g_weight);

        for (vertex_t d_u = 0; d_u < n_dense; ++d_u) {
            vertex_t s_u = dense_to_sparse[d_u];
            csr_g.v_weights[d_u] = g.v_weights[s_u];

            size_t edge_count = 0;
            for (const auto &nbr: g.neighbors[s_u]) {
                if (nbr.u < g.n && sparse_to_dense[nbr.u] != invalid_v) {
                    edge_count++;
                }
            }
            csr_g.neighborhoods[d_u + 1] = csr_g.neighborhoods[d_u] + edge_count;

            size_t edge_idx = 0;
            for (const auto &nbr: g.neighbors[s_u]) {
                if (nbr.u < g.n && sparse_to_dense[nbr.u] != invalid_v) {
                    csr_g.edges_v[csr_g.neighborhoods[d_u] + edge_idx] = sparse_to_dense[nbr.u];
                    csr_g.edges_w[csr_g.neighborhoods[d_u] + edge_idx] = nbr.w;
                    edge_idx++;
                }
            }
        }
        csr_g.m = (vertex_t) csr_g.neighborhoods[n_dense];

        AlgorithmConfiguration ac;
        ac.hierarchy = hierarchy;
        ac.distance = distance;
        ac.imbalance = imbalance;
        ac.seed = seed;
        ac.threads = n_threads;
        ac.k = 1;
        for (auto h: hierarchy) ac.k *= h;

        if (mode_str == "fast") ac.set_fast();
        else if (mode_str == "eco") ac.set_eco();
        else if (mode_str == "strong") ac.set_strong();
        else if (mode_str == "super-strong") ac.set_super_strong();
        else ac.set_strong();

        ac.mapping_out = ""; // Avoid writing to file unless explicitly requested

        std::vector<vertex_t> result;
        if (ac.use_binary_oracle()) {
            HeiProMapSolver<BinaryDistanceOracle> solver(std::move(csr_g), ac);
            result = solver.solve();
        } else {
            HeiProMapSolver<DistanceOracle> solver(std::move(csr_g), ac);
            result = solver.solve();
        }

        // 3. Translate dense result back to sparse partition
        partition.assign(g.n, 0);
        for (vertex_t d_u = 0; d_u < n_dense; ++d_u) {
            partition[dense_to_sparse[d_u]] = result[d_u];
        }
    }
}

#endif //HEIPROMAP_DYN_PARTITIONER_H
