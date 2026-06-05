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
#include "../datastructures/HeiProMap_solver.h"

namespace HeiProMap {
    inline void heipromap_partition(const DynGraph &g,
                                    const std::vector<partition_t> &hierarchy,
                                    const std::vector<weight_t> &distance,
                                    f64 imbalance,
                                    u64 seed,
                                    u64 n_threads,
                                    const std::string &mode_str,
                                    std::vector<partition_t> &partition) {
        CSRGraph csr_g(g.n, g.m, g.g_weight);
        for (vertex_t u = 0; u < g.n; ++u) {
            csr_g.v_weights[u] = g.v_weights[u];
            csr_g.neighborhoods[u + 1] = csr_g.neighborhoods[u] + g.neighbors[u].size();
            for (size_t i = 0; i < g.neighbors[u].size(); ++i) {
                csr_g.edges_v[csr_g.neighborhoods[u] + i] = g.neighbors[u][i].u;
                csr_g.edges_w[csr_g.neighborhoods[u] + i] = g.neighbors[u][i].w;
            }
        }

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

        HeiProMapSolver solver(std::move(csr_g), ac);
        std::vector<vertex_t> result = solver.solve();

        partition.resize(g.n);
        for (vertex_t u = 0; u < g.n; ++u) {
            partition[u] = result[u];
        }
    }
}

#endif //HEIPROMAP_DYN_PARTITIONER_H
