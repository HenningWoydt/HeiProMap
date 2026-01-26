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

#ifndef HEIPROMAP_MTKAHYPAR_PARTITION_H
#define HEIPROMAP_MTKAHYPAR_PARTITION_H

#include "../definitions.h"
#include "../utility/assert_state.h"
#include "../../extern/local/mt-kahypar/include/mtkahypar.h"

namespace HeiProMap {
    enum MtKahyparPartitionMode {
        MTKAHYPAR_DEFAULT         = 0,
        MTKAHYPAR_QUALITY         = 1,
        MTKAHYPAR_HIGHEST_QUALITY = 2,
    };

    // NOTE: your signature does not have n_threads, so we use 1 thread by default.
    // If you want threads, add a parameter or read from a global config.
    inline void mtkahypar_partition(graph_t &g,
                                    partition_t k,
                                    f64 imb,
                                    MtKahyparPartitionMode mode,
                                    u64 seed,
                                    AlignedArray<partition_t> &partition) {
        ASSERT(assert_graph(g));

        imb = std::max(imb, 0.0);

        mt_kahypar_error_t error{};
        error.status = SUCCESS;
        error.msg = nullptr;
        error.msg_len = 0;

        // Init library (thread pool, etc.)
        const u64 n_threads = 1;
        mt_kahypar_initialize((size_t) n_threads,
                              false /* interleaved NUMA allocation policy */);

        // Preset -> context
        mt_kahypar_preset_type_t preset;
        switch (mode) {
            case MTKAHYPAR_DEFAULT: preset = DEFAULT;
                break;
            case MTKAHYPAR_QUALITY: preset = QUALITY;
                break;
            case MTKAHYPAR_HIGHEST_QUALITY: preset = HIGHEST_QUALITY;
                break;
            default:
                std::cerr << "Mt-KaHyPar mode " << (u64) mode << " not known!\n";
                std::abort();
        }

        mt_kahypar_context_t *context = mt_kahypar_context_from_preset(preset);
        if (!context) {
            std::cerr << "mt_kahypar_context_from_preset returned null\n";
            std::abort();
        }

        mt_kahypar_set_partitioning_parameters(context,
                                               (mt_kahypar_partition_id_t) k,
                                               (double) imb,
                                               CUT /* objective for graphs */);

        // NOTE: in your other code this is mt_kahypar_set_seed(seed) (no context).
        // The C API differs by version; this matches your working snippet.
        mt_kahypar_set_seed((uint64_t) seed);

        // Optional: silence logging
        mt_kahypar_set_context_parameter(context, VERBOSE, "0", &error);
        if (error.status != SUCCESS) {
            std::cerr << "mt_kahypar_set_context_parameter failed: " << (error.msg ? error.msg : "<null>") << "\n";
            std::abort();
        }

        const auto n = (mt_kahypar_hypernode_id_t) g.n;

        // Vertex weights
        std::vector<mt_kahypar_hypernode_weight_t> v_weights((size_t) n);
        for (u64 i = 0; i < (u64) n; ++i) {
            v_weights[i] = (mt_kahypar_hypernode_weight_t) g.v_weights[i];
        }

        // Build undirected edge list: each undirected edge becomes a 2-pin "hyperedge"
        // Only keep u < v to avoid duplicating symmetric CSR.
        std::vector<mt_kahypar_hypernode_id_t> edges;
        std::vector<mt_kahypar_hyperedge_weight_t> e_weights;
        edges.reserve((size_t) 2 * (size_t) g.m);
        e_weights.reserve((size_t) g.m);

        for (u64 u = 0; u < g.n; ++u) {
            for (u64 j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                const u64 v = (u64) g.edges_v[j];
                const u64 w = (u64) g.edges_w[j];
                if (u >= v) continue;

                edges.push_back((mt_kahypar_hypernode_id_t) u);
                edges.push_back((mt_kahypar_hypernode_id_t) v);
                e_weights.push_back((mt_kahypar_hyperedge_weight_t) w);
            }
        }

        const auto num_edges = (mt_kahypar_hyperedge_id_t) e_weights.size();
        if (num_edges == 0) {
            // Degenerate graph: all isolated / no u<v edges
            for (u64 i = 0; i < g.n; ++i) partition[i] = 0;
            mt_kahypar_free_context(context);
            return;
        }

        // Build KaHyPar graph (special helper for 2-pin hyperedges)
        mt_kahypar_hypergraph_t hg = mt_kahypar_create_graph(context,
                                                             n,
                                                             num_edges,
                                                             edges.data(),
                                                             e_weights.data(),
                                                             v_weights.data(),
                                                             &error);

        if (error.status != SUCCESS) {
            std::cerr << "mt_kahypar_create_graph failed: " << (error.msg ? error.msg : "<null>") << "\n";
            std::abort();
        }

        mt_kahypar_partitioned_hypergraph_t phg =
                mt_kahypar_partition(hg, context, &error);

        if (error.status != SUCCESS) {
            std::cerr << "mt_kahypar_partition failed: " << (error.msg ? error.msg : "<null>") << "\n";
            std::abort();
        }

        auto part = std::make_unique<mt_kahypar_partition_id_t[]>((size_t) mt_kahypar_num_hypernodes(hg));
        mt_kahypar_get_partition(phg, part.get());

        for (u64 i = 0; i < (u64) n; ++i) {
            partition[i] = (partition_t) part[i];
            ASSERT((u64)partition[i] < (u64)k);
        }

        mt_kahypar_free_partitioned_hypergraph(phg);
        mt_kahypar_free_hypergraph(hg);
        mt_kahypar_free_context(context);
    }
}

#endif //HEIPROMAP_MTKAHYPAR_PARTITION_H
