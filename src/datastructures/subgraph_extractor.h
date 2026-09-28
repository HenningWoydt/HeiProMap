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

#ifndef HEIPROMAP_SUBGRAPH_EXTRACTOR_H
#define HEIPROMAP_SUBGRAPH_EXTRACTOR_H

#include <vector>
#include <limits>

#include "csr_graph.h"
#include "../utility/translation_table.h"
#include "../utility/small_translation_table.h"

namespace HeiProMap {
    class SubgraphExtractor {
    public:
        /**
         * Extracts a subgraph from g based on the side array.
         * Only vertices u with side[u] == target_side are included.
         *
         * @param g The original graph.
         * @param side Side array of size g.n.
         * @param target_side The side to extract (0 or 1).
         * @param sub_g The output subgraph.
         * @param tt The output translation table (mapping sub_g vertices back to g vertices).
         */
        static void extract(const CSRGraph &g,
                            const std::vector<u8> &side,
                            u8 target_side,
                            CSRGraph &sub_g,
                            TranslationTable<vertex_t> &tt) {
            vertex_t sub_n = 0;
            for (vertex_t u = 0; u < g.n; ++u) {
                if (side[u] == target_side) sub_n++;
            }

            sub_g.n = sub_n;
            sub_g.v_weights.initialize(sub_n);
            sub_g.neighborhoods.initialize(sub_n + 1);
            sub_g.neighborhoods[0] = 0;
            sub_g.g_weight = 0;
            sub_g.uniform_v_weights = g.uniform_v_weights;
            sub_g.uniform_e_weights = g.uniform_e_weights;

            tt.reserve(sub_n, g.n);
            std::vector<vertex_t> old_to_new(g.n, std::numeric_limits<vertex_t>::max());

            vertex_t sub_u = 0;
            for (vertex_t u = 0; u < g.n; ++u) {
                if (side[u] == target_side) {
                    tt.add(u, sub_u);
                    old_to_new[u] = sub_u;
                    sub_g.v_weights[sub_u] = g.v_weights[u];
                    sub_g.g_weight += g.v_weights[u];
                    sub_u++;
                }
            }

            // First pass: count edges
            size_t sub_m = 0;
            for (vertex_t i = 0; i < sub_n; ++i) {
                vertex_t u = tt.get_o(i);
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    vertex_t v = g.edges_v[j];
                    if (old_to_new[v] != std::numeric_limits<vertex_t>::max()) {
                        sub_m++;
                    }
                }
                sub_g.neighborhoods[i + 1] = sub_m;
            }

            sub_g.m = sub_m;
            sub_g.edges_v.initialize(sub_m);
            sub_g.edges_w.initialize(sub_m);

            // Second pass: fill edges
            size_t edge_cursor = 0;
            for (vertex_t i = 0; i < sub_n; ++i) {
                vertex_t u = tt.get_o(i);
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    vertex_t v = g.edges_v[j];
                    vertex_t sub_v = old_to_new[v];
                    if (sub_v != std::numeric_limits<vertex_t>::max()) {
                        sub_g.edges_v[edge_cursor] = sub_v;
                        sub_g.edges_w[edge_cursor] = g.edges_w[i];
                        edge_cursor++;
                    }
                }
            }
        }

        // Multi-block batch extraction
        partition_t k = 0;
        u64 threads = 1;
        std::vector<CSRGraph> graphs;
        std::vector<SmallTranslationTable<vertex_t> > tts;
        std::vector<u64> idxs;

        void initialize(partition_t t_k, u64 t_threads) {
            k = t_k;
            threads = t_threads;
            graphs.resize(t_k);
            tts.resize(t_k);
            idxs.resize(t_k);
        }

        template<typename PartitionManagerT>
        void extract(const CSRGraph &g,
                     PartitionManagerT &p_manager,
                     const std::vector<u32> &one_hot) {
            // Reset
            for (size_t i = 0; i < k; ++i) {
                if (one_hot[i] == 1) {
                    graphs[i].clear();
                    tts[i].clear();
                }
            }
            std::fill(idxs.begin(), idxs.end(), 0);

            // Count
            for (vertex_t u = 0; u < g.n; ++u) {
                partition_t u_id = p_manager[u];
                if (one_hot[u_id] == 0) continue;

                tts[u_id].add(u, graphs[u_id].n);
                graphs[u_id].n += 1;
                graphs[u_id].g_weight += g.v_weights[u];

                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    vertex_t v = g.edges_v[i];
                    if (p_manager[v] == u_id) {
                        graphs[u_id].m += 1;
                    }
                }
            }

            // Allocate
            for (size_t i = 0; i < k; ++i) {
                if (one_hot[i] == 1) {
                    graphs[i].resize(graphs[i].n, graphs[i].m, graphs[i].g_weight);
                }
            }

            // Build
            for (vertex_t u = 0; u < g.n; ++u) {
                partition_t u_id = p_manager[u];
                if (one_hot[u_id] == 0) continue;

                vertex_t sub_u = tts[u_id].get_n(u);
                graphs[u_id].v_weights[sub_u] = g.v_weights[u];

                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    vertex_t v = g.edges_v[i];
                    if (p_manager[v] == u_id) {
                        vertex_t sub_v = tts[u_id].get_n(v);
                        graphs[u_id].edges_v[idxs[u_id]] = sub_v;
                        graphs[u_id].edges_w[idxs[u_id]] = g.edges_w[i];
                        idxs[u_id] += 1;
                    }
                }
                graphs[u_id].neighborhoods[sub_u + 1] = idxs[u_id];
            }
        }
    };
}

#endif // HEIPROMAP_SUBGRAPH_EXTRACTOR_H
