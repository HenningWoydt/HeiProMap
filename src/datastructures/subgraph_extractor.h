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
#include <omp.h>

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
        // Reusable scratch buffers
        std::vector<std::vector<vertex_t> > block_vertices;
        std::vector<weight_t> block_weights;
        std::vector<size_t> scratch_degs;

        // Reusable buffers for parallel extraction
        std::vector<std::vector<std::vector<vertex_t> > > par_thread_vertices;
        std::vector<std::vector<weight_t> > par_thread_weights;

        void initialize(partition_t t_k, u64 t_threads) {
            k = t_k;
            threads = t_threads;
            graphs.resize(t_k);
            tts.resize(t_k);
            block_vertices.resize(t_k);
            block_weights.resize(t_k);
            par_thread_vertices.resize(t_threads);
            for (u64 tid = 0; tid < t_threads; ++tid) {
                par_thread_vertices[tid].resize(t_k);
            }
            par_thread_weights.assign(t_threads, std::vector<weight_t>(t_k, 0));
        }

        template<typename PartitionManagerT>
        void extract(const CSRGraph &g,
                     PartitionManagerT &p_manager,
                     const std::vector<partition_t> &ids,
                     const std::vector<partition_t> &id_to_dense) {
            if (ids.empty()) return;

            if (threads > 1 && g.n >= 1024) {
                extract_parallel(g, p_manager, ids, id_to_dense);
            } else {
                extract_serial(g, p_manager, ids, id_to_dense);
            }
        }

        template<typename PartitionManagerT>
        void extract_serial(const CSRGraph &g,
                            PartitionManagerT &p_manager,
                            const std::vector<partition_t> &ids,
                            const std::vector<partition_t> &id_to_dense) {
            const size_t num_active = ids.size();
            if (num_active == 0) return;

            for (size_t dense_idx = 0; dense_idx < num_active; ++dense_idx) {
                tts[dense_idx].clear();
                block_vertices[dense_idx].clear();
                block_weights[dense_idx] = 0;
            }

            // Single pass over g.n
            for (vertex_t u = 0; u < g.n; ++u) {
                partition_t u_id = p_manager[u];
                partition_t dense_u = id_to_dense[u_id];
                if (dense_u == std::numeric_limits<partition_t>::max()) continue;

                block_vertices[dense_u].push_back(u);
                block_weights[dense_u] += g.v_weights[u];
            }

            // Process each active block independently
            for (size_t dense_idx = 0; dense_idx < num_active; ++dense_idx) {
                const partition_t orig_block_id = ids[dense_idx];
                const auto &verts = block_vertices[dense_idx];
                const vertex_t sub_n = static_cast<vertex_t>(verts.size());

                for (vertex_t sub_u = 0; sub_u < sub_n; ++sub_u) {
                    tts[dense_idx].add(verts[sub_u], sub_u);
                }

                if (scratch_degs.size() < sub_n) {
                    scratch_degs.resize(sub_n);
                }
                size_t block_m = 0;
                for (vertex_t sub_u = 0; sub_u < sub_n; ++sub_u) {
                    const vertex_t orig_u = verts[sub_u];
                    size_t deg = 0;
                    for (size_t i = g.neighborhoods[orig_u]; i < g.neighborhoods[orig_u + 1]; ++i) {
                        const vertex_t orig_v = g.edges_v[i];
                        if (p_manager[orig_v] == orig_block_id) {
                            deg++;
                        }
                    }
                    scratch_degs[sub_u] = deg;
                    block_m += deg;
                }

                // Allocate CSRGraph
                graphs[dense_idx].resize(sub_n, static_cast<vertex_t>(block_m), block_weights[dense_idx]);

                // Prefix-sum neighborhoods
                graphs[dense_idx].neighborhoods[0] = 0;
                for (vertex_t sub_u = 0; sub_u < sub_n; ++sub_u) {
                    graphs[dense_idx].neighborhoods[sub_u + 1] = graphs[dense_idx].neighborhoods[sub_u] + scratch_degs[sub_u];
                }

                // Fill edges and vertex weights
                for (vertex_t sub_u = 0; sub_u < sub_n; ++sub_u) {
                    const vertex_t orig_u = verts[sub_u];
                    graphs[dense_idx].v_weights[sub_u] = g.v_weights[orig_u];

                    size_t edge_cursor = graphs[dense_idx].neighborhoods[sub_u];
                    for (size_t i = g.neighborhoods[orig_u]; i < g.neighborhoods[orig_u + 1]; ++i) {
                        const vertex_t orig_v = g.edges_v[i];
                        if (p_manager[orig_v] == orig_block_id) {
                            graphs[dense_idx].edges_v[edge_cursor] = tts[dense_idx].get_n(orig_v);
                            graphs[dense_idx].edges_w[edge_cursor] = g.edges_w[i];
                            edge_cursor++;
                        }
                    }
                }
            }
        }

        template<typename PartitionManagerT>
        void extract_parallel(const CSRGraph &g,
                              PartitionManagerT &p_manager,
                              const std::vector<partition_t> &ids,
                              const std::vector<partition_t> &id_to_dense) {
            const size_t num_active = ids.size();
            const u64 num_threads = threads;

            #pragma omp parallel for schedule(static) num_threads(num_threads)
            for (size_t dense_idx = 0; dense_idx < num_active; ++dense_idx) {
                tts[dense_idx].clear();
            }

            #pragma omp parallel num_threads(num_threads)
            {
                u64 tid = omp_get_thread_num();
                for (size_t dense_idx = 0; dense_idx < num_active; ++dense_idx) {
                    par_thread_vertices[tid][dense_idx].clear();
                    par_thread_weights[tid][dense_idx] = 0;
                }

                vertex_t chunk = (g.n + num_threads - 1) / num_threads;
                vertex_t start_u = std::min(g.n, tid * chunk);
                vertex_t end_u = std::min(g.n, start_u + chunk);

                for (vertex_t u = start_u; u < end_u; ++u) {
                    partition_t u_id = p_manager[u];
                    partition_t dense_u = id_to_dense[u_id];
                    if (dense_u == std::numeric_limits<partition_t>::max()) continue;

                    par_thread_vertices[tid][dense_u].push_back(u);
                    par_thread_weights[tid][dense_u] += g.v_weights[u];
                }
            }

            // Parallel assembly per active block
            #pragma omp parallel for schedule(dynamic) num_threads(num_threads)
            for (size_t dense_idx = 0; dense_idx < num_active; ++dense_idx) {
                partition_t orig_block_id = ids[dense_idx];

                // Merge thread-local vertices into single list
                std::vector<vertex_t> block_vertices;
                weight_t total_w = 0;
                size_t total_n = 0;
                for (u64 tid = 0; tid < num_threads; ++tid) {
                    total_w += par_thread_weights[tid][dense_idx];
                    total_n += par_thread_vertices[tid][dense_idx].size();
                }

                block_vertices.reserve(total_n);
                for (u64 tid = 0; tid < num_threads; ++tid) {
                    block_vertices.insert(block_vertices.end(), par_thread_vertices[tid][dense_idx].begin(), par_thread_vertices[tid][dense_idx].end());
                }

                // Populate translation table
                for (vertex_t sub_u = 0; sub_u < (vertex_t) block_vertices.size(); ++sub_u) {
                    vertex_t orig_u = block_vertices[sub_u];
                    tts[dense_idx].add(orig_u, sub_u);
                }

                // Count edges
                size_t block_m = 0;
                std::vector<size_t> degs(block_vertices.size(), 0);
                for (vertex_t sub_u = 0; sub_u < (vertex_t) block_vertices.size(); ++sub_u) {
                    vertex_t orig_u = block_vertices[sub_u];
                    for (size_t i = g.neighborhoods[orig_u]; i < g.neighborhoods[orig_u + 1]; ++i) {
                        vertex_t orig_v = g.edges_v[i];
                        if (p_manager[orig_v] == orig_block_id) {
                            degs[sub_u]++;
                            block_m++;
                        }
                    }
                }

                // Allocate CSRGraph
                graphs[dense_idx].resize(static_cast<vertex_t>(block_vertices.size()), static_cast<vertex_t>(block_m), total_w);

                // Build neighborhoods prefix sum
                graphs[dense_idx].neighborhoods[0] = 0;
                for (vertex_t sub_u = 0; sub_u < (vertex_t) block_vertices.size(); ++sub_u) {
                    graphs[dense_idx].neighborhoods[sub_u + 1] = graphs[dense_idx].neighborhoods[sub_u] + degs[sub_u];
                }

                // Fill edges and vertex weights
                for (vertex_t sub_u = 0; sub_u < (vertex_t) block_vertices.size(); ++sub_u) {
                    vertex_t orig_u = block_vertices[sub_u];
                    graphs[dense_idx].v_weights[sub_u] = g.v_weights[orig_u];

                    size_t edge_cursor = graphs[dense_idx].neighborhoods[sub_u];
                    for (size_t i = g.neighborhoods[orig_u]; i < g.neighborhoods[orig_u + 1]; ++i) {
                        vertex_t orig_v = g.edges_v[i];
                        if (p_manager[orig_v] == orig_block_id) {
                            vertex_t sub_v = tts[dense_idx].get_n(orig_v);
                            graphs[dense_idx].edges_v[edge_cursor] = sub_v;
                            graphs[dense_idx].edges_w[edge_cursor] = g.edges_w[i];
                            edge_cursor++;
                        }
                    }
                }
            }
        }
    };
}

#endif // HEIPROMAP_SUBGRAPH_EXTRACTOR_H
