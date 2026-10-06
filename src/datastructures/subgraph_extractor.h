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

#include <algorithm>
#include <limits>
#include <vector>
#include <omp.h>

#include "csr_graph.h"
#include "../utility/aligned_array.h"
#include "../utility/translation_table.h"
#include "../utility/small_translation_table.h"

namespace HeiProMap {
    class SubgraphExtractor {
    public:
        // Multi-block batch extraction
        partition_t k = 0;
        u64 threads = 1;
        std::vector<CSRGraph> graphs;
        std::vector<SmallTranslationTable<vertex_t>> tts;

        // Reusable scratch buffers — flat CSR layout for block vertices
        std::vector<vertex_t> bv_data;
        std::vector<size_t> bv_offsets;
        AlignedArray<weight_t> block_weights;

        // Reusable buffers for parallel extraction
        struct VertexEntry { vertex_t u; partition_t dense_id; };
        std::vector<std::vector<VertexEntry>> par_thread_entries;
        std::vector<AlignedArray<weight_t>> par_thread_weights;

        size_t heap_bytes() const {
            size_t bytes = block_weights.heap_bytes()
                         + bv_data.capacity() * sizeof(vertex_t)
                         + bv_offsets.capacity() * sizeof(size_t);
            for (size_t i = 0; i < graphs.size(); ++i) {
                bytes += graphs[i].heap_bytes();
            }
            for (size_t i = 0; i < par_thread_weights.size(); ++i) {
                bytes += par_thread_weights[i].heap_bytes();
                bytes += par_thread_entries[i].capacity() * sizeof(VertexEntry);
            }
            return bytes;
        }

        MemoryPool *m_pool = nullptr;

        void set_pool(MemoryPool *pool) {
            m_pool = pool;
            block_weights.set_pool(pool);
            for (auto &tw : par_thread_weights) {
                tw.set_pool(pool);
            }
            for (auto &g : graphs) {
                g.set_pool(pool);
            }
        }

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
                            const u8 target_side,
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
            AlignedArray<vertex_t> old_to_new;
            old_to_new.initialize(g.n, std::numeric_limits<vertex_t>::max());

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
                const vertex_t u = tt.get_o(i);
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const vertex_t v = g.edges_v[j];
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
                const vertex_t u = tt.get_o(i);
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const vertex_t v = g.edges_v[j];
                    const vertex_t sub_v = old_to_new[v];
                    if (sub_v != std::numeric_limits<vertex_t>::max()) {
                        sub_g.edges_v[edge_cursor] = sub_v;
                        sub_g.edges_w[edge_cursor] = g.edges_w[j];
                        edge_cursor++;
                    }
                }
            }
        }

    private:
        template<typename PartitionManagerT, typename VerticesT>
        void build_subgraph_from_vertices(const CSRGraph &g,
                                          const PartitionManagerT &p_manager,
                                          const partition_t orig_block_id,
                                          const size_t dense_idx,
                                          const VerticesT &verts,
                                          const weight_t total_w) {
            const vertex_t sub_n = static_cast<vertex_t>(verts.size());
            auto &cur_tt = tts[dense_idx];

            for (vertex_t sub_u = 0; sub_u < sub_n; ++sub_u) {
                cur_tt.add(verts[sub_u], sub_u);
            }

            CSRGraph &sub_g = graphs[dense_idx];
            sub_g.resize(sub_n, 0, total_w);
            sub_g.uniform_v_weights = g.uniform_v_weights;
            sub_g.uniform_e_weights = g.uniform_e_weights;

            // Pass 1: count induced degrees directly into neighborhoods
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
                sub_g.neighborhoods[sub_u + 1] = deg;
                block_m += deg;
            }

            // Prefix sum in-place
            sub_g.m = static_cast<vertex_t>(block_m);
            sub_g.edges_v.initialize(block_m);
            sub_g.edges_w.initialize(block_m);

            size_t sum = 0;
            for (vertex_t sub_u = 0; sub_u < sub_n; ++sub_u) {
                const size_t deg = sub_g.neighborhoods[sub_u + 1];
                sub_g.neighborhoods[sub_u] = sum;
                sum += deg;
            }
            sub_g.neighborhoods[sub_n] = sum;

            // Pass 2: fill edges and vertex weights
            for (vertex_t sub_u = 0; sub_u < sub_n; ++sub_u) {
                const vertex_t orig_u = verts[sub_u];
                sub_g.v_weights[sub_u] = g.v_weights[orig_u];

                size_t edge_cursor = sub_g.neighborhoods[sub_u];
                for (size_t i = g.neighborhoods[orig_u]; i < g.neighborhoods[orig_u + 1]; ++i) {
                    const vertex_t orig_v = g.edges_v[i];
                    if (p_manager[orig_v] == orig_block_id) {
                        sub_g.edges_v[edge_cursor] = cur_tt.get_n(orig_v);
                        sub_g.edges_w[edge_cursor] = g.edges_w[i];
                        edge_cursor++;
                    }
                }
            }
        }

    public:
        void initialize(const partition_t t_k, const u64 t_threads, MemoryPool *pool = nullptr) {
            if (pool) set_pool(pool);
            k = t_k;
            threads = t_threads;
            block_weights.initialize(t_k, 0);

            par_thread_entries.resize(t_threads);
            par_thread_weights.resize(t_threads);
            for (u64 tid = 0; tid < t_threads; ++tid) {
                par_thread_weights[tid].initialize(t_k, 0);
            }
        }

        template<typename PartitionManagerT>
        void extract(const CSRGraph &g,
                     PartitionManagerT &p_manager,
                     const std::vector<partition_t> &ids,
                     const std::vector<partition_t> &id_to_dense) {
            if (ids.empty()) return;

            const size_t num_active = ids.size();
            if (graphs.size() < num_active) {
                size_t old_size = graphs.size();
                graphs.resize(num_active);
                tts.resize(num_active);
                if (m_pool) {
                    for (size_t i = old_size; i < num_active; ++i) {
                        graphs[i].set_pool(m_pool);
                    }
                }
            }

            if (threads > 1 && g.n >= 1024) {
                extract_parallel(g, p_manager, ids, id_to_dense);
            } else {
                extract_serial(g, p_manager, ids, id_to_dense);
            }
        }

        struct VertexSpan {
            const vertex_t *m_begin;
            size_t m_size;
            const vertex_t *begin() const { return m_begin; }
            const vertex_t *end() const { return m_begin + m_size; }
            size_t size() const { return m_size; }
            const vertex_t &operator[](size_t i) const { return m_begin[i]; }
        };

        template<typename PartitionManagerT>
        void build_bv_flat(const CSRGraph &g,
                           const PartitionManagerT &p_manager,
                           const std::vector<partition_t> &ids,
                           const std::vector<partition_t> &id_to_dense,
                           const size_t num_active) {
            bv_offsets.resize(num_active + 1);
            bv_offsets[0] = 0;
            for (size_t i = 0; i < num_active; ++i) {
                bv_offsets[i + 1] = bv_offsets[i] + p_manager.n_vertices[ids[i]];
                block_weights[i] = 0;
            }
            const size_t total_verts = bv_offsets[num_active];
            bv_data.resize(total_verts);

            std::vector<size_t> counters(num_active, 0);

            for (vertex_t u = 0; u < g.n; ++u) {
                const partition_t u_id = p_manager[u];
                const partition_t dense_u = id_to_dense[u_id];
                if (dense_u == std::numeric_limits<partition_t>::max()) continue;

                bv_data[bv_offsets[dense_u] + counters[dense_u]] = u;
                counters[dense_u]++;
                block_weights[dense_u] += g.v_weights[u];
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
            }

            build_bv_flat(g, p_manager, ids, id_to_dense, num_active);

            for (size_t dense_idx = 0; dense_idx < num_active; ++dense_idx) {
                VertexSpan verts{bv_data.data() + bv_offsets[dense_idx],
                                bv_offsets[dense_idx + 1] - bv_offsets[dense_idx]};
                build_subgraph_from_vertices(g, p_manager, ids[dense_idx], dense_idx,
                                            verts, block_weights[dense_idx]);
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
                const u64 tid = omp_get_thread_num();
                auto &entries = par_thread_entries[tid];
                entries.clear();
                for (size_t dense_idx = 0; dense_idx < num_active; ++dense_idx) {
                    par_thread_weights[tid][dense_idx] = 0;
                }

                const vertex_t chunk = (g.n + num_threads - 1) / num_threads;
                const vertex_t start_u = std::min(g.n, static_cast<vertex_t>(tid * chunk));
                const vertex_t end_u = std::min(g.n, static_cast<vertex_t>(start_u + chunk));

                for (vertex_t u = start_u; u < end_u; ++u) {
                    const partition_t u_id = p_manager[u];
                    const partition_t dense_u = id_to_dense[u_id];
                    if (dense_u == std::numeric_limits<partition_t>::max()) continue;

                    entries.push_back({u, dense_u});
                    par_thread_weights[tid][dense_u] += g.v_weights[u];
                }
            }

            // Build flat block_vertices from p_manager sizes
            bv_offsets.resize(num_active + 1);
            bv_offsets[0] = 0;
            for (size_t i = 0; i < num_active; ++i) {
                bv_offsets[i + 1] = bv_offsets[i] + p_manager.n_vertices[ids[i]];
            }
            const size_t total_verts = bv_offsets[num_active];
            bv_data.resize(total_verts);

            std::vector<size_t> counters(num_active, 0);
            for (u64 tid = 0; tid < num_threads; ++tid) {
                for (const auto &e : par_thread_entries[tid]) {
                    bv_data[bv_offsets[e.dense_id] + counters[e.dense_id]] = e.u;
                    counters[e.dense_id]++;
                }
            }

            // Parallel assembly per active block
            #pragma omp parallel for schedule(dynamic) num_threads(num_threads)
            for (size_t dense_idx = 0; dense_idx < num_active; ++dense_idx) {
                weight_t total_w = 0;
                for (u64 tid = 0; tid < num_threads; ++tid) {
                    total_w += par_thread_weights[tid][dense_idx];
                }

                VertexSpan verts{bv_data.data() + bv_offsets[dense_idx],
                                bv_offsets[dense_idx + 1] - bv_offsets[dense_idx]};
                build_subgraph_from_vertices(g, p_manager, ids[dense_idx], dense_idx,
                                            verts, total_w);
            }
        }
    };
}

#endif // HEIPROMAP_SUBGRAPH_EXTRACTOR_H
