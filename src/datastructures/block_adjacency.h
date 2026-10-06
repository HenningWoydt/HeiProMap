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

#ifndef HEIPROMAP_BLOCK_ADJACENCY_H
#define HEIPROMAP_BLOCK_ADJACENCY_H

#include <vector>

#include <omp.h>

#include "../definitions.h"
#include "../utility/aligned_array.h"
#include "../utility/profiler.h"

namespace HeiProMap {
    class BlockAdjacency {
    public:
        struct Entry {
            partition_t target = 0;
            mutable u32 used_epoch = 0;
        };

        struct Span {
            const Entry *m_begin;
            const Entry *m_end;
            const Entry *begin() const { return m_begin; }
            const Entry *end() const { return m_end; }
            size_t size() const { return m_end - m_begin; }
        };

    private:
        partition_t m_k = 0;
        MemoryPool *m_pool = nullptr;
        AlignedArray<size_t> m_offsets;
        AlignedArray<size_t> m_degrees;
        AlignedArray<Entry> m_data;
        size_t m_data_capacity = 0;
        std::vector<AlignedArray<u32>> m_seen;
        std::vector<u32> m_epoch;

    public:
        size_t heap_bytes() const {
            size_t bytes = m_offsets.heap_bytes() + m_degrees.heap_bytes() + m_data.heap_bytes();
            for (size_t t = 0; t < m_seen.size(); ++t) {
                bytes += m_seen[t].heap_bytes();
            }
            return bytes;
        }

        void set_pool(MemoryPool *pool) {
            m_pool = pool;
            m_offsets.set_pool(pool);
            m_degrees.set_pool(pool);
            m_data.set_pool(pool);
            for (auto &s : m_seen) {
                s.set_pool(pool);
            }
        }

        void initialize(const partition_t k, MemoryPool *pool = nullptr) {
            if (pool) set_pool(pool);
            m_k = k;
            m_offsets.initialize(m_k + 1, 0);
            m_degrees.initialize(m_k, 0);
        }

        void initialize(const partition_t k, const u64 num_threads, MemoryPool *pool = nullptr) {
            if (pool) set_pool(pool);
            m_k = k;
            m_offsets.initialize(m_k + 1, 0);
            m_degrees.initialize(m_k, 0);
            m_seen.resize(num_threads);
            m_epoch.assign(num_threads, 0);
            for (u64 t = 0; t < num_threads; ++t) {
                m_seen[t].set_pool(m_pool);
                m_seen[t].initialize(m_k, 0);
            }
        }

        template<typename GraphT, typename PartitionManagerT, typename BoundaryVertexManagerT>
        void compute(const GraphT &g,
                     const PartitionManagerT &p_manager,
                     const BoundaryVertexManagerT &bv_manager,
                     const u64 num_threads) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "BlockAdjacency", "compute");

            if (m_seen.size() < num_threads) {
                m_seen.resize(num_threads);
                m_epoch.resize(num_threads, 0);
                for (u64 t = 0; t < num_threads; ++t) {
                    if (m_seen[t].size() < m_k) {
                        m_seen[t].set_pool(m_pool);
                        m_seen[t].initialize(m_k, 0);
                    }
                }
            }

            // Pass 1: count degrees
            #pragma omp parallel num_threads(num_threads)
            {
                const u64 tid = omp_get_thread_num();
                auto &seen = m_seen[tid];
                u32 &epoch = m_epoch[tid];

                #pragma omp for schedule(static)
                for (partition_t id = 0; id < m_k; ++id) {
                    if (bv_manager.size(id) == 0) {
                        m_degrees[id] = 0;
                        continue;
                    }
                    epoch++;
                    if (epoch == 0) { seen.fill(0); epoch = 1; }

                    size_t deg = 0;
                    for (const vertex_t u : bv_manager.boundary(id)) {
                        for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                            const partition_t v_id = p_manager[g.edges_v[i]];
                            if (id != v_id && seen[v_id] != epoch) {
                                seen[v_id] = epoch;
                                deg++;
                            }
                        }
                    }
                    m_degrees[id] = deg;
                }
            }

            // Prefix sum for offsets
            m_offsets[0] = 0;
            for (partition_t id = 0; id < m_k; ++id) {
                m_offsets[id + 1] = m_offsets[id] + m_degrees[id];
            }

            // Allocate flat array
            const size_t total = m_offsets[m_k];
            if (total > m_data_capacity) {
                m_data.initialize(total);
                m_data_capacity = total;
            }

            // Pass 2: fill
            #pragma omp parallel num_threads(num_threads)
            {
                const u64 tid = omp_get_thread_num();
                auto &seen = m_seen[tid];
                u32 &epoch = m_epoch[tid];

                #pragma omp for schedule(static)
                for (partition_t id = 0; id < m_k; ++id) {
                    if (m_degrees[id] == 0) continue;
                    epoch++;
                    if (epoch == 0) { seen.fill(0); epoch = 1; }

                    size_t pos = m_offsets[id];
                    for (const vertex_t u : bv_manager.boundary(id)) {
                        for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                            const partition_t v_id = p_manager[g.edges_v[i]];
                            if (id != v_id && seen[v_id] != epoch) {
                                seen[v_id] = epoch;
                                m_data[pos++] = Entry{v_id};
                            }
                        }
                    }
                }
            }
        }

        Span neighbors(const partition_t x) const {
            return {m_data.get_ptr() + m_offsets[x], m_data.get_ptr() + m_offsets[x] + m_degrees[x]};
        }
        size_t degree(const partition_t x) const { return m_degrees[x]; }

        template<typename F>
        void for_each_neighbor(const partition_t x, F &&f) const {
            const Entry *ptr = m_data.get_ptr() + m_offsets[x];
            for (size_t i = 0; i < m_degrees[x]; ++i) {
                f(ptr[i].target, weight_t(1));
            }
        }

        size_t edge_index(const partition_t u_id, const partition_t v_id) const {
            const partition_t min_part = std::min(u_id, v_id);
            const partition_t max_part = std::max(u_id, v_id);
            return static_cast<size_t>(min_part) * static_cast<size_t>(m_k) + static_cast<size_t>(max_part);
        }
    };
}

#endif // HEIPROMAP_BLOCK_ADJACENCY_H
