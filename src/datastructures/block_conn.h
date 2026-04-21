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

#ifndef HEIPROMAP_BLOCK_CONN_H
#define HEIPROMAP_BLOCK_CONN_H

#include "csr_graph.h"
#include "distance_oracle.h"
#include "partition_manager.h"
#include "../definitions.h"
#include "../utility/macros.h"

namespace HeiProMap {
    class BlockConn {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;

        AlignedArray<size_t> m_sizes; // number of conns for each vertex
        AlignedArray<size_t> m_start; // start idx in total array
        AlignedArray<partition_t> m_arr_ids;
        AlignedArray<weight_t> m_arr_weights;
        size_t total_size = 0;

    public:
        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k) {
            ScopedTimer _t("misc", "BlockConn", "initialize");
            m_n = t_n;
            m_m = t_m;
            m_k = t_k;

            m_sizes.initialize(m_n);
            m_start.initialize(m_n);
            m_arr_ids.initialize(m_m);
            m_arr_weights.initialize(m_m);
            std::fill_n(m_sizes.get_ptr(), m_n, 0);
            total_size = 0;
        }

        size_t size(const vertex_t u) const { return m_sizes[u]; }
        size_t start(const vertex_t u) const { return m_start[u]; }
        size_t end(const vertex_t u) const { return m_start[u] + m_sizes[u]; }

        partition_t get_id(const size_t i) const { return m_arr_ids[i]; }
        weight_t get_w(const size_t i) const { return m_arr_weights[i]; }

        /**
         * O(max_deg*max_conn)
         *
         * @param g
         * @param u
         * @param old_id
         * @param new_id
         */
        void move(const graph_t &g,
                  const vertex_t u,
                  const partition_t old_id,
                  const partition_t new_id) {
            if (old_id == new_id) return;

            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) { const vertex_t v = g.edges_v[i]; const weight_t w = g.edges_w[i];
                {
                    update(v, old_id, new_id, w);
                }
            }
        }

        void compute_from_scratch(const graph_t &g, const p_manager_t &p_manager) {
            ScopedTimer _t("uncontraction", "BlockConn", "compute_from_scratch");

            m_sizes.initialize(g.n);
            m_start.initialize(g.n);
            m_arr_ids.initialize(g.m);
            m_arr_weights.initialize(g.m);
            std::fill_n(m_sizes.get_ptr(), g.n, 0);

            total_size = 0;
            for (vertex_t u = 0; u < g.n; ++u) {
                {
                    m_start[u] = total_size;
                    total_size += std::min(m_k, g.deg(u));
                    for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) { const vertex_t v = g.edges_v[i]; const weight_t w = g.edges_w[i];
                        {
                            partition_t v_id = p_manager[v];
                            add(u, v_id, w);
                        }
                    }
                }
            }
        }

        void copy_from(const BlockConn &bm) {
            m_n = bm.m_n;
            m_k = bm.m_k;
            total_size = bm.total_size;

            m_sizes.initialize(m_n);
            m_start.initialize(m_n);
            m_arr_ids.initialize(total_size);
            m_arr_weights.initialize(total_size);

            for (vertex_t u = 0; u < m_n; u++) {
                m_sizes[u] = bm.m_sizes[u];
                m_start[u] = bm.m_start[u];
            }
            for (size_t i = 0; i < total_size; i++) {
                m_arr_ids[i] = bm.m_arr_ids[i];
                m_arr_weights[i] = bm.m_arr_weights[i];
            }
        }

        void reset_build() {
            std::fill_n(m_sizes.get_ptr(), m_n, 0);
            total_size = 0;
        }

        void begin_vertex(const graph_t &g, const vertex_t u) {
            m_start[u] = total_size;
            total_size += std::min(m_k, g.deg(u));
        }

        void add_connection(const vertex_t u, const partition_t id, const weight_t w) {
            add(u, id, w);
        }

    private:
        void add(const vertex_t u, const partition_t id, const weight_t w) {
            size_t start = m_start[u];
            size_t end = start + m_sizes[u];

            bool found = false;
            for (size_t l = start; l < end; ++l) {
                if (m_arr_ids[l] == id) {
                    m_arr_weights[l] += w;
                    found = true;
                    break;
                }
            }
            if (!found) {
                m_arr_ids[end] = id;
                m_arr_weights[end] = w;
                m_sizes[u] += 1;
            }
        }

        void remove(const vertex_t u, const partition_t id, const weight_t w) {
            size_t start = m_start[u];
            size_t end = start + m_sizes[u];

            for (size_t l = start; l < end; ++l) {
                if (m_arr_ids[l] == id) {
                    m_arr_weights[l] -= w;

                    if (m_arr_weights[l] == 0) {
                        std::swap(m_arr_ids[l], m_arr_ids[end - 1]);
                        std::swap(m_arr_weights[l], m_arr_weights[end - 1]);
                        m_sizes[u] -= 1;
                    }

                    break;
                }
            }
        }

        void update(const vertex_t u,
                    const partition_t old_id,
                    const partition_t new_id,
                    const weight_t w) {
            size_t start = m_start[u];
            size_t end = start + m_sizes[u];

            size_t old_pos = end;
            size_t new_pos = end;

            for (size_t l = start; l < end; ++l) {
                if (m_arr_ids[l] == old_id) {
                    old_pos = l;
                }
                if (m_arr_ids[l] == new_id) {
                    new_pos = l;
                }
            }

            // old_id should normally exist
            if (old_pos == end) {
                return;
            }

            // Case 1: old_id and new_id are both present
            if (new_pos != end) {
                m_arr_weights[old_pos] -= w;
                m_arr_weights[new_pos] += w;

                if (m_arr_weights[old_pos] == 0) {
                    // If new_pos was the last live entry and old_pos is removed,
                    // the swap moves new_id from last -> old_pos.
                    size_t last = end - 1;
                    if (old_pos != last) {
                        std::swap(m_arr_ids[old_pos], m_arr_ids[last]);
                        std::swap(m_arr_weights[old_pos], m_arr_weights[last]);
                    }
                    m_sizes[u] -= 1;
                }
                return;
            }

            // Case 2: new_id not present
            // We can reuse the old slot if old_id disappears completely.
            m_arr_weights[old_pos] -= w;

            if (m_arr_weights[old_pos] == 0) {
                m_arr_ids[old_pos] = new_id;
                m_arr_weights[old_pos] = w;
            } else {
                m_arr_ids[end] = new_id;
                m_arr_weights[end] = w;
                m_sizes[u] += 1;
            }
        }
    };
}

#endif //HEIPROMAP_BLOCK_CONN_H
