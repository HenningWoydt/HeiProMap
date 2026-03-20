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

#ifndef HEIPROMAP_BOUNDARY_VERTEX_MANGER_H
#define HEIPROMAP_BOUNDARY_VERTEX_MANGER_H

#include "../definitions_1.h"
#include "../definitions_2.h"
#include "../definitions.h"
#include "../utility/macros.h"

namespace HeiProMap {
    class BoundaryVertexManager {
        vertex_t m_n = 0;
        partition_t m_k = 0;

        AlignedArray<vertex_t> m_n_boundary_edges; // number of boundary edges for each vertex
        std::vector<std::vector<vertex_t> > m_boundaries; // boundary vertices per partition
        AlignedArray<size_t> m_vertex_idx; // index of vertex inside its partition boundary vector

    public:
        void initialize(const vertex_t t_n,
                        const partition_t t_k) {
            ScopedTimer _t("io", "BoundaryVertexManager", "initialize");

            m_n = t_n;
            m_k = t_k;

            m_n_boundary_edges.initialize(m_n, 0);
            m_vertex_idx.initialize(m_n, 0);

            m_boundaries.clear();
            m_boundaries.resize(m_k);
        }

        size_t size(const partition_t id) const {
            return m_boundaries[id].size();
        }

        size_t size() const {
            size_t n = 0;
            for (partition_t id = 0; id < m_k; ++id) {
                n += m_boundaries[id].size();
            }
            return n;
        }

        vertex_t get(const partition_t id, const size_t i) const {
            return m_boundaries[id][i];
        }

        partition_t get_k() const {
            return m_k;
        }

        bool is_boundary(const vertex_t u) const {
            return m_n_boundary_edges[u] > 0;
        }

        void add(const vertex_t u,
                 const partition_t id) {
            if (m_n_boundary_edges[u] == 0) {
                m_vertex_idx[u] = m_boundaries[id].size();
                m_boundaries[id].push_back(u);
            }
            m_n_boundary_edges[u] += 1;
        }

        void add_new(const graph_t &g,
                     const p_manager_t &p_manager,
                     const vertex_t u,
                     const partition_t id) {
            u64 n_neighbors = 0;
            forall_guiv(g, u, i, v)
                {
                    partition_t v_id = p_manager[v];
                    if (v_id != id) {
                        n_neighbors += 1;
                        add(v, v_id);
                    }
                }
            endfor

            if (n_neighbors > 0) {
                m_n_boundary_edges[u] = n_neighbors;
                m_vertex_idx[u] = m_boundaries[id].size();
                m_boundaries[id].push_back(u);
            }
        }

        /**
         * O(max_deg)
         *
         * @param g
         * @param p_manager
         * @param u
         * @param old_id
         * @param new_id
         */
        void move(const graph_t &g,
                  const p_manager_t &p_manager,
                  const vertex_t u,
                  const partition_t old_id,
                  const partition_t new_id) {
            bool u_was_boundary = is_boundary(u);

            if (u_was_boundary) {
                remove(u, old_id);
            }

            forall_guiv(g, u, i, v)
                {
                    partition_t v_id = p_manager[v];

                    if (v_id == new_id) {
                        ASSERT(m_n_boundary_edges[u] > 0);
                        ASSERT(m_n_boundary_edges[v] > 0);
                        m_n_boundary_edges[u] -= 1;
                        m_n_boundary_edges[v] -= 1;
                        if (m_n_boundary_edges[v] == 0) {
                            remove(v, v_id);
                        }
                    } else if (v_id == old_id) {
                        m_n_boundary_edges[u] += 1;
                        m_n_boundary_edges[v] += 1;
                        if (m_n_boundary_edges[v] == 1) {
                            emplace(v, v_id);
                        }
                    }
                }
            endfor

            if (m_n_boundary_edges[u] > 0) {
                emplace(u, new_id);
            }
        }

        void compute_from_scratch(const graph_t &g,
                                  const p_manager_t &p_manager) {
            ScopedTimer _t("uncontraction", "BoundaryVertexManager", "compute_from_scratch");

            m_n_boundary_edges.initialize(m_n, 0);

            for (partition_t id = 0; id < m_k; ++id) {
                m_boundaries[id].clear();
            }

            forall_gu(g, u)
                {
                    size_t n_different = 0;
                    partition_t u_id = p_manager[u];

                    forall_guiv(g, u, i, v)
                        {
                            n_different += (u_id != p_manager[v]);
                        }
                    endfor

                    if (n_different > 0) {
                        m_n_boundary_edges[u] = n_different;
                        m_vertex_idx[u] = m_boundaries[u_id].size();
                        m_boundaries[u_id].push_back(u);
                    }
                }
            endfor
        }

        void add_boundary_vertex_from_count(const vertex_t u,
                                    const partition_t id,
                                    const size_t n_boundary_edges) {
            if (n_boundary_edges == 0) return;

            m_n_boundary_edges[u] = n_boundary_edges;
            m_vertex_idx[u] = m_boundaries[id].size();
            m_boundaries[id].push_back(u);
        }

        void reset() {
            ScopedTimer _t("misc", "BoundaryVertexManager", "reset");

            m_n_boundary_edges.initialize(m_n, 0);
            for (partition_t id = 0; id < m_k; ++id) {
                m_boundaries[id].clear();
            }
        }

        void copy_from(const BoundaryVertexManager &bm) {
            for (vertex_t u = 0; u < m_n; u++) {
                m_n_boundary_edges[u] = bm.m_n_boundary_edges[u];
                m_vertex_idx[u] = bm.m_vertex_idx[u];
            }
            m_boundaries = bm.m_boundaries;
        }

    private:
        void remove(const vertex_t u,
                    const partition_t id) {
            auto &boundary = m_boundaries[id];
            size_t u_idx = m_vertex_idx[u];
            vertex_t last_vertex = boundary.back();

            boundary[u_idx] = last_vertex;
            m_vertex_idx[last_vertex] = u_idx;
            boundary.pop_back();
        }

        void emplace(const vertex_t u,
                     const partition_t id) {
            m_vertex_idx[u] = m_boundaries[id].size();
            m_boundaries[id].push_back(u);
        }
    };
}

#endif //HEIPROMAP_BOUNDARY_VERTEX_MANGER_H
