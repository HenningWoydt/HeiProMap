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

#ifndef HEIPROMAP_DEEP_BOUNDARY_VERTEX_MANAGER_H
#define HEIPROMAP_DEEP_BOUNDARY_VERTEX_MANAGER_H

#include "../distance_oracle.h"
#include "../../serial_definitions_1.h"
#include "../../serial_definitions_2.h"
#include "../../../commons/definitions.h"
#include "../../../commons/macros.h"

namespace HeiProMap {
    class DeepBoundaryVertexManager {
        vertex_t    m_n = 0;                                    // number of vertices in the graph
        partition_t m_k = 0;                                    // number of partitions

        AlignedArray<vertex_t> m_boundary;               // complete boundary
        AlignedArray<vertex_t> m_vertex_idx;             // index of the vertex in the complete boundary array for each vertex
        AlignedArray<vertex_t> m_n_boundary_edges;       // number of boundary edges for each vertex
        size_t                 m_boundary_size = 0;      // size of the complete boundary array

        std::vector<std::vector<vertex_t>> m_sub_boundary;  // boundary for each block
        AlignedArray<vertex_t>             m_sub_vertex_idx;            // index of the vertex in the sub boundary array

    public:
        void initialize(const vertex_t t_n,
                        const partition_t t_k) {
            m_n = t_n;
            m_k = t_k;

            m_boundary.initialize(m_n);
            m_vertex_idx.initialize(m_n);
            m_n_boundary_edges.initialize(m_n, 0);
            m_boundary_size = 0;

            m_sub_boundary.resize(m_k);
            m_sub_vertex_idx.initialize(m_n);
        }

        size_t size() const { return m_boundary_size; }

        size_t size(const partition_t id) const {
            return m_sub_boundary[id].size();
        }

        vertex_t get(const size_t i) const { return m_boundary[i]; }

        vertex_t get(const partition_t id, const size_t i) const {
            return m_sub_boundary[id][i];
        }

        bool is_boundary(const vertex_t u) const { return m_n_boundary_edges[u] > 0; }

        void add(const vertex_t u,
                 const partition_t id) {}

        void add_new(const graph_t &g,
                     const deep_p_manager_t &p_manager,
                     const vertex_t u,
                     const partition_t u_id) {}

        void move(const graph_t &g,
                  const deep_p_manager_t &p_manager,
                  const vertex_t u,
                  const partition_t old_id,
                  const partition_t new_id) {
            if (is_boundary(u)) {
                // remove from the old sub boundary
                size_t idx = m_sub_vertex_idx[u];
                ASSERT(idx < m_sub_boundary[old_id].size());
                m_sub_boundary[old_id][idx]                   = m_sub_boundary[old_id].back();
                m_sub_vertex_idx[m_sub_boundary[old_id][idx]] = idx;
                m_sub_boundary[old_id].pop_back();
            }

            forall_guiv(g, u, i, v)
                {
                    partition_t v_id = p_manager[v];
                    if (v_id == new_id) { remove_edge(u, old_id, new_id, v, v_id); }
                    else if (v_id == old_id) { add_edge(u, old_id, new_id, v, v_id); }
                }
            endfor

            if (is_boundary(u)) {
                // add to the new sub boundary
                m_sub_boundary[new_id].emplace_back(u);
                m_sub_vertex_idx[u] = m_sub_boundary[new_id].size() - 1;
            }
        }

        void compute_from_scratch(const graph_t &g,
                                  const deep_p_manager_t &p_manager) {
            // compute all from scratch
            m_n_boundary_edges.initialize(m_n, 0);
            m_boundary_size = 0;
            for (auto &vec: m_sub_boundary) { vec.clear(); }

            forall_gu(g, u)
                {
                    partition_t u_id = p_manager[u];
                    forall_guiv(g, u, i, v)
                        {
                            partition_t v_id = p_manager[v];
                            if (u_id != v_id) {
                                if (m_n_boundary_edges[u] == 0) {
                                    // add to the complete boundary
                                    m_boundary[m_boundary_size] = u;
                                    m_vertex_idx[u]             = m_boundary_size;
                                    m_boundary_size += 1;
                                }
                                m_n_boundary_edges[u] += 1;
                            }
                        }
                    endfor
                    if (is_boundary(u)) {
                        // add to the new sub boundary
                        m_sub_boundary[u_id].emplace_back(u);
                        m_sub_vertex_idx[u] = m_sub_boundary[u_id].size() - 1;
                    }
                }
            endfor
        }

        void reset() {
            m_n_boundary_edges.initialize(m_n, 0);
            m_boundary_size = 0;
            for (auto &vec: m_sub_boundary) { vec.clear(); }
        }

    private:

        void add_edge(const vertex_t u,
                      const partition_t u_old_id,
                      const partition_t u_new_id,
                      const vertex_t v,
                      const partition_t v_id) {
            ASSERT(u != v);
            ASSERT(u_old_id < m_k);
            ASSERT(u_new_id < m_k);
            ASSERT(v_id < m_k);
            ASSERT(u_old_id == v_id);
            ASSERT(u_new_id != v_id);

            if (m_n_boundary_edges[u] == 0) {
                // add to the complete boundary
                m_boundary[m_boundary_size] = u;
                m_vertex_idx[u]             = m_boundary_size;
                m_boundary_size += 1;
            }
            m_n_boundary_edges[u] += 1;

            if (m_n_boundary_edges[v] == 0) {
                // add to the complete boundary
                m_boundary[m_boundary_size] = v;
                m_vertex_idx[v]             = m_boundary_size;
                m_boundary_size += 1;

                // add to the sub boundary
                m_sub_boundary[v_id].emplace_back(v);
                m_sub_vertex_idx[v] = m_sub_boundary[v_id].size() - 1;
            }
            m_n_boundary_edges[v] += 1;
        }

        void remove_edge(const vertex_t u,
                         const partition_t u_old_id,
                         const partition_t u_new_id,
                         const vertex_t v,
                         const partition_t v_id) {
            ASSERT(u != v);
            ASSERT(u_old_id < m_k);
            ASSERT(u_new_id < m_k);
            ASSERT(v_id < m_k);
            ASSERT(u_old_id != v_id);
            ASSERT(u_new_id == v_id);
            ASSERT(m_n_boundary_edges[u] > 0);
            ASSERT(m_n_boundary_edges[v] > 0);

            if (m_n_boundary_edges[u] == 1) {
                // remove from the complete boundary
                size_t idx = m_vertex_idx[u];
                m_boundary[idx]               = m_boundary[m_boundary_size - 1];
                m_vertex_idx[m_boundary[idx]] = idx;
                m_boundary_size -= 1;
            }
            m_n_boundary_edges[u] -= 1;

            if (m_n_boundary_edges[v] == 1) {
                // remove from the complete boundary
                size_t idx = m_vertex_idx[v];
                m_boundary[idx]               = m_boundary[m_boundary_size - 1];
                m_vertex_idx[m_boundary[idx]] = idx;
                m_boundary_size -= 1;

                // remove from the sub boundary
                idx = m_sub_vertex_idx[v];
                ASSERT(idx < m_sub_boundary[v_id].size());
                m_sub_boundary[v_id][idx]                   = m_sub_boundary[v_id].back();
                m_sub_vertex_idx[m_sub_boundary[v_id][idx]] = idx;
                m_sub_boundary[v_id].pop_back();
            }
            m_n_boundary_edges[v] -= 1;
        }
    };
}

#endif //HEIPROMAP_DEEP_BOUNDARY_VERTEX_MANAGER_H
