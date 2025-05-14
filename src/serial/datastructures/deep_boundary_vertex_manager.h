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

#include "distance_oracle.h"
#include "../serial_definitions_1.h"
#include "../serial_definitions_2.h"
#include "../../commons/definitions.h"
#include "../../commons/macros.h"
#include "../interfaces/ISerialBoundaryVertexManager.h"

namespace HeiProMap {
    class DeepBoundaryVertexManager final : public ISerialBoundaryVertexManager {
        vertex_t    m_n = 0;                                    // number of vertices in the graph
        partition_t m_k = 0;                                    // number of partitions

        AlignedArray<vertex_t>    m_boundary;             // complete boundary, i.e. all vertices that are in a partition
        AlignedArray<vertex_t>    m_vertex_idx;  // index of the vertex in the complete boundary array for each vertex
        AlignedArray<vertex_t>    m_n_boundary_edges;              // number of boundary edges for each vertex
        AlignedArray<partition_t> m_partition;
        size_t                    m_boundary_size = 0;    // size of the complete boundary array

    public:
        void initialize(const vertex_t t_n,
                        const partition_t t_k) override {
            m_n = t_n;
            m_k = t_k;

            m_boundary.initialize(m_n);
            m_vertex_idx.initialize(m_n);
            m_n_boundary_edges.initialize(m_n, 0);
            m_partition.initialize(m_n);
            m_boundary_size = 0;
        }

        size_t size() const override { return m_boundary_size; }

        size_t size(const partition_t id) const override {
            size_t      size = 0;
            for (size_t i    = 0; i < m_boundary_size; ++i) {
                size += m_partition[i] == id;
            }
            return size;
        }

        vertex_t get(const size_t i) const override { return m_boundary[i]; }

        vertex_t get(const partition_t id, const size_t i) const override {
            size_t      size = 0;
            for (size_t j    = 0; j < m_boundary_size; ++j) {
                if (m_partition[j] == id) {
                    if (size == i) { return m_partition[j]; }
                    size += 1;
                }
            }
            abort();
        }

        bool is_boundary(const vertex_t u) const override { return m_n_boundary_edges[u] > 0; }

        void add(const vertex_t u,
                 const partition_t id) override {
            if (m_n_boundary_edges[u] == 0) {
                m_boundary[m_boundary_size] = u;
                m_vertex_idx[u]             = m_boundary_size;
                m_partition[u]              = id;
                m_boundary_size += 1;
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
                m_boundary[m_boundary_size] = u;
                m_vertex_idx[u]             = m_boundary_size;
                m_n_boundary_edges[u]       = n_neighbors;
                m_partition[u]              = id;
                m_boundary_size += 1;
            }
        }

        void move(const graph_t &g,
                  const p_manager_t &p_manager,
                  const vertex_t u,
                  const partition_t old_id,
                  const partition_t new_id) override {

            bool u_was_boundary = is_boundary(u);

            m_partition[u] = new_id;

            // check how many connections u still has and if the neighbor are still boundary
            forall_guiv(g, u, i, v)
                {
                    partition_t v_id = p_manager[v];

                    if (v_id == new_id) {
                        // u was moved to the same block as v, both loose 1 edge
                        ASSERT(m_n_boundary_edges[u] > 0);
                        ASSERT(m_n_boundary_edges[v] > 0);
                        m_n_boundary_edges[u] -= 1;
                        m_n_boundary_edges[v] -= 1;
                        if (m_n_boundary_edges[v] == 0) {
                            vertex_t last_v             = m_boundary[m_boundary_size - 1];
                            m_vertex_idx[last_v]        = m_vertex_idx[v];
                            m_boundary[m_vertex_idx[v]] = last_v;
                            m_boundary_size - 1;
                        }
                    } else if (v_id == old_id) {
                        // u was moved to a different block, both gain 1 edge
                        m_n_boundary_edges[u] += 1;
                        m_n_boundary_edges[v] += 1;
                        if (m_n_boundary_edges[v] == 1) {
                            m_boundary[m_boundary_size] = v;
                            m_vertex_idx[v]             = m_boundary_size;
                            m_n_boundary_edges[v]       = 1;
                            m_partition[v]              = v_id;
                            m_boundary_size += 1;
                        }
                    }
                    // else, v and u are in different blocks and still connected, nothing changes
                }
            endfor

            if (u_was_boundary && m_n_boundary_edges[u] == 0) {
                // remove u
                vertex_t last_v = m_boundary[m_boundary_size - 1];
                m_vertex_idx[last_v]        = m_vertex_idx[u];
                m_boundary[m_vertex_idx[u]] = last_v;
                m_boundary_size - 1;
            } else if (!u_was_boundary && m_n_boundary_edges[u] > 0) {
                // add u
                m_boundary[m_boundary_size] = u;
                m_vertex_idx[u]             = m_boundary_size;
                m_n_boundary_edges[u]       = 1;
                m_partition[u]              = new_id;
                m_boundary_size += 1;
            }
        }

        void compute_from_scratch(const graph_t &g,
                                  const p_manager_t &p_manager) override {
            // compute all from scratch
            m_n_boundary_edges.initialize(m_n, 0);
            m_boundary_size = 0;

            forall_gu(g, u)
                {
                    size_t      n_different = 0;
                    partition_t u_id        = p_manager[u];

                    forall_guiv(g, u, i, v)
                        {
                            n_different += u_id != p_manager[v];
                        }
                    endfor

                    if (n_different > 0) {
                        m_n_boundary_edges[u] = n_different;

                        m_boundary[m_boundary_size] = u;
                        m_vertex_idx[u]             = m_boundary_size;
                        m_partition[u]              = u_id;
                        m_boundary_size += 1;
                    }
                }
            endfor
        }

        void reset() {
            m_n_boundary_edges.initialize(m_n, 0);
            m_boundary_size = 0;
        }
    };
}

#endif //HEIPROMAP_DEEP_BOUNDARY_VERTEX_MANAGER_H
