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

#include "distance_oracle.h"
#include "../serial_definitions_1.h"
#include "../serial_definitions_2.h"
#include "../../commons/definitions.h"
#include "../../commons/macros.h"

namespace HeiProMap {
    class BoundaryVertexManager {
        vertex_t    m_n = 0;                                    // number of vertices in the graph
        partition_t m_k = 0;                                    // number of partitions

        // vertex_t               m_n_boundary = 0;                // number of boundary vertices
        AlignedArray<vertex_t> m_n_boundary_edges;              // number of boundary edges for each vertex

        AlignedArray<vertex_t> m_boundaries;
        AlignedArray<size_t>   m_boundaries_size;               // number of boundary vertices for each partition
        AlignedArray<size_t>   m_vertex_idx;                    // index of the vertex in the boundary array for each vertex

        // AlignedArray<vertex_t> m_complete_boundary;             // complete boundary, i.e. all vertices that are in a partition
        // AlignedArray<vertex_t> m_complete_boundary_vertex_idx;  // index of the vertex in the complete boundary array for each vertex
        // size_t                 m_complete_boundary_size = 0;    // size of the complete boundary array

    public:
        void initialize(const vertex_t t_n,
                        const partition_t t_k) {
            m_n = t_n;
            m_k = t_k;

            // m_n_boundary = 0;
            m_n_boundary_edges.initialize(m_n, 0);

            size_t size = (size_t) m_k * (size_t) m_n;
            m_boundaries.initialize(size);
            m_boundaries_size.initialize(m_k, 0);
            m_vertex_idx.initialize(m_n);

            // m_complete_boundary.initialize(m_n);
            // m_complete_boundary_vertex_idx.initialize(m_n);
            // m_complete_boundary_size = 0;
        }

        size_t size(const partition_t id) const { return m_boundaries_size[id]; }

        vertex_t get(const partition_t id, const size_t i) const { return m_boundaries[id * m_n + i]; }

        partition_t get_k() const { return m_k; }

        bool is_boundary(const vertex_t u) const { return m_n_boundary_edges[u] > 0; }

        void add(const vertex_t u,
                 const partition_t id) {
            if (m_n_boundary_edges[u] == 0) {
                // m_n_boundary += 1;
                m_boundaries[id * m_n + (m_boundaries_size[id]++)] = u;
                m_vertex_idx[u]                                    = m_boundaries_size[id] - 1;

                // m_complete_boundary[m_complete_boundary_size] = u;
                // m_complete_boundary_vertex_idx[u]             = m_complete_boundary_size;
                // m_complete_boundary_size += 1;
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
                // m_n_boundary += 1;
                m_n_boundary_edges[u]                              = n_neighbors;
                m_boundaries[id * m_n + (m_boundaries_size[id]++)] = u;
                m_vertex_idx[u]                                    = m_boundaries_size[id] - 1;

                // m_complete_boundary[m_complete_boundary_size] = u;
                // m_complete_boundary_vertex_idx[u]             = m_complete_boundary_size;
                // m_complete_boundary_size += 1;
            }
        }

        void move(const graph_t &g,
                  const p_manager_t &p_manager,
                  const vertex_t u,
                  const partition_t old_id,
                  const partition_t new_id) {
            bool u_was_boundary = is_boundary(u);

            // remove u from its old id
            if (u_was_boundary) {
                remove(u, old_id);
            }

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
                        // m_n_boundary -= m_n_boundary_edges[u] == 0; // decrease by 1 if it is 0
                        // m_n_boundary -= m_n_boundary_edges[v] == 0; // decrease by 1 if it is 0
                        if (m_n_boundary_edges[v] == 0) {
                            remove(v, v_id);
                            remove_from_complete(v);
                        }
                    } else if (v_id == old_id) {
                        // u was moved to a different block, both gain 1 edge
                        m_n_boundary_edges[u] += 1;
                        m_n_boundary_edges[v] += 1;
                        // m_n_boundary += m_n_boundary_edges[u] == 1; // add by 1 if 0
                        // m_n_boundary += m_n_boundary_edges[v] == 1; // add by 1 if 0
                        if (m_n_boundary_edges[v] == 1) {
                            emplace(v, v_id);
                            emplace_in_complete(v);
                        }
                    }
                    // else, v and u are in different blocks and still connected, nothing changes
                }
            endfor

            if (m_n_boundary_edges[u] > 0) { emplace(u, new_id); } // emplace u into the sub-boundary
            if (u_was_boundary && m_n_boundary_edges[u] == 0) { remove_from_complete(u); } // if u has no more edges remove it from the complete-boundary
            if (!u_was_boundary && m_n_boundary_edges[u] > 0) { emplace_in_complete(u); }
        }

        void compute_from_scratch(const graph_t &g,
                                  const p_manager_t &p_manager) {
            // compute all from scratch
            m_n_boundary_edges.initialize(m_n, 0);
            m_boundaries_size.initialize(m_k, 0);
            // m_complete_boundary_size = 0;
            // m_n_boundary             = 0;

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
                        m_n_boundary_edges[u]                                  = n_different;
                        m_boundaries[u_id * m_n + (m_boundaries_size[u_id]++)] = u;
                        m_vertex_idx[u]                                        = m_boundaries_size[u_id] - 1;

                        // m_complete_boundary[m_complete_boundary_size] = u;
                        // m_complete_boundary_vertex_idx[u]             = m_complete_boundary_size;
                        // m_complete_boundary_size += 1;

                        // m_n_boundary += 1;
                    }
                }
            endfor
        }

        void reset() {
            m_n_boundary_edges.initialize(m_n, 0);
            m_boundaries_size.initialize(m_k, 0);
            // m_complete_boundary_size = 0;
            // m_n_boundary             = 0;
        }

    private:
        void remove_from_complete(const vertex_t u) {
            // vertex_t last_vertex = m_complete_boundary[m_complete_boundary_size - 1];
            // m_complete_boundary_size -= 1;
            // size_t u_idx = m_complete_boundary_vertex_idx[u];

            // m_complete_boundary[u_idx]                  = last_vertex;
            // m_complete_boundary_vertex_idx[last_vertex] = u_idx;
        }

        void emplace_in_complete(const vertex_t u) {
            // m_complete_boundary[m_complete_boundary_size] = u;
            // m_complete_boundary_vertex_idx[u]             = m_complete_boundary_size;
            // m_complete_boundary_size += 1;
        }

        void remove(const vertex_t u,
                    const partition_t id) {
            vertex_t last_vertex = m_boundaries[id * m_n + (--m_boundaries_size[id])];
            size_t   u_idx       = m_vertex_idx[u];

            m_boundaries[id * m_n + u_idx] = last_vertex;
            m_vertex_idx[last_vertex]      = u_idx;
        }

        void emplace(const vertex_t u,
                     const partition_t id) {
            m_boundaries[id * m_n + (m_boundaries_size[id]++)] = u;
            m_vertex_idx[u]                                    = m_boundaries_size[id] - 1;
        }
    };
}

#endif //HEIPROMAP_BOUNDARY_VERTEX_MANGER_H
