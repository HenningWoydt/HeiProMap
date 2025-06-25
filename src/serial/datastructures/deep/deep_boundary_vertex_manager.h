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
        vertex_t m_n    = 0; // number of vertices in the graph
        partition_t m_k = 0; // number of partitions

        AlignedArray<vertex_t> m_n_boundary_edges; // number of boundary edges for each vertex

        std::vector<std::vector<vertex_t>> m_sub_boundary; // boundary for each block
        AlignedArray<vertex_t> m_sub_vertex_idx; // index of the vertex in the sub boundary array

    public:
        void initialize(const vertex_t t_n,
                        const partition_t t_k) {
            m_n = t_n;
            m_k = t_k;

            m_n_boundary_edges.initialize(m_n, 0);

            m_sub_boundary.resize(m_k);
            m_sub_vertex_idx.initialize(m_n);
        }

        size_t size(const partition_t id) const { return m_sub_boundary[id].size(); }
        partition_t get_k() const { return m_k; }
        vertex_t get(const partition_t id, const size_t i) const { return m_sub_boundary[id][i]; }
        bool is_boundary(const vertex_t u) const { return m_n_boundary_edges[u] > 0; }

        void add_new(const graph_t& g,
                     const deep_p_manager_t& p_manager,
                     const vertex_t u,
                     const partition_t u_id) {
            m_n_boundary_edges[u] = 0;
            forall_guiv(g, u, i, v)
                {
                    partition_t v_id = p_manager[v];
                    if (v_id != u_id) {
                        m_n_boundary_edges[u] += 1;

                        if (m_n_boundary_edges[v] == 0) {
                            // v is now also boundary since it is connected to u
                            m_sub_boundary[v_id].emplace_back(v);
                            m_sub_vertex_idx[v] = m_sub_boundary[v_id].size() - 1;
                        }
                        m_n_boundary_edges[v] += 1;
                    }
                }
            endfor
            if (m_n_boundary_edges[u] > 0) {
                // add to the sub boundary
                m_sub_boundary[u_id].emplace_back(u);
                m_sub_vertex_idx[u] = m_sub_boundary[u_id].size() - 1;
            }
        }

        void move(const graph_t& g,
                  const deep_p_manager_t& p_manager,
                  const vertex_t u,
                  const partition_t old_id,
                  const partition_t new_id) {
            ASSERT(old_id != new_id);

            if (!is_boundary(u)) {
                add_new(g, p_manager, u, new_id);
                return;
            }

            // remove u from its old id
            remove_from_sub(u, old_id);

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
                            remove_from_sub(v, v_id);
                        }
                    } else if (v_id == old_id) {
                        // u was moved to a different block, both gain 1 edge
                        m_n_boundary_edges[u] += 1;
                        m_n_boundary_edges[v] += 1;
                        if (m_n_boundary_edges[v] == 1) {
                            emplace_in_sub(v, v_id);
                        }
                    }
                    // else, v and u are in different blocks and still connected, nothing changes
                }
            endfor

            if (m_n_boundary_edges[u] > 0) { emplace_in_sub(u, new_id); } // emplace u into the sub-boundary
        }

        void compute_from_scratch(const graph_t& g,
                                  const deep_p_manager_t& p_manager) {
            // compute all from scratch
            for (auto& vec : m_sub_boundary) { vec.clear(); }

            forall_gu(g, u)
                {
                    partition_t u_id = p_manager[u];
                    vertex_t n_edges = 0;
                    forall_guiv(g, u, i, v)
                        {
                            partition_t v_id = p_manager[v];
                            n_edges += u_id != v_id;
                        }
                    endfor
                    m_n_boundary_edges[u] = n_edges;
                    if(n_edges != 0){
                        emplace_in_sub(u, u_id);
                    }
                }
            endfor
        }

    private:
        void remove_from_sub(const vertex_t u,
                             const partition_t id) {
            vertex_t last_vertex = m_sub_boundary[id].back();
            size_t u_idx         = m_sub_vertex_idx[u];

            m_sub_boundary[id][u_idx]     = last_vertex;
            m_sub_vertex_idx[last_vertex] = u_idx;
            m_sub_boundary[id].pop_back();
        }

        void emplace_in_sub(const vertex_t u,
                            const partition_t id) {
            m_sub_boundary[id].push_back(u);
            m_sub_vertex_idx[u] = m_sub_boundary[id].size() - 1;
        }
    };
}

#endif //HEIPROMAP_DEEP_BOUNDARY_VERTEX_MANAGER_H
