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
#include "../interfaces/ISerialBoundaryVertexManager.h"

namespace HeiProMap {
    class BoundaryVertexManager final : public ISerialBoundaryVertexManager {
        vertex_t    m_n = 0;
        partition_t m_k = 0;

        vertex_t m_n_boundary        = 0;
        vertex_t *m_n_boundary_edges = nullptr;

        vertex_t **m_boundaries     = nullptr;
        size_t   *m_boundaries_size = nullptr;
        size_t   *m_vertex_idx      = nullptr;

        vertex_t *m_complete_boundary            = nullptr;
        size_t   *m_complete_boundary_vertex_idx = nullptr;
        size_t   m_complete_boundary_size        = 0;

    public:
        ~BoundaryVertexManager() override {
            free(m_n_boundary_edges);
            free(m_vertex_idx);
            free(m_complete_boundary);
            free(m_complete_boundary_vertex_idx);
            for (partition_t i = 0; i < m_k; ++i) {
                free(m_boundaries[i]);
            }
            free(m_boundaries);
            free(m_boundaries_size);
        }

        void initialize(const vertex_t t_n,
                        const partition_t t_k) override {
            m_n = t_n;
            m_k = t_k;

            vertex_t m_n_64 = round_up_64(m_n);
            m_n_boundary_edges = (vertex_t *) aligned_alloc(64, m_n_64 * sizeof(vertex_t));
            std::fill_n(m_n_boundary_edges, m_n_64, 0);

            vertex_t m_k_64 = round_up_64(m_k);
            m_boundaries = (vertex_t **) aligned_alloc(64, m_k_64 * sizeof(vertex_t *));
            std::fill_n(m_boundaries, m_k_64, nullptr);
            for (partition_t i = 0; i < m_k; i++) {
                m_boundaries[i] = (vertex_t *) aligned_alloc(64, m_n_64 * sizeof(vertex_t));
            }
            m_boundaries_size = (size_t *) aligned_alloc(64, m_k_64 * sizeof(size_t));
            std::fill_n(m_boundaries_size, m_k_64, 0);
            m_vertex_idx = (size_t *) aligned_alloc(64, m_n_64 * sizeof(size_t));

            m_complete_boundary            = (vertex_t *) aligned_alloc(64, m_n_64 * sizeof(vertex_t));
            m_complete_boundary_vertex_idx = (size_t *) aligned_alloc(64, m_n_64 * sizeof(size_t));
            m_complete_boundary_size       = 0;
        }

        size_t size() const override { return m_n_boundary; }

        size_t size(const partition_t id) const override { return m_boundaries_size[id]; }

        vertex_t get(const size_t i) const override { return m_complete_boundary[i]; }

        vertex_t get(const partition_t id, const size_t i) const override { return m_boundaries[id][i]; }

        bool is_boundary(const vertex_t u) const override { return m_n_boundary_edges[u] > 0; }

        void add(const vertex_t u,
                 const partition_t id) override {
            if (m_n_boundary_edges[u] == 0) {
                m_n_boundary += 1;
                m_boundaries[id][m_boundaries_size[id]++] = u;
                m_vertex_idx[u]                           = m_boundaries_size[id] - 1;

                m_complete_boundary[m_complete_boundary_size] = u;
                m_complete_boundary_vertex_idx[u]             = m_complete_boundary_size;
                m_complete_boundary_size += 1;
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
                    if(v_id != id){
                        n_neighbors += 1;
                        add(v, v_id);
                    }
                }
            endfor
            if(n_neighbors > 0) {
                m_n_boundary += 1;
                m_n_boundary_edges[u] = n_neighbors;
                m_boundaries[id][m_boundaries_size[id]++] = u;
                m_vertex_idx[u]                           = m_boundaries_size[id] - 1;

                m_complete_boundary[m_complete_boundary_size] = u;
                m_complete_boundary_vertex_idx[u]             = m_complete_boundary_size;
                m_complete_boundary_size += 1;
            }
        }

        void move(const graph_t &g,
                  const p_manager_t &p_manager,
                  const vertex_t u,
                  const partition_t old_id,
                  const partition_t new_id) override {
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
                        m_n_boundary -= m_n_boundary_edges[u] == 0; // decrease by 1 if it is 0
                        m_n_boundary -= m_n_boundary_edges[v] == 0; // decrease by 1 if it is 0
                        if (m_n_boundary_edges[v] == 0) {
                            remove(v, v_id);
                            remove_from_complete(v);
                        }
                    } else if (v_id == old_id) {
                        // u was moved to a different block, both gain 1 edge
                        m_n_boundary_edges[u] += 1;
                        m_n_boundary_edges[v] += 1;
                        m_n_boundary += m_n_boundary_edges[u] == 1; // add by 1 if 0
                        m_n_boundary += m_n_boundary_edges[v] == 1; // add by 1 if 0
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
                                  const p_manager_t &p_manager) override {
            // compute all from scratch
            std::fill_n(m_n_boundary_edges, m_n, 0);
            std::fill_n(m_boundaries_size, m_k, 0);
            m_complete_boundary_size = 0;
            m_n_boundary             = 0;

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
                        m_n_boundary_edges[u]                         = n_different;
                        m_boundaries[u_id][m_boundaries_size[u_id]++] = u;
                        m_vertex_idx[u]                               = m_boundaries_size[u_id] - 1;

                        m_complete_boundary[m_complete_boundary_size] = u;
                        m_complete_boundary_vertex_idx[u]             = m_complete_boundary_size;
                        m_complete_boundary_size += 1;

                        m_n_boundary += 1;
                    }
                }
            endfor
        }

    private:
        void remove_from_complete(const vertex_t u) {
            vertex_t last_vertex = m_complete_boundary[m_complete_boundary_size - 1];
            m_complete_boundary_size -= 1;
            size_t u_idx = m_complete_boundary_vertex_idx[u];

            m_complete_boundary[u_idx]                  = last_vertex;
            m_complete_boundary_vertex_idx[last_vertex] = u_idx;
        }

        void emplace_in_complete(const vertex_t u) {
            m_complete_boundary[m_complete_boundary_size] = u;
            m_complete_boundary_vertex_idx[u]             = m_complete_boundary_size;
            m_complete_boundary_size += 1;
        }

        void remove(const vertex_t u,
                    const partition_t id) const {
            vertex_t last_vertex = m_boundaries[id][--m_boundaries_size[id]];
            size_t   u_idx       = m_vertex_idx[u];

            m_boundaries[id][u_idx]   = last_vertex;
            m_vertex_idx[last_vertex] = u_idx;
        }

        void emplace(const vertex_t u,
                     const partition_t id) const {
            m_boundaries[id][m_boundaries_size[id]++] = u;
            m_vertex_idx[u]                           = m_boundaries_size[id] - 1;
        }
    };
}

#endif //HEIPROMAP_BOUNDARY_VERTEX_MANGER_H
