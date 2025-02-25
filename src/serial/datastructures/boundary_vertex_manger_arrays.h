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

#ifndef HEIPROMAP_BOUNDARY_VERTEX_MANGER_ARRAYS_H
#define HEIPROMAP_BOUNDARY_VERTEX_MANGER_ARRAYS_H

#include "distance_oracle.h"
#include "../../definitions.h"
#include "../../macros.h"
#include "../interfaces/ISerialBoundaryVertexManager.h"

namespace HeiProMap {
    class BoundaryVertexManagerArrays final : public ISerialBoundaryVertexManager {
        vertex_t    m_n = 0;
        partition_t m_k = 0;

        vertex_t m_n_boundary        = 0;
        vertex_t *m_n_boundary_edges = nullptr;

        std::vector<std::vector<vertex_t>> m_boundaries;
        size_t                             *m_vertex_idx = nullptr;

        vertex_t *m_complete_boundary            = nullptr;
        size_t   *m_complete_boundary_vertex_idx = nullptr;
        size_t   m_complete_boundary_size        = 0;

    public:
        ~BoundaryVertexManagerArrays() override {
            free(m_n_boundary_edges);
            free(m_vertex_idx);
            free(m_complete_boundary);
            free(m_complete_boundary_vertex_idx);
        }

        class SubBoundaryVertexManager {
            std::vector<vertex_t> &m_sub_boundaries;

        public:
            explicit SubBoundaryVertexManager(std::vector<vertex_t> &t_sub_boundaries) : m_sub_boundaries(t_sub_boundaries) {}

            class SubIterator {
                std::vector<vertex_t> &m_sub_boundaries;
                size_t                m_idx;

            public:
                // Constructor
                SubIterator(std::vector<vertex_t> &t_sub_boundaries,
                            size_t t_idx) : m_sub_boundaries(t_sub_boundaries) {
                    m_idx = t_idx;
                }

                // Dereference operator
                vertex_t operator*() const { return m_sub_boundaries[m_idx]; }

                // Pre-increment operator
                SubIterator &operator++() {
                    m_idx++;
                    return *this;
                }

                bool operator!=(const SubIterator &other) const { return m_idx != other.m_idx; }
            };

            SubIterator begin() const { return {m_sub_boundaries, 0}; }

            SubIterator end() const { return {m_sub_boundaries, m_sub_boundaries.size()}; }
        };

        SubBoundaryVertexManager operator[](const vertex_t u) {
            return SubBoundaryVertexManager(m_boundaries[u]);
        }

        void initialize(const vertex_t t_n,
                        const partition_t t_k) override {
            m_n = t_n;
            m_k = t_k;

            vertex_t m_n_64 = round_up_64(m_n);
            m_n_boundary_edges = (vertex_t *) aligned_alloc(64, m_n_64 * sizeof(vertex_t));

            m_boundaries.resize(m_k);
            m_vertex_idx = (size_t *) aligned_alloc(64, m_n_64 * sizeof(size_t));

            m_complete_boundary            = (vertex_t *) aligned_alloc(64, m_n_64 * sizeof(vertex_t));
            m_complete_boundary_vertex_idx = (size_t *) aligned_alloc(64, m_n_64 * sizeof(size_t));
            m_complete_boundary_size       = 0;
        }

        vertex_t get_n_boundary() const override { return m_n_boundary; }

        bool is_boundary(const vertex_t u) const override { return m_n_boundary_edges[u] > 0; }

        void add(const vertex_t u, const partition_t id) override {
            if (m_n_boundary_edges[u] == 0) {
                m_n_boundary += 1;
                m_boundaries[id].emplace_back(u);
                m_vertex_idx[u] = m_boundaries.size() - 1;

                m_complete_boundary[m_complete_boundary_size] = u;
                m_complete_boundary_vertex_idx[u]             = m_complete_boundary_size;
                m_complete_boundary_size += 1;
            }
            m_n_boundary_edges[u] += 1;
        }

        template<typename TSerialGraph, typename TSerialPartitionManager>
        void move(TSerialGraph &g, TSerialPartitionManager &p_manager, vertex_t u, partition_t old_id, partition_t new_id) {
            ASSERT(new_id < m_k);
            ASSERT(new_id != old_id);

            // remove u from its old id
            remove(u, old_id);

            // check how many connections u still has an if the neighbor are still boundary
            for (const auto [v, w]: g[u]) {
                partition_t v_id = p_manager[v];

                if (v_id == new_id) {
                    // u was moved to the same block as v, both loose 1 edge
                    ASSERT(m_n_boundary_edges[u] > 0);
                    ASSERT(m_n_boundary_edges[v] > 0);
                    m_n_boundary_edges[u] -= 1;
                    m_n_boundary_edges[v] -= 1;
                    m_n_boundary -= m_n_boundary_edges[u] == 0; // decrease by 1 if 0
                    m_n_boundary -= m_n_boundary_edges[v] == 0; // decrease by 1 if 0
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
            if (m_n_boundary_edges[u] > 0) { emplace(u, new_id); }
            if (m_n_boundary_edges[u] == 0) { remove_from_complete(u); }
        }

        void remove_from_complete(vertex_t u) {
            vertex_t last_vertex = m_complete_boundary[m_complete_boundary_size - 1];
            m_complete_boundary_size -= 1;
            size_t u_idx = m_complete_boundary_vertex_idx[u];

            m_complete_boundary[u_idx]                  = last_vertex;
            m_complete_boundary_vertex_idx[last_vertex] = u_idx;
        }

        void emplace_in_complete(vertex_t u) {
            m_complete_boundary[m_complete_boundary_size] = u;
            m_complete_boundary_vertex_idx[u]             = m_complete_boundary_size;
            m_complete_boundary_size += 1;
        }

        void remove(vertex_t u, partition_t id) {
            vertex_t last_vertex = m_boundaries[id].back();
            m_boundaries[id].pop_back();
            size_t u_idx = m_vertex_idx[u];

            m_boundaries[id][u_idx]   = last_vertex;
            m_vertex_idx[last_vertex] = u_idx;
        }

        void emplace(vertex_t u, partition_t id) {
            m_boundaries[id].emplace_back(u);
            m_vertex_idx[u] = m_boundaries[id].size() - 1;
        }


        template<typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialPartitionManager>
        void uncontract(const EdgeUV *matches,
                        size_t &matches_size,
                        TSerialGraph &new_g, // the larger uncontracted graph
                        [[maybe_unused]] TSerialGraph &old_g, // the smaller not contracted graph
                        [[maybe_unused]] TSerialActiveVertexManager &av_manager,
                        TSerialPartitionManager &p_manager) {
            matches = ASSUME_ALIGNED(EdgeUV*, matches, 64);

            // compute all from scratch
            std::fill_n(m_n_boundary_edges, m_n, 0);
            std::fill_n(m_vertex_idx, m_n, 0);
            for (auto &vec: m_boundaries) { vec.clear(); }
            m_n_boundary = 0;

            for (vertex_t u: av_manager) {
                size_t      n_different = 0;
                partition_t u_id        = p_manager[u];

                for (size_t i = 0; i < new_g.size(u); ++i) {
                    const vertex_t v    = new_g.neighbor(u, i);
                    partition_t    v_id = p_manager[v];

                    n_different += u_id != v_id;
                }

                if (n_different > 0) {
                    m_n_boundary_edges[u] = n_different;
                    m_boundaries[u_id].push_back(u);
                    m_vertex_idx[u] = m_boundaries[u_id].size() - 1;

                    m_complete_boundary[m_complete_boundary_size] = u;
                    m_complete_boundary_vertex_idx[u]             = m_complete_boundary_size;
                    m_complete_boundary_size += 1;

                    m_n_boundary += 1;
                }
            }
        }

        class Iterator {
            vertex_t *m_complete_boundary;
            size_t   &m_complete_boundary_size;
            size_t   m_idx;

        public:
            // Constructor

            Iterator(vertex_t *t_complete_boundary,
                     size_t &t_complete_boundary_size) : m_complete_boundary_size(t_complete_boundary_size) {
                m_complete_boundary = ASSUME_ALIGNED(vertex_t*, t_complete_boundary, 64);
                m_idx               = 0;
            }

            // Dereference operator
            vertex_t operator*() const { return m_complete_boundary[m_idx]; }

            // Pre-increment operator
            Iterator &operator++() {
                m_idx++;
                return *this;
            }

            bool operator!=(const Iterator &other) const { return m_idx != other.m_complete_boundary_size; }
        };

        Iterator begin() { return {m_complete_boundary, m_complete_boundary_size}; }

        Iterator end() { return {m_complete_boundary, m_complete_boundary_size}; }
    };
}

#endif //HEIPROMAP_BOUNDARY_VERTEX_MANGER_ARRAYS_H
