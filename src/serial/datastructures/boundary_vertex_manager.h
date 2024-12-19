#ifndef HEIDELBERGPROCESSMAPPING_BOUNDARY_VERTEX_MANAGER_H
#define HEIDELBERGPROCESSMAPPING_BOUNDARY_VERTEX_MANAGER_H

#include "distance_oracle.h"
#include "../../definitions.h"
#include "../../macros.h"
#include "../interfaces/ISerialBoundaryVertexManager.h"

namespace HeiProMap {
    template <typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialPartitionManager>
    class BoundaryVertexManager final : public ISerialBoundaryVertexManager<TSerialGraph, TSerialActiveVertexManager, TSerialPartitionManager> {
        vertex_t m_n    = 0;
        partition_t m_k = 0;

        vertex_t m_n_boundary = 0;
        std::vector<vertex_t> m_n_boundary_edges;
        std::vector<std::vector<vertex_t>> m_boundaries;

    public:
        class SubBoundaryVertexManager {
            std::vector<vertex_t>& m_sub_boundaries;
            std::vector<vertex_t>& m_n_boundary_edges;

        public:
            SubBoundaryVertexManager(std::vector<vertex_t>& m_sub_boundaries,
                                     std::vector<vertex_t>& m_n_boundary_edges) : m_sub_boundaries(m_sub_boundaries), m_n_boundary_edges(m_n_boundary_edges) {}

            class SubIterator {
                std::vector<vertex_t>& m_sub_boundaries;
                std::vector<vertex_t>& m_n_boundary_edges;
                size_t m_idx;

            public:
                // Constructor
                SubIterator(std::vector<vertex_t>& t_sub_boundaries,
                            std::vector<vertex_t>& t_n_boundary_edges): m_sub_boundaries(t_sub_boundaries), m_n_boundary_edges(t_n_boundary_edges) {
                    m_idx = 0;
                    advance_to_next_valid();
                }

                // Dereference operator
                vertex_t operator*() const { return m_sub_boundaries[m_idx]; }

                // Pre-increment operator
                SubIterator& operator++() {
                    // Remove the current inactive vertex from the vector
                    if (m_n_boundary_edges[m_sub_boundaries[m_idx]] == 0) {
                        m_sub_boundaries[m_idx] = m_sub_boundaries.back();
                        m_sub_boundaries.pop_back();
                    } else {
                        ++m_idx;
                    }

                    // Advance to the next valid vertex
                    advance_to_next_valid();
                    return *this;
                }

                bool operator!=(const SubIterator& other) const { return m_idx != other.m_sub_boundaries.size(); }

            private:
                void advance_to_next_valid() {
                    while (m_idx < m_sub_boundaries.size() && m_n_boundary_edges[m_sub_boundaries[m_idx]] == 0) {
                        m_sub_boundaries[m_idx] = m_sub_boundaries.back();
                        m_sub_boundaries.pop_back();
                    }
                }
            };

            SubIterator begin() { return {m_sub_boundaries, m_n_boundary_edges}; }
            SubIterator end() { return {m_sub_boundaries, m_n_boundary_edges}; }
        };

        void initialize(const vertex_t t_n,
                        const partition_t t_k) override {
            m_n = t_n;
            m_k = t_k;

            m_n_boundary_edges.resize(m_n, 0);
            m_boundaries.resize(m_k);
        }

        SubBoundaryVertexManager operator[](const vertex_t u) {
            return SubBoundaryVertexManager(m_boundaries[u], m_n_boundary_edges);
        }

        SubBoundaryVertexManager& operator[](const vertex_t u) const {
            return SubBoundaryVertexManager(m_boundaries[u], m_n_boundary_edges);
        }

        vertex_t get_n_boundary() const override { return m_n_boundary; }
        bool is_boundary(const vertex_t u) const override { return m_n_boundary_edges[u] > 0; }

        void add(const vertex_t u, const partition_t id) override {
            if (m_n_boundary_edges[u] == 0) {
                m_n_boundary += 1;
                m_boundaries[id].emplace_back(u);
            }
            m_n_boundary_edges[u] += 1;
        }

        void move(TSerialGraph &g, TSerialPartitionManager &p_manager, vertex_t u, partition_t old_id, partition_t new_id) {
            ASSERT(new_id < m_k);
            ASSERT(new_id != old_id);

            remove_if_exists(old_id, u);

            // new boundary vertices could be discovered and other could be removed
            for (size_t i = 0; i < g.size(u); ++i) {
                vertex_t v       = g.neighbor(u, i);
                partition_t v_id = p_manager[v];

                if (v_id == new_id) {
                    // u was moved to the same block as v, both loose 1 boundary edge
                    ASSERT(m_n_boundary_edges[u] > 0);
                    ASSERT(m_n_boundary_edges[v] > 0);
                    m_n_boundary_edges[u] -= 1;
                    m_n_boundary_edges[v] -= 1;
                    if (m_n_boundary_edges[u] == 0) { m_n_boundary -= 1; }
                    if (m_n_boundary_edges[v] == 0) { m_n_boundary -= 1; }
                }
                if (v_id == old_id) {
                    // u was moved to a different block as v, both gain 1 boundary edge
                    m_n_boundary_edges[u] += 1;
                    m_n_boundary_edges[v] += 1;
                    if (m_n_boundary_edges[u] == 1) {
                        m_n_boundary += 1;
                    }
                    if (m_n_boundary_edges[v] == 1) {
                        m_n_boundary += 1;
                        emplace_if_not_exists(v_id, v);
                    }
                }
                // else, v and u are in different blocks and still connected, nothing changes
            }

            if (m_n_boundary_edges[u] > 0) {
                emplace_if_not_exists(new_id, u);
            }
        }

        void emplace_if_not_exists(partition_t b, vertex_t u) {
            for (size_t i = 0; i < m_boundaries[b].size(); ++i) {
                if (m_boundaries[b][i] == u) {
                    return;
                }
            }
            m_boundaries[b].emplace_back(u);
        }

        void remove_if_exists(partition_t b, vertex_t u) {
            for (size_t i = 0; i < m_boundaries[b].size(); ++i) {
                if (m_boundaries[b][i] == u) {
                    m_boundaries[b][i] = m_boundaries[b].back();
                    m_boundaries[b].pop_back();
                    return;
                }
            }
        }

        void uncontract(std::vector<EdgeUV>& matches,
                        TSerialGraph& new_g, // the larger uncontracted graph
                        TSerialGraph& old_g, // the smaller not contracted graph
                        TSerialActiveVertexManager& av_manager,
                        TSerialPartitionManager& p_manager) override {
            // get current boundary vertices
            std::vector<vertex_t> curr_boundary;
            for (vertex_t u : *this) {
                curr_boundary.emplace_back(u);
                m_n_boundary_edges[u] = 0;
            }

            // remove all
            for (auto& vec : m_boundaries) { vec.clear(); }
            m_n_boundary = 0;

            // add all second matched vertices
            for (auto& [u, v] : matches) { curr_boundary.emplace_back(v); }

            // check all active vertices
            for (vertex_t u : curr_boundary) {
                for (size_t i = 0; i < new_g.size(u); ++i) {
                    const vertex_t v = new_g.neighbor(u, i);
                    if (p_manager[u] != p_manager[v]) {
                        add(u, p_manager[u]);
                    }
                }
            }
        }

        class Iterator {
            std::vector<std::vector<vertex_t>>& m_boundaries;
            std::vector<vertex_t>& m_n_boundary_edges;
            size_t m_partition_idx;
            size_t m_idx;

        public:
            // Constructor

            Iterator(std::vector<std::vector<vertex_t>>& t_boundaries,
                     std::vector<vertex_t>& t_n_boundary_edges): m_boundaries(t_boundaries), m_n_boundary_edges(t_n_boundary_edges) {
                m_partition_idx = 0;
                m_idx           = 0;
                advance_to_next_valid();
            }

            // Dereference operator
            vertex_t operator*() const { return m_boundaries[m_partition_idx][m_idx]; }

            // Pre-increment operator
            Iterator& operator++() {
                ++m_idx;
                advance_to_next_valid();
                return *this;
            }

            bool operator!=(const Iterator& other) const { return m_partition_idx != other.m_boundaries.size(); }

        private:
            void advance_to_next_valid() {
                while (m_partition_idx < m_boundaries.size()) {
                    while (m_idx < m_boundaries[m_partition_idx].size() && m_n_boundary_edges[m_boundaries[m_partition_idx][m_idx]] == 0) {
                        m_boundaries[m_partition_idx][m_idx] = m_boundaries[m_partition_idx].back();
                        m_boundaries[m_partition_idx].pop_back();
                    }

                    if (m_idx >= m_boundaries[m_partition_idx].size()) {
                        m_partition_idx += 1;
                        m_idx = 0;
                    } else {
                        break;
                    }
                }
            }
        };

        Iterator begin() { return {m_boundaries, m_n_boundary_edges}; }
        Iterator end() { return {m_boundaries, m_n_boundary_edges}; }
    };
}

#endif //HEIDELBERGPROCESSMAPPING_BOUNDARY_VERTEX_MANAGER_H
