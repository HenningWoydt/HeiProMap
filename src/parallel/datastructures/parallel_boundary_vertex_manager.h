#ifndef HEIDELBERGPROCESSMAPPING_PARALLEL_BOUNDARY_VERTEX_MANAGER_H
#define HEIDELBERGPROCESSMAPPING_PARALLEL_BOUNDARY_VERTEX_MANAGER_H

#include <atomic>

#include "../../definitions.h"
#include "../../macros.h"
#include "../../interfaces/IBoundaryVertexManager.h"
#include "../interfaces/IParallelBoundaryVertexManager.h"

namespace HeiProMap {

    template<typename TParallelGraph, typename TParallelActiveVertexManager, typename TParallelPartitionManager>
    class ParallelBoundaryVertexManager : public IParallelBoundaryVertexManager<TParallelGraph, TParallelActiveVertexManager, TParallelPartitionManager> {
    private:
        TParallelGraph *m_p_g = nullptr;
        TParallelActiveVertexManager *m_p_av_manager = nullptr;
        TParallelPartitionManager *m_p_p_manger = nullptr;

        partition_t m_k = 0;

        std::atomic<vertex_t> m_n_boundary = 0;

        // Current number of edges to other boundary vertices
        std::vector<vertex_t> m_n_boundary_edges;

        // Current boundary for each block. Can be inconsistent!
        std::vector<std::vector<vertex_t>> m_boundaries;

        u64 m_n_threads = 1;

        // iteration
        size_t iterator_idx = 0;
        std::vector<size_t> iterator_indices;

        // block iteration
        std::vector<size_t> b_iterator_indices;

    public:
        // initialization
        void initialize(TParallelGraph *t_p_g,
                        TParallelActiveVertexManager *t_p_av_manager,
                        TParallelPartitionManager *t_p_p_manger,
                        partition_t t_k,
                        u64 n_threads) final {
            ASSERT(t_p_g != nullptr);
            ASSERT(t_p_av_manager != nullptr);
            ASSERT(t_p_p_manger != nullptr);

            m_p_g = t_p_g;
            m_p_av_manager = t_p_av_manager;
            m_p_p_manger = t_p_p_manger;
            m_k = t_k;

            m_n_boundary_edges.resize(m_p_g->get_n(), 0);
            m_boundaries.resize(m_k);

            m_n_threads = n_threads;

            // iteration
            iterator_idx = 0;
            iterator_indices.resize(m_k, 0);

            // block iteration
            b_iterator_indices.resize(m_k, 0);
        }

        // add
        vertex_t get_n_boundary() const final { return m_n_boundary; }

        vertex_t get_vertex(size_t idx) const final {
            for(size_t i = 0; i < m_k; ++i){
                if(idx < m_boundaries[i].size()){
                    return m_boundaries[i][idx];
                }
                idx -= m_boundaries[i].size();
            }
            ASSERT(false);
            abort();
        }

        void insert(vertex_t u, partition_t id) final {
            ASSERT(u < m_p_g->get_n());
            ASSERT(id < m_k);

            // add connections to other boundary vertices
            bool is_boundary = false;
            for (size_t i = 0; i < m_p_g->size(u); ++i) {
                vertex_t v = m_p_g->neighbor(u, i);
                partition_t v_id = (*m_p_p_manger)[v];
                if (v_id != id) {
                    // u and v are boundary vertices in different blocks
                    m_n_boundary_edges[u] += 1;
                    is_boundary = true;
                }
            }

            if (is_boundary) {
                // put u in its respective boundary queue
                m_n_boundary += 1;
                m_boundaries[id].emplace_back(u);
            }
        }

        void move(vertex_t u, partition_t old_id, partition_t new_id) final {
            ASSERT(u < m_p_g->get_n());
            ASSERT(new_id < m_k);
            ASSERT(new_id != old_id);

            remove_if_exists(old_id, u);

            // new boundary vertices could be discovered and other could be removed
            for (size_t i = 0; i < m_p_g->size(u); ++i) {
                vertex_t v = m_p_g->neighbor(u, i);
                partition_t v_id = (*m_p_p_manger)[v];

                if (v_id == new_id) {
                    // u was moved to the same block as v, both loose 1 boundary edge
                    ASSERT(m_n_boundary_edges[u] > 0);
                    ASSERT(m_n_boundary_edges[v] > 0);
                    m_n_boundary_edges[u] -= 1;
                    m_n_boundary_edges[v] -= 1;
                    if (m_n_boundary_edges[u] == 0) { m_n_boundary -= 1; }
                    if (m_n_boundary_edges[v] == 0) { m_n_boundary -= 1; remove_if_exists(new_id, v);}
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
                        m_boundaries[v_id].push_back(v);
                        //emplace_if_not_exists(v_id, v);
                    }
                }
                // else, v and b are in different blocks and still connected, nothing changes
            }

            if (m_n_boundary_edges[u] > 0) {
                emplace_if_not_exists(new_id, u);
            }
        }

        // check
        bool is_boundary(vertex_t u) const final {
            return m_n_boundary_edges[u] > 0;
        }

        // uncontract
        void uncontract(vertex_t u, vertex_t v) final {
            ASSERT(u != v);
            ASSERT((*m_p_p_manger)[u] == (*m_p_p_manger)[v]);
            ASSERT(m_p_av_manager->is_active(u));
            ASSERT(m_p_av_manager->is_active(v));
            ASSERT(!is_boundary(v));

            if (!is_boundary(u)) {
                return;
            }
            ASSERT(m_n_boundary_edges[u] > 0);

            partition_t v_id = (*m_p_p_manger)[v];

            for (size_t i = 0; i < m_p_g->size(v); ++i) {
                vertex_t ev = m_p_g->neighbor(v, i);
                ASSERT(m_p_av_manager->is_active(ev));

                partition_t ev_id = (*m_p_p_manger)[ev];

                if (ev_id != v_id) {
                    m_n_boundary_edges[v] += 1;
                    if (m_n_boundary_edges[v] == 1) {
                        m_n_boundary += 1;
                        m_boundaries[v_id].push_back(v);
                    }
                    if (m_p_g->edge_exists(u, ev)) {
                        m_n_boundary_edges[ev] += 1;
                        // both u and v are connected to ev, so ev gains
                    } else {
                        ASSERT(m_n_boundary_edges[u] > 0);
                        m_n_boundary_edges[u] -= 1;
                        if (m_n_boundary_edges[u] == 0) {
                            m_n_boundary -= 1;
                        }
                        // only v is connected to ev, so u looses
                    }
                }
            }
        }

        // iteration
        void reset_iterator() final {
            iterator_idx = 0;
            std::fill(iterator_indices.begin(), iterator_indices.end(), 0);
        }

        vertex_t get() final { return m_boundaries[iterator_idx][iterator_indices[iterator_idx]]; }

        void next() final { iterator_indices[iterator_idx] += 1; }

        bool available() final {
            for (partition_t i = 0; i < m_k; ++i) {
                while (iterator_indices[iterator_idx] < m_boundaries[iterator_idx].size() && !is_boundary(m_boundaries[iterator_idx][iterator_indices[iterator_idx]])) {
                    m_boundaries[iterator_idx][iterator_indices[iterator_idx]] = m_boundaries[iterator_idx].back();
                    m_boundaries[iterator_idx].pop_back();
                }

                if (iterator_indices[iterator_idx] < m_boundaries[iterator_idx].size()) {
                    return true;
                }

                iterator_idx = (iterator_idx + 1) % m_k;
            }
            return false;
        }

        // block iteration
        void reset_iterator(partition_t id) final { b_iterator_indices[id] = 0; }

        vertex_t get(partition_t id) final { return m_boundaries[id][b_iterator_indices[id]]; }

        void next(partition_t id) final { b_iterator_indices[id] += 1; }

        bool available(partition_t id) final {
            while (b_iterator_indices[id] < m_boundaries[id].size() && !is_boundary(m_boundaries[id][b_iterator_indices[id]])) {
                m_boundaries[id][b_iterator_indices[id]] = m_boundaries[id].back();
                m_boundaries[id].pop_back();
            }

            return b_iterator_indices[id] < m_boundaries[id].size();
        }

    private:
        void emplace_if_not_exists(partition_t b, vertex_t u) {
            for (size_t i = 0; i < m_boundaries[b].size(); ++i) {
                if (m_boundaries[b][i] == u) {
                    return;
                }
                /*
                if ((*m_p_p_manger)[m_boundaries[b][i]] != b) {
                    m_boundaries[b][i] = m_boundaries[b].back();
                    m_boundaries[b].pop_back();
                    i -= 1;
                    continue;
                }
                 */
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
                /*
                if ((*m_p_p_manger)[m_boundaries[b][i]] != b) {
                    m_boundaries[b][i] = m_boundaries[b].back();
                    m_boundaries[b].pop_back();
                    i -= 1;
                    continue;
                }
                 */
            }
        }
    };

}

#endif //HEIDELBERGPROCESSMAPPING_PARALLEL_BOUNDARY_VERTEX_MANAGER_H
