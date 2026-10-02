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

#ifndef HEIPROMAP_LARGE_QUOTIENT_GRAPH_H
#define HEIPROMAP_LARGE_QUOTIENT_GRAPH_H

#include <algorithm>
#include <tuple>
#include <utility>
#include <vector>

#include "../definitions.h"
#include "../utility/aligned_array.h"
#include "../utility/macros.h"
#include "../utility/profiler.h"

namespace HeiProMap {
    /**
     * LargeQuotientGraph represents the quotient graph using sparse adjacency lists
     * per block, achieving O(k + |E_Q|) memory instead of O(k^2) dense matrix storage.
     * Suitable for very large numbers of blocks k.
     */
    class LargeQuotientGraph {
    public:
        struct HalfEdge {
            partition_t target = 0;
            weight_t weight = 0;
            mutable u32 used_epoch = 0;
        };

    private:
        partition_t m_k = 0;
        std::vector<std::vector<HalfEdge> > m_adj;
        AlignedArray<weight_t> m_self_weights;
        u64 m_total_half_edges = 0;

        static auto find_edge(std::vector<HalfEdge> &vec, const partition_t target) {
            return std::find_if(vec.begin(), vec.end(), [target](const HalfEdge &e) {
                return e.target == target;
            });
        }

        static auto find_edge(const std::vector<HalfEdge> &vec, const partition_t target) {
            return std::find_if(vec.begin(), vec.end(), [target](const HalfEdge &e) {
                return e.target == target;
            });
        }

        void add_half_edge(const partition_t u_id, const partition_t v_id, const weight_t w) {
            auto &vec = m_adj[u_id];
            auto it = find_edge(vec, v_id);
            if (it != vec.end()) {
                it->weight += w;
            } else {
                vec.push_back(HalfEdge{v_id, w});
                #pragma omp atomic
                m_total_half_edges++;
            }
        }

        void remove_half_edge(const partition_t u_id, const partition_t v_id, const weight_t w) {
            auto &vec = m_adj[u_id];
            auto it = find_edge(vec, v_id);
            ASSERT(it != vec.end());
            ASSERT(it->weight >= w);

            it->weight -= w;
            if (it->weight == 0) {
                std::swap(*it, vec.back());
                vec.pop_back();
                #pragma omp atomic
                m_total_half_edges--;
            }
        }

    public:
        LargeQuotientGraph() = default;

        void initialize(const partition_t t_k) {
            HEIPROMAP_PROFILE_SCOPE("recompute_datastructures", "qgraph", "initialize");
            m_k = t_k;
            m_adj.clear();
            m_adj.resize(m_k);
            m_self_weights.initialize(m_k, 0);
            m_total_half_edges = 0;
        }

        template<typename GraphT, typename PartitionManagerT>
        void compute_from_scratch(const GraphT &g, const PartitionManagerT &p_manager) {
            HEIPROMAP_PROFILE_SCOPE("recompute_datastructures", "qgraph", "compute_from_scratch");
            initialize(p_manager.get_k());

            for (vertex_t u = 0; u < g.n; ++u) {
                const partition_t u_id = p_manager[u];
                ASSERT(u_id < m_k);

                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    const vertex_t v = g.edges_v[i];
                    if (u < v) {
                        const partition_t v_id = p_manager[v];
                        ASSERT(v_id < m_k);
                        const weight_t w = g.edges_w[i];
                        add_edge(u_id, v_id, w);
                    }
                }
            }
        }

        template<typename GraphT, typename PartitionManagerT, typename BoundaryVertexManagerT>
        void compute_from_scratch(const GraphT &g,
                                  const PartitionManagerT &p_manager,
                                  const BoundaryVertexManagerT &bv_manager,
                                  const u64 num_threads) {
            HEIPROMAP_PROFILE_SCOPE("recompute_datastructures", "LargeQuotientGraph", "compute_from_scratch_bv");
            initialize(p_manager.get_k());

            if (num_threads <= 1 || bv_manager.size() < 1024) {
                for (partition_t id = 0; id < m_k; ++id) {
                    for (const vertex_t u: bv_manager.boundary(id)) {
                        for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                            const vertex_t v = g.edges_v[i];
                            if (u < v) {
                                const partition_t v_id = p_manager[v];
                                if (id != v_id) {
                                    add_edge(id, v_id, g.edges_w[i]);
                                }
                            }
                        }
                    }
                }
            } else {
                #pragma omp parallel for schedule(dynamic) num_threads(num_threads)
                for (partition_t id = 0; id < m_k; ++id) {
                    auto &adj = m_adj[id];
                    for (const vertex_t u: bv_manager.boundary(id)) {
                        for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                            const vertex_t v = g.edges_v[i];
                            const partition_t v_id = p_manager[v];
                            if (id != v_id) {
                                const weight_t w = g.edges_w[i];
                                auto it = find_edge(adj, v_id);
                                if (it != adj.end()) {
                                    it->weight += w;
                                } else {
                                    adj.push_back(HalfEdge{v_id, w});
                                }
                            }
                        }
                    }
                }

                u64 total = 0;
                for (partition_t id = 0; id < m_k; ++id) {
                    total += m_adj[id].size();
                }
                m_total_half_edges = total;
            }
        }

        void add_edge(const partition_t u_id, const partition_t v_id, const weight_t w) {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);
            if (w == 0) return;

            if (u_id == v_id) {
                m_self_weights[u_id] += w;
            } else {
                add_half_edge(u_id, v_id, w);
                add_half_edge(v_id, u_id, w);
            }
        }

        void remove_edge(const partition_t u_id, const partition_t v_id, const weight_t w) {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);
            if (w == 0) return;

            if (u_id == v_id) {
                ASSERT(m_self_weights[u_id] >= w);
                m_self_weights[u_id] -= w;
            } else {
                remove_half_edge(u_id, v_id, w);
                remove_half_edge(v_id, u_id, w);
            }
        }

        bool has_edge(const partition_t u_id, const partition_t v_id) const {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);
            if (u_id == v_id) {
                return m_self_weights[u_id] > 0;
            }
            const auto &vec = m_adj[u_id];
            auto it = find_edge(vec, v_id);
            return (it != vec.end() && it->weight > 0);
        }

        weight_t get_weight(const partition_t u_id, const partition_t v_id) const {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);
            if (u_id == v_id) {
                return m_self_weights[u_id];
            }
            const auto &vec = m_adj[u_id];
            auto it = find_edge(vec, v_id);
            if (it != vec.end()) {
                return it->weight;
            }
            return 0;
        }

        template<typename F>
        void for_each_neighbor(const partition_t x, F &&f) const {
            ASSERT(x < m_k);
            for (const auto &edge: m_adj[x]) {
                if (edge.weight > 0) {
                    f(edge.target, edge.weight);
                }
            }
        }

        const std::vector<HalfEdge> &neighbors(const partition_t x) const {
            ASSERT(x < m_k);
            return m_adj[x];
        }

        size_t degree(const partition_t x) const {
            ASSERT(x < m_k);
            return m_adj[x].size();
        }

        partition_t get_k() const {
            return m_k;
        }

        u64 total_half_edges() const {
            return m_total_half_edges;
        }

        template<typename GraphT, typename PartitionManagerT>
        void move(const GraphT &g,
                  const PartitionManagerT &p_manager,
                  const vertex_t u,
                  const partition_t old_id,
                  const partition_t new_id) {
            ASSERT(new_id < m_k);
            ASSERT(old_id < m_k);
            ASSERT(new_id != old_id);

            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t v = g.edges_v[i];
                const weight_t w = g.edges_w[i];
                const partition_t v_id = p_manager[v];

                if (old_id != v_id) {
                    remove_edge(old_id, v_id, w);
                }
                if (new_id != v_id) {
                    add_edge(new_id, v_id, w);
                }
            }
        }

        size_t edge_index(const partition_t u_id, const partition_t v_id) const {
            const partition_t min_part = std::min(u_id, v_id);
            const partition_t max_part = std::max(u_id, v_id);
            return static_cast<size_t>(min_part) * static_cast<size_t>(m_k) + static_cast<size_t>(max_part);
        }
    };
}

#endif //HEIPROMAP_LARGE_QUOTIENT_GRAPH_H
