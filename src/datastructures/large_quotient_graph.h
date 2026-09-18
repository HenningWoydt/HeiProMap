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
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

#include "../definitions.h"
#include "../utility/aligned_array.h"
#include "../utility/macros.h"
#include "../utility/profiler.h"

namespace HeiProMap {
    /**
     * LargeQuotientGraph represents the quotient graph using sparse adjacency lists
     * per block with sorted neighbor targets, achieving O(k + |E_Q|) memory instead
     * of O(k^2) dense matrix storage. Suitable for very large numbers of blocks k.
     */
    class LargeQuotientGraph {
    public:
        struct HalfEdge {
            partition_t target = 0;
            weight_t weight = 0;

            bool operator<(const HalfEdge &other) const {
                return target < other.target;
            }
        };

    private:
        partition_t m_k = 0;
        std::vector<std::vector<HalfEdge>> m_adj;
        std::vector<weight_t> m_self_weights;

        static auto find_edge(std::vector<HalfEdge> &vec, const partition_t target) {
            return std::lower_bound(
                vec.begin(), vec.end(), target,
                [](const HalfEdge &e, const partition_t t) { return e.target < t; }
            );
        }

        static auto find_edge(const std::vector<HalfEdge> &vec, const partition_t target) {
            return std::lower_bound(
                vec.begin(), vec.end(), target,
                [](const HalfEdge &e, const partition_t t) { return e.target < t; }
            );
        }

        void add_half_edge(const partition_t u_id, const partition_t v_id, const weight_t w) {
            auto &vec = m_adj[u_id];
            auto it = find_edge(vec, v_id);
            if (it != vec.end() && it->target == v_id) {
                it->weight += w;
            } else {
                vec.insert(it, HalfEdge{v_id, w});
            }
        }

        void remove_half_edge(const partition_t u_id, const partition_t v_id, const weight_t w) {
            auto &vec = m_adj[u_id];
            auto it = find_edge(vec, v_id);
            ASSERT(it != vec.end() && it->target == v_id);
            ASSERT(it->weight >= w);

            it->weight -= w;
            if (it->weight == 0) {
                vec.erase(it);
            }
        }

    public:
        LargeQuotientGraph() = default;

        void initialize(const partition_t t_k) {
            HEIPROMAP_PROFILE_SCOPE("misc", "LargeQuotientGraph", "initialize");
            m_k = t_k;
            m_adj.clear();
            m_adj.resize(m_k);
            m_self_weights.assign(m_k, 0);
        }

        template<typename GraphT, typename PartitionManagerT>
        void compute_from_scratch(const GraphT &g, const PartitionManagerT &p_manager) {
            HEIPROMAP_PROFILE_SCOPE("misc", "LargeQuotientGraph", "compute_from_scratch");
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
            return (it != vec.end() && it->target == v_id && it->weight > 0);
        }

        weight_t get_weight(const partition_t u_id, const partition_t v_id) const {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);
            if (u_id == v_id) {
                return m_self_weights[u_id];
            }
            const auto &vec = m_adj[u_id];
            auto it = find_edge(vec, v_id);
            if (it != vec.end() && it->target == v_id) {
                return it->weight;
            }
            return 0;
        }

        template<typename F>
        void for_each_neighbor(const partition_t x, F &&f) const {
            ASSERT(x < m_k);
            for (const auto &edge : m_adj[x]) {
                if (edge.weight > 0) {
                    f(edge.target, edge.weight);
                }
            }
        }

        size_t degree(const partition_t x) const {
            ASSERT(x < m_k);
            return m_adj[x].size();
        }

        partition_t get_k() const {
            return m_k;
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

        bool find_distance_3_matching(AlignedArray<u8> &active_this_round,
                                      AlignedArray<u8> &used_edges_this_round,
                                      std::vector<std::pair<partition_t, partition_t>> &matching) {
            matching.clear();

            std::vector<u8> vertex_frozen(m_k, 0);

            auto freeze_distance_2 = [&](const partition_t x) {
                vertex_frozen[x] = 1;

                for_each_neighbor(x, [&](const partition_t n1, const weight_t) {
                    vertex_frozen[n1] = 1;

                    for_each_neighbor(n1, [&](const partition_t n2, const weight_t) {
                        vertex_frozen[n2] = 1;
                    });
                });
            };

            for (partition_t u_id = 0; u_id < m_k; ++u_id) {
                if (vertex_frozen[u_id] == 1) {
                    continue;
                }

                bool matched_u = false;

                for_each_neighbor(u_id, [&](const partition_t v_id, const weight_t) {
                    if (matched_u) {
                        return;
                    }

                    if (v_id <= u_id) {
                        return;
                    }
                    if (vertex_frozen[v_id] == 1) {
                        return;
                    }
                    if (active_this_round[u_id] == 0 && active_this_round[v_id] == 0) {
                        return;
                    }

                    const size_t eidx = edge_index(u_id, v_id);
                    if (eidx < used_edges_this_round.size() && used_edges_this_round[eidx] == 1) {
                        return;
                    }

                    matching.emplace_back(u_id, v_id);
                    if (eidx < used_edges_this_round.size()) {
                        used_edges_this_round[eidx] = 1;
                    }

                    freeze_distance_2(u_id);
                    freeze_distance_2(v_id);
                    matched_u = true;
                });
            }

            return !matching.empty();
        }

        bool find_all_pairs(AlignedArray<u8> &active_this_round,
                            AlignedArray<u8> &used_edges_this_round,
                            std::vector<std::pair<partition_t, partition_t>> &matching) {
            matching.clear();

            for (partition_t u_id = 0; u_id < m_k; ++u_id) {
                for_each_neighbor(u_id, [&](const partition_t v_id, const weight_t) {
                    if (v_id <= u_id) {
                        return;
                    }
                    if (active_this_round[u_id] == 0 && active_this_round[v_id] == 0) {
                        return;
                    }

                    const size_t eidx = edge_index(u_id, v_id);
                    if (eidx < used_edges_this_round.size() && used_edges_this_round[eidx] == 1) {
                        return;
                    }

                    matching.emplace_back(u_id, v_id);
                    if (eidx < used_edges_this_round.size()) {
                        used_edges_this_round[eidx] = 1;
                    }
                });
            }

            return !matching.empty();
        }

        size_t edge_index(const partition_t u_id, const partition_t v_id) const {
            const partition_t min_part = std::min(u_id, v_id);
            const partition_t max_part = std::max(u_id, v_id);
            return static_cast<size_t>(min_part) * static_cast<size_t>(m_k) + static_cast<size_t>(max_part);
        }
    };
}

#endif //HEIPROMAP_LARGE_QUOTIENT_GRAPH_H
