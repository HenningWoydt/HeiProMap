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

#ifndef HEIPROMAP_QUOTIENT_GRAPH_H
#define HEIPROMAP_QUOTIENT_GRAPH_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../definitions.h"
#include "../utility/aligned_array.h"
#include "../utility/profiler.h"

namespace HeiProMap {
    class QuotientGraph {
        static constexpr size_t INVALID_INDEX = std::numeric_limits<size_t>::max();

        struct EdgeRecord {
            weight_t weight = 0;

            // canonical endpoints: a <= b
            partition_t a = 0;
            partition_t b = 0;

            // intrusive links in adjacency list of endpoint a
            size_t prev_a = INVALID_INDEX;
            size_t next_a = INVALID_INDEX;

            // intrusive links in adjacency list of endpoint b
            size_t prev_b = INVALID_INDEX;
            size_t next_b = INVALID_INDEX;
        };

        partition_t m_k = 0;

        // Canonical edge storage addressed via edge_index(min(u,v), max(u,v)).
        // We keep the same k*k indexing scheme as before for compatibility.
        AlignedArray<EdgeRecord> m_edges;

        // Head of adjacency list for each partition. Only non-self-loop edges
        // with positive weight are linked into these lists.
        AlignedArray<size_t> m_head;

        static partition_t min_id(const partition_t u, const partition_t v) {
            return std::min(u, v);
        }

        static partition_t max_id(const partition_t u, const partition_t v) {
            return std::max(u, v);
        }

        const EdgeRecord &edge_ref(const partition_t u_id, const partition_t v_id) const {
            return m_edges[edge_index(u_id, v_id)];
        }

        static partition_t other_endpoint(const EdgeRecord &e, const partition_t x) {
            return (e.a == x) ? e.b : e.a;
        }

        static size_t &prev_for(EdgeRecord &e, const partition_t x) {
            return (e.a == x) ? e.prev_a : e.prev_b;
        }

        static size_t &next_for(EdgeRecord &e, const partition_t x) {
            return (e.a == x) ? e.next_a : e.next_b;
        }

        static size_t prev_for(const EdgeRecord &e, const partition_t x) {
            return (e.a == x) ? e.prev_a : e.prev_b;
        }

        static size_t next_for(const EdgeRecord &e, const partition_t x) {
            return (e.a == x) ? e.next_a : e.next_b;
        }

        static bool is_self_loop(const EdgeRecord &e) {
            return e.a == e.b;
        }

        void link_edge(const size_t idx) {
            EdgeRecord &e = m_edges[idx];

            if (is_self_loop(e) || e.weight <= 0) {
                return;
            }

            // Insert into adjacency list of a
            e.prev_a = INVALID_INDEX;
            e.next_a = m_head[e.a];
            if (m_head[e.a] != INVALID_INDEX) {
                EdgeRecord &old_head = m_edges[m_head[e.a]];
                prev_for(old_head, e.a) = idx;
            }
            m_head[e.a] = idx;

            // Insert into adjacency list of b
            e.prev_b = INVALID_INDEX;
            e.next_b = m_head[e.b];
            if (m_head[e.b] != INVALID_INDEX) {
                EdgeRecord &old_head = m_edges[m_head[e.b]];
                prev_for(old_head, e.b) = idx;
            }
            m_head[e.b] = idx;
        }

        void unlink_one_endpoint(const size_t idx, const partition_t x) {
            EdgeRecord &e = m_edges[idx];

            size_t &prev = prev_for(e, x);
            size_t &next = next_for(e, x);

            if (prev != INVALID_INDEX) {
                EdgeRecord &p = m_edges[prev];
                next_for(p, x) = next;
            } else {
                m_head[x] = next;
            }

            if (next != INVALID_INDEX) {
                EdgeRecord &n = m_edges[next];
                prev_for(n, x) = prev;
            }

            prev = INVALID_INDEX;
            next = INVALID_INDEX;
        }

        void unlink_edge(const size_t idx) {
            EdgeRecord &e = m_edges[idx];

            if (is_self_loop(e)) {
                return;
            }

            unlink_one_endpoint(idx, e.a);
            unlink_one_endpoint(idx, e.b);
        }

    public:
        void initialize(const partition_t t_k) {
            HEIPROMAP_PROFILE_SCOPE("misc", "QuotientGraph", "initialize");

            m_k = t_k;

            const size_t size = static_cast<size_t>(m_k) * static_cast<size_t>(m_k);
            m_edges.initialize(size);
            m_head.initialize(m_k, INVALID_INDEX);

            for (partition_t u = 0; u < m_k; ++u) {
                for (partition_t v = 0; v < m_k; ++v) {
                    const size_t idx = static_cast<size_t>(u) * static_cast<size_t>(m_k) + static_cast<size_t>(v);
                    m_edges[idx].weight = 0;
                    m_edges[idx].a = min_id(u, v);
                    m_edges[idx].b = max_id(u, v);
                    m_edges[idx].prev_a = INVALID_INDEX;
                    m_edges[idx].next_a = INVALID_INDEX;
                    m_edges[idx].prev_b = INVALID_INDEX;
                    m_edges[idx].next_b = INVALID_INDEX;
                }
            }
        }

        template<typename F>
        void for_each_neighbor(const partition_t x, F &&f) const {
            size_t idx = m_head[x];
            while (idx != INVALID_INDEX) {
                const EdgeRecord &e = m_edges[idx];
                const partition_t y = other_endpoint(e, x);
                f(y, e.weight);
                idx = next_for(e, x);
            }
        }

        void add_edge(const partition_t u_id, const partition_t v_id, const weight_t w) {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);

            const size_t idx = edge_index(u_id, v_id);
            EdgeRecord &e = m_edges[idx];

            const bool was_inactive = (e.weight <= 0);
            e.weight += w;

            if (u_id != v_id && was_inactive && e.weight > 0) {
                link_edge(idx);
            }
        }

        void remove_edge(const partition_t u_id, const partition_t v_id, const weight_t w) {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);

            const size_t idx = edge_index(u_id, v_id);
            EdgeRecord &e = m_edges[idx];

            ASSERT(e.weight >= w);

            e.weight -= w;

            if (u_id != v_id && e.weight == 0) {
                unlink_edge(idx);
            }
        }

        bool has_edge(const partition_t u_id, const partition_t v_id) const {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);
            return edge_ref(u_id, v_id).weight > 0;
        }

        weight_t get_weight(const partition_t u_id, const partition_t v_id) const {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);
            return edge_ref(u_id, v_id).weight;
        }

        /**
         * O(max_deg)
         *
         * @param g
         * @param p_manager
         * @param u
         * @param old_id
         * @param new_id
         */
        template<typename GraphT, typename PartitionManagerT>
        void move(const GraphT &g,
                  const PartitionManagerT &p_manager,
                  const vertex_t u,
                  const partition_t old_id,
                  const partition_t new_id) {
            ASSERT(new_id < m_k);
            ASSERT(old_id < m_k);
            ASSERT(new_id != old_id);

            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) { const vertex_t v = g.edges_v[i]; const weight_t w = g.edges_w[i];
                {
                    const partition_t v_id = p_manager[v];

                    // remove old edge, if existed
                    if (old_id != v_id) {
                        remove_edge(old_id, v_id, w);
                    }

                    // add new edge, if has to exist
                    if (new_id != v_id) {
                        add_edge(new_id, v_id, w);
                    }
                }
            }
        }


        bool find_distance_3_matching(AlignedArray<u8> &active_this_round,
                                      AlignedArray<u8> &used_edges_this_round,
                                      std::vector<std::pair<partition_t, partition_t> > &matching) {
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
                    if (used_edges_this_round[eidx] == 1) {
                        return;
                    }

                    matching.emplace_back(u_id, v_id);
                    used_edges_this_round[eidx] = 1;

                    freeze_distance_2(u_id);
                    freeze_distance_2(v_id);
                    matched_u = true;
                });
            }

            return !matching.empty();
        }

        bool find_all_pairs(AlignedArray<u8> &active_this_round,
                            AlignedArray<u8> &used_edges_this_round,
                            std::vector<std::pair<partition_t, partition_t> > &matching) {
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
                    if (used_edges_this_round[eidx] == 1) {
                        return;
                    }

                    matching.emplace_back(u_id, v_id);
                    used_edges_this_round[eidx] = 1;
                });
            }

            return !matching.empty();
        }

        size_t edge_index(const partition_t u_id, const partition_t v_id) const {
            const partition_t min_part = min_id(u_id, v_id);
            const partition_t max_part = max_id(u_id, v_id);
            return static_cast<size_t>(min_part) * static_cast<size_t>(m_k) + static_cast<size_t>(max_part);
        }
    };
}

#endif //HEIPROMAP_QUOTIENT_GRAPH_H
