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

#ifndef HEIPROMAP_DEEP_QUOTIENT_GRAPH_H
#define HEIPROMAP_DEEP_QUOTIENT_GRAPH_H


#include "../../../commons/definitions.h"

namespace HeiProMap {
    class DeepQuotientGraph {
        partition_t m_k = 0;

        std::vector<std::vector<std::pair<partition_t, weight_t>>> edges;

        std::vector<std::pair<partition_t, partition_t>> pairs;

    public:
        void initialize(const partition_t t_k) {
            m_k = t_k;

            edges.resize(m_k);
        }

        void add_edge(const partition_t u_id, const partition_t v_id, const weight_t w) {
            ASSERT(u_id != v_id);

            add(u_id, v_id, w);
            add(v_id, u_id, w);
        }

        void remove_edge(const partition_t u_id, const partition_t v_id, const weight_t w) {
            ASSERT(u_id != v_id);

            remove(u_id, v_id, w);
            remove(v_id, u_id, w);
        }

        bool has_edge(const partition_t u_id, const partition_t v_id) const {
            ASSERT(u_id != v_id);

            for (auto &[id, w]: edges[u_id]) {
                if (id == v_id) {
                    return true;
                }
            }

            return false;
        }

        std::vector<std::pair<partition_t, weight_t>> &neighborhood(const partition_t id){
            return edges[id];
        }

        weight_t get_weight(const partition_t u_id, const partition_t v_id) const {
            ASSERT(u_id != v_id);

            for (const auto &[id, w]: edges[u_id]) {
                if (id == v_id) {
                    return w;
                }
            }

            return 0;
        }

        void move(const graph_t &g,
                  const deep_p_manager_t &p_manager,
                  const vertex_t u,
                  const partition_t old_id,
                  const partition_t new_id) {
            ASSERT(new_id < m_k);
            ASSERT(new_id != old_id);

            forall_guivw(g, u, i, v, w)
                {
                    partition_t v_id = p_manager[v];

                    // remove old edge, if existed
                    if (old_id != v_id) {
                        remove_edge(old_id, v_id, w);
                    }
                    // add new edge, if has to exist
                    if (new_id != v_id) {
                        add_edge(new_id, v_id, w);
                    }
                }
            endfor
        }

        void reset() {
            for (auto &edge: edges) {
                edge.clear();
            }
        }

        size_t n_pairs(AlignedArray<u8> &active) {
            size_t n = 0;
            pairs.clear();

            for (partition_t id1 = 0; id1 < m_k; ++id1) {
                for (auto &[id2, w]: edges[id1]) {
                    if (active[id1] == 0 && active[id2] == 0) { continue; }
                    if (w == 0) { continue; }
                    n += 1;
                    pairs.emplace_back(id1, id2);
                }
            }
            return n;
        }

        std::pair<partition_t, partition_t> get_pair(size_t idx) { return pairs[idx]; }

        std::vector<std::vector<std::pair<partition_t, partition_t>>> get_distance_2_matchings(AlignedArray<u8> &active) {
            std::vector<std::vector<std::pair<partition_t, partition_t>>> matchings;

            size_t m = n_pairs(active);

            std::vector<u8> edge_included(m, 0);
            std::vector<u8> vertex_frozen(m_k, 0);

            while (true) {
                matchings.emplace_back();

                for (size_t i = 0; i < m; ++i) {
                    if (edge_included[i] == 1) { continue; }
                    auto [u_id, v_id] = get_pair(i);

                    u64 neighbors_frozen = 0;
                    for (auto &[id, w]: edges[u_id]) {
                        neighbors_frozen += vertex_frozen[id];
                    }
                    for (auto &[id, w]: edges[v_id]) {
                        neighbors_frozen += vertex_frozen[id];
                    }

                    if (neighbors_frozen > 0) { continue; }

                    matchings.back().emplace_back(u_id, v_id);
                    edge_included[i] = 1;

                    for (auto &[id, w]: edges[u_id]) {
                        vertex_frozen[id] = 1;
                    }
                    for (auto &[id, w]: edges[v_id]) {
                        vertex_frozen[id] = 1;
                    }
                }

                if (matchings.back().empty()) {
                    matchings.pop_back();
                    break;
                }

                std::fill_n(vertex_frozen.begin(), m_k, 0);
            }

            return matchings;
        }

    private:
        void add(const partition_t id1, const partition_t id2, const weight_t w) {
            for (auto &[id, id_w]: edges[id1]) {
                if (id == id2) {
                    id_w += w;
                    return;
                }
            }
            edges[id1].emplace_back(id2, w);
        }

        void remove(const partition_t id1, const partition_t id2, const weight_t w) {
            for (size_t i = 0; i < edges[id1].size(); ++i) {
                if (edges[id1][i].first == id2) {
                    edges[id1][i].second -= w;
                    if (edges[id1][i].second == 0) {
                        std::swap(edges[id1][i], edges[id1].back());
                        edges[id1].pop_back();
                    }
                    return;
                }
            }
        }
    };
}

#endif //HEIPROMAP_DEEP_QUOTIENT_GRAPH_H
