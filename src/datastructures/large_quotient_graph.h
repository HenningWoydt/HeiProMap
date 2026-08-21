#ifndef HEIPROMAP_LARGE_QUOTIENT_GRAPH_H
#define HEIPROMAP_LARGE_QUOTIENT_GRAPH_H

#include <vector>
#include <algorithm>
#include "../definitions.h"
#include "../utility/aligned_array.h"

namespace HeiProMap {
    class LargeQuotientGraph {
    private:
        partition_t m_k = 0;

        // Dynamic graph: vector of vectors
        std::vector<std::vector<std::pair<partition_t, weight_t>>> m_dynamic_adj;

        // CSR graph
        mutable bool m_dirty = false;
        mutable std::vector<size_t> m_neighborhoods;
        mutable std::vector<partition_t> m_edges_v;
        mutable std::vector<weight_t> m_edges_w;

        void ensure_csr() const {
            if (!m_dirty) return;

            m_neighborhoods.assign(m_k + 1, 0);
            size_t num_edges = 0;
            for (partition_t u = 0; u < m_k; ++u) {
                num_edges += m_dynamic_adj[u].size();
                m_neighborhoods[u + 1] = num_edges;
            }

            m_edges_v.resize(num_edges);
            m_edges_w.resize(num_edges);

            size_t idx = 0;
            for (partition_t u = 0; u < m_k; ++u) {
                for (const auto& edge : m_dynamic_adj[u]) {
                    m_edges_v[idx] = edge.first;
                    m_edges_w[idx] = edge.second;
                    idx++;
                }
            }
            m_dirty = false;
        }

    public:
        void initialize(const partition_t t_k) {
            HEIPROMAP_PROFILE_SCOPE("misc", "LargeQuotientGraph", "initialize");
            m_k = t_k;
            m_dynamic_adj.assign(m_k, std::vector<std::pair<partition_t, weight_t>>());
            
            // Clean up CSR
            m_neighborhoods.clear();
            m_edges_v.clear();
            m_edges_w.clear();
            m_dirty = true;
        }

        template<typename F>
        void for_each_neighbor(const partition_t x, F &&f) const {
            ensure_csr();
            for (size_t i = m_neighborhoods[x]; i < m_neighborhoods[x + 1]; ++i) {
                f(m_edges_v[i], m_edges_w[i]);
            }
        }

        void add_edge(const partition_t u_id, const partition_t v_id, const weight_t w) {
            if (u_id == v_id) return;
            m_dirty = true;

            auto add_half = [&](partition_t from, partition_t to) {
                for (auto& edge : m_dynamic_adj[from]) {
                    if (edge.first == to) {
                        edge.second += w;
                        return;
                    }
                }
                m_dynamic_adj[from].push_back({to, w});
            };

            add_half(u_id, v_id);
            add_half(v_id, u_id);
        }

        void remove_edge(const partition_t u_id, const partition_t v_id, const weight_t w) {
            if (u_id == v_id) return;
            m_dirty = true;

            auto remove_half = [&](partition_t from, partition_t to) {
                auto& adj = m_dynamic_adj[from];
                for (size_t i = 0; i < adj.size(); ++i) {
                    if (adj[i].first == to) {
                        ASSERT(adj[i].second >= w);
                        adj[i].second -= w;
                        if (adj[i].second <= 0) {
                            adj[i] = adj.back();
                            adj.pop_back();
                        }
                        return;
                    }
                }
            };

            remove_half(u_id, v_id);
            remove_half(v_id, u_id);
        }

        bool has_edge(const partition_t u_id, const partition_t v_id) const {
            ensure_csr();
            for (size_t i = m_neighborhoods[u_id]; i < m_neighborhoods[u_id + 1]; ++i) {
                if (m_edges_v[i] == v_id) return true;
            }
            return false;
        }

        weight_t get_weight(const partition_t u_id, const partition_t v_id) const {
            ensure_csr();
            for (size_t i = m_neighborhoods[u_id]; i < m_neighborhoods[u_id + 1]; ++i) {
                if (m_edges_v[i] == v_id) return m_edges_w[i];
            }
            return 0;
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
                                      std::vector<std::pair<partition_t, partition_t> > &matching) {
            ensure_csr();
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
            ensure_csr();
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
            const partition_t min_part = std::min(u_id, v_id);
            const partition_t max_part = std::max(u_id, v_id);
            return static_cast<size_t>(min_part) * static_cast<size_t>(m_k) + static_cast<size_t>(max_part);
        }
    };
}

#endif //HEIPROMAP_LARGE_QUOTIENT_GRAPH_H
