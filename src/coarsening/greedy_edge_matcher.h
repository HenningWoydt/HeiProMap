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

#ifndef HEIPROMAP_GREEDY_EDGE_MATCHER_H
#define HEIPROMAP_GREEDY_EDGE_MATCHER_H

#include <algorithm>
#include <parallel/algorithm>
#include <execution>

#include "../definitions.h"
#include "../utility/random_engine.h"
#include "../utility/aligned_array.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/distance_oracle.h"
#include "../datastructures/partition_manager.h"

namespace HeiProMap {
    class GreedyEdgeMatcherConfiguration {
    public:
    };

    class GreedyEdgeMatcher {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
        u64 m_threads = 1;

        GreedyEdgeMatcherConfiguration config;
        RandomEngine random_engine = RandomEngine(0);

        AlignedArray<EdgeUVW> edges;
        size_t edges_size = 0;

        std::vector<std::vector<EdgeUVW> > thread_edges;

    public:
        GreedyEdgeMatcher() = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const u64 t_threads,
                        const u64 t_seed,
                        const GreedyEdgeMatcherConfiguration &t_config) {
            ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "initialize");

            m_n = t_n;
            m_m = t_m;
            m_k = t_k;
            m_threads = t_threads;

            config = t_config;
            random_engine = RandomEngine(t_seed);

            edges.initialize(m_m);

            thread_edges.resize(m_threads);
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void match([[maybe_unused]] const size_t level,
                   const graph_t &g,
                   const p_manager_t &p_manager,
                   Mapping &mapping,
                   f64 imbalance) {
            if (m_threads == 1) {
                match_serial<t_uniform_v_weights, t_uniform_e_weights>(level, g, p_manager, mapping, imbalance);
                return;
            }

            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            Matching matching;

            {
                ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "allocate");

                matching.initialize(g.n);
            }

            edges_size = 0;

            std::vector<f32> min_ratings(m_threads, std::numeric_limits<f32>::max());
            std::vector<f32> max_ratings(m_threads, std::numeric_limits<f32>::min());

            // handle all other vertices
            {
                ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "rate_edges");

                for (u64 t_id = 0; t_id < m_threads; ++t_id) {
                    thread_edges[t_id].clear();
                }

                #pragma omp parallel num_threads(m_threads)
                {
                    u64 t_id = omp_get_thread_num();

                    f32 min_rating = std::numeric_limits<f32>::max();
                    f32 max_rating = std::numeric_limits<f32>::min();

                    #pragma omp for
                    for (vertex_t u = 0; u < g.n; ++u) {
                        {
                            weight_t u_w = t_uniform_v_weights ? 1 : g.v_weights[u];

                            for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) { const vertex_t v = g.edges_v[j]; const weight_t w = g.edges_w[j];
                                {
                                    if (u >= v) { continue; }
                                    weight_t v_w = t_uniform_v_weights ? 1 : g.v_weights[v];
                                    if (u_w + v_w > lmax) { continue; }
                                    if (p_manager[u] != p_manager[v]) { continue; }

                                    weight_t ew = t_uniform_e_weights ? 1 : w;

                                    f32 edge_rating;
                                    if constexpr (t_uniform_v_weights && t_uniform_e_weights) {
                                        edge_rating = 1.0f;
                                    } else if constexpr (t_uniform_v_weights) {
                                        edge_rating = (f32) ew;
                                    } else {
                                        edge_rating = (f32) ew / std::sqrt((f32) u_w * (f32) v_w);
                                    }

                                    min_rating = std::min(min_rating, edge_rating);
                                    max_rating = std::max(max_rating, edge_rating);

                                    thread_edges[t_id].emplace_back(u, v, edge_rating);
                                }
                            }
                        }
                    }

                    min_ratings[t_id] = min_rating;
                    max_ratings[t_id] = max_rating;
                }
            }

            {
                ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "sort");

                #pragma omp parallel for num_threads(m_threads)
                for (u64 t_id = 0; t_id < m_threads; ++t_id) {
                    if (min_ratings[t_id] != max_ratings[t_id]) {
                        std::sort(thread_edges[t_id].begin(), thread_edges[t_id].end(), std::greater<>());
                    }
                }
            }

            f32 global_min_rating = std::numeric_limits<f32>::max();
            f32 global_max_rating = std::numeric_limits<f32>::min();

            {
                ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "global_min_max");

                for (u64 t_id = 0; t_id < m_threads; ++t_id) {
                    global_min_rating = std::min(global_min_rating, min_ratings[t_id]);
                    global_max_rating = std::max(global_max_rating, max_ratings[t_id]);
                }
            }

            {
                ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "choose");

                if (global_min_rating == global_max_rating) {
                    std::vector<size_t> pos(m_threads, 0);
                    u64 remaining_threads = 0;

                    for (u64 t_id = 0; t_id < m_threads; ++t_id) {
                        if (!thread_edges[t_id].empty()) {
                            ++remaining_threads;
                        }
                    }

                    while (remaining_threads > 0) {
                        for (u64 t_id = 0; t_id < m_threads; ++t_id) {
                            if (pos[t_id] >= thread_edges[t_id].size()) {
                                continue;
                            }

                            const auto &[u, v, w] = thread_edges[t_id][pos[t_id]];
                            ++pos[t_id];

                            if (pos[t_id] == thread_edges[t_id].size()) {
                                --remaining_threads;
                            }

                            if (!matching.is_matched(u) && !matching.is_matched(v)) {
                                matching.add(u, v);
                            }
                        }
                    }
                } else {
                    struct HeapEntry {
                        vertex_t u;
                        vertex_t v;
                        f32 w;
                        u64 tid;
                        size_t idx;

                        bool operator<(const HeapEntry &other) const {
                            return w < other.w; // max-heap
                        }
                    };

                    std::priority_queue<HeapEntry> pq;

                    for (u64 t_id = 0; t_id < m_threads; ++t_id) {
                        if (!thread_edges[t_id].empty()) {
                            const auto &[u, v, w] = thread_edges[t_id][0];
                            pq.push({u, v, w, t_id, 0});
                        }
                    }

                    while (!pq.empty()) {
                        auto top = pq.top();
                        pq.pop();

                        if (!matching.is_matched(top.u) && !matching.is_matched(top.v)) {
                            matching.add(top.u, top.v);
                        }

                        size_t next_idx = top.idx + 1;
                        if (next_idx < thread_edges[top.tid].size()) {
                            const auto &[u, v, w] = thread_edges[top.tid][next_idx];
                            pq.push({u, v, w, top.tid, next_idx});
                        }
                    }
                }
            }

            {
                ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "mapping");

                matching.set_translation();
                mapping.set_coarse_n(matching.get_n_coarse_nodes());
                for (vertex_t u = 0; u < matching.get_n(); ++u) {
                    mapping.set(u, matching.get_n(u));
                }
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void match_serial([[maybe_unused]] const size_t level,
                          const graph_t &g,
                          const p_manager_t &p_manager,
                          Mapping &mapping,
                          f64 imbalance) {
            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            Matching matching;

            {
                ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "allocate");

                matching.initialize(g.n);
            }

            if constexpr (t_uniform_v_weights && t_uniform_e_weights) {
                // Fast path: all edges equally rated, skip rating/sorting.
                // For each unmatched vertex, pick a random eligible neighbor.
                {
                    ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "choose");

                    for (vertex_t u = 0; u < g.n; ++u) {
                        if (matching.is_matched(u)) { continue; }

                        for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) { const vertex_t v = g.edges_v[j];
                            {
                                if (matching.is_matched(v)) { continue; }
                                if (p_manager[u] != p_manager[v]) { continue; }
                                matching.add(u, v);
                                break;
                            }
                        }
                    }
                }
            } else {
                edges_size = 0;

                f32 min_rating = std::numeric_limits<f32>::max();
                f32 max_rating = std::numeric_limits<f32>::min();

                {
                    ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "rate_edges");

                    for (vertex_t u = 0; u < g.n; ++u) {
                        {
                            weight_t u_w = t_uniform_v_weights ? 1 : g.v_weights[u];

                            for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) { const vertex_t v = g.edges_v[j]; const weight_t w = g.edges_w[j];
                                {
                                    if (u >= v) { continue; }
                                    weight_t v_w = t_uniform_v_weights ? 1 : g.v_weights[v];
                                    if (u_w + v_w > lmax) { continue; }
                                    if (p_manager[u] != p_manager[v]) { continue; }

                                    weight_t ew = t_uniform_e_weights ? 1 : w;

                                    f32 edge_rating;
                                    if constexpr (t_uniform_v_weights) {
                                        edge_rating = (f32) ew;
                                    } else {
                                        edge_rating = (f32) ew / std::sqrt((f32) u_w * (f32) v_w);
                                    }

                                    min_rating = std::min(min_rating, edge_rating);
                                    max_rating = std::max(max_rating, edge_rating);

                                    edges[edges_size++] = {u, v, edge_rating};
                                }
                            }
                        }
                    }
                }

                if (min_rating != max_rating) {
                    ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "sort");

                    std::sort(edges.get_ptr(), edges.get_ptr() + edges_size, std::greater<>());
                }

                {
                    ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "choose");

                    for (size_t i = 0; i < edges_size; ++i) {
                        const auto &[u, v, w] = edges[i];
                        if (!matching.is_matched(u) && !matching.is_matched(v)) {
                            matching.add(u, v);
                        }
                    }
                }
            }

            {
                ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "mapping");

                matching.set_translation();
                mapping.set_coarse_n(matching.get_n_coarse_nodes());
                for (vertex_t u = 0; u < matching.get_n(); ++u) {
                    mapping.set(u, matching.get_n(u));
                }
            }
        }
    };
}

#endif //HEIPROMAP_GREEDY_EDGE_MATCHER_H
