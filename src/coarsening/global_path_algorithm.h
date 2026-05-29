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

#ifndef HEIPROMAP_GLOBAL_PATH_ALGORITHM_H
#define HEIPROMAP_GLOBAL_PATH_ALGORITHM_H

#include <algorithm>
#include <chrono>
#include <queue>
#include <vector>
#include <cmath>
#include <omp.h>

#include "../definitions.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../utility/aligned_array.h"
#include "../utility/JSON_utils.h"
#include "../utility/mapping.h"
#include "../utility/matching.h"
#include "../utility/profiler.h"
#include "../utility/random_engine.h"

namespace HeiProMap {
    struct Neighbors {
        vertex_t n1;
        vertex_t n2;
        f32 w1;
        f32 w2;
    };

    static inline bool is_not_endpoint(const Neighbors &n, vertex_t u) { return n.n2 != u; }

    static inline bool is_one_endpoint(const Neighbors &n, vertex_t u) { return n.n1 != u && n.n2 == u; }

    static inline bool is_unmatched(const Neighbors &n, vertex_t u) { return n.n1 == u; }

    class GlobalPathAlgorithmConfiguration {
    public:
        size_t random_level = 0;
        EdgeRatingFunction rating_function = EdgeRatingFunction::HEAVY_EDGE;
        bool use_adaptive_max_vertex_weight = false;
        bool use_edge_rating_tiebreaking = true;
        f64 two_hop_threshold = 0.75;
    };

    /**
     * Computes a matching based on the Global Path Algorithm from
     * > Jens Maue and Peter Sanders.
     * > Engineering Algorithms for Approximate Weighted Matching.
     * > Experimental Algorithms, 6th International Workshop, WEA 2007, Rome, Italy, June 6-8, 2007, Proceedings.
     */
    class GlobalPathAlgorithmMatcher {
        struct HeapEntry {
            size_t thread_idx;
            size_t edge_idx;
            EdgeUVW edge;

            bool operator<(const HeapEntry &other) const {
                return edge < other.edge;
            }
        };

        struct ThreadInfo {
            struct DPState {
                f32 w;
                s64 m;
                u8 take;
                vertex_t edge;
            };

            std::vector<EdgeUVW> local_edges;
            size_t edge_idx = 0;
            f32 min_rating = std::numeric_limits<f32>::max();
            f32 max_rating = std::numeric_limits<f32>::min();

            // for DP
            AlignedArray<DPState> dp;

            std::vector<std::pair<vertex_t, vertex_t> > dp_cycle_matches1;
            std::vector<std::pair<vertex_t, vertex_t> > dp_cycle_matches2;
        };

        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
        u64 m_threads = 1;

        const GlobalPathAlgorithmConfiguration *config = nullptr;
        RandomEngine *random_engine = nullptr;

        AlignedArray<Neighbors> m_neighbors;
        AlignedArray<u32> path_id;
        AlignedArray<u32> path_length;

        std::vector<ThreadInfo> m_thread_infos;
        std::vector<vertex_t> cycles;

    public:
        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const u64 t_threads,
                        RandomEngine &t_random_engine,
                        const GlobalPathAlgorithmConfiguration &i_config) {
            HEIPROMAP_PROFILE_SCOPE("coarsening", "GlobalPathAlgorithmMatcher", "initialize");

            m_n = t_n;
            m_m = t_m;
            m_k = t_k;
            m_threads = t_threads;

            config = &i_config;
            random_engine = &t_random_engine;

            m_neighbors.initialize(m_n);
            path_id.initialize(m_n);
            path_length.initialize(m_n);

            m_thread_infos.resize(m_threads);
            size_t edges_per_thread = (m_m / m_threads) + 1;
            for (u64 i = 0; i < m_threads; ++i) {
                m_thread_infos[i].local_edges.reserve(edges_per_thread);
                m_thread_infos[i].dp.initialize(m_n);
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights, EdgeRatingFunction t_rating_function>
        void match_templated(const size_t level,
                             const graph_t &g,
                             const p_manager_t &p_manager,
                             Mapping &mapping,
                             f64 imbalance) {
            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            if (config->use_adaptive_max_vertex_weight) {
                // KaHIP's fast mode logic: 1.5 * W_total / num_stop
                // num_stop = max(N / (2 * 60 * k), 60 * k)
                f64 x = 60.0;
                f64 num_stop = std::max((f64) g.n / (2.0 * x * (f64) p_manager.k), x * (f64) p_manager.k);
                weight_t adaptive_lmax = (weight_t) (1.5 * (f64) g.g_weight / num_stop);
                lmax = std::min(lmax, std::max((weight_t) 2, adaptive_lmax));
            }

            Matching matching;
            HEIPROMAP_PROFILE_SCOPE("coarsening", "GlobalPathAlgorithmMatcher", "allocate_matching");
            matching.initialize(g.n);

            if (m_threads == 1) {
                if constexpr (t_uniform_v_weights && t_uniform_e_weights) {
                    HEIPROMAP_PROFILE_SCOPE("coarsening", "GlobalPathAlgorithmMatcher", "simple_loop");
                    for (vertex_t u = 0; u < g.n; ++u) {
                        if (matching.is_matched(u)) { continue; }
                        partition_t u_id = p_manager[u];
                        for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                            const vertex_t v = g.edges_v[j];
                            if (matching.is_matched(v)) { continue; }
                            if (u_id != p_manager[v]) { continue; }
                            matching.add(u, v);
                            break;
                        }
                    }
                    finalize_matching(g, matching, mapping);
                    return;
                }
            }

            if (level < config->random_level) {
                random_matching<t_uniform_v_weights>(level, g, matching, lmax);
                finalize_matching(g, matching, mapping);
                return;
            }

            compute_ratings<t_uniform_v_weights, t_uniform_e_weights, t_rating_function>(g, p_manager, lmax);

            f32 global_min_rating = std::numeric_limits<f32>::max();
            f32 global_max_rating = std::numeric_limits<f32>::min();
            for (u64 i = 0; i < m_threads; ++i) {
                global_min_rating = std::min(global_min_rating, m_thread_infos[i].min_rating);
                global_max_rating = std::max(global_max_rating, m_thread_infos[i].max_rating);
            }

            HEIPROMAP_PROFILE_SCOPE("coarsening", "GlobalPathAlgorithmMatcher", "sort_ratings");
            std::vector<u64> seeds(m_threads);
            for (u64 i = 0; i < m_threads; ++i) { seeds[i] = random_engine->get_u64(); }

            #pragma omp parallel for num_threads(m_threads) schedule(static, 1)
            for (u64 i = 0; i < m_threads; ++i) {
                if (config->use_edge_rating_tiebreaking) {
                    std::mt19937 g(seeds[i]);
                    std::shuffle(m_thread_infos[i].local_edges.begin(), m_thread_infos[i].local_edges.end(), g);
                }
                if (m_thread_infos[i].min_rating != m_thread_infos[i].max_rating) {
                    std::sort(m_thread_infos[i].local_edges.begin(), m_thread_infos[i].local_edges.end(), std::greater<>());
                }
                m_thread_infos[i].edge_idx = 0;
            }

            HEIPROMAP_PROFILE_SCOPE("coarsening", "GlobalPathAlgorithmMatcher", "init_paths");
            #pragma omp parallel for num_threads(m_threads)
            for (vertex_t u = 0; u < g.n; ++u) {
                m_neighbors[u].n1 = u;
                m_neighbors[u].n2 = u;
                path_id[u] = u;
                path_length[u] = 0;
            }

            cycles.clear();
            HEIPROMAP_PROFILE_SCOPE("coarsening", "GlobalPathAlgorithmMatcher", "extract_paths");

            auto process_edge = [&](vertex_t u, vertex_t v, f32 w) {
                if (is_not_endpoint(m_neighbors[u], u) || is_not_endpoint(m_neighbors[v], v)) {
                    return;
                }

                bool u_unmatched = is_unmatched(m_neighbors[u], u);
                bool v_unmatched = is_unmatched(m_neighbors[v], v);

                if (u_unmatched && v_unmatched) {
                    m_neighbors[u].n1 = v;
                    m_neighbors[u].w1 = w;
                    m_neighbors[v].n1 = u;
                    m_neighbors[v].w1 = w;
                    path_id[u] = u;
                    path_id[v] = u;
                    path_length[u] = 1;
                    return;
                }

                u32 u_id = path_id[u];
                u32 v_id = path_id[v];

                if (u_unmatched || v_unmatched) {
                    if (v_unmatched) {
                        std::swap(u, v);
                        std::swap(u_id, v_id);
                    }
                    m_neighbors[u].n1 = v;
                    m_neighbors[u].w1 = w;
                    m_neighbors[v].n2 = u;
                    m_neighbors[v].w2 = w;
                    path_id[u] = v_id;
                    path_length[v_id] += 1;
                    return;
                }

                if (u_id == v_id) {
                    if (path_length[u_id] & 1) {
                        path_length[u_id] += 1;
                        m_neighbors[u].n2 = v;
                        m_neighbors[u].w2 = w;
                        m_neighbors[v].n2 = u;
                        m_neighbors[v].w2 = w;
                        cycles.push_back(u);
                    }
                    return;
                }

                m_neighbors[u].n2 = v;
                m_neighbors[u].w2 = w;
                m_neighbors[v].n2 = u;
                m_neighbors[v].w2 = w;

                vertex_t v1 = v;
                vertex_t v2 = u;
                u32 id1 = v_id;
                u32 id2 = u_id;
                if (path_length[u_id] > path_length[v_id]) {
                    std::swap(v1, v2);
                    std::swap(id1, id2);
                }

                path_length[id1] += 1 + path_length[id2];
                while (m_neighbors[v2].n2 != v2) {
                    path_id[v2] = id1;
                    vertex_t temp_last_vertex = v1;
                    v1 = v2;
                    v2 = m_neighbors[v2].n1 == temp_last_vertex ? m_neighbors[v2].n2 : m_neighbors[v2].n1;
                }
                path_id[v2] = id1;
            };

            if (global_min_rating == global_max_rating) {
                // Fast path: Process edges in arbitrary order since all ratings are equal
                for (u64 t_idx = 0; t_idx < m_threads; ++t_idx) {
                    for (const auto &e: m_thread_infos[t_idx].local_edges) {
                        process_edge(e.u, e.v, e.w);
                    }
                }
            } else {
                std::priority_queue<HeapEntry> edge_queue;
                for (u64 i = 0; i < m_threads; ++i) {
                    if (!m_thread_infos[i].local_edges.empty()) {
                        edge_queue.push({i, 0, m_thread_infos[i].local_edges[0]});
                    }
                }

                while (!edge_queue.empty()) {
                    HeapEntry top = edge_queue.top();
                    edge_queue.pop();

                    u64 t_idx = top.thread_idx;
                    size_t next_idx = top.edge_idx + 1;

                    while (next_idx < m_thread_infos[t_idx].local_edges.size()) {
                        const auto &e = m_thread_infos[t_idx].local_edges[next_idx];
                        if (is_not_endpoint(m_neighbors[e.u], e.u) || is_not_endpoint(m_neighbors[e.v], e.v)) {
                            next_idx++;
                            continue;
                        }
                        edge_queue.push({t_idx, next_idx, e});
                        break;
                    }

                    process_edge(top.edge.u, top.edge.v, top.edge.w);
                }
            }

            HEIPROMAP_PROFILE_SCOPE("coarsening", "GlobalPathAlgorithmMatcher", "solve_paths");
            #pragma omp parallel for num_threads(m_threads) schedule(static, 32768)
            for (vertex_t u = 0; u < g.n; ++u) {
                u64 thread_id = omp_get_thread_num();
                if (is_one_endpoint(m_neighbors[u], u)) {
                    vertex_t v1 = u;
                    vertex_t v2 = m_neighbors[u].n1;
                    while (m_neighbors[v2].n2 != v2) {
                        vertex_t temp_last_vertex = v1;
                        v1 = v2;
                        v2 = m_neighbors[v2].n1 == temp_last_vertex ? m_neighbors[v2].n2 : m_neighbors[v2].n1;
                    }
                    if (u < v2) {
                        solve_path(g, u, path_length[path_id[u]], [&](vertex_t uu, vertex_t vv) { matching.add(uu, vv); }, thread_id);
                    }
                }
            }

            HEIPROMAP_PROFILE_SCOPE("coarsening", "GlobalPathAlgorithmMatcher", "solve_cycles");
            #pragma omp parallel for num_threads(m_threads) schedule(static)
            for (size_t i = 0; i < cycles.size(); ++i) {
                u64 thread_id = omp_get_thread_num();
                vertex_t u = cycles[i];
                solve_cycle(g, u, path_length[path_id[u]], matching, thread_id);
            }

            if ((f64) matching.size() * 2 < config->two_hop_threshold * (f64) g.n) {
                two_hop_degree_one<t_uniform_v_weights, t_uniform_e_weights>(level, g, p_manager, matching, imbalance);
            }

            if ((f64) matching.size() * 2 < config->two_hop_threshold * (f64) g.n) {
                two_hop_twins<t_uniform_v_weights, t_uniform_e_weights>(level, g, p_manager, matching, imbalance);
            }

            if ((f64) matching.size() * 2 < config->two_hop_threshold * (f64) g.n) {
                two_hop_matchmaker<t_uniform_v_weights, t_uniform_e_weights>(level, g, p_manager, matching, imbalance);
            }

            finalize_matching(g, matching, mapping);
        }

        void match(const size_t level,
                   const graph_t &g,
                   const p_manager_t &p_manager,
                   Mapping &mapping,
                   f64 imbalance) {
            auto dispatch_with_rating = [&](auto rating_func_const) {
                constexpr EdgeRatingFunction rating_func = rating_func_const;
                if (g.uniform_v_weights && g.uniform_e_weights) {
                    match_templated<true, true, rating_func>(level, g, p_manager, mapping, imbalance);
                } else if (g.uniform_v_weights) {
                    match_templated<true, false, rating_func>(level, g, p_manager, mapping, imbalance);
                } else if (g.uniform_e_weights) {
                    match_templated<false, true, rating_func>(level, g, p_manager, mapping, imbalance);
                } else {
                    match_templated<false, false, rating_func>(level, g, p_manager, mapping, imbalance);
                }
            };

            switch (config->rating_function) {
                case EdgeRatingFunction::WEIGHT:
                    dispatch_with_rating(std::integral_constant<EdgeRatingFunction, EdgeRatingFunction::WEIGHT>{});
                    break;
                case EdgeRatingFunction::EXPANSION:
                    dispatch_with_rating(std::integral_constant<EdgeRatingFunction, EdgeRatingFunction::EXPANSION>{});
                    break;
                case EdgeRatingFunction::HEAVY_EDGE:
                    dispatch_with_rating(std::integral_constant<EdgeRatingFunction, EdgeRatingFunction::HEAVY_EDGE>{});
                    break;
                case EdgeRatingFunction::GREEDY:
                    dispatch_with_rating(std::integral_constant<EdgeRatingFunction, EdgeRatingFunction::GREEDY>{});
                    break;
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights, EdgeRatingFunction t_rating_function>
        void compute_ratings(const graph_t &g, const p_manager_t &p_manager, weight_t lmax) {
            HEIPROMAP_PROFILE_SCOPE("coarsening", "GlobalPathAlgorithmMatcher", "compute_ratings");
            for (u64 i = 0; i < m_threads; ++i) {
                m_thread_infos[i].local_edges.clear();
                m_thread_infos[i].min_rating = std::numeric_limits<f32>::max();
                m_thread_infos[i].max_rating = std::numeric_limits<f32>::min();
            }

            #pragma omp parallel for num_threads(m_threads) schedule(guided)
            for (vertex_t u = 0; u < g.n; ++u) {
                u64 thread_id = omp_get_thread_num();
                weight_t u_w = t_uniform_v_weights ? 1 : g.v_weights[u];
                partition_t u_id = p_manager[u];

                f32 local_min = m_thread_infos[thread_id].min_rating;
                f32 local_max = m_thread_infos[thread_id].max_rating;

                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const vertex_t v = g.edges_v[j];
                    const weight_t w = g.edges_w[j];
                    if (u >= v) { continue; }
                    if (u_id != p_manager[v]) { continue; }
                    weight_t v_w = t_uniform_v_weights ? 1 : g.v_weights[v];

                    if (u_w + v_w > lmax) { continue; }

                    weight_t ew = t_uniform_e_weights ? 1 : w;

                    f32 edge_rating;
                    if constexpr (t_uniform_v_weights && t_uniform_e_weights) {
                        edge_rating = 1.0f;
                    } else {
                        if constexpr (t_rating_function == EdgeRatingFunction::WEIGHT) {
                            edge_rating = (f32) ew;
                        } else if constexpr (t_rating_function == EdgeRatingFunction::EXPANSION) {
                            edge_rating = (f32) ew / (f32) (u_w * v_w);
                        } else if constexpr (t_rating_function == EdgeRatingFunction::HEAVY_EDGE) {
                            edge_rating = (f32) (ew * ew) / (f32) (u_w * v_w);
                        } else if constexpr (t_rating_function == EdgeRatingFunction::GREEDY) {
                            edge_rating = (f32) ew / std::sqrt((f32) u_w * (f32) v_w);
                        }
                    }

                    local_min = std::min(local_min, edge_rating);
                    local_max = std::max(local_max, edge_rating);
                    m_thread_infos[thread_id].local_edges.emplace_back(u, v, edge_rating);
                }

                m_thread_infos[thread_id].min_rating = local_min;
                m_thread_infos[thread_id].max_rating = local_max;
            }
        }

        template<typename AddMatchFunc>
        f32 solve_path(const graph_t &g, const vertex_t u, const u32 length, AddMatchFunc add_match, u64 thread_id) {
            if (length == 1) {
                vertex_t v = m_neighbors[u].n1;
                add_match(u, v);
                return m_neighbors[u].w1;
            }
            if (length == 2) {
                vertex_t v1 = u;
                vertex_t v2 = m_neighbors[u].n1;
                vertex_t v3 = (m_neighbors[v2].n1 == v1) ? m_neighbors[v2].n2 : m_neighbors[v2].n1;
                f32 w1 = m_neighbors[u].w1;
                f32 w2 = (m_neighbors[v2].n1 == v1) ? m_neighbors[v2].w2 : m_neighbors[v2].w1;

                if (w1 > w2) {
                    add_match(v1, v2);
                    return w1;
                } else {
                    add_match(v2, v3);
                    return w2;
                }
            }

            auto &ti = m_thread_infos[thread_id];
            vertex_t v1 = u, v2, v3;
            f32 w;
            s64 i = 0;

            ti.dp[i].edge = v1;
            v2 = (m_neighbors[v1].n1 == v1) ? m_neighbors[v1].n2 : m_neighbors[v1].n1;
            w = (m_neighbors[v1].n1 == v1) ? m_neighbors[v1].w2 : m_neighbors[v1].w1;
            ti.dp[i + 1].edge = v2;
            ti.dp[i].w = w;
            ti.dp[i].m = -1;
            ti.dp[i].take = 1;
            i++;

            v3 = (m_neighbors[v2].n1 == v1) ? m_neighbors[v2].n2 : m_neighbors[v2].n1;
            w = (m_neighbors[v2].n1 == v1) ? m_neighbors[v2].w2 : m_neighbors[v2].w1;
            ti.dp[i + 1].edge = v3;
            if (w > ti.dp[i - 1].w) {
                ti.dp[i].w = w;
                ti.dp[i].m = -1;
                ti.dp[i].take = 1;
            } else {
                ti.dp[i].w = ti.dp[i - 1].w;
                ti.dp[i].m = 0;
                ti.dp[i].take = 0;
            }
            i++;

            v1 = v2;
            v2 = v3;
            while (m_neighbors[v2].n1 != v2 && m_neighbors[v2].n2 != v2) {
                v3 = (m_neighbors[v2].n1 == v1) ? m_neighbors[v2].n2 : m_neighbors[v2].n1;
                w = (m_neighbors[v2].n1 == v1) ? m_neighbors[v2].w2 : m_neighbors[v2].w1;
                ti.dp[i + 1].edge = v3;

                if (w + ti.dp[i - 2].w > ti.dp[i - 1].w) {
                    ti.dp[i].w = w + ti.dp[i - 2].w;
                    ti.dp[i].m = i - 2;
                    ti.dp[i].take = 1;
                } else {
                    ti.dp[i].w = ti.dp[i - 1].w;
                    ti.dp[i].m = i - 1;
                    ti.dp[i].take = 0;
                }
                v1 = v2;
                v2 = v3;
                i++;
            }

            s64 idx = i - 1;
            while (idx != -1) {
                if (ti.dp[idx].take) {
                    add_match(ti.dp[idx].edge, ti.dp[idx + 1].edge);
                }
                idx = ti.dp[idx].m;
            }
            return ti.dp[i - 1].w;
        }

        void solve_cycle(const graph_t &g, const vertex_t u, const u32 length, Matching &matching, u64 thread_id) {
            auto &ti = m_thread_infos[thread_id];
            vertex_t n1 = m_neighbors[u].n1, n2 = m_neighbors[u].n2;
            Neighbors original_u = m_neighbors[u], original_n1 = m_neighbors[n1], original_n2 = m_neighbors[n2];

            ti.dp_cycle_matches1.clear();
            ti.dp_cycle_matches2.clear();

            m_neighbors[u].n1 = original_u.n2;
            m_neighbors[u].w1 = original_u.w2;
            m_neighbors[u].n2 = u;
            if (m_neighbors[n1].n1 == u) {
                m_neighbors[n1].n1 = original_n1.n2;
                m_neighbors[n1].w1 = original_n1.w2;
                m_neighbors[n1].n2 = n1;
            } else { m_neighbors[n1].n2 = n1; }
            f32 w1 = solve_path(g, u, length - 1, [&](vertex_t uu, vertex_t vv) { ti.dp_cycle_matches1.emplace_back(uu, vv); }, thread_id);
            m_neighbors[u] = original_u;
            m_neighbors[n1] = original_n1;

            m_neighbors[u].n2 = u;
            if (m_neighbors[n2].n1 == u) {
                m_neighbors[n2].n1 = original_n2.n2;
                m_neighbors[n2].w1 = original_n2.w2;
                m_neighbors[n2].n2 = n2;
            } else { m_neighbors[n2].n2 = n2; }
            f32 w2 = solve_path(g, u, length - 1, [&](vertex_t uu, vertex_t vv) { ti.dp_cycle_matches2.emplace_back(uu, vv); }, thread_id);
            m_neighbors[u] = original_u;
            m_neighbors[n2] = original_n2;

            if (w1 > w2) {
                for (auto &[uu, vv]: ti.dp_cycle_matches1) matching.add(uu, vv);
            } else {
                for (auto &[uu, vv]: ti.dp_cycle_matches2) matching.add(uu, vv);
            }
        }

        template<bool t_uniform_v_weights>
        void random_matching(const size_t, const graph_t &g, Matching &matching, weight_t lmax) {
            HEIPROMAP_PROFILE_SCOPE("coarsening", "GlobalPathAlgorithmMatcher", "random_matching");
            for (vertex_t u = 0; u < g.n; ++u) {
                if (matching.is_matched(u)) { continue; }
                weight_t u_w = t_uniform_v_weights ? 1 : g.v_weights[u];
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const vertex_t v = g.edges_v[j];
                    if (matching.is_matched(v)) { continue; }
                    weight_t v_w = t_uniform_v_weights ? 1 : g.v_weights[v];
                    if (u_w + v_w > lmax) { continue; }
                    matching.add(u, v);
                    break;
                }
            }
        }

        void finalize_matching(const graph_t &g, Matching &matching, Mapping &mapping) {
            HEIPROMAP_PROFILE_SCOPE("coarsening", "GlobalPathAlgorithmMatcher", "get_mapping");
            matching.set_translation();
            mapping.set_coarse_n(matching.get_n_coarse_nodes());
            for (vertex_t u = 0; u < matching.get_n(); ++u) {
                mapping.set(u, matching.get_n(u));
            }
        }

        static inline u64 splitmix64(u64 x) {
            x += 0x9e3779b97f4a7c15ULL;
            x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
            x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
            return x ^ (x >> 31);
        }

        static f64 small_noise(vertex_t u, vertex_t v) {
            u64 a = static_cast<u64>(std::min(u, v)), b = static_cast<u64>(std::max(u, v));
            u64 key = a * 0x9e3779b97f4a7c15ULL + b;
            u64 h = splitmix64(key);
            return 1e-12 * (static_cast<f64>(h) / static_cast<f64>(std::numeric_limits<u64>::max()));
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void two_hop_degree_one(const size_t, const graph_t &g, const p_manager_t &p_manager, Matching &matching, f64 imbalance) {
            HEIPROMAP_PROFILE_SCOPE("coarsening", "GlobalPathAlgorithmMatcher", "two_hop_degree_one");
            std::vector<vertex_t> preferred(g.n);
            std::iota(preferred.begin(), preferred.end(), 0);
            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            for (vertex_t u = 0; u < g.n; ++u) {
                if (g.deg(u) != 1 || matching.is_matched(u)) { continue; }
                vertex_t mid = g.edges_v[g.neighborhoods[u]];
                weight_t mid_w = t_uniform_e_weights ? 1 : g.edges_w[g.neighborhoods[u]];
                weight_t u_w = t_uniform_v_weights ? 1 : g.v_weights[u];
                f64 best_rating = -1e18;

                for (size_t i = g.neighborhoods[mid]; i < g.neighborhoods[mid + 1]; ++i) {
                    const vertex_t v = g.edges_v[i];
                    const weight_t w = g.edges_w[i];
                    if (u == v || g.deg(v) != 1 || matching.is_matched(v)) { continue; }
                    weight_t v_w = t_uniform_v_weights ? 1 : g.v_weights[v];
                    if (u_w + v_w > lmax) { continue; }

                    weight_t mw = t_uniform_e_weights ? 1 : w;
                    f64 rating = (f64) (mw + mid_w) + small_noise(u, v);
                    if (rating > best_rating) {
                        best_rating = rating;
                        preferred[u] = v;
                    }
                }
            }
            for (vertex_t u = 0; u < g.n; ++u) { if (g.deg(u) == 1 && !matching.is_matched(u) && preferred[u] != u && preferred[preferred[u]] == u && u < preferred[u]) matching.add(u, preferred[u]); }
        }

        static inline uint64_t hash_combine_u64(uint64_t a, uint64_t b) { return splitmix64(a ^ (b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2))); }

        template<bool t_uniform_e_weights>
        static inline uint64_t hash_edge(vertex_t v, weight_t w) {
            uint64_t hv = splitmix64(v);
            if constexpr (t_uniform_e_weights) {
                return hv;
            } else {
                return hash_combine_u64(hv, splitmix64(w));
            }
        }

        template<bool t_uniform_e_weights>
        static inline uint64_t neighborhood_hash(const graph_t &g, vertex_t u) {
            uint64_t x = splitmix64(g.deg(u)), s1 = 0, s2 = 0;
            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t v = g.edges_v[i];
                const weight_t w = g.edges_w[i];
                uint64_t he = hash_edge<t_uniform_e_weights>(v, w);
                x ^= he;
                s1 += he;
                s2 += splitmix64(he);
            }
            return hash_combine_u64(x, hash_combine_u64(s1, s2));
        }

        template<bool t_uniform_e_weights>
        static inline bool same_neighborhood(const graph_t &g, vertex_t u, vertex_t v) {
            if (g.deg(u) != g.deg(v)) return false;
            std::vector<std::pair<vertex_t, weight_t> > nu, nv;
            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t x = g.edges_v[i];
                const weight_t w = g.edges_w[i];
                nu.emplace_back(x, t_uniform_e_weights ? 1 : w);
            }
            for (size_t i = g.neighborhoods[v]; i < g.neighborhoods[v + 1]; ++i) {
                const vertex_t x = g.edges_v[i];
                const weight_t w = g.edges_w[i];
                nv.emplace_back(x, t_uniform_e_weights ? 1 : w);
            }
            std::sort(nu.begin(), nu.end());
            std::sort(nv.begin(), nv.end());
            return nu == nv;
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void two_hop_twins(const size_t, const graph_t &g, const p_manager_t &p_manager, Matching &matching, f64 imbalance) {
            HEIPROMAP_PROFILE_SCOPE("coarsening", "GlobalPathAlgorithmMatcher", "two_hop_twins");
            struct Candidate {
                uint64_t hash;
                vertex_t u;
            };
            std::vector<Candidate> candidates(g.n);
            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            size_t n_unmatched = 0;
            //
            {
                std::vector<std::vector<Candidate> > local_candidates(m_threads);
                for (auto &v: local_candidates) v.reserve(g.n / m_threads);

                #pragma omp parallel num_threads(m_threads)
                {
                    u64 tid = omp_get_thread_num();
                    #pragma omp for schedule(static)
                    for (vertex_t u = 0; u < g.n; ++u) {
                        if (!matching.is_matched(u)) {
                            local_candidates[tid].push_back({neighborhood_hash<t_uniform_e_weights>(g, u), u});
                        }
                    }
                }

                std::vector<size_t> offsets(m_threads + 1, 0);
                for (u64 i = 0; i < m_threads; ++i) offsets[i + 1] = offsets[i] + local_candidates[i].size();
                n_unmatched = offsets[m_threads];
                candidates.resize(n_unmatched);

                #pragma omp parallel for num_threads(m_threads) schedule(static)
                for (u64 i = 0; i < m_threads; ++i) {
                    std::copy(local_candidates[i].begin(), local_candidates[i].end(), candidates.begin() + offsets[i]);
                }
            }

            std::sort(candidates.begin(), candidates.end(), [](const Candidate &a, const Candidate &b) { return (a.hash != b.hash) ? a.hash < b.hash : a.u < b.u; });

            for (size_t begin = 0, end; begin < candidates.size(); begin = end) {
                for (end = begin + 1; end < candidates.size() && candidates[end].hash == candidates[begin].hash; ++end) {}
                std::vector<std::vector<vertex_t> > groups;
                for (size_t i = begin; i < end; ++i) {
                    vertex_t u = candidates[i].u;
                    if (matching.is_matched(u)) continue;
                    bool placed = false;
                    for (auto &group: groups)
                        if (same_neighborhood<t_uniform_e_weights>(g, u, group.front())) {
                            group.push_back(u);
                            placed = true;
                            break;
                        }
                    if (!placed) groups.push_back({u});
                }
                for (auto &group: groups) {
                    if constexpr (!t_uniform_v_weights) {
                        std::sort(group.begin(), group.end(), [&](vertex_t a, vertex_t b) { return (g.v_weights[a] != g.v_weights[b]) ? g.v_weights[a] < g.v_weights[b] : a < b; });
                    }
                    for (size_t i = 0; i + 1 < group.size(); ++i) {
                        if (matching.is_matched(group[i])) continue;
                        weight_t u_w = t_uniform_v_weights ? 1 : g.v_weights[group[i]];
                        for (size_t j = i + 1; j < group.size(); ++j) {
                            weight_t v_w = t_uniform_v_weights ? 1 : g.v_weights[group[j]];
                            if (!matching.is_matched(group[j]) && u_w + v_w <= lmax) {
                                matching.add(group[i], group[j]);
                                break;
                            }
                        }
                    }
                }
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void two_hop_matchmaker(const size_t, const graph_t &g, const p_manager_t &p_manager, Matching &matching, f64 imbalance) {
            HEIPROMAP_PROFILE_SCOPE("coarsening", "GlobalPathAlgorithmMatcher", "two_hop_matchmaker");
            std::vector<vertex_t> preferred(g.n);
            std::iota(preferred.begin(), preferred.end(), 0);
            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));
            for (vertex_t u = 0; u < g.n; ++u) {
                if (matching.is_matched(u)) continue;
                weight_t u_w = t_uniform_v_weights ? 1 : g.v_weights[u];
                f64 best_rating = -1e18;

                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const vertex_t mid = g.edges_v[j];
                    const weight_t mw_orig = g.edges_w[j];
                    weight_t mid_w = t_uniform_e_weights ? 1 : mw_orig;
                    for (size_t i = g.neighborhoods[mid]; i < g.neighborhoods[mid + 1]; ++i) {
                        const vertex_t v = g.edges_v[i];
                        const weight_t w = g.edges_w[i];
                        if (u != v && !matching.is_matched(v)) {
                            weight_t v_w = t_uniform_v_weights ? 1 : g.v_weights[v];
                            if (u_w + v_w <= lmax) {
                                weight_t mw = t_uniform_e_weights ? 1 : w;
                                f64 rating = (f64) (mw + mid_w) + small_noise(u, v);
                                if (rating > best_rating) {
                                    best_rating = rating;
                                    preferred[u] = v;
                                }
                            }
                        }
                    }
                }
            }
            for (vertex_t u = 0; u < g.n; ++u) { if (!matching.is_matched(u) && preferred[u] != u && preferred[preferred[u]] == u && u < preferred[u]) matching.add(u, preferred[u]); }
        }
    };
}

#endif //HEIPROMAP_GLOBAL_PATH_ALGORITHM_H
