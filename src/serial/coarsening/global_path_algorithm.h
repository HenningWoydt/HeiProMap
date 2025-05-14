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

#include "../../commons/definitions.h"
#include "../../commons/JSON_utils.h"
#include "../../commons/random_engine.h"
#include "../../commons/statistic_collector.h"
#include "../../commons/utils.h"
#include "../interfaces/ISerialMatcher.h"

namespace HeiProMap {
    struct Neighbors {
        vertex_t n1;
        vertex_t n2;
        f32      w1;
        f32      w2;
    };

#define IS_ENDPOINT(n, u) n.n2 == u
#define IS_NOT_ENDPOINT(n, u) n.n2 != u
#define IS_ONE_ENDPOINT(n, u) n.n1 != u && n.n2 == u
#define IS_UNMATCHED(n, u) n.n1 == u

    class GlobalPathAlgorithmConfiguration final : public ISerialMatcherConfiguration {
    public:
        size_t random_level = 4;
    };

    /**
     * Computes a matching based on the Global Path Algorithm from
     * > Jens Maue and Peter Sanders.
     * > Engineering Algorithms for Approximate Weighted Matching.
     * > Experimental Algorithms, 6th International Workshop, WEA 2007, Rome, Italy, June 6-8, 2007, Proceedings.
     */
    class GlobalPathAlgorithmMatcher final : public ISerialMatcher {
        vertex_t    m_n     = 0;
        vertex_t    m_m     = 0;
        partition_t m_k     = 0;
        weight_t    m_l_max = 0;

        const GlobalPathAlgorithmConfiguration *config           = nullptr;
        RandomEngine                           *random_engine    = nullptr;
        StatisticCollector                     *m_stat_collector = nullptr;

        AlignedArray<Neighbors> m_neighbors;
        AlignedArray<u32> path_id;
        AlignedArray<u32> path_length;

        AlignedArray<EdgeUVW> edges;
        size_t edges_size = 0;

        // for DP
        AlignedArray<f32> dp_w;
        AlignedArray<s64> dp_m;
        AlignedArray<u8> dp_take;
        AlignedArray<vertex_t> dp_edges;

        Matching dp_cycle_matches1;
        Matching dp_cycle_matches2;

        METRICS(std::vector<f64> level_time_compute_ratings;)
        METRICS(std::vector<f64> level_time_sorting;)
        METRICS(std::vector<f64> level_time_build_paths;)
        METRICS(std::vector<f64> level_time_solve_paths;)
        METRICS(std::vector<f64> level_time_solve_cycles;)
        METRICS(std::vector<f64> level_time_random;)

        METRICS(std::vector<u64> level_edges;)
        METRICS(std::vector<u64> level_edges_skipped;)
        METRICS(std::vector<u64> level_edges_form_cycle;)
        METRICS(std::vector<u64> level_edges_new_paths;)
        METRICS(std::vector<u64> level_edges_enlarge_path;)
        METRICS(std::vector<u64> level_edges_combine_paths;)

    public:
        GlobalPathAlgorithmMatcher() = default;

        ~GlobalPathAlgorithmMatcher() override = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const weight_t t_l_max,
                        RandomEngine &t_random_engine,
                        const ISerialMatcherConfiguration &i_config,
                        StatisticCollector &t_stat_collect) override {
            m_n     = t_n;
            m_m     = t_m;
            m_k     = t_k;
            m_l_max = t_l_max;

            config           = dynamic_cast<const GlobalPathAlgorithmConfiguration *>(&i_config);
            random_engine    = &t_random_engine;
            m_stat_collector = &t_stat_collect;

            vertex_t t_n_64 = round_up_64(t_n);
            vertex_t t_m_64 = round_up_64(t_m);

            m_neighbors.initialize(m_n);
            path_id.initialize(m_n);
            path_length.initialize(m_n);

            edges.initialize(m_m);

            dp_w.initialize(m_n);
            dp_m.initialize(m_n);
            dp_take.initialize(m_n);
            dp_edges.initialize(m_n);

            dp_cycle_matches1.initialize(t_n_64);
            dp_cycle_matches2.initialize(t_n_64);
        }

        void match(const size_t level,
                   const graph_t &g,
                   p_manager_t &p_manager,
                   Matching &matching) override {
            METRICS(level_time_compute_ratings.emplace_back(0.0);)
            METRICS(level_time_sorting.emplace_back(0.0);)
            METRICS(level_time_build_paths.emplace_back(0.0);)
            METRICS(level_time_solve_paths.emplace_back(0.0);)
            METRICS(level_time_solve_cycles.emplace_back(0.0);)
            METRICS(level_time_random.emplace_back(0.0);)

            METRICS(level_edges.emplace_back(0);)
            METRICS(level_edges_skipped.emplace_back(0);)
            METRICS(level_edges_form_cycle.emplace_back(0);)
            METRICS(level_edges_new_paths.emplace_back(0);)
            METRICS(level_edges_enlarge_path.emplace_back(0);)
            METRICS(level_edges_combine_paths.emplace_back(0);)

            if (level < config->random_level) {
                // use a random matching
                random_matching(level, g, matching);
                return;
            }

            compute_ratings(g, p_manager);

            METRICS_TIME(sp_sorting)
            std::sort(edges.get_ptr(), edges.get_ptr() + edges_size, std::greater<>());
            METRICS_TIME(ep_sorting)
            METRICS(level_time_sorting.back() += get_seconds(sp_sorting, ep_sorting);)

            METRICS_TIME(sp_build_paths)
            for (vertex_t u = 0; u < g.get_n(); ++u) {
                m_neighbors[u].n1 = u;
                m_neighbors[u].n2 = u;
            }

            for (size_t i = 0; i < edges_size; ++i) {
                auto [u, v, w] = edges[i];

                if (IS_NOT_ENDPOINT(m_neighbors[u], u) || IS_NOT_ENDPOINT(m_neighbors[v], v)) {
                    // u or v is not an endpoint
                    METRICS(level_edges_skipped.back() += 1;)
                    continue;
                }

                bool u_unmatched = IS_UNMATCHED(m_neighbors[u], u);
                bool v_unmatched = IS_UNMATCHED(m_neighbors[v], v);

                if (u_unmatched && v_unmatched) {
                    // both are unmatched, only one new path of length 1
                    m_neighbors[u].n1 = v;
                    m_neighbors[u].w1 = w;
                    m_neighbors[v].n1 = u;
                    m_neighbors[v].w1 = w;
                    path_id[u]     = u;
                    path_id[v]     = u;
                    path_length[u] = 1;
                    METRICS(level_edges_new_paths.back() += 1;)
                    continue;
                }

                u32 u_id = path_id[u];
                u32 v_id = path_id[v];

                if (u_unmatched || v_unmatched) {
                    if (v_unmatched) {
                        std::swap(u, v);
                        std::swap(u_id, v_id);
                    }
                    // only one unmatched, enlarge path
                    m_neighbors[u].n1 = v;
                    m_neighbors[u].w1 = w;
                    m_neighbors[v].n2 = u;
                    m_neighbors[v].w2 = w;
                    path_id[u] = v_id;
                    path_length[v_id] += 1;
                    METRICS(level_edges_enlarge_path.back() += 1;)
                    continue;
                }

                // cycle
                if (u_id == v_id) {
                    if (path_length[u_id] & 1) {
                        // same path and odd length size, close the cycle
                        path_length[u_id] += 1; // increase path length

                        // for u set v as a neighbor
                        m_neighbors[u].n2 = v;
                        m_neighbors[u].w2 = w;

                        // for v set u as a neighbor
                        m_neighbors[v].n2 = u;
                        m_neighbors[v].w2 = w;

                        // solve the cycle
                        METRICS_TIME(sp_solve_cycle)
                        solve_cycle(g, u, path_length[u_id], matching);
                        METRICS_TIME(ep_solve_cycle)
                        path_length[u_id] = 0;

                        METRICS(level_time_solve_cycles.back() += get_seconds(sp_solve_cycle, ep_solve_cycle);)
                        METRICS(level_time_build_paths.back() -= get_seconds(sp_solve_cycle, ep_solve_cycle);)
                        METRICS(level_edges_form_cycle.back() += 1;)
                    }
                    continue;
                }

                // two paths, both u and v connect larger paths

                // for u set v as a neighbor
                m_neighbors[u].n2 = v;
                m_neighbors[u].w2 = w;

                // for v set u as a neighbor
                m_neighbors[v].n2 = u;
                m_neighbors[v].w2 = w;

                vertex_t v1  = v;
                vertex_t v2  = u;
                u32      id1 = v_id;
                u32      id2 = u_id;
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

                METRICS(level_edges_combine_paths.back() += 1;)
            }
            METRICS_TIME(ep_build_paths)
            METRICS(level_time_build_paths.back() += get_seconds(sp_build_paths, ep_build_paths);)

            // process all paths
            METRICS_TIME(sp_solve_paths)
            forall_gu(g, u)
                {
                    if (IS_ONE_ENDPOINT(m_neighbors[u], u) && path_length[path_id[u]] > 0) {
                        solve_path(g, u, path_length[path_id[u]], matching);
                        path_length[path_id[u]] = 0;
                    }
                }
            endfor
            METRICS_TIME(ep_solve_paths)
            METRICS(level_time_solve_paths.back() += get_seconds(sp_solve_paths, ep_solve_paths);)

#if ASSERT_ENABLED
            for (size_t i = 0; i < matching.size(); ++i) {
                const auto& [u, v] = matching[i];
                ASSERT(u != v);
            }
#endif

#if ASSERT_ENABLED
            std::vector<u8> hit(g.get_n(), 0);
            for (size_t i = 0; i < matching.size(); ++i) {
                const auto& [u, v] = matching[i];
                hit[u] += 1;
                hit[v] += 1;

                ASSERT(hit[u] == 1);
                ASSERT(hit[v] == 1);
            }
#endif
        }

        void compute_ratings_4(const graph_t &g) {
            METRICS_TIME(sp_compute_ratings)

            edges_size = 0;
            forall_gu(g, u)
                {
                    weight_t u_w        = g.weight(u);
                    weight_t max_weight = m_l_max - u_w;

                    size_t j = 0;
                    for (; j + 4 < g.size(u); j += 4) {
                        const vertex_t v0 = g.neighbor(u, j + 0);
                        const vertex_t v1 = g.neighbor(u, j + 1);
                        const vertex_t v2 = g.neighbor(u, j + 2);
                        const vertex_t v3 = g.neighbor(u, j + 3);

                        const weight_t w0 = g.weight(u, j + 0);
                        const weight_t w1 = g.weight(u, j + 1);
                        const weight_t w2 = g.weight(u, j + 2);
                        const weight_t w3 = g.weight(u, j + 3);

                        weight_t v0_w = g.weight(v0);
                        weight_t v1_w = g.weight(v1);
                        weight_t v2_w = g.weight(v2);
                        weight_t v3_w = g.weight(v3);

                        f32 edge_rating_0 = ((f32) (w0 * w0)) / ((f32) (u_w * v0_w));
                        f32 edge_rating_1 = ((f32) (w1 * w1)) / ((f32) (u_w * v1_w));
                        f32 edge_rating_2 = ((f32) (w2 * w2)) / ((f32) (u_w * v2_w));
                        f32 edge_rating_3 = ((f32) (w3 * w3)) / ((f32) (u_w * v3_w));

                        edges[edges_size] = {u, v0, edge_rating_0};
                        edges_size += (u < v0) && (v0_w <= max_weight);

                        edges[edges_size] = {u, v1, edge_rating_1};
                        edges_size += (u < v1) && (v1_w <= max_weight);

                        edges[edges_size] = {u, v2, edge_rating_2};
                        edges_size += (u < v2) && (v2_w <= max_weight);

                        edges[edges_size] = {u, v3, edge_rating_3};
                        edges_size += (u < v3) && (v3_w <= max_weight);
                    }

                    for (; j < g.size(u); ++j) {
                        const vertex_t v = g.neighbor(u, j);
                        const weight_t w = g.weight(u, j);

                        if (u > v) { continue; }
                        weight_t v_w = g.weight(v);

                        // if (u_w > 1.5*av_manager.get_n_active() / 20.0 * m_k) { continue; }
                        // if (v_w > 1.5*av_manager.get_n_active() / 20.0 * m_k) { continue; }

                        if (v_w > max_weight) { continue; }

                        f32 edge_rating;

                        // edge_rating = ((f32) w) / (g.size(u) * g.size(v));
                        // edge_rating = (f32) w / (f32) (u_w * v_w);
                        // edge_rating = (f32) w / (f32) (u_w * v_w);
                        edge_rating = ((f32) (w * w)) / ((f32) (u_w * v_w));
                        // edge_rating = ((f32) (w * w)) / ((f32) (u_w + v_w));
                        // edge_rating = ((f32) (w * w * w)) / ((f32) (u_w * v_w));
                        // edge_rating = (f32) w / (f32) (u_w * v_w * u_w * v_w);
                        // edge_rating = (f32)w;

                        edges[edges_size++] = {u, v, edge_rating};
                    }
                }
            endfor
            METRICS_TIME(ep_compute_ratings)
            METRICS(level_time_compute_ratings.back() += get_seconds(sp_compute_ratings, ep_compute_ratings);)
            METRICS(level_edges.back() += edges_size;)
        }

        void compute_ratings(const graph_t &g, const p_manager_t &p_manager) {
            // compute_ratings_4(g);
            // return;
            METRICS_TIME(sp_compute_ratings)

            edges_size = 0;
            forall_gu(g, u)
                {
                    weight_t u_w = g.weight(u);

                    forall_guivw(g, u, j, v, w)
                        {
                            if (u > v) { continue; }
                            if (p_manager[u] != p_manager[v]) { continue; }
                            weight_t v_w = g.weight(v);

                            // if (u_w > 1.5*av_manager.get_n_active() / 20.0 * m_k) { continue; }
                            // if (v_w > 1.5*av_manager.get_n_active() / 20.0 * m_k) { continue; }

                            if (u_w + v_w > m_l_max) { continue; }

                            f32 edge_rating;

                            // edge_rating = ((f32) w) / (g.size(u) * g.size(v));
                            // edge_rating = (f32) w / (f32) (u_w * v_w);
                            edge_rating = ((f32) (w * w)) / ((f32) (u_w * v_w));
                            // edge_rating = ((f32) (w * w)) / ((f32) (u_w + v_w));
                            // edge_rating = ((f32) (w * w * w)) / ((f32) (u_w * v_w));
                            // edge_rating = (f32) w / (f32) (u_w * v_w * u_w * v_w);
                            // edge_rating = (f32)w;

                            edges[edges_size++] = {u, v, edge_rating};
                        }
                    endfor
                }
            endfor
            METRICS_TIME(ep_compute_ratings)
            METRICS(level_time_compute_ratings.back() += get_seconds(sp_compute_ratings, ep_compute_ratings);)
            METRICS(level_edges.back() += edges_size;)
        }

        f32 solve_path_length_1(const graph_t &g,
                                const vertex_t u,
                                Matching &matching) {
            vertex_t uu = u;
            vertex_t vv = m_neighbors[u].n1;
            f32      w  = m_neighbors[u].w1;

            matching.add(uu, vv);

            return w;
        }

        f32 solve_path_length_2(const graph_t &g,
                                const vertex_t u,
                                Matching &matching) {
            vertex_t v1 = u;
            vertex_t v2 = m_neighbors[u].n1;
            vertex_t v3;
            f32      w1 = m_neighbors[u].w1;
            f32      w2;

            if (m_neighbors[v2].n1 == v1) {
                v3 = m_neighbors[v2].n2;
                w2 = m_neighbors[v2].w2;
            } else {
                v3 = m_neighbors[v2].n1;
                w2 = m_neighbors[v2].w1;
            }

            vertex_t uu, vv;
            f32      w;
            if (w1 > w2) {
                uu = v1;
                vv = v2;
                w  = w1;
            } else {
                uu = v2;
                vv = v3;
                w  = w2;
            }

            matching.add(uu, vv);
            return w;
        }

        f32 solve_path(const graph_t &g,
                       const vertex_t u,
                       const u32 length,
                       Matching &matching) {
            // special case of length 1
            if (length == 1) {
                return solve_path_length_1(g, u, matching);
            }

            // special case of length 2
            if (length == 2) {
                return solve_path_length_2(g, u, matching);
            }

            vertex_t v1, v2, v3;
            f32      w;

            s64 i = 0;

            // first edge of the path
            v1 = u;
            dp_edges[i] = v1;
            if (m_neighbors[v1].n1 == v1) {
                v2 = m_neighbors[u].n2;
                w  = m_neighbors[u].w2;
            } else {
                v2 = m_neighbors[u].n1;
                w  = m_neighbors[u].w1;
            }
            dp_edges[i + 1] = v2; // save edge

            // init dp
            dp_w[i]    = w;
            dp_m[i]    = -1;
            dp_take[i] = 1;
            i += 1;

            // second edge of the path
            if (m_neighbors[v2].n1 == v1) {
                v3 = m_neighbors[v2].n2;
                w  = m_neighbors[v2].w2;
            } else {
                v3 = m_neighbors[v2].n1;
                w  = m_neighbors[v2].w1;
            }
            dp_edges[i + 1] = v3; // save edge

            // init dp
            if (w > dp_w[i - 1]) {
                dp_w[i]    = w;
                dp_m[i]    = -1;
                dp_take[i] = 1;
            } else {
                dp_w[i]    = dp_w[i - 1];
                dp_m[i]    = 0;
                dp_take[i] = 0;
            }
            i += 1;

            // all other edges of the path
            v1 = v2;
            v2 = v3;
            while (m_neighbors[v2].n1 != v2 && m_neighbors[v2].n2 != v2) {
                if (m_neighbors[v2].n1 == v1) {
                    v3 = m_neighbors[v2].n2;
                    w  = m_neighbors[v2].w2;
                } else {
                    v3 = m_neighbors[v2].n1;
                    w  = m_neighbors[v2].w1;
                }
                dp_edges[i + 1] = v3; // save edge

                // dp
                if (w + dp_w[i - 2] > dp_w[i - 1]) {
                    dp_w[i]    = w + dp_w[i - 2];
                    dp_m[i]    = i - 2;
                    dp_take[i] = 1;
                } else {
                    dp_w[i]    = dp_w[i - 1];
                    dp_m[i]    = i - 1;
                    dp_take[i] = 0;
                }

                v1 = v2;
                v2 = v3;
                i += 1;
            }

            s64 idx = i - 1;
            while (idx != -1) {
                if (dp_take[idx]) {
                    vertex_t uu = dp_edges[idx];
                    vertex_t vv = dp_edges[idx + 1];

                    matching.add(uu, vv);
                }
                idx = dp_m[idx];
            }
            return dp_w[i - 1];
        }

        void solve_cycle(const graph_t &g,
                         const vertex_t u,
                         const u32 length,
                         Matching &matching) {
            vertex_t  n1          = m_neighbors[u].n1;
            vertex_t  n2          = m_neighbors[u].n2;
            Neighbors original_u  = m_neighbors[u];
            Neighbors original_n1 = m_neighbors[original_u.n1];
            Neighbors original_n2 = m_neighbors[original_u.n2];

            f32 matching_weight1 = 0.0;
            f32 matching_weight2 = 0.0;
            dp_cycle_matches1.clear();
            dp_cycle_matches2.clear();

            // cut connection between u and n1, n2 should point to self,
            m_neighbors[u].n1 = original_u.n2;
            m_neighbors[u].w1 = original_u.w2;
            m_neighbors[u].n2 = u;

            if (m_neighbors[n1].n1 == u) {
                m_neighbors[n1].n1 = original_n1.n2;
                m_neighbors[n1].w1 = original_n1.w2;
                m_neighbors[n1].n2 = n1;
            } else {
                m_neighbors[n1].n2 = n1;
            }
            matching_weight1 = solve_path(g, u, length - 1, dp_cycle_matches1);
            m_neighbors[u]  = original_u;
            m_neighbors[n1] = original_n1;

            // cut connection between u and n2
            m_neighbors[u].n2 = u;

            if (m_neighbors[n2].n1 == u) {
                m_neighbors[n2].n1 = original_n2.n2;
                m_neighbors[n2].w1 = original_n2.w2;
                m_neighbors[n2].n2 = n2;
            } else {
                m_neighbors[n2].n2 = n2;
            }
            matching_weight2 = solve_path(g, u, length - 1, dp_cycle_matches2);
            m_neighbors[u]  = original_u;
            m_neighbors[n2] = original_n2;

            Matching *dp_cycle_matches = &dp_cycle_matches1;
            if (matching_weight2 > matching_weight1) {
                dp_cycle_matches = &dp_cycle_matches2;
            }

            for (size_t i = 0; i < dp_cycle_matches->size(); ++i) {
                matching.add((*dp_cycle_matches)[i].u, (*dp_cycle_matches)[i].v);
            }
        }

        void random_matching(const size_t level,
                             const graph_t &g,
                             Matching &matching) {
            METRICS_TIME(sp)

            std::vector<u8> is_matched(g.get_n(), 0);

            forall_gu(g, u)
                {
                    if (is_matched[u]) { continue; }
                    weight_t u_w = g.weight(u);

                    forall_guiv(g, u, j, v)
                        {
                            if (is_matched[v]) { continue; }
                            weight_t v_w = g.weight(v);

                            if (u_w + v_w > m_l_max) { continue; }

                            is_matched[u] = 1;
                            is_matched[v] = 1;

                            matching.add(u, v);
                            break;
                        }
                    endfor
                }
            endfor

            METRICS_TIME(ep);
            METRICS(level_time_random.back() += get_seconds(sp, ep);)

#if ASSERT_ENABLED
            for (size_t i = 0; i < matching.size(); ++i) {
                const auto& [u, v] = matching[i];
                ASSERT(u != v);
            }
#endif

#if ASSERT_ENABLED
            std::vector<u8> hit(g.get_n(), 0);
            for (size_t i = 0; i < matching.size(); ++i) {
                const auto& [u, v] = matching[i];
                hit[u] += 1;
                hit[v] += 1;

                ASSERT(hit[u] == 1);
                ASSERT(hit[v] == 1);
            }
#endif
        }

        JSONString get_stats() override {
            std::string stats = "{ \n";
#if COLLECT_METRICS
            std::vector<f64> level_time(level_time_compute_ratings.size(), 0.0);
            for (size_t i = 0; i < level_time_compute_ratings.size(); ++i) {
                level_time[i] = level_time_compute_ratings[i] + level_time_sorting[i] + level_time_build_paths[i] + level_time_solve_paths[i] + level_time_solve_cycles[i] + level_time_random[i];
            }

            f64 global_time                 = sum<f64>(level_time);
            f64 global_time_compute_ratings = sum<f64>(level_time_compute_ratings);
            f64 global_time_sorting         = sum<f64>(level_time_sorting);
            f64 global_time_build_paths     = sum<f64>(level_time_build_paths);
            f64 global_time_solve_paths     = sum<f64>(level_time_solve_paths);
            f64 global_time_solve_cycles    = sum<f64>(level_time_solve_cycles);
            f64 global_time_random          = sum<f64>(level_time_random);

            u64 global_edges               = sum<u64>(level_edges);
            u64 global_edges_skipped       = sum<u64>(level_edges_skipped);
            u64 global_edges_form_cycle    = sum<u64>(level_edges_form_cycle);
            u64 global_edges_new_paths     = sum<u64>(level_edges_new_paths);
            u64 global_edges_enlarge_path  = sum<u64>(level_edges_enlarge_path);
            u64 global_edges_combine_paths = sum<u64>(level_edges_combine_paths);

            stats += to_JSON_MACRO(global_time);
            stats += to_JSON_MACRO(global_time_compute_ratings);
            stats += to_JSON_MACRO(global_time_sorting);
            stats += to_JSON_MACRO(global_time_build_paths);
            stats += to_JSON_MACRO(global_time_solve_paths);
            stats += to_JSON_MACRO(global_time_solve_cycles);
            stats += to_JSON_MACRO(global_time_random);
            stats += to_JSON_MACRO(global_edges);
            stats += to_JSON_MACRO(global_edges_skipped);
            stats += to_JSON_MACRO(global_edges_form_cycle);
            stats += to_JSON_MACRO(global_edges_new_paths);
            stats += to_JSON_MACRO(global_edges_enlarge_path);
            stats += to_JSON_MACRO(global_edges_combine_paths);
            stats += to_JSON_MACRO(level_time);
            stats += to_JSON_MACRO(level_time_compute_ratings);
            stats += to_JSON_MACRO(level_time_sorting);
            stats += to_JSON_MACRO(level_time_build_paths);
            stats += to_JSON_MACRO(level_time_solve_paths);
            stats += to_JSON_MACRO(level_time_solve_cycles);
            stats += to_JSON_MACRO(level_time_random);
            stats += to_JSON_MACRO(level_edges);
            stats += to_JSON_MACRO(level_edges_skipped);
            stats += to_JSON_MACRO(level_edges_form_cycle);
            stats += to_JSON_MACRO(level_edges_new_paths);
            stats += to_JSON_MACRO(level_edges_enlarge_path);
            stats += to_JSON_MACRO(level_edges_combine_paths);
#endif
            stats.pop_back();
            stats.pop_back();
            stats += "\n}";

            JSONString json_stats;
            json_stats.s = stats;
            return json_stats;
        }
    };
}

#endif //HEIPROMAP_GLOBAL_PATH_ALGORITHM_H
