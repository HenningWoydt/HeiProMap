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
#include <iomanip>
#include <queue>
#include <vector>

#include "../../definitions.h"
#include "../../macros.h"
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

    inline bool is_endpoint_fast(const Neighbors &n, const vertex_t u) { return n.n2 == u; }

    inline bool is_endpoint(const Neighbors &n, const vertex_t u) { return is_endpoint_fast(n, u); }

    inline bool is_one_endpoint_fast(const Neighbors &n, const vertex_t u) { return n.n1 != u && n.n2 == u; }

    inline bool is_one_endpoint(const Neighbors &n, const vertex_t u) { return is_one_endpoint_fast(n, u); }

    inline bool is_unmatched_fast(const Neighbors &n, const vertex_t u) { return n.n1 == u; }

    inline bool is_unmatched(const Neighbors &n, const vertex_t u) { return is_unmatched_fast(n, u); }

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

        Neighbors *m_neighbors = nullptr;
        u32       *path_id     = nullptr;
        u32       *path_length = nullptr;

        EdgeUVW *edges     = nullptr;
        size_t  edges_size = 0;

        // for DP
        f32      *dp_w     = nullptr;
        s64      *dp_m     = nullptr;
        u8       *dp_take  = nullptr;
        vertex_t *dp_edges = nullptr;

        Matching dp_cycle_matches1;
        Matching dp_cycle_matches2;

        METRICS(f64 global_time_compute_ratings = 0.0;)
        METRICS(f64 global_time_sorting         = 0.0;)
        METRICS(f64 global_time_build_paths     = 0.0;)
        METRICS(f64 global_time_solve_paths     = 0.0;)
        METRICS(f64 global_time_solve_cycle     = 0.0;)
        METRICS(u64 global_edges_skipped        = 0;)
        METRICS(u64 global_edges_form_cycle     = 0;)
        METRICS(u64 global_edges_new_paths      = 0;)
        METRICS(u64 global_edges_enlarge_path   = 0;)
        METRICS(u64 global_edges_combine_paths  = 0;)

    public:
        GlobalPathAlgorithmMatcher() = default;

        ~GlobalPathAlgorithmMatcher() override {
            free(m_neighbors);
            free(path_id);
            free(path_length);

            free(edges);

            free(dp_w);
            free(dp_m);
            free(dp_take);
            free(dp_edges);
        }

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

            m_neighbors = (Neighbors *) aligned_alloc(64, t_n_64 * sizeof(Neighbors));
            path_id     = (u32 *) aligned_alloc(64, t_n_64 * sizeof(u32));
            path_length = (u32 *) aligned_alloc(64, t_n_64 * sizeof(u32));

            edges = (EdgeUVW *) aligned_alloc(64, t_m_64 * sizeof(EdgeUVW));

            dp_w     = (f32 *) aligned_alloc(64, t_n_64 * sizeof(f32));
            dp_m     = (s64 *) aligned_alloc(64, t_n_64 * sizeof(s64));
            dp_take  = (u8 *) aligned_alloc(64, t_n_64 * sizeof(u8));
            dp_edges = (vertex_t *) aligned_alloc(64, t_n_64 * sizeof(vertex_t));

            dp_cycle_matches1.initialize(t_n_64);
            dp_cycle_matches2.initialize(t_n_64);
        }

        void match(const size_t level,
                   const graph_t &g,
                   Matching &matching) override {
            if (level < config->random_level) {
                // use a random matching
                random_matching(level, g, matching);
                return;
            }
            METRICS(f64 time_compute_ratings = 0.0;)
            METRICS(f64 time_sorting         = 0.0;)
            METRICS(f64 time_build_paths     = 0.0;)
            METRICS(f64 time_solve_paths     = 0.0;)
            METRICS(f64 time_solve_cycle     = 0.0;)
            METRICS(u64 edges_skipped        = 0;)
            METRICS(u64 edges_form_cycle     = 0;)
            METRICS(u64 edges_new_paths      = 0;)
            METRICS(u64 edges_enlarge_path   = 0;)
            METRICS(u64 edges_combine_paths  = 0;)

            METRICS_TIME(sp_compute_ratings);
            for (vertex_t u = 0; u < m_n; ++u) {
                m_neighbors[u].n1 = u;
                m_neighbors[u].n2 = u;
            }

            edges_size = 0;
            f32 max_rating = -std::numeric_limits<f32>::max();
            f32 min_rating = std::numeric_limits<f32>::max();
            forall_gu(g, u)
                {
                    weight_t u_w = g.weight(u);

                    forall_guivw(g, u, j, v, w)
                        {
                            if (u > v) { continue; }
                            weight_t v_w = g.weight(v);

                            // if (u_w > 1.5*av_manager.get_n_active() / 20.0 * m_k) { continue; }
                            // if (v_w > 1.5*av_manager.get_n_active() / 20.0 * m_k) { continue; }

                            if (u_w + v_w > m_l_max) { continue; }

                            f32 edge_rating;

                            // edge_rating = ((f32) w) / (g.size(u) * g.size(v));
                            // edge_rating = (f32) w / (f32) (u_w * v_w);
                            // edge_rating = (f32) w / (f32) (u_w * v_w);
                            edge_rating = ((f32) (w * w)) / ((f32) (u_w * v_w));
                            // edge_rating = ((f32) (w * w)) / ((f32) (u_w + v_w));
                            // edge_rating = ((f32) (w * w * w)) / ((f32) (u_w * v_w));
                            // edge_rating = (f32) w / (f32) (u_w * v_w * u_w * v_w);
                            // edge_rating = (f32)w;

                            max_rating = std::max(max_rating, edge_rating);
                            min_rating = std::min(min_rating, edge_rating);

                            edges[edges_size++] = {u, v, edge_rating};
                        }
                    endfor
                }
            endfor
            METRICS_TIME(ep_compute_ratings);

            METRICS_TIME(sp_sorting);
            if (max_rating != min_rating) {
                std::sort(edges, edges + edges_size, std::greater<>());
            }
            METRICS_TIME(ep_sorting);

            METRICS_TIME(sp_build_paths);
            for (size_t i = 0; i < edges_size; ++i) {
                auto [u, v, w] = edges[i];

                if (!is_endpoint(m_neighbors[u], u) || !is_endpoint(m_neighbors[v], v)) {
                    // u or v is not an endpoint
                    METRICS(edges_skipped += 1);
                    continue;
                }

                bool u_unmatched = is_unmatched(m_neighbors[u], u);
                bool v_unmatched = is_unmatched(m_neighbors[v], v);

                if (u_unmatched && v_unmatched) {
                    // both are unmatched, only one new path of length 1
                    m_neighbors[u].n1 = v;
                    m_neighbors[u].w1 = w;
                    m_neighbors[v].n1 = u;
                    m_neighbors[v].w1 = w;
                    path_id[u]     = u;
                    path_id[v]     = u;
                    path_length[u] = 1;
                    METRICS(edges_new_paths += 1);
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
                    METRICS(edges_enlarge_path += 1);
                    continue;
                }

                // cycle
                if (u_id == v_id) {
                    if (path_length[u_id] % 2 == 1) {
                        // same path and odd length size, close the cycle
                        path_length[u_id] += 1; // increase path length

                        // for u set v as a neighbor
                        m_neighbors[u].n2 = v;
                        m_neighbors[u].w2 = w;

                        // for v set u as a neighbor
                        m_neighbors[v].n2 = u;
                        m_neighbors[v].w2 = w;

                        // solve the cycle
                        METRICS_TIME(sp_solve_cycle);
                        solve_cycle(g, u, path_length[u_id], matching);
                        METRICS_TIME(ep_solve_cycle);
                        path_length[u_id] = 0;

                        METRICS(f64 t_solve = get_seconds(sp_solve_cycle, ep_solve_cycle);)
                        METRICS(time_solve_cycle += t_solve;)
                        METRICS(time_build_paths -= t_solve;)
                        METRICS(edges_form_cycle += 1);
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

                METRICS(edges_combine_paths += 1);
            }
            METRICS_TIME(ep_build_paths)

            // process all paths
            METRICS_TIME(sp_solve_paths)
            forall_gu(g, u)
                {
                    if (is_one_endpoint(m_neighbors[u], u) && path_length[path_id[u]] > 0) {
                        solve_path(g, u, path_length[path_id[u]], matching);
                        path_length[path_id[u]] = 0;
                    }
                }
            endfor
            METRICS_TIME(ep_solve_paths)

#if COLLECT_METRICS
            time_compute_ratings += get_seconds(sp_compute_ratings, ep_compute_ratings);
            time_build_paths += get_seconds(sp_build_paths, ep_build_paths);
            time_sorting += get_seconds(sp_sorting, ep_sorting);
            time_solve_paths += get_seconds(sp_solve_paths, ep_solve_paths);
            f64 time_total = time_compute_ratings + time_sorting + time_build_paths + time_solve_paths + time_solve_cycle;

            global_time_compute_ratings += time_compute_ratings;
            global_time_sorting += time_sorting;
            global_time_build_paths += time_build_paths;
            global_time_solve_paths += time_solve_paths;
            global_time_solve_cycle += time_solve_cycle;

            global_edges_skipped += edges_skipped;
            global_edges_form_cycle += edges_form_cycle;
            global_edges_new_paths += edges_new_paths;
            global_edges_enlarge_path += edges_enlarge_path;
            global_edges_combine_paths += edges_combine_paths;
            u64 total_edges = edges_skipped + edges_form_cycle + edges_new_paths + edges_enlarge_path + edges_combine_paths;

            std::string stats = "{ \n";
            stats += to_JSON_MACRO(time_total);
            stats += to_JSON_MACRO(time_compute_ratings);
            stats += to_JSON_MACRO(time_sorting);
            stats += to_JSON_MACRO(time_build_paths);
            stats += to_JSON_MACRO(time_solve_cycle);
            stats += to_JSON_MACRO(time_solve_paths);
            stats += to_JSON_MACRO(total_edges);
            stats += to_JSON_MACRO(edges_skipped);
            stats += to_JSON_MACRO(edges_new_paths);
            stats += to_JSON_MACRO(edges_enlarge_path);
            stats += to_JSON_MACRO(edges_combine_paths);
            stats += to_JSON_MACRO(edges_form_cycle);

            stats.pop_back();
            stats.pop_back();
            stats += "\n}";
            JSONString json_stats;
            json_stats.s = stats;
            m_stat_collector->add_matching_method_stats(level, json_stats);
#endif

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

#if COLLECT_METRICS
            f64 time          = get_seconds(sp, ep);
            std::string stats = "{ \n";
            stats += to_JSON_MACRO(time);

            stats.pop_back();
            stats.pop_back();
            stats += "\n}";
            JSONString json_stats;
            json_stats.s = stats;
            m_stat_collector->add_matching_method_stats(level, json_stats);
#endif

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
    };
}

#endif //HEIPROMAP_GLOBAL_PATH_ALGORITHM_H
