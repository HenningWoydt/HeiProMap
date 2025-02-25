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

#ifndef HEIPROMAP_GLOBAL_PATH_ALGORITHM_ARRAYS_H
#define HEIPROMAP_GLOBAL_PATH_ALGORITHM_ARRAYS_H

#include <algorithm>
#include <iomanip>
#include <numeric>
#include <queue>
#include <vector>
#include <chrono>

#include "../../definitions.h"
#include "../../macros.h"
#include "../interfaces/ISerialMatcher.h"
#include "global_path_algorithm.h"
#include "../utility/utils.h"

namespace HeiProMap {
#ifndef COLLECT_METRICS
#define COLLECT_METRICS false
#endif

#if (COLLECT_METRICS)
#define TIME_POINT(x) auto x = std::chrono::high_resolution_clock::now()
#else
#define TIME_POINT(x) ((void)0)
#endif

    /**
     * Computes a matching based on the Global Path Algorithm from
     * > Jens Maue and Peter Sanders.
     * > Engineering Algorithms for Approximate Weighted Matching.
     * > Experimental Algorithms, 6th International Workshop, WEA 2007, Rome, Italy, June 6-8, 2007, Proceedings.
     */
    class GlobalPathAlgorithmArraysMatcher final : public ISerialMatcher {
        vertex_t    m_n     = 0;
        vertex_t    m_m     = 0;
        partition_t m_k     = 0;
        weight_t    m_l_max = 0;

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

        std::vector<EdgeUV> dp_cycle_matches1;
        std::vector<EdgeUV> dp_cycle_matches2;

#if COLLECT_METRICS
        f64 global_time_compute_ratings = 0.0;
        f64 global_time_sorting         = 0.0;
        f64 global_time_build_paths     = 0.0;
        f64 global_time_solve_paths     = 0.0;
        f64 global_time_solve_cycle     = 0.0;
        u64 global_edges_skipped        = 0;
        u64 global_edges_form_cycle     = 0;
        u64 global_edges_new_paths      = 0;
        u64 global_edges_enlarge_path   = 0;
        u64 global_edges_combine_paths  = 0;
#endif

    public:
        GlobalPathAlgorithmArraysMatcher() = default;

        ~GlobalPathAlgorithmArraysMatcher() override {
            free(m_neighbors);
            free(path_id);
            free(path_length);

            free(edges);

            free(dp_w);
            free(dp_m);
            free(dp_take);
            free(dp_edges);
        }

        void initialize(const vertex_t t_n, const vertex_t t_m, const partition_t t_k, const weight_t t_l_max) override {
            vertex_t t_n_64 = round_up_64(t_n);
            vertex_t t_m_64 = round_up_64(t_m);

            m_n     = t_n;
            m_m     = t_m;
            m_k     = t_k;
            m_l_max = t_l_max;

            m_neighbors = (Neighbors *) aligned_alloc(64, t_n_64 * sizeof(Neighbors));
            path_id     = (u32 *) aligned_alloc(64, t_n_64 * sizeof(u32));
            path_length = (u32 *) aligned_alloc(64, t_n_64 * sizeof(u32));

            edges = (EdgeUVW *) aligned_alloc(64, t_m_64 * sizeof(EdgeUVW));

            dp_w     = (f32 *) aligned_alloc(64, t_n_64 * sizeof(f32));
            dp_m     = (s64 *) aligned_alloc(64, t_n_64 * sizeof(s64));
            dp_take  = (u8 *) aligned_alloc(64, t_n_64 * sizeof(u8));
            dp_edges = (vertex_t *) aligned_alloc(64, t_n_64 * sizeof(vertex_t));
        }

        template<typename TSerialGraph, typename TSerialActiveVertexManager>
        void match(GlobalPathAlgorithmConfiguration &config,
                   TSerialGraph &g,
                   TSerialActiveVertexManager &av_manager,
                   std::vector<EdgeUV> &matches) {
            for (vertex_t u = 0; u < m_n; ++u) {
                m_neighbors[u].n1 = u;
                m_neighbors[u].n2 = u;
            }

#if COLLECT_METRICS
            f64 time_compute_ratings = 0.0;
            f64 time_sorting         = 0.0;
            f64 time_build_paths     = 0.0;
            f64 time_solve_paths     = 0.0;
            f64 time_solve_cycle     = 0.0;
            u64 edges_skipped        = 0;
            u64 edges_form_cycle     = 0;
            u64 edges_new_paths      = 0;
            u64 edges_enlarge_path   = 0;
            u64 edges_combine_paths  = 0;
#endif

            TIME_POINT(sp_compute_ratings);

            edges_size = 0;

            for (vertex_t u: av_manager) {
                weight_t u_w = g.get_weight(u);

                for (const auto &[v, w]: g[u]) {
                    weight_t v_w = g.get_weight(v);

                    if (u_w + v_w > m_l_max) { continue; }
                    if (v < u) { continue; }

                    // f32 edge_rating = (((f64) w) * ((f64) w)) / (g.get_weight(u) * g.get_weight(v));
                    // f32 edge_rating = ((f64) w) / (g.size(u) * g.size(v));
                    f32 edge_rating = (f32) w / (f32) (u_w * v_w);

                    if (edge_rating == 0.0) { continue; }

                    edges[edges_size] = {u, v, edge_rating};
                    edges_size += 1;
                }
            }
            TIME_POINT(ep_compute_ratings);

            TIME_POINT(sp_sorting);
            std::sort(edges, edges + edges_size, std::greater<>());
            TIME_POINT(ep_sorting);

            TIME_POINT(sp_build_paths);
            for (size_t   i = 0; i < edges_size; ++i) {
                auto &[u, v, w] = edges[i];

                if (!is_endpoint(m_neighbors[u], u) || !is_endpoint(m_neighbors[v], v)) {
                    // u or v is not an endpoint
                    // edges_skipped++;
                    continue;
                }

                bool u_unmatched = is_unmatched(m_neighbors[u], u);
                bool v_unmatched = is_unmatched(m_neighbors[v], v);
                u32  u_id        = path_id[u];
                u32  v_id        = path_id[v];

                if (u_unmatched && v_unmatched) {
                    // both unmatched, one new path of length 1
                    m_neighbors[u].n1 = v;
                    m_neighbors[u].w1 = w;
                    m_neighbors[v].n1 = u;
                    m_neighbors[v].w1 = w;
                    path_id[u]     = u;
                    path_id[v]     = u;
                    path_length[u] = 1;
                    // edges_new_paths++;
                    continue;
                }
                if (u_unmatched && !v_unmatched) {
                    // only one unmatched, enlarge path
                    m_neighbors[u].n1 = v;
                    m_neighbors[u].w1 = w;
                    m_neighbors[v].n2 = u;
                    m_neighbors[v].w2 = w;
                    path_id[u] = v_id;
                    path_length[v_id] += 1;
                    // edges_enlarge_path++;
                    continue;
                }
                if (!u_unmatched && v_unmatched) {
                    // only one unmatched, enlarge path
                    m_neighbors[u].n2 = v;
                    m_neighbors[u].w2 = w;
                    m_neighbors[v].n1 = u;
                    m_neighbors[v].w1 = w;
                    path_id[v] = u_id;
                    path_length[u_id] += 1;
                    // edges_enlarge_path++;
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
                        TIME_POINT(sp_solve_cycle);
                        solve_cycle(g, u, matches);
                        TIME_POINT(ep_solve_cycle);

#if COLLECT_METRICS
                        f64 t_solve = get_seconds(sp_solve_cycle, ep_solve_cycle);
                        time_solve_cycle += t_solve;
                        time_build_paths -= t_solve;
                        edges_form_cycle++;
#endif
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

                vertex_t last_vertex;
                vertex_t curr_vertex;
                u32      id1;
                u32      id2;
                if (path_length[u_id] <= path_length[v_id]) {
                    last_vertex = v;
                    curr_vertex = u;
                    id1         = v_id;
                    id2         = u_id;
                } else {
                    last_vertex = u;
                    curr_vertex = v;
                    id1         = u_id;
                    id2         = v_id;
                }

                path_length[id1] += 1 + path_length[id2];
                while (m_neighbors[curr_vertex].n1 != curr_vertex && m_neighbors[curr_vertex].n2 != curr_vertex) {
                    vertex_t temp_last_vertex = last_vertex;
                    last_vertex = curr_vertex;
                    curr_vertex = m_neighbors[curr_vertex].n1 == temp_last_vertex ? m_neighbors[curr_vertex].n2 : m_neighbors[curr_vertex].n1;
                }
                path_id[curr_vertex] = id1;

                // edges_combine_paths++;

            }
            TIME_POINT(ep_build_paths);

            TIME_POINT(sp_solve_paths);
            // process all paths
            for (vertex_t u: av_manager) {
                if (is_one_endpoint(m_neighbors[u], u)) {
                    solve_path(g, u, matches);
                }
            }
            TIME_POINT(ep_solve_paths);

#if COLLECT_METRICS
            time_compute_ratings += get_seconds(sp_compute_ratings, ep_compute_ratings);
            time_sorting += get_seconds(sp_sorting, ep_sorting);
            time_build_paths += get_seconds(sp_build_paths, ep_build_paths);

            time_solve_paths += get_seconds(sp_solve_paths, ep_solve_paths);

            global_time_compute_ratings += time_compute_ratings;
            global_time_sorting += time_sorting;
            global_time_build_paths += time_build_paths;
            global_time_solve_paths += time_solve_paths;
            global_time_solve_cycle += time_solve_cycle;
            f64 global_time_total = global_time_compute_ratings + global_time_sorting + global_time_build_paths + global_time_solve_paths + global_time_solve_cycle;

            f64 time_total = time_compute_ratings + time_sorting + time_build_paths + time_solve_paths + time_solve_cycle;
            std::cout << std::fixed << std::setprecision(4);
            std::cout << "time compute ratings: " << std::setprecision(4) << time_compute_ratings << " - " << time_compute_ratings / time_total << " -- global time compute ratings: " << global_time_compute_ratings << " - " << global_time_compute_ratings / global_time_total << std::endl;
            std::cout << "time sorting        : " << std::setprecision(4) << time_sorting << " - " << time_sorting / time_total << " -- global time sorting        : " << global_time_sorting << " - " << global_time_sorting / global_time_total << std::endl;
            std::cout << "time build paths    : " << std::setprecision(4) << time_build_paths << " - " << time_build_paths / time_total << " -- global time build paths    : " << global_time_build_paths << " - " << global_time_build_paths / global_time_total << std::endl;
            std::cout << "time solve paths    : " << std::setprecision(4) << time_solve_paths << " - " << time_solve_paths / time_total << " -- global time solve paths    : " << global_time_solve_paths << " - " << global_time_solve_paths / global_time_total << std::endl;
            std::cout << "time solve cycle    : " << std::setprecision(4) << time_solve_paths << " - " << time_solve_cycle / time_total << " -- global time solve cycle    : " << global_time_solve_cycle << " - " << global_time_solve_cycle / global_time_total << std::endl;
            std::cout << "time total          : " << std::setprecision(4) << time_total << " - " << time_total / time_total << " -- global time total          : " << global_time_total << " - " << global_time_total / global_time_total << std::endl;

            global_edges_skipped += edges_skipped;
            global_edges_form_cycle += edges_form_cycle;
            global_edges_new_paths += edges_new_paths;
            global_edges_enlarge_path += edges_enlarge_path;
            global_edges_combine_paths += edges_combine_paths;
            u64 global_total_edges = global_edges_skipped + global_edges_form_cycle + global_edges_new_paths + global_edges_enlarge_path + global_edges_combine_paths;
            u64 local_total_edges  = edges_skipped + edges_form_cycle + edges_new_paths + edges_enlarge_path + edges_combine_paths;
            std::cout << "edges skipped       : " << edges_skipped << "\t - " << (f64) edges_skipped / (f64) local_total_edges << " -- global edges skipped       : " << global_edges_skipped << " - " << (f64) global_edges_skipped / (f64) global_total_edges << std::endl;
            std::cout << "edges form cycle    : " << edges_form_cycle << "\t - " << (f64) edges_form_cycle / (f64) local_total_edges << " -- global edges form cycle    : " << global_edges_form_cycle << " - " << (f64) global_edges_form_cycle / (f64) global_total_edges << std::endl;
            std::cout << "edges new paths     : " << edges_new_paths << "\t - " << (f64) edges_new_paths / (f64) local_total_edges << " -- global edges new paths     : " << global_edges_new_paths << " - " << (f64) global_edges_new_paths / (f64) global_total_edges << std::endl;
            std::cout << "edges enlarge path  : " << edges_enlarge_path << "\t - " << (f64) edges_enlarge_path / (f64) local_total_edges << " -- global edges enlarge path  : " << global_edges_enlarge_path << " - " << (f64) global_edges_enlarge_path / (f64) global_total_edges << std::endl;
            std::cout << "edges combine paths : " << edges_combine_paths << "\t - " << (f64) edges_combine_paths / (f64) local_total_edges << " -- global edges combine paths : " << global_edges_combine_paths << " - " << (f64) global_edges_combine_paths / (f64) global_total_edges << std::endl;
#endif

#if ASSERT_ENABLED
            for (const EdgeUV &e: matches) {
                ASSERT(e.u != e.v);
                ASSERT(av_manager.is_active(e.u));
                ASSERT(av_manager.is_active(e.v));
            }
#endif

#if ASSERT_ENABLED
            std::vector<u8> hit(g.get_n(), 0);
            for (const auto &e: matches) {
                hit[e.u] += 1;
                hit[e.v] += 1;

                ASSERT(hit[e.u] != 2);
                ASSERT(hit[e.v] != 2);
            }
#endif
        }

        template<typename TSerialGraph>
        f32 solve_path_length_1(TSerialGraph &g,
                                const vertex_t u,
                                std::vector<EdgeUV> &matches) {
            vertex_t uu = u;
            vertex_t vv = m_neighbors[u].n1;
            f32      w  = m_neighbors[u].w1;

            if (g.size(uu) < g.size(vv)) {
                matches.emplace_back(uu, vv);
            } else {
                matches.emplace_back(vv, uu);
            }

            // destroy the endpoints
            m_neighbors[uu].n1 = -1;
            m_neighbors[uu].n2 = -1;
            m_neighbors[vv].n1 = -1;
            m_neighbors[vv].n2 = -1;
            return w;
        }

        template<typename TSerialGraph>
        f32 solve_path_length_2(TSerialGraph &g,
                                const vertex_t u,
                                std::vector<EdgeUV> &matches) {
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

            if (g.size(uu) < g.size(vv)) {
                matches.emplace_back(uu, vv);
            } else {
                matches.emplace_back(vv, uu);
            }

            // destroy the endpoints
            m_neighbors[v1].n1 = -1;
            m_neighbors[v1].n2 = -1;
            m_neighbors[v3].n1 = -1;
            m_neighbors[v3].n2 = -1;
            return w;
        }

        template<typename TSerialGraph>
        f32 solve_path(TSerialGraph &g,
                       const vertex_t u,
                       std::vector<EdgeUV> &matches) {
            // return solve_path_new(g, u, matches);

            // special case of length 1
            if (path_length[path_id[u]] == 1) {
                return solve_path_length_1(g, u, matches);
            }

            // special case of length 2
            if (path_length[path_id[u]] == 2) {
                return solve_path_length_2(g, u, matches);
            }

            vertex_t endpoint_1 = u;
            vertex_t endpoint_2;
            vertex_t v1, v2, v3;
            f32      w1, w2;

            s64 i = 0;

            // first edge of the path
            v1 = u;
            dp_edges[i] = v1;
            if (m_neighbors[v1].n1 == v1) {
                v2 = m_neighbors[u].n2;
                w1 = m_neighbors[u].w2;
            } else {
                v2 = m_neighbors[u].n1;
                w1 = m_neighbors[u].w1;
            }
            dp_edges[i + 1] = v2; // save edge

            // init dp
            dp_w[i]    = w1;
            dp_m[i]    = -1;
            dp_take[i] = 1;
            i += 1;

            // second edge of the path
            if (m_neighbors[v2].n1 == v1) {
                v3 = m_neighbors[v2].n2;
                w2 = m_neighbors[v2].w2;
            } else {
                v3 = m_neighbors[v2].n1;
                w2 = m_neighbors[v2].w1;
            }
            dp_edges[i + 1] = v3; // save edge

            // init dp
            if (w2 > dp_w[i - 1]) {
                dp_w[i]    = w2;
                dp_m[i]    = -1;
                dp_take[i] = 1;
            } else {
                dp_w[i]    = dp_w[i - 1];
                dp_m[i]    = 0;
                dp_take[i] = 0;
            }
            i += 1;

            // all other edges of the path
            v1         = v2;
            v2         = v3;
            w1         = w2;
            while (m_neighbors[v2].n1 != v2 && m_neighbors[v2].n2 != v2) {
                if (m_neighbors[v2].n1 == v1) {
                    v3 = m_neighbors[v2].n2;
                    w2 = m_neighbors[v2].w2;
                } else {
                    v3 = m_neighbors[v2].n1;
                    w2 = m_neighbors[v2].w1;
                }
                dp_edges[i + 1] = v3; // save edge

                // dp
                if (w2 + dp_w[i - 2] > dp_w[i - 1]) {
                    dp_w[i]    = w2 + dp_w[i - 2];
                    dp_m[i]    = i - 2;
                    dp_take[i] = 1;
                } else {
                    dp_w[i]    = dp_w[i - 1];
                    dp_m[i]    = i - 1;
                    dp_take[i] = 0;
                }

                v1 = v2;
                v2 = v3;
                w1 = w2;
                i += 1;
            }
            endpoint_2 = v2;

            s64 idx = i - 1;
            while (idx != -1) {
                if (dp_take[idx]) {

                    vertex_t uu = dp_edges[idx];
                    vertex_t vv = dp_edges[idx + 1];

                    if (g.size(uu) < g.size(vv)) {
                        matches.emplace_back(uu, vv);
                    } else {
                        matches.emplace_back(vv, uu);
                    }
                }
                idx = dp_m[idx];
            }

            // destroy the endpoints
            m_neighbors[endpoint_1].n1 = -1;
            m_neighbors[endpoint_1].n2 = -1;
            m_neighbors[endpoint_2].n1 = -1;
            m_neighbors[endpoint_2].n2 = -1;
            return dp_w[i - 1];
        }

        template<typename TSerialGraph>
        void solve_cycle(TSerialGraph &g,
                         vertex_t u,
                         std::vector<EdgeUV> &matches) {
            vertex_t left  = m_neighbors[u].n1;
            vertex_t right = m_neighbors[u].n2;

            vertex_t left_left   = m_neighbors[left].n1;
            vertex_t left_right  = m_neighbors[left].n2;
            vertex_t right_left  = m_neighbors[right].n1;
            vertex_t right_right = m_neighbors[right].n2;

            f32 matching_weight1 = 0.0;
            f32 matching_weight2 = 0.0;
            dp_cycle_matches1.clear();
            dp_cycle_matches2.clear();

            m_neighbors[u].n1 = u;
            if (m_neighbors[left].n1 == u) { m_neighbors[left].n1 = left; } else { m_neighbors[left].n2 = left; }
            path_length[path_id[u]] -= 1;
            matching_weight1 = solve_path(g, u, dp_cycle_matches1);
            path_length[path_id[u]] += 1;
            m_neighbors[u].n1    = left;
            m_neighbors[u].n2    = right;
            m_neighbors[left].n1 = left_left;
            m_neighbors[left].n2 = left_right;

            m_neighbors[u].n1 = u;
            if (m_neighbors[right].n2 == u) { m_neighbors[right].n1 = right; } else { m_neighbors[right].n2 = right; }
            path_length[path_id[u]] -= 1;
            matching_weight2 = solve_path(g, u, dp_cycle_matches2);
            path_length[path_id[u]] += 1;
            m_neighbors[u].n1     = left;
            m_neighbors[u].n2     = right;
            m_neighbors[right].n1 = right_left;
            m_neighbors[right].n2 = right_right;

            if (matching_weight1 >= matching_weight2) {
                for (auto &e: dp_cycle_matches1) {
                    matches.emplace_back(e);
                }
            } else {
                for (auto &e: dp_cycle_matches2) {
                    matches.emplace_back(e);
                }
            }
        }
    };
}

#endif //HEIPROMAP_GLOBAL_PATH_ALGORITHM_ARRAYS_H
