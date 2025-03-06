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
#include <chrono>
#include <iomanip>
#include <numeric>
#include <queue>
#include <random>
#include <vector>

#include <boost/sort/sort.hpp>

#include "global_path_algorithm.h"
#include "../../definitions.h"
#include "../../macros.h"
#include "../interfaces/ISerialMatcher.h"
#include "../utility/utils.h"

namespace HeiProMap {
    /**
     * Computes a matching based on the Global Path Algorithm from
     * > Jens Maue and Peter Sanders.
     * > Engineering Algorithms for Approximate Weighted Matching.
     * > Experimental Algorithms, 6th International Workshop, WEA 2007, Rome, Italy, June 6-8, 2007, Proceedings.
     */
    class GlobalPathAlgorithmArraysMatcher final : public ISerialMatcher {
        vertex_t m_n     = 0;
        vertex_t m_m     = 0;
        partition_t m_k  = 0;
        weight_t m_l_max = 0;

        Neighbors* m_neighbors = nullptr;
        u32* path_id           = nullptr;
        u32* path_length       = nullptr;

        EdgeUVW* edges    = nullptr;
        size_t edges_size = 0;

        // for DP
        f32* dp_w          = nullptr;
        s64* dp_m          = nullptr;
        u8* dp_take        = nullptr;
        vertex_t* dp_edges = nullptr;

        EdgeUV* dp_cycle_matches1 = nullptr;
        EdgeUV* dp_cycle_matches2 = nullptr;

        std::mt19937 gen;
        std::uniform_real_distribution<float> dis;

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

            free(dp_cycle_matches1);
            free(dp_cycle_matches2);
        }

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const weight_t t_l_max,
                        const u64 t_seed) override {
            vertex_t t_n_64 = round_up_64(t_n);
            vertex_t t_m_64 = round_up_64(t_m);

            m_n     = t_n;
            m_m     = t_m;
            m_k     = t_k;
            m_l_max = t_l_max;

            m_neighbors = (Neighbors*)aligned_alloc(64, t_n_64 * sizeof(Neighbors));
            path_id     = (u32*)aligned_alloc(64, t_n_64 * sizeof(u32));
            path_length = (u32*)aligned_alloc(64, t_n_64 * sizeof(u32));

            edges = (EdgeUVW*)aligned_alloc(64, t_m_64 * sizeof(EdgeUVW));

            dp_w     = (f32*)aligned_alloc(64, t_n_64 * sizeof(f32));
            dp_m     = (s64*)aligned_alloc(64, t_n_64 * sizeof(s64));
            dp_take  = (u8*)aligned_alloc(64, t_n_64 * sizeof(u8));
            dp_edges = (vertex_t*)aligned_alloc(64, t_n_64 * sizeof(vertex_t));

            dp_cycle_matches1 = (EdgeUV*)aligned_alloc(64, t_n_64 * sizeof(EdgeUV));
            dp_cycle_matches2 = (EdgeUV*)aligned_alloc(64, t_n_64 * sizeof(EdgeUV));

            gen.seed(t_seed);
            dis = std::uniform_real_distribution<float>(0.0f, 1.0f);
        }

        template <typename TSerialGraph, typename TSerialActiveVertexManager>
        void match(size_t level,
                   GlobalPathAlgorithmConfiguration& config,
                   TSerialGraph& g,
                   TSerialActiveVertexManager& av_manager,
                   EdgeUV* matches,
                   size_t& matches_size) {
            if (level < config.random_level) {
                // use a random matching
                random_matching(g, av_manager, matches, matches_size);
                return;
            }

            matches      = ASSUME_ALIGNED(EdgeUV*, matches, 64);
            matches_size = 0;

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

            for (vertex_t u = 0; u < m_n; ++u) {
                m_neighbors[u].n1 = u;
                m_neighbors[u].n2 = u;
            }

            edges_size     = 0;
            f32 max_rating = -std::numeric_limits<f32>::max();
            f32 min_rating = std::numeric_limits<f32>::max();
            for (vertex_t u : av_manager) {
                weight_t u_w = g.get_weight(u);

                for (const auto [v, w] : g[u]) {
                    if (u > v) { continue; }
                    weight_t v_w = g.get_weight(v);

                    // if (u_w > 1.5*av_manager.get_n_active() / 20.0 * m_k) { continue; }
                    // if (v_w > 1.5*av_manager.get_n_active() / 20.0 * m_k) { continue; }

                    if (u_w + v_w > m_l_max) { continue; }

                    f32 edge_rating;

                    // edge_rating = ((f32) w) / (g.size(u) * g.size(v));
                    // edge_rating = (f32) w / (f32) (u_w * v_w);
                    // edge_rating = (f32) w / (f32) (u_w * v_w);
                    edge_rating = ((f32)(w * w)) / ((f32)(u_w * v_w));
                    // edge_rating = ((f32) (w * w)) / ((f32) (u_w + v_w));
                    // edge_rating = ((f32) (w * w * w)) / ((f32) (u_w * v_w));
                    // edge_rating = (f32) w / (f32) (u_w * v_w * u_w * v_w);
                    // edge_rating = (f32)w;

                    max_rating = std::max(max_rating, edge_rating);
                    min_rating = std::min(min_rating, edge_rating);

                    edges[edges_size++] = {u, v, edge_rating};
                }
            }
            TIME_POINT(ep_compute_ratings);

            TIME_POINT(sp_sorting);
            if (max_rating != min_rating) {
                // boost::sort::pdqsort(edges, edges + edges_size, std::greater<>());
                // boost::sort::spinsort(edges, edges + edges_size, std::greater<>());
                // boost::sort::flat_stable_sort(edges, edges + edges_size, std::greater<>());
                std::sort(edges, edges + edges_size, std::greater<>());
            }
            TIME_POINT(ep_sorting);

            TIME_POINT(sp_build_paths);

            for (size_t i = 0; i < edges_size; ++i) {
                auto [u, v, w] = edges[i];

                if (!is_endpoint(m_neighbors[u], u) || !is_endpoint(m_neighbors[v], v)) {
                    // u or v is not an endpoint
                    INCREASE_COUNTER(edges_skipped);
                    continue;
                }

                bool u_unmatched = is_unmatched(m_neighbors[u], u);
                bool v_unmatched = is_unmatched(m_neighbors[v], v);

                if (u_unmatched && v_unmatched) {
                    // both unmatched, one new path of length 1
                    m_neighbors[u].n1 = v;
                    m_neighbors[u].w1 = w;
                    m_neighbors[v].n1 = u;
                    m_neighbors[v].w1 = w;
                    path_id[u]        = u;
                    path_id[v]        = u;
                    path_length[u]    = 1;
                    INCREASE_COUNTER(edges_new_paths);
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
                    path_id[u]        = v_id;
                    path_length[v_id] += 1;
                    INCREASE_COUNTER(edges_enlarge_path);
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
                        TIME_POINT(sp_solve_cycle);
                        solve_cycle(g, u, path_length[u_id], matches, matches_size);
                        TIME_POINT(ep_solve_cycle);
                        path_length[u_id] = 0;

#if COLLECT_METRICS
                        f64 t_solve = get_seconds(sp_solve_cycle, ep_solve_cycle);
                        time_solve_cycle += t_solve;
                        time_build_paths -= t_solve;
#endif
                        INCREASE_COUNTER(edges_form_cycle);
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

                vertex_t v1 = v;
                vertex_t v2 = u;
                u32 id1     = v_id;
                u32 id2     = u_id;
                if (path_length[u_id] > path_length[v_id]) {
                    std::swap(v1, v2);
                    std::swap(id1, id2);
                }

                path_length[id1] += 1 + path_length[id2];
                while (m_neighbors[v2].n2 != v2) {
                    path_id[v2]               = id1;
                    vertex_t temp_last_vertex = v1;
                    v1                        = v2;
                    v2                        = m_neighbors[v2].n1 == temp_last_vertex ? m_neighbors[v2].n2 : m_neighbors[v2].n1;
                }
                path_id[v2] = id1;

                INCREASE_COUNTER(edges_combine_paths);
            }
            TIME_POINT(ep_build_paths);

            TIME_POINT(sp_solve_paths);
            // process all paths
            for (vertex_t u : av_manager) {
                if (is_one_endpoint(m_neighbors[u], u) && path_length[path_id[u]] > 0) {
                    solve_path(g, u, path_length[path_id[u]], matches, matches_size);
                    path_length[path_id[u]] = 0;
                }
            }
            TIME_POINT(ep_solve_paths);

#if COLLECT_METRICS
            time_compute_ratings += get_seconds(sp_compute_ratings, ep_compute_ratings);
            time_build_paths += get_seconds(sp_build_paths, ep_build_paths);
            time_sorting += get_seconds(sp_sorting, ep_sorting);
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
            for (size_t i = 0; i < matches_size; ++i) {
                const auto& [u, v] = matches[i];
                ASSERT(u != v);
                ASSERT(av_manager.is_active(u));
                ASSERT(av_manager.is_active(v));
            }
#endif

#if ASSERT_ENABLED
            std::vector<u8> hit(g.get_n(), 0);
            for (size_t i = 0; i < matches_size; ++i) {
                const auto& [u, v] = matches[i];
                hit[u] += 1;
                hit[v] += 1;

                ASSERT(hit[u] == 1);
                ASSERT(hit[v] == 1);
            }
#endif
        }

        template <typename TSerialGraph>
        f32 solve_path_length_1(TSerialGraph& g,
                                const vertex_t u,
                                EdgeUV* matches,
                                size_t& matches_size) {
            matches = ASSUME_ALIGNED(EdgeUV*, matches, 64);

            vertex_t uu = u;
            vertex_t vv = m_neighbors[u].n1;
            f32 w       = m_neighbors[u].w1;

            if (g.size(uu) >= g.size(vv)) {
                std::swap(uu, vv);
            }
            matches[matches_size++] = {uu, vv};

            return w;
        }

        template <typename TSerialGraph>
        f32 solve_path_length_2(TSerialGraph& g,
                                const vertex_t u,
                                EdgeUV* matches,
                                size_t& matches_size) {
            matches = ASSUME_ALIGNED(EdgeUV*, matches, 64);

            vertex_t v1 = u;
            vertex_t v2 = m_neighbors[u].n1;
            vertex_t v3;
            f32 w1 = m_neighbors[u].w1;
            f32 w2;

            if (m_neighbors[v2].n1 == v1) {
                v3 = m_neighbors[v2].n2;
                w2 = m_neighbors[v2].w2;
            } else {
                v3 = m_neighbors[v2].n1;
                w2 = m_neighbors[v2].w1;
            }

            vertex_t uu, vv;
            f32 w;
            if (w1 > w2) {
                uu = v1;
                vv = v2;
                w  = w1;
            } else {
                uu = v2;
                vv = v3;
                w  = w2;
            }

            if (g.size(uu) >= g.size(vv)) {
                std::swap(uu, vv);
            }
            matches[matches_size++] = {uu, vv};
            return w;
        }

        template <typename TSerialGraph>
        f32 solve_path(TSerialGraph& g, const vertex_t u, const u32 length, EdgeUV* matches, size_t& matches_size) {
            matches = ASSUME_ALIGNED(EdgeUV*, matches, 64);

            // special case of length 1
            if (length == 1) {
                return solve_path_length_1(g, u, matches, matches_size);
            }

            // special case of length 2
            if (length == 2) {
                return solve_path_length_2(g, u, matches, matches_size);
            }

            vertex_t v1, v2, v3;
            f32 w;

            s64 i = 0;

            // first edge of the path
            v1          = u;
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

                    if (g.size(uu) >= g.size(vv)) {
                        std::swap(uu, vv);
                    }
                    matches[matches_size++] = {uu, vv};
                }
                idx = dp_m[idx];
            }
            return dp_w[i - 1];
        }

        template <typename TSerialGraph>
        void solve_cycle(TSerialGraph& g, vertex_t u, const u32 length, EdgeUV* matches, size_t& matches_size) {
            matches = ASSUME_ALIGNED(EdgeUV*, matches, 64);

            vertex_t n1           = m_neighbors[u].n1;
            vertex_t n2           = m_neighbors[u].n2;
            Neighbors original_u  = m_neighbors[u];
            Neighbors original_n1 = m_neighbors[original_u.n1];
            Neighbors original_n2 = m_neighbors[original_u.n2];

            f32 matching_weight1          = 0.0;
            f32 matching_weight2          = 0.0;
            size_t dp_cycle_matches1_size = 0;
            size_t dp_cycle_matches2_size = 0;

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
            matching_weight1 = solve_path(g, u, length - 1, dp_cycle_matches1, dp_cycle_matches1_size);
            m_neighbors[u]   = original_u;
            m_neighbors[n1]  = original_n1;

            // cut connection between u and n2
            m_neighbors[u].n2 = u;

            if (m_neighbors[n2].n1 == u) {
                m_neighbors[n2].n1 = original_n2.n2;
                m_neighbors[n2].w1 = original_n2.w2;
                m_neighbors[n2].n2 = n2;
            } else {
                m_neighbors[n2].n2 = n2;
            }
            matching_weight2 = solve_path(g, u, length - 1, dp_cycle_matches2, dp_cycle_matches2_size);
            m_neighbors[u]   = original_u;
            m_neighbors[n2]  = original_n2;

            EdgeUV* dp_cycle_matches     = dp_cycle_matches1;
            size_t dp_cycle_matches_size = dp_cycle_matches1_size;
            if (matching_weight2 > matching_weight1) {
                dp_cycle_matches      = dp_cycle_matches2;
                dp_cycle_matches_size = dp_cycle_matches2_size;
            }
            dp_cycle_matches = ASSUME_ALIGNED(EdgeUV*, dp_cycle_matches, 64);

            for (size_t i = 0; i < dp_cycle_matches_size; ++i) {
                matches[matches_size++] = dp_cycle_matches[i];
            }
        }

        template <typename TSerialGraph, typename TSerialActiveVertexManager>
        void random_matching(TSerialGraph& g,
                             TSerialActiveVertexManager& av_manager,
                             EdgeUV* matches,
                             size_t& matches_size) {
            matches      = ASSUME_ALIGNED(EdgeUV*, matches, 64);
            matches_size = 0;

            std::vector<u8> is_matched(g.get_n(), 0);

            for (vertex_t u : av_manager) {
                if (is_matched[u]) { continue; }
                weight_t u_w = g.get_weight(u);

                for (auto[v, w] : g[u]) {
                    if (is_matched[v]) { continue; }
                    weight_t v_w = g.get_weight(v);

                    if (u_w + v_w > m_l_max) { continue; }

                    is_matched[u] = 1;
                    is_matched[v] = 1;

                    matches[matches_size++] = {u, v};
                    break;
                }
            }
        }
    };
}

#endif //HEIPROMAP_GLOBAL_PATH_ALGORITHM_ARRAYS_H
