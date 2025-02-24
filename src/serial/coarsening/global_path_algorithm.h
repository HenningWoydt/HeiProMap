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
#include <iomanip>
#include <numeric>
#include <queue>
#include <vector>

#include "../../definitions.h"
#include "../../macros.h"
#include "../interfaces/ISerialMatcher.h"

namespace HeiProMap {
    struct Neighbors {
        vertex_t n1;
        vertex_t n2;
        f32      w1;
        f32      w2;
    };

    inline bool is_endpoint_fast(const Neighbors &n, const vertex_t u) {
        return n.n2 == u;
    }

    inline bool is_endpoint(const Neighbors &n, const vertex_t u) {
        return is_endpoint_fast(n, u);
        bool b = n.n1 == u || n.n2 == u;
        ASSERT(b == is_endpoint_fast(n, u));
        return b;
    }

    inline bool is_one_endpoint_fast(const Neighbors &n, const vertex_t u) {
        return n.n1 != u && n.n2 == u;
    }

    inline bool is_one_endpoint(const Neighbors &n, const vertex_t u) {
        return is_one_endpoint_fast(n, u);
        bool b = (n.n1 == u && n.n2 != u) || (n.n1 != u && n.n2 == u);
        ASSERT(b == is_one_endpoint_fast(n, u));
        return b;
    }

    inline bool is_unmatched_fast(const Neighbors &n, const vertex_t u) {
        return n.n1 == u;
    }

    inline bool is_unmatched(const Neighbors &n, const vertex_t u) {
        return is_unmatched_fast(n, u);
        bool b = n.n1 == u && n.n2 == u;
        ASSERT(b == is_unmatched_fast(n, u));
        return b;
    }

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

        std::vector<Neighbors> m_neighbors;
        std::vector<u32>       path_id;
        std::vector<u32>       path_length;

        std::vector<EdgeUVW> edges_greater_1;
        std::vector<EdgeUVW> edges_equal_1;
        std::vector<EdgeUVW> edges_smaller_1;

        // for DP
        std::vector<f32>     w;
        std::vector<s64>     m;
        std::vector<u8>      take;
        std::vector<EdgeUVW> dp_edges;

        std::vector<EdgeUV> dp_cycle_matches1;
        std::vector<EdgeUV> dp_cycle_matches2;

        f64 global_time_compute_ratings = 0.0;
        f64 global_time_sorting         = 0.0;
        f64 global_time_build_paths     = 0.0;
        f64 global_time_solve_paths     = 0.0;
        u64 global_edges_skipped        = 0;
        u64 global_edges_form_cycle     = 0;
        u64 global_edges_new_paths      = 0;
        u64 global_edges_enlarge_path   = 0;
        u64 global_edges_combine_paths  = 0;

    public:
        GlobalPathAlgorithmMatcher() = default;

        void initialize(const vertex_t t_n, const vertex_t t_m, const partition_t t_k, const weight_t t_l_max) override {
            m_n     = t_n;
            m_m     = t_m;
            m_k     = t_k;
            m_l_max = t_l_max;

            m_neighbors.resize(m_n);
            path_id.resize(m_n);
            path_length.resize(m_n);

            edges_greater_1.reserve(m_n);
            edges_equal_1.reserve(m_n);
            edges_smaller_1.reserve(m_n);
        }

        template<typename TSerialGraph, typename TSerialActiveVertexManager>
        void match(TSerialGraph &g,
                   TSerialActiveVertexManager &av_manager,
                   std::vector<EdgeUV> &matches) {
            for (vertex_t u = 0; u < m_neighbors.size(); ++u) {
                m_neighbors[u].n1 = u;
                m_neighbors[u].n2 = u;
            }

            f64 time_compute_ratings = 0.0;
            f64 time_sorting         = 0.0;
            f64 time_build_paths     = 0.0;
            f64 time_solve_paths     = 0.0;
            u64 edges_skipped        = 0;
            u64 edges_form_cycle     = 0;
            u64 edges_new_paths      = 0;
            u64 edges_enlarge_path   = 0;
            u64 edges_combine_paths  = 0;

            auto sp_compute_ratings = std::chrono::high_resolution_clock::now();

            edges_greater_1.clear();
            edges_equal_1.clear();
            edges_smaller_1.clear();
            for (vertex_t u: av_manager) {
                weight_t u_w = g.get_weight(u);

                for (const auto &[v, w]: g[u]) {
                    weight_t v_w = g.get_weight(v);

                    if (u_w + v_w > m_l_max) {
                        continue;
                    }

                    if (v < u) {
                        continue;
                    }
                    // f32 edge_rating = (((f64) w) * ((f64) w)) / (g.get_weight(u) * g.get_weight(v));
                    // f32 edge_rating = ((f64) w) / (g.size(u) * g.size(v));
                    f32 edge_rating = (f32) w / (f32) (u_w * v_w);

                    if (edge_rating == 0.0) {
                        continue;
                    }

                    if (edge_rating > 1.0) {
                        edges_greater_1.emplace_back(u, v, edge_rating);
                    } else if (edge_rating == 1.0) {
                        edges_equal_1.emplace_back(u, v, edge_rating);
                    } else {
                        edges_smaller_1.emplace_back(u, v, edge_rating);
                    }

                }
            }
            auto ep_compute_ratings = std::chrono::high_resolution_clock::now();
            time_compute_ratings += get_seconds(sp_compute_ratings, ep_compute_ratings);

            auto                              sp_build_paths    = std::chrono::high_resolution_clock::now();
            std::vector<std::vector<EdgeUVW>> buckets           = {edges_greater_1, edges_equal_1, edges_smaller_1};
            std::vector<bool>                 all_ratings_equal = {false, true, false};

            for (size_t j = 0; j < buckets.size(); ++j) {
                std::vector<EdgeUVW> &edge_list    = buckets[j];
                bool                 equal_ratings = all_ratings_equal[j];

                auto sp_sorting = std::chrono::high_resolution_clock::now();
                // filter out all edges that cannot be used anymore, before sorting the bucket
                if (!equal_ratings) {
                    auto is_not_endpoint = [&](const EdgeUVW &edge) {
                        return !is_endpoint(m_neighbors[edge.u], edge.u) || !is_endpoint(m_neighbors[edge.v], edge.v);
                    };

                    edge_list.erase(std::remove_if(edge_list.begin(), edge_list.end(), is_not_endpoint), edge_list.end());
                    std::sort(edge_list.begin(), edge_list.end(), std::greater<>());
                }
                auto ep_sorting = std::chrono::high_resolution_clock::now();
                f64  t          = get_seconds(sp_sorting, ep_sorting);
                time_sorting += t;
                time_build_paths -= t;


                for (const auto &[u, v, w]: edge_list) {
                    if (!is_endpoint(m_neighbors[u], u) || !is_endpoint(m_neighbors[v], v)) {
                        // u is not an endpoint
                        edges_skipped++;
                        continue;
                    }

                    bool u_unmatched = is_unmatched(m_neighbors[u], u);
                    bool v_unmatched = is_unmatched(m_neighbors[v], v);
                    u32  u_id        = path_id[u];
                    u32  v_id        = path_id[v];

                    if (u_unmatched && v_unmatched) {
                        m_neighbors[u].n1 = v;
                        m_neighbors[u].w1 = w;
                        m_neighbors[v].n1 = u;
                        m_neighbors[v].w1 = w;
                        path_id[u]     = u;
                        path_id[v]     = u;
                        path_length[u] = 1;
                        edges_new_paths++;
                        continue;
                    }
                    if (u_unmatched && !v_unmatched) {
                        m_neighbors[u].n1 = v;
                        m_neighbors[u].w1 = w;
                        m_neighbors[v].n2 = u;
                        m_neighbors[v].w2 = w;
                        path_id[u] = v_id;
                        path_length[v_id] += 1;
                        edges_enlarge_path++;
                        continue;
                    }
                    if (!u_unmatched && v_unmatched) {
                        m_neighbors[u].n2 = v;
                        m_neighbors[u].w2 = w;
                        m_neighbors[v].n1 = u;
                        m_neighbors[v].w1 = w;
                        path_id[v] = u_id;
                        path_length[u_id] += 1;
                        edges_enlarge_path++;
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
                            auto sp_solve_paths = std::chrono::high_resolution_clock::now();
                            solve_cycle(g, u, matches);
                            auto ep_solve_paths = std::chrono::high_resolution_clock::now();
                            f64  t              = get_seconds(sp_solve_paths, ep_solve_paths);
                            time_solve_paths += t;
                            time_build_paths -= t;
                            edges_form_cycle++;
                        }
                        continue;
                    }

                    // two paths

                    // both u and v connect larger paths

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
                    if (path_length[u_id] <= path_length[u_id]) {
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

                    edges_combine_paths++;
                }
            }
            auto ep_build_paths = std::chrono::high_resolution_clock::now();
            time_build_paths += get_seconds(sp_build_paths, ep_build_paths);

            auto          sp_solve_paths = std::chrono::high_resolution_clock::now();
            // process all paths
            for (vertex_t u: av_manager) {
                if (is_one_endpoint(m_neighbors[u], u)) {
                    solve_path(g, u, matches);
                }
            }
            auto          ep_solve_paths = std::chrono::high_resolution_clock::now();
            time_solve_paths += get_seconds(sp_solve_paths, ep_solve_paths);

            /*
            global_time_compute_ratings += time_compute_ratings;
            global_time_sorting += time_sorting;
            global_time_build_paths += time_build_paths;
            global_time_solve_paths += time_solve_paths;
            f64 global_time_total = global_time_compute_ratings + global_time_sorting + global_time_build_paths + global_time_solve_paths;

            f64 time_total = time_compute_ratings + time_sorting + time_build_paths + time_solve_paths;
            std::cout << std::fixed << std::setprecision(4);
            std::cout << "time compute ratings: " << std::setprecision(4) << time_compute_ratings << " - " << time_compute_ratings / time_total << " -- global time compute ratings: " << global_time_compute_ratings << " - " << global_time_compute_ratings / global_time_total << std::endl;
            std::cout << "time sorting        : " << std::setprecision(4) << time_sorting << " - " << time_sorting / time_total << " -- global time sorting        : " << global_time_sorting << " - " << global_time_sorting / global_time_total << std::endl;
            std::cout << "time build paths    : " << std::setprecision(4) << time_build_paths << " - " << time_build_paths / time_total << " -- global time build paths    : " << global_time_build_paths << " - " << global_time_build_paths / global_time_total << std::endl;
            std::cout << "time solve paths    : " << std::setprecision(4) << time_solve_paths << " - " << time_solve_paths / time_total << " -- global time solve paths    : " << global_time_solve_paths << " - " << global_time_solve_paths / global_time_total << std::endl;
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
            */

#if ASSERT_ENABLED
            for (const EdgeUV& e : matches) {
                ASSERT(e.u != e.v);
                ASSERT(av_manager.is_active(e.u));
                ASSERT(av_manager.is_active(e.v));
            }
#endif

#if ASSERT_ENABLED
            std::vector<u8> hit(g.get_n(), 0);
            for (const auto& e : matches) {
                hit[e.u] += 1;
                hit[e.v] += 1;

                ASSERT(hit[e.u] != 2);
                ASSERT(hit[e.v] != 2);
            }
#endif
        }

        template<typename TSerialGraph>
        f64 solve_path(TSerialGraph &g,
                       const vertex_t u,
                       std::vector<EdgeUV> &matches) {
            // save the path into better format
            dp_edges.clear();
            dp_edges.reserve(path_length[path_id[u]]);

            if (m_neighbors[u].n1 == u) {
                dp_edges.emplace_back(u, m_neighbors[u].n2, m_neighbors[u].w2);
            } else {
                dp_edges.emplace_back(u, m_neighbors[u].n1, m_neighbors[u].w1);
            }

            vertex_t last_vertex = u;
            vertex_t curr_vertex = m_neighbors[u].n1 == u ? m_neighbors[u].n2 : m_neighbors[u].n1;
            while (m_neighbors[curr_vertex].n1 != curr_vertex && m_neighbors[curr_vertex].n2 != curr_vertex) {
                if (m_neighbors[curr_vertex].n1 == last_vertex) {
                    dp_edges.emplace_back(curr_vertex, m_neighbors[curr_vertex].n2, m_neighbors[curr_vertex].w2);

                    last_vertex = curr_vertex;
                    curr_vertex = m_neighbors[curr_vertex].n2;
                } else {
                    dp_edges.emplace_back(curr_vertex, m_neighbors[curr_vertex].n1, m_neighbors[curr_vertex].w1);

                    last_vertex = curr_vertex;
                    curr_vertex = m_neighbors[curr_vertex].n1;
                }
            }

            if (dp_edges.size() == 1) {
                size_t idx = 0;

                vertex_t uu = dp_edges[idx].u;
                vertex_t vv = dp_edges[idx].v;
                if (g.size(uu) < g.size(vv)) {
                    matches.emplace_back(uu, vv);
                } else {
                    matches.emplace_back(vv, uu);
                }

                // destroy the endpoints
                vertex_t v1 = dp_edges[0].u;
                vertex_t v2 = dp_edges.back().v;
                m_neighbors[v1].n1 = -1;
                m_neighbors[v1].n2 = -1;
                m_neighbors[v2].n1 = -1;
                m_neighbors[v2].n2 = -1;
                return dp_edges[idx].w;
            }
            if (dp_edges.size() == 2) {
                size_t idx = dp_edges[0].w > dp_edges[1].w ? 0 : 1;

                vertex_t uu = dp_edges[idx].u;
                vertex_t vv = dp_edges[idx].v;
                if (g.size(uu) < g.size(vv)) {
                    matches.emplace_back(uu, vv);
                } else {
                    matches.emplace_back(vv, uu);
                }

                // destroy the endpoints
                vertex_t v1 = dp_edges[0].u;
                vertex_t v2 = dp_edges.back().v;
                m_neighbors[v1].n1 = -1;
                m_neighbors[v1].n2 = -1;
                m_neighbors[v2].n1 = -1;
                m_neighbors[v2].n2 = -1;
                return dp_edges[idx].w;
            }

            w.resize(dp_edges.size());
            m.resize(dp_edges.size());
            take.resize(dp_edges.size());

            w[0]    = dp_edges[0].w;
            m[0]    = -1;
            take[0] = 1;

            if (dp_edges[1].w > w[0]) {
                w[1]    = dp_edges[1].w;
                m[1]    = -1;
                take[1] = 1;
            } else {
                w[1]    = w[0];
                m[1]    = 0;
                take[1] = 0;
            }

            for (size_t i = 2; i < dp_edges.size(); ++i) {
                if (dp_edges[i].w + w[i - 2] > w[i - 1]) {
                    w[i]    = dp_edges[i].w + w[i - 2];
                    m[i]    = i - 2;
                    take[i] = 1;
                } else {
                    w[i]    = w[i - 1];
                    m[i]    = i - 1;
                    take[i] = 0;
                }
            }

            s64 idx             = m[dp_edges.size() - 1];
            f64 matching_weight = 0.0;
            while (idx != -1) {
                if (take[idx]) {

                    vertex_t uu = dp_edges[idx].u;
                    vertex_t vv = dp_edges[idx].v;
                    if (g.size(uu) < g.size(vv)) {
                        matches.emplace_back(uu, vv);
                    } else {
                        matches.emplace_back(vv, uu);
                    }

                    matching_weight += dp_edges[idx].w;
                }
                idx = m[idx];
            }

            // destroy the endpoints
            vertex_t v1 = dp_edges[0].u;
            vertex_t v2 = dp_edges.back().v;
            m_neighbors[v1].n1 = -1;
            m_neighbors[v1].n2 = -1;
            m_neighbors[v2].n1 = -1;
            m_neighbors[v2].n2 = -1;
            return matching_weight;
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

            f64 matching_weight1 = 0.0;
            f64 matching_weight2 = 0.0;
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

#endif //HEIPROMAP_GLOBAL_PATH_ALGORITHM_H
