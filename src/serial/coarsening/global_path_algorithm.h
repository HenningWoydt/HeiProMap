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
#include <execution>
#include <queue>
#include <vector>

#include "../../commons/definitions.h"
#include "../../commons/JSON_utils.h"
#include "../../commons/random_engine.h"

namespace HeiProMap {
    struct Neighbors {
        vertex_t n1;
        vertex_t n2;
        f32 w1;
        f32 w2;
    };

#define IS_ENDPOINT(n, u) n.n2 == u
#define IS_NOT_ENDPOINT(n, u) n.n2 != u
#define IS_ONE_ENDPOINT(n, u) n.n1 != u && n.n2 == u
#define IS_UNMATCHED(n, u) n.n1 == u

    class GlobalPathAlgorithmConfiguration {
    public:
        size_t random_level = 4;
    };

    /**
     * Computes a matching based on the Global Path Algorithm from
     * > Jens Maue and Peter Sanders.
     * > Engineering Algorithms for Approximate Weighted Matching.
     * > Experimental Algorithms, 6th International Workshop, WEA 2007, Rome, Italy, June 6-8, 2007, Proceedings.
     */
    class GlobalPathAlgorithmMatcher {
        vertex_t m_n     = 0;
        vertex_t m_m     = 0;
        partition_t m_k  = 0;
        weight_t m_l_max = 0;
        u64 m_threads    = 1;

        const GlobalPathAlgorithmConfiguration* config = nullptr;
        RandomEngine* random_engine                    = nullptr;

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

    public:
        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const weight_t t_l_max,
                        const u64 t_threads,
                        RandomEngine& t_random_engine,
                        const GlobalPathAlgorithmConfiguration& i_config) {
            m_n       = t_n;
            m_m       = t_m;
            m_k       = t_k;
            m_l_max   = t_l_max;
            m_threads = t_threads;

            config           = dynamic_cast<const GlobalPathAlgorithmConfiguration*>(&i_config);
            random_engine    = &t_random_engine;

            m_neighbors.initialize(m_n);
            path_id.initialize(m_n);
            path_length.initialize(m_n);

            edges.initialize(m_m);

            dp_w.initialize(m_n);
            dp_m.initialize(m_n);
            dp_take.initialize(m_n);
            dp_edges.initialize(m_n);

            dp_cycle_matches1.initialize(m_n);
            dp_cycle_matches2.initialize(m_n);
        }

        template <typename PartitionManagerT>
        void match(const size_t level,
                   const graph_t& g,
                   const PartitionManagerT& p_manager,
                   Matching& matching) {
            if (level < config->random_level) {
                // use a random matching
                random_matching(level, g, matching);
                return;
            }

            compute_ratings(g, p_manager);

            std::sort(std::execution::par, edges.get_ptr(), edges.get_ptr() + edges_size, std::greater<>());

            for (vertex_t u = 0; u < g.get_n(); ++u) {
                m_neighbors[u].n1 = u;
                m_neighbors[u].n2 = u;
            }

            for (size_t i = 0; i < edges_size; ++i) {
                auto [u, v, w] = edges[i];

                if (IS_NOT_ENDPOINT(m_neighbors[u], u) || IS_NOT_ENDPOINT(m_neighbors[v], v)) {
                    // u or v is not an endpoint
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
                    path_id[u]        = u;
                    path_id[v]        = u;
                    path_length[u]    = 1;
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
                        solve_cycle(g, u, path_length[u_id], matching);
                        path_length[u_id] = 0;
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
            }

            // process all paths
            forall_gu(g, u)
                {
                    if (IS_ONE_ENDPOINT(m_neighbors[u], u) && path_length[path_id[u]] > 0) {
                        solve_path(g, u, path_length[path_id[u]], matching);
                        path_length[path_id[u]] = 0;
                    }
                }
            endfor

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

        template <typename PartitionManagerT>
        void compute_ratings(const graph_t& g, const PartitionManagerT& p_manager) {

            edges_size = 0;
#pragma omp parallel num_threads(m_threads)
            {
                std::vector<EdgeUVW> local_edges;

#pragma omp for
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
                                edge_rating = ((f32)(w * w)) / ((f32)(u_w * v_w));
                                // edge_rating = ((f32) (w * w)) / ((f32) (u_w + v_w));
                                // edge_rating = ((f32) (w * w * w)) / ((f32) (u_w * v_w));
                                // edge_rating = (f32) w / (f32) (u_w * v_w * u_w * v_w);
                                // edge_rating = (f32)w;

                                // edges[edges_size++] = {u, v, edge_rating};
                                local_edges.push_back({u, v, edge_rating});
                            }
                        endfor
                    }
                endfor

                // Critical section or lock-free append
#pragma omp critical
                {
                    for (const auto& e : local_edges) {
                        edges[edges_size++] = e;
                    }
                }
            }
        }

        f32 solve_path_length_1(const graph_t& g,
                                const vertex_t u,
                                Matching& matching) {
            vertex_t uu = u;
            vertex_t vv = m_neighbors[u].n1;
            f32 w       = m_neighbors[u].w1;

            matching.add(uu, vv);

            return w;
        }

        f32 solve_path_length_2(const graph_t& g,
                                const vertex_t u,
                                Matching& matching) {
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

            matching.add(uu, vv);
            return w;
        }

        f32 solve_path(const graph_t& g,
                       const vertex_t u,
                       const u32 length,
                       Matching& matching) {
            // special case of length 1
            if (length == 1) {
                return solve_path_length_1(g, u, matching);
            }

            // special case of length 2
            if (length == 2) {
                return solve_path_length_2(g, u, matching);
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

                    matching.add(uu, vv);
                }
                idx = dp_m[idx];
            }
            return dp_w[i - 1];
        }

        void solve_cycle(const graph_t& g,
                         const vertex_t u,
                         const u32 length,
                         Matching& matching) {
            vertex_t n1           = m_neighbors[u].n1;
            vertex_t n2           = m_neighbors[u].n2;
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
            matching_weight2 = solve_path(g, u, length - 1, dp_cycle_matches2);
            m_neighbors[u]   = original_u;
            m_neighbors[n2]  = original_n2;

            Matching* dp_cycle_matches = &dp_cycle_matches1;
            if (matching_weight2 > matching_weight1) {
                dp_cycle_matches = &dp_cycle_matches2;
            }

            for (size_t i = 0; i < dp_cycle_matches->size(); ++i) {
                matching.add((*dp_cycle_matches)[i].u, (*dp_cycle_matches)[i].v);
            }
        }

        void random_matching(const size_t level,
                             const graph_t& g,
                             Matching& matching) {
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
