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
#include <vector>

#include "../definitions.h"
#include "../utility/random_engine.h"
#include "../utility/aligned_array.h"
#include "../definitions_1.h"
#include "../definitions_2.h"

namespace HeiProMap {
    class GreedyEdgeMatcherConfiguration {
    public:
        bool match_pendant_vertices_first = false; // Vertices with only one neighbor should be handled first.
        bool match_triangles = false; // triangle where one vertex is heavily connected and the other two are not.
    };

    class GreedyEdgeMatcher {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;

        const GreedyEdgeMatcherConfiguration *config = nullptr;
        RandomEngine *random_engine = nullptr;

        u32 mark = 0;
        AlignedArray<u32> used;

    public:
        GreedyEdgeMatcher() = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        RandomEngine &t_random_engine,
                        const GreedyEdgeMatcherConfiguration &i_config) {
            ScopedTimer _t("io", "GreedyEdgeMatcher", "initialize");

            m_n = t_n;
            m_m = t_m;
            m_k = t_k;

            config = dynamic_cast<const GreedyEdgeMatcherConfiguration *>(&i_config);
            random_engine = &t_random_engine;

            mark = 0;
            used.initialize(m_n, 0);
        }

        void match([[maybe_unused]] const size_t level,
                   const graph_t &g,
                   [[maybe_unused]] p_manager_t &p_manager,
                   Mapping &mapping,
                   f64 imbalance,
                   u64 threads = 1) {
            if (threads == 1) {
                match_serial(level, g, p_manager, mapping, imbalance);
                return;
            }

            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            std::vector<EdgeUVW> edges;

            Matching matching;
            {
                ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "initialize");

                edges.reserve(g.m);
                matching.initialize(g.n);
            }

            mark += 1;

            // first handle pendant vertices
            if (config->match_pendant_vertices_first) {
                edges.clear();
                {
                    ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "pendant_rate");

                    forall_gu(g, u)
                        {
                            if (g.deg(u) != 1) { continue; }

                            vertex_t v = g.edges_v[g.neighborhoods[u]];
                            weight_t ew = g.edges_w[g.neighborhoods[u]];

                            if (g.v_weights[u] + g.v_weights[v] > lmax) { continue; }

                            f32 rating = (f32) ew / (f32) (g.deg(u) * g.deg(v));
                            edges.emplace_back(u, v, rating);
                        }
                    endfor
                }
                {
                    ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "pendant_sort");
                    std::sort(edges.begin(), edges.end(), std::greater<>());
                }
                {
                    ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "pendant_choose");

                    for (size_t i = 0; i < edges.size(); ++i) {
                        auto &[u, v, w] = edges[i];
                        if (used[u] == mark || used[v] == mark) { continue; }

                        used[u] = mark;
                        used[v] = mark;
                        matching.add(u, v);
                    }
                }
            }

            if (config->match_triangles) {
                edges.clear();
                {
                    ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "triangles_rate");

                    forall_gu(g, u)
                        {
                            if (used[u] == mark) { continue; }
                            if (g.deg(u) != 2) { continue; }

                            vertex_t u0 = g.neighborhoods[u];
                            vertex_t u1 = u0 + 1;

                            vertex_t a = g.edges_v[u0];
                            vertex_t b = g.edges_v[u1];
                            weight_t wa = g.edges_w[u0];
                            weight_t wb = g.edges_w[u1];

                            if (g.deg(a) == 2 && used[a] != mark) {
                                vertex_t a0 = g.neighborhoods[a];
                                vertex_t a1 = a0 + 1;

                                vertex_t a_n0 = g.edges_v[a0];
                                vertex_t a_n1 = g.edges_v[a1];

                                if ((a_n0 == u && a_n1 == b) || (a_n1 == u && a_n0 == b)) {
                                    if (u < a && g.v_weights[u] + g.v_weights[a] <= lmax) {
                                        f32 rating = (f32) wa / (f32) (g.deg(u) * g.deg(a));
                                        edges.emplace_back(u, a, rating);
                                    }
                                }
                            }

                            if (g.deg(b) == 2 && used[b] != mark) {
                                vertex_t b0 = g.neighborhoods[b];
                                vertex_t b1 = b0 + 1;

                                vertex_t b_n0 = g.edges_v[b0];
                                vertex_t b_n1 = g.edges_v[b1];

                                if ((b_n0 == u && b_n1 == a) || (b_n1 == u && b_n0 == a)) {
                                    if (u < b && g.v_weights[u] + g.v_weights[b] <= lmax) {
                                        f32 rating = (f32) wb / (f32) (g.deg(u) * g.deg(b));
                                        edges.emplace_back(u, b, rating);
                                    }
                                }
                            }
                        }
                    endfor
                }
                {
                    ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "triangles_sort");
                    std::sort(edges.begin(), edges.end(), std::greater<>());
                }
                {
                    ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "triangles_choose");

                    for (size_t i = 0; i < edges.size(); ++i) {
                        auto &[u, v, w] = edges[i];
                        if (used[u] == mark || used[v] == mark) { continue; }

                        used[u] = mark;
                        used[v] = mark;
                        matching.add(u, v);
                    }
                }
            }

            edges.clear();

            std::vector<std::vector<EdgeUVW> > thread_edges(threads);
            //
            {
                ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "rate");

                #pragma omp parallel num_threads(threads)
                {
                    int tid = omp_get_thread_num();
                    auto &local_edges = thread_edges[tid];
                    local_edges.reserve(g.m / (size_t) threads + 1);

                    #pragma omp for
                    forall_gu(g, u)
                        {
                            if (used[u] == mark) { continue; }

                            weight_t u_w = g.v_weights[u];

                            forall_guivw(g, u, j, v, w)
                                {
                                    if (u >= v) { continue; } // avoid duplicates
                                    if (used[v] == mark) { continue; }

                                    weight_t v_w = g.v_weights[v];
                                    if (u_w + v_w > lmax) { continue; }

                                    f32 edge_rating = (f32) w / (f32) (u_w * v_w);
                                    local_edges.emplace_back(u, v, edge_rating);
                                }
                            endfor
                        }
                    endfor
                }

                {
                    ScopedTimer _t2("coarsening", "GreedyEdgeMatcher", "sort");

                    // sort each thread-local vector in parallel
                    #pragma omp parallel for num_threads(threads)
                    for (u64 t = 0; t < threads; ++t) {
                        std::sort(thread_edges[t].begin(), thread_edges[t].end(), std::greater<>());
                    }
                }
            }

            {
                ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "choose");

                std::vector<size_t> pos(threads, 0);

                while (true) {
                    bool found = false;
                    u64 best_t = 0;

                    // find best front edge among all thread-local vectors
                    for (u64 t = 0; t < threads; ++t) {
                        if (pos[t] >= thread_edges[t].size()) { continue; }

                        if (!found || thread_edges[t][pos[t]] > thread_edges[best_t][pos[best_t]]) {
                            found = true;
                            best_t = t;
                        }
                    }

                    if (!found) { break; }

                    auto &[u, v, w] = thread_edges[best_t][pos[best_t]];
                    ++pos[best_t];

                    if (used[u] != mark && used[v] != mark) {
                        used[u] = mark;
                        used[v] = mark;
                        matching.add(u, v);
                    }
                }
            }

            {
                ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "map");

                matching.set_translation();
                mapping.set_coarse_n(matching.get_n_coarse_nodes());
                for (vertex_t u = 0; u < matching.get_n(); ++u) {
                    mapping.set(u, matching.get_n(u));
                }
            }
        }

        void match_serial([[maybe_unused]] const size_t level,
                          const graph_t &g,
                          [[maybe_unused]] p_manager_t &p_manager,
                          Mapping &mapping,
                          f64 imbalance) {
            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));
            std::vector<EdgeUVW> edges;
            Matching matching;
            //
            {
                ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "initialize");
                edges.reserve(g.m);
                matching.initialize(g.n);
            }
            mark += 1;
            // first handle pendant vertices
            if (config->match_pendant_vertices_first) {
                edges.clear();
                //
                {
                    ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "pendant_rate");
                    forall_gu(g, u)
                        {
                            if (g.deg(u) != 1) { continue; }
                            vertex_t v = g.edges_v[g.neighborhoods[u]];
                            weight_t ew = g.edges_w[g.neighborhoods[u]];
                            if (g.v_weights[u] + g.v_weights[v] > lmax) { continue; }
                            f32 rating = (f32) ew / (f32) (g.deg(u) * g.deg(v));
                            edges.emplace_back(u, v, rating);
                        }
                    endfor
                }
                //
                {
                    ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "pendant_sort");
                    std::sort(edges.begin(), edges.end(), std::greater<>());
                }
                //
                {
                    ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "pendant_choose");
                    for (size_t i = 0; i < edges.size(); ++i) {
                        auto &[u, v, w] = edges[i];
                        if (used[u] == mark || used[v] == mark) { continue; }
                        // use this edge
                        used[u] = mark;
                        used[v] = mark;
                        matching.add(u, v);
                    }
                }
            }
            if (config->match_triangles) {
                edges.clear();
                //
                {
                    ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "triangles_rate");
                    //
                    forall_gu(g, u)
                        {
                            if (used[u] == mark) { continue; }
                            if (g.deg(u) != 2) { continue; }
                            vertex_t u0 = g.neighborhoods[u];
                            vertex_t u1 = u0 + 1;
                            vertex_t a = g.edges_v[u0];
                            vertex_t b = g.edges_v[u1];
                            weight_t wa = g.edges_w[u0];
                            weight_t wb = g.edges_w[u1];
                            // Try matching u with a if a is also degree 2 and shares b.
                            if (g.deg(a) == 2 && used[a] != mark) {
                                vertex_t a0 = g.neighborhoods[a];
                                vertex_t a1 = a0 + 1;
                                vertex_t a_n0 = g.edges_v[a0];
                                vertex_t a_n1 = g.edges_v[a1];
                                if ((a_n0 == u && a_n1 == b) || (a_n1 == u && a_n0 == b)) {
                                    if (u < a && g.v_weights[u] + g.v_weights[a] <= lmax) {
                                        // weight of edge (u, a) is wa
                                        f32 rating = (f32) wa / (f32) (g.deg(u) * g.deg(a));
                                        edges.emplace_back(u, a, rating);
                                    }
                                }
                            } // Try matching u with b if b is also degree 2 and shares a.
                            if (g.deg(b) == 2 && used[b] != mark) {
                                vertex_t b0 = g.neighborhoods[b];
                                vertex_t b1 = b0 + 1;
                                vertex_t b_n0 = g.edges_v[b0];
                                vertex_t b_n1 = g.edges_v[b1];
                                if ((b_n0 == u && b_n1 == a) || (b_n1 == u && b_n0 == a)) {
                                    if (u < b && g.v_weights[u] + g.v_weights[b] <= lmax) {
                                        // weight of edge (u, b) is wb
                                        f32 rating = (f32) wb / (f32) (g.deg(u) * g.deg(b));
                                        edges.emplace_back(u, b, rating);
                                    }
                                }
                            }
                        }
                    endfor
                } //
                {
                    ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "triangles_sort");
                    std::sort(edges.begin(), edges.end(), std::greater<>());
                } //
                {
                    ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "triangles_choose");
                    for (size_t i = 0; i < edges.size(); ++i) {
                        auto &[u, v, w] = edges[i];
                        if (used[u] == mark || used[v] == mark) { continue; }
                        used[u] = mark;
                        used[v] = mark;
                        matching.add(u, v);
                    }
                }
            }
            edges.clear();
            // handle all other vertices
            {
                ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "rate");
                forall_gu(g, u)
                    {
                        weight_t u_w = g.v_weights[u];
                        if (used[u] == mark) { continue; }
                        forall_guivw(g, u, j, v, w)
                            {
                                if (u >= v) { continue; }
                                if (used[v] == mark) { continue; }
                                weight_t v_w = g.v_weights[v];
                                if (u_w + v_w > lmax) { continue; }
                                f32 edge_rating = (f32) w / (f32) (u_w * v_w);
                                edges.emplace_back(u, v, edge_rating);
                            }
                        endfor
                    }
                endfor
            } //
            {
                ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "sort");
                std::sort(edges.begin(), edges.end(), std::greater<>());
            } //
            {
                ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "choose");
                for (size_t i = 0; i < edges.size(); ++i) {
                    auto &[u, v, w] = edges[i];
                    if (used[u] != mark && used[v] != mark) {
                        // use this edge
                        used[u] = mark;
                        used[v] = mark;
                        matching.add(u, v);
                    }
                }
            } //
            {
                ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "map");
                matching.set_translation();
                mapping.set_coarse_n(matching.get_n_coarse_nodes());
                for (vertex_t u = 0; u < matching.get_n(); ++u) { mapping.set(u, matching.get_n(u)); }
            }
        }
    };
}


#endif //HEIPROMAP_GREEDY_EDGE_MATCHER_H
