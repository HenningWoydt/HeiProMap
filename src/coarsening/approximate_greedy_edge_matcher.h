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

#ifndef HEIPROMAP_APPROXIMATE_GREEDY_EDGE_MATCHER_H
#define HEIPROMAP_APPROXIMATE_GREEDY_EDGE_MATCHER_H

#include <algorithm>
#include <vector>
#include <limits>
#include <cmath>

#include "../definitions.h"
#include "../utility/random_engine.h"
#include "../utility/aligned_array.h"
#include "../definitions_1.h"
#include "../definitions_2.h"

namespace HeiProMap {
    class ApproximateGreedyEdgeMatcherConfiguration {
    public:
        bool match_pendant_vertices_first = false;
        bool match_triangles = false;
        size_t num_buckets = 2048;

        bool sample_edges = true;
        size_t max_samples_per_vertex = 16;
        size_t sample_degree_threshold = 32;
    };

    class ApproximateGreedyEdgeMatcher {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;

        ApproximateGreedyEdgeMatcherConfiguration config;
        RandomEngine random_engine = RandomEngine(0);

        AlignedArray<EdgeUVW> edges;
        AlignedArray<EdgeUVW> tmp_edges;
        size_t edges_size = 0;

        u32 mark = 0;
        AlignedArray<u32> used;

        std::vector<size_t> bucket_sizes;
        std::vector<size_t> bucket_offsets;
        std::vector<size_t> bucket_next;

    public:
        ApproximateGreedyEdgeMatcher() = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const u64 t_seed,
                        const ApproximateGreedyEdgeMatcherConfiguration &t_config) {
            ScopedTimer _t("coarsening", "ApproximateGreedyEdgeMatcher", "initialize");

            m_n = t_n;
            m_m = t_m;
            m_k = t_k;

            config = t_config;
            random_engine = RandomEngine(t_seed);

            edges.initialize(m_m);
            tmp_edges.initialize(m_m);

            mark = 0;
            used.initialize(m_n, 0);

            if (config.num_buckets == 0) { config.num_buckets = 2048; }
            bucket_sizes.resize(config.num_buckets, 0);
            bucket_offsets.resize(config.num_buckets, 0);
            bucket_next.resize(config.num_buckets, 0);
        }

    private:
        void approximate_sort_edges(const size_t n_edges) {
            if (n_edges <= 1) { return; }

            f32 min_rating = edges[0].w;
            f32 max_rating = edges[0].w;

            for (size_t i = 1; i < n_edges; ++i) {
                min_rating = std::min(min_rating, edges[i].w);
                max_rating = std::max(max_rating, edges[i].w);
            }

            if (min_rating == max_rating) { return; }

            std::fill(bucket_sizes.begin(), bucket_sizes.end(), 0);

            const f32 range = max_rating - min_rating;
            const size_t last_bucket = config.num_buckets - 1;

            auto bucket_of = [&](const f32 rating) -> size_t {
                const f32 x = (rating - min_rating) / range;
                const size_t b = static_cast<size_t>(x * static_cast<f32>(last_bucket));
                return std::min(b, last_bucket);
            };

            for (size_t i = 0; i < n_edges; ++i) {
                ++bucket_sizes[bucket_of(edges[i].w)];
            }

            size_t sum = 0;
            for (size_t b = 0; b < config.num_buckets; ++b) {
                bucket_offsets[b] = sum;
                sum += bucket_sizes[b];
            }

            for (size_t b = 0; b < config.num_buckets; ++b) {
                bucket_next[b] = bucket_offsets[b];
            }

            for (size_t i = 0; i < n_edges; ++i) {
                const size_t b = bucket_of(edges[i].w);
                tmp_edges[bucket_next[b]++] = edges[i];
            }

            size_t out = 0;
            for (size_t b = config.num_buckets; b-- > 0;) {
                const size_t begin = bucket_offsets[b];
                const size_t end = begin + bucket_sizes[b];
                for (size_t i = begin; i < end; ++i) {
                    edges[out++] = tmp_edges[i];
                }
            }
        }

    public:
        void match([[maybe_unused]] const size_t level,
                   const graph_t &g,
                   const p_manager_t &p_manager,
                   Mapping &mapping,
                   f64 imbalance) {
            const weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            Matching matching;

            {
                ScopedTimer _t("coarsening", "ApproximateGreedyEdgeMatcher", "allocate");
                matching.initialize(g.n);
            }

            mark += 1;
            edges_size = 0;

            // first handle pendant vertices
            if (config.match_pendant_vertices_first) {
                ScopedTimer _t("coarsening", "ApproximateGreedyEdgeMatcher", "match_pendant");

                forall_gu(g, u)
                    {
                        if (g.deg(u) != 1) { continue; }

                        const vertex_t v = g.edges_v[g.neighborhoods[u]];
                        const weight_t ew = g.edges_w[g.neighborhoods[u]];

                        if (u >= v) { continue; }
                        if (p_manager[u] != p_manager[v]) { continue; }
                        if (g.v_weights[u] + g.v_weights[v] > lmax) { continue; }

                        const f32 rating = (f32) ew / (f32) (g.deg(u) * g.deg(v));
                        edges[edges_size++] = {u, v, rating};
                    }
                endfor

                std::sort(edges.get_ptr(), edges.get_ptr() + edges_size, std::greater<>());

                for (size_t i = 0; i < edges_size; ++i) {
                    const auto &[u, v, w] = edges[i];
                    if (used[u] == mark || used[v] == mark) { continue; }

                    used[u] = mark;
                    used[v] = mark;
                    matching.add(u, v);
                }

                edges_size = 0;
            }

            if (config.match_triangles) {
                ScopedTimer _t("coarsening", "ApproximateGreedyEdgeMatcher", "match_triangles");

                forall_gu(g, u)
                    {
                        if (used[u] == mark) { continue; }
                        if (g.deg(u) != 2) { continue; }

                        const vertex_t u0 = g.neighborhoods[u];
                        const vertex_t u1 = u0 + 1;

                        const vertex_t a = g.edges_v[u0];
                        const vertex_t b = g.edges_v[u1];
                        const weight_t wa = g.edges_w[u0];
                        const weight_t wb = g.edges_w[u1];

                        if (p_manager[u] != p_manager[a] || p_manager[u] != p_manager[b]) { continue; }

                        // Try matching u with a if a is also degree 2 and shares b.
                        if (g.deg(a) == 2 && used[a] != mark) {
                            const vertex_t a0 = g.neighborhoods[a];
                            const vertex_t a1 = a0 + 1;

                            const vertex_t a_n0 = g.edges_v[a0];
                            const vertex_t a_n1 = g.edges_v[a1];

                            if ((a_n0 == u && a_n1 == b) || (a_n1 == u && a_n0 == b)) {
                                if (u < a && g.v_weights[u] + g.v_weights[a] <= lmax) {
                                    const f32 rating = (f32) wa / (f32) (g.deg(u) * g.deg(a));
                                    edges[edges_size++] = {u, a, rating};
                                }
                            }
                        }

                        // Try matching u with b if b is also degree 2 and shares a.
                        if (g.deg(b) == 2 && used[b] != mark) {
                            const vertex_t b0 = g.neighborhoods[b];
                            const vertex_t b1 = b0 + 1;

                            const vertex_t b_n0 = g.edges_v[b0];
                            const vertex_t b_n1 = g.edges_v[b1];

                            if ((b_n0 == u && b_n1 == a) || (b_n1 == u && b_n0 == a)) {
                                if (u < b && g.v_weights[u] + g.v_weights[b] <= lmax) {
                                    const f32 rating = (f32) wb / (f32) (g.deg(u) * g.deg(b));
                                    edges[edges_size++] = {u, b, rating};
                                }
                            }
                        }
                    }
                endfor

                std::sort(edges.get_ptr(), edges.get_ptr() + edges_size, std::greater<>());

                for (size_t i = 0; i < edges_size; ++i) {
                    const auto &[u, v, w] = edges[i];
                    if (used[u] == mark || used[v] == mark) { continue; }

                    used[u] = mark;
                    used[v] = mark;
                    matching.add(u, v);
                }

                edges_size = 0;
            }

            // handle all other vertices
            {
                ScopedTimer _t("coarsening", "ApproximateGreedyEdgeMatcher", "rate_edges");

                forall_gu(g, u)
                    {
                        if (used[u] == mark) { continue; }

                        const weight_t u_w = g.v_weights[u];
                        const size_t deg_u = g.deg(u);

                        const size_t begin = g.neighborhoods[u];
                        const size_t end = g.neighborhoods[u + 1];

                        if (!config.sample_edges || deg_u <= config.sample_degree_threshold ||
                            config.max_samples_per_vertex >= deg_u) {
                            // full scan
                            for (size_t j = begin; j < end; ++j) {
                                const vertex_t v = g.edges_v[j];
                                const weight_t w = g.edges_w[j];

                                if (u >= v) { continue; }
                                if (used[v] == mark) { continue; }
                                if (p_manager[u] != p_manager[v]) { continue; }
                                if (u_w + g.v_weights[v] > lmax) { continue; }

                                const f32 edge_rating = (f32) w / (f32) (u_w * g.v_weights[v]);
                                edges[edges_size++] = {u, v, edge_rating};
                            }
                        } else {
                            // sample a fixed number of neighbors
                            const size_t samples = config.max_samples_per_vertex;
                            for (size_t s = 0; s < samples; ++s) {
                                const size_t j = begin + (random_engine.get_u64() % deg_u);

                                const vertex_t v = g.edges_v[j];
                                const weight_t w = g.edges_w[j];

                                if (u >= v) { continue; }
                                if (used[v] == mark) { continue; }
                                if (p_manager[u] != p_manager[v]) { continue; }
                                if (u_w + g.v_weights[v] > lmax) { continue; }

                                const f32 edge_rating = (f32) w / (f32) (u_w * g.v_weights[v]);
                                edges[edges_size++] = {u, v, edge_rating};
                            }
                        }
                    }
                endfor
            }

            {
                ScopedTimer _t("coarsening", "ApproximateGreedyEdgeMatcher", "approx_sort");

                approximate_sort_edges(edges_size);
            }

            {
                ScopedTimer _t("coarsening", "ApproximateGreedyEdgeMatcher", "choose");

                for (size_t i = 0; i < edges_size; ++i) {
                    const auto &[u, v, w] = edges[i];
                    if (used[u] != mark && used[v] != mark) {
                        used[u] = mark;
                        used[v] = mark;
                        matching.add(u, v);
                    }
                }
            }

            {
                ScopedTimer _t("coarsening", "ApproximateGreedyEdgeMatcher", "mapping");

                matching.set_translation();
                mapping.set_coarse_n(matching.get_n_coarse_nodes());
                for (vertex_t u = 0; u < matching.get_n(); ++u) {
                    mapping.set(u, matching.get_n(u));
                }
            }
        }
    };
}

#endif //HEIPROMAP_APPROXIMATE_GREEDY_EDGE_MATCHER_H
