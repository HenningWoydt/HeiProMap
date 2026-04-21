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

#ifndef HEIPROMAP_HEAVY_EDGE_MATCHER_H
#define HEIPROMAP_HEAVY_EDGE_MATCHER_H

#include <vector>

#include "../definitions.h"
#include "../utility/random_engine.h"

namespace HeiProMap {
    f32 edge_noise(vertex_t u, vertex_t v, u32 round) {
        // order independent
        const u32 a = u < v ? u : v;
        const u32 b = u < v ? v : u;

        // simple 32-bit mix (Wang-style)
        u32 x = a * 0x9e3779b1u;
        x ^= b * 0x85ebca77u;
        x ^= round * 0xc2b2ae3du;

        x ^= x >> 16;
        x *= 0x7feb352du;
        x ^= x >> 15;

        // map to tiny float in [0, 1e-6)
        return (f32) (x & 0x00ffffffu) * (1.0f / 16777216.0f) * 1e-6f;
    }

    class HeavyEdgeMatcherConfiguration {
    public:
    };

    class HeavyEdgeMatcher {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;

        const HeavyEdgeMatcherConfiguration *config = nullptr;
        RandomEngine *random_engine = nullptr;

        AlignedArray<vertex_t> preferred;
        AlignedArray<vertex_t> vertex_list;
        size_t vertex_list_size = 0;
        AlignedArray<vertex_t> next_vertex_list;
        size_t next_vertex_list_size = 0;

    public:
        HeavyEdgeMatcher() = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        RandomEngine &t_random_engine,
                        const HeavyEdgeMatcherConfiguration &i_config) {
            ScopedTimer _t("coarsening", "HeavyEdgeMatcher", "initialize");

            m_n = t_n;
            m_m = t_m;
            m_k = t_k;

            config = dynamic_cast<const HeavyEdgeMatcherConfiguration *>(&i_config);
            random_engine = &t_random_engine;

            preferred.initialize(m_n);
            vertex_list.initialize(m_n);
            next_vertex_list.initialize(m_n);
        }

        void match([[maybe_unused]] const size_t level,
                   const graph_t &g,
                   [[maybe_unused]] p_manager_t &p_manager,
                   Mapping &mapping,
                   f64 imbalance) {
            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            Matching matching;
            //
            {
                ScopedTimer _t("coarsening", "HeavyEdgeMatcher", "initialize");

                matching.initialize(g.n);

                std::iota(preferred.get_ptr(), preferred.get_ptr() + g.n, 0);
                std::iota(vertex_list.get_ptr(), vertex_list.get_ptr() + g.n, 0);
                vertex_list_size = g.n;
            }

            size_t round = 0;
            size_t n_matched = 0;
            size_t n_new_matched = 1;
            while (n_new_matched > 0) {
                //
                {
                    ScopedTimer _t("coarsening", "HeavyEdgeMatcher", "preferred_neighbor");

                    for (size_t i = 0; i < vertex_list_size; ++i) {
                        vertex_t u = vertex_list[i];

                        weight_t u_w = g.v_weights[u];
                        vertex_t best_v = u;
                        f32 best_rating = 0;

                        for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) { const vertex_t v = g.edges_v[j]; const weight_t w = g.edges_w[j];
                            {
                                if (matching.is_matched(v)) { continue; }

                                weight_t v_w = g.v_weights[v];

                                if (u_w + v_w > lmax) { continue; }

                                // f32 edge_rating = ((f32) (w * w)) / ((f32) (g.deg(u) * g.deg(v)));
                                f32 edge_rating = ((f32) (w * w)) / ((f32) (u_w * v_w));
                                // std::cout << edge_rating << " " << edge_noise(u, v, round) << std::endl;
                                edge_rating += edge_noise(u, v, round);

                                if (edge_rating > best_rating) {
                                    best_v = v;
                                    best_rating = edge_rating;
                                }
                            }
                        }

                        preferred[u] = best_v;
                    }
                }

                n_new_matched = 0;
                next_vertex_list_size = 0;
                //
                {
                    ScopedTimer _t("coarsening", "HeavyEdgeMatcher", "match_preferred");

                    for (size_t i = 0; i < vertex_list_size; ++i) {
                        vertex_t u = vertex_list[i];
                        vertex_t v = preferred[u];
                        vertex_t pref_v = preferred[v];

                        if (u < v && u == pref_v) {
                            matching.add(u, v);
                            n_new_matched += 2;
                        } else if (u == v) {
                            next_vertex_list[next_vertex_list_size++] = u;
                        }
                    }
                }
                n_matched += n_new_matched;

                std::swap(next_vertex_list, vertex_list);
                std::swap(next_vertex_list_size, vertex_list_size);
                round += 1;
            }

            {
                ScopedTimer _t("coarsening", "HeavyEdgeMatcher", "mapping");

                matching.set_translation();
                mapping.set_coarse_n(matching.get_n_coarse_nodes());
                for (vertex_t u = 0; u < matching.get_n(); ++u) {
                    mapping.set(u, matching.get_n(u));
                }
            }
        }
    };
}

#endif //HEIPROMAP_HEAVY_EDGE_MATCHER_H
