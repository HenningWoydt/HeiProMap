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
    };

    class GreedyEdgeMatcher {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;

        const GreedyEdgeMatcherConfiguration *config = nullptr;
        RandomEngine *random_engine = nullptr;

        AlignedArray<EdgeUVW> edges;
        size_t edges_size = 0;

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

            edges.initialize(m_m);

            mark = 0;
            used.initialize(m_n, 0);
        }

        void match([[maybe_unused]] const size_t level,
                   const graph_t &g,
                   [[maybe_unused]] p_manager_t &p_manager,
                   Mapping &mapping,
                   f64 imbalance) {
            ScopedTimer _t("coarsening", "GreedyEdgeMatcher", "match");

            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            Matching matching;
            matching.initialize(g.n);

            mark += 1;
            edges_size = 0;

            // first handle pendant vertices
            if (config->match_pendant_vertices_first) {
                forall_gu(g, u)
                    {
                        if (g.deg(u) != 1) {
                            continue;
                        }

                        const vertex_t v = g.edges_v[g.neighborhoods[u]];
                        const weight_t ew = g.edges_w[g.neighborhoods[u]];
                        const f32 rating = (f32) ew / (f32) (g.deg(u) * g.deg(v));
                        edges[edges_size++] = {u, v, rating};
                    }
                endfor
                std::sort(edges.get_ptr(), edges.get_ptr() + edges_size, std::greater<>());

                for (size_t i = 0; i < edges_size; ++i) {
                    const auto &[u, v, w] = edges[i];
                    if (used[u] == mark || used[v] == mark) {
                        continue;
                    }

                    if (g.v_weights[u] + g.v_weights[v] > lmax) {
                        continue;
                    }

                    if (used[u] != mark && used[v] != mark) {
                        // use this edge
                        used[u] = mark;
                        used[v] = mark;

                        matching.add(u, v);
                    }
                }
            }

            // handle all other vertices
            forall_gu(g, u)
                {
                    weight_t u_w = g.v_weights[u];

                    if (used[u] == mark) {
                        continue;
                    }

                    forall_guivw(g, u, j, v, w)
                        {
                            weight_t v_w = g.v_weights[v];

                            if (used[v] == mark) {
                                continue;
                            }

                            if (g.v_weights[u] + g.v_weights[v] > lmax) {
                                continue;
                            }

                            const f32 edge_rating = (f32) w / (f32) (u_w * v_w);
                            edges[edges_size++] = {u, v, edge_rating};
                        }
                    endfor
                }
            endfor
            std::sort(edges.get_ptr(), edges.get_ptr() + edges_size, std::greater<>());

            for (size_t i = 0; i < edges_size; ++i) {
                const auto &[u, v, w] = edges[i];
                if (used[u] != mark && used[v] != mark) {
                    // use this edge
                    used[u] = mark;
                    used[v] = mark;

                    matching.add(u, v);
                }
            }

            matching.set_translation();
            mapping.set_coarse_n(matching.get_n_coarse_nodes());
            for (vertex_t u = 0; u < matching.get_n(); ++u) {
                mapping.set(u, matching.get_n(u));
            }
        }
    };
}

#endif //HEIPROMAP_GREEDY_EDGE_MATCHER_H
