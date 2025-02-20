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

#include "../../definitions.h"
#include "../../macros.h"
#include "../interfaces/ISerialMatcher.h"

namespace HeiProMap {
    struct GreedyEdgeMatcherConfiguration {
        bool match_pendant_vertices_first = false; // Vertices with only one neighbor should be handled first.
        bool no_overload                  = false; // Matching an edge, should not create a vertex with a weight greater l_max.
    };

    class GreedyEdgeMatcher final : public ISerialMatcher {
        u32 mark = 0;
        std::vector<u32> used;

        weight_t l_max = 0;

    public:
        GreedyEdgeMatcher() = default;

        void initialize(const vertex_t n, const vertex_t m, const partition_t k, const weight_t t_l_max) override {
            mark = 0;
            used.resize(n, 0);

            this->l_max = l_max;
        }

        template <typename TSerialGraph, typename TSerialActiveVertexManager>
        void match(GreedyEdgeMatcherConfiguration& config,
                   TSerialGraph& g,
                   TSerialActiveVertexManager& av_manager,
                   std::vector<EdgeUV>& matches) {
            mark += 1;
            matches.clear();

            std::vector<EdgeUVW> edges;
            edges.reserve(g.get_m());

            // first handle pendant vertices
            if (config.match_pendant_vertices_first) {
                for (vertex_t u : av_manager) {
                    ASSERT(av_manager.is_active(u));

                    if (g.size(u) != 1) {
                        continue;
                    }

                    const vertex_t v  = g.neighbor(u, 0);
                    const weight_t ew = g.get_weight(u, 0);
                    const f64 rating  = (f64)ew / (f64)(g.size(u) * g.size(v));
                    edges.emplace_back(u, v, rating);
                }
                std::sort(edges.begin(), edges.end(), std::greater<>());

                for (const auto& [u, v, w] : edges) {
                    if (used[u] == mark || used[v] == mark) {
                        continue;
                    }

                    if (config.no_overload && g.get_weight(u) + g.get_weight(v) > l_max) {
                        continue;
                    }

                    if (used[u] != mark && used[v] != mark) {
                        // use this edge
                        used[u] = mark;
                        used[v] = mark;
                        if (g.size(u) > g.size(v)) {
                            matches.emplace_back(u, v);
                        } else {
                            matches.emplace_back(v, u);
                        }
                    }
                }
            }

            // handle all other vertices
            for (vertex_t u : av_manager) {
                ASSERT(av_manager.is_active(u));

                if (used[u] == mark) {
                    continue;
                }

                for (const auto& [v, w] : g[u]) {
                    if (used[v] == mark) {
                        continue;
                    }

                    if (config.no_overload && g.get_weight(u) + g.get_weight(v) > l_max) {
                        continue;
                    }

                    const f64 rating = (f64)w / (f64)(g.size(u) * g.size(v));
                    edges.emplace_back(u, v, rating);
                }
            }
            std::sort(edges.begin(), edges.end(), std::greater<>());

            for (const auto& [u, v, w] : edges) {
                if (used[u] != mark && used[v] != mark) {
                    // use this edge
                    used[u] = mark;
                    used[v] = mark;
                    if (g.size(u) > g.size(v)) {
                        matches.emplace_back(u, v);
                    } else {
                        matches.emplace_back(v, u);
                    }
                }
            }

#if ASSERT_ENABLED
            for (const auto& [u, v] : matches) {
                ASSERT(u != v);
                ASSERT(av_manager.is_active(u));
                ASSERT(av_manager.is_active(v));
            }
#endif

#if ASSERT_ENABLED
            std::vector<u8> hit(g.get_n(), 0);
            for (auto& [u, v] : matches) {
                hit[u] += 1;
                hit[v] += 1;

                if (hit[u] == 2) {
                    ASSERT(false);
                }
                if (hit[v] == 2) {
                    ASSERT(false);
                }
            }
#endif
        }
    };
}

#endif //HEIPROMAP_GREEDY_EDGE_MATCHER_H
