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

#include "../../commons/definitions.h"
#include "../../commons/random_engine.h"
#include "../../commons/statistic_collector.h"
#include "../interfaces/ISerialMatcher.h"

namespace HeiProMap {
    class GreedyEdgeMatcherConfiguration final : public ISerialMatcherConfiguration {
    public:
        bool match_pendant_vertices_first = false; // Vertices with only one neighbor should be handled first.
    };

    class GreedyEdgeMatcher final : public ISerialMatcher {
        vertex_t    m_n     = 0;
        vertex_t    m_m     = 0;
        partition_t m_k     = 0;
        weight_t    m_l_max = 0;

        const GreedyEdgeMatcherConfiguration *config           = nullptr;
        RandomEngine                         *random_engine    = nullptr;
        StatisticCollector                   *m_stat_collector = nullptr;

        std::vector<EdgeUVW> edges;

        u32              mark = 0;
        std::vector<u32> used;

    public:
        GreedyEdgeMatcher() = default;

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

            config           = dynamic_cast<const GreedyEdgeMatcherConfiguration *>(&i_config);
            random_engine    = &t_random_engine;
            m_stat_collector = &t_stat_collect;

            edges.reserve(m_m);

            mark = 0;
            used.resize(m_n, 0);
        }

        void match(const size_t level,
                   const graph_t &g,
                   p_manager_t& p_manager,
                   Matching &matching) override {
            mark += 1;
            edges.clear();

            // first handle pendant vertices
            if (config->match_pendant_vertices_first) {
                forall_gu(g, u)
                    {
                        if (g.size(u) != 1) {
                            continue;
                        }

                        const vertex_t v      = g.neighbor(u, 0);
                        const weight_t ew     = g.weight(u, 0);
                        const f64      rating = (f64) ew / (f64) (g.size(u) * g.size(v));
                        edges.emplace_back(u, v, rating);
                    }
                endfor
                std::sort(edges.begin(), edges.end(), std::greater<>());

                for (const auto &[u, v, w]: edges) {
                    if (used[u] == mark || used[v] == mark) {
                        continue;
                    }

                    if (g.weight(u) + g.weight(v) > m_l_max) {
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
                    weight_t u_w = g.weight(u);

                    if (used[u] == mark) {
                        continue;
                    }

                    forall_guivw(g, u, j, v, w)
                        {
                            weight_t v_w = g.weight(v);

                            if (used[v] == mark) {
                                continue;
                            }

                            if (g.weight(u) + g.weight(v) > m_l_max) {
                                continue;
                            }

                            const f32 edge_rating = (f32) w / (f32) (u_w * v_w);
                            edges.emplace_back(u, v, edge_rating);
                        }
                    endfor
                }
            endfor
            std::sort(edges.begin(), edges.end(), std::greater<>());

            for (const auto &[u, v, w]: edges) {
                if (used[u] != mark && used[v] != mark) {
                    // use this edge
                    used[u] = mark;
                    used[v] = mark;

                    matching.add(u, v);
                }
            }

#if ASSERT_ENABLED
            for (size_t i = 0; i < matching.size(); ++i) {
                const auto &[u, v] = matching[i];
                ASSERT(u != v);
            }
#endif

#if ASSERT_ENABLED
            std::vector<u8> hit(g.get_n(), 0);
            for (size_t     i = 0; i < matching.size(); ++i) {
                const auto &[u, v] = matching[i];
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

        JSONString get_stats() override {
            return {};
        }
    };
}

#endif //HEIPROMAP_GREEDY_EDGE_MATCHER_H
