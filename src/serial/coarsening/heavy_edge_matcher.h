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

#include "../../definitions.h"
#include "../../macros.h"
#include "../../commons/JSON_utils.h"
#include "../../commons/random_engine.h"
#include "../../commons/statistic_collector.h"
#include "../interfaces/ISerialMatcher.h"

namespace HeiProMap {
    class HeavyEdgeMatcherConfiguration final : public ISerialMatcherConfiguration {
    public:
        bool match_pendant_vertices_first = false; // Vertices with only one neighbor should be handled first.
    };

    class HeavyEdgeMatcher final : public ISerialMatcher {
        vertex_t m_n     = 0;
        vertex_t m_m     = 0;
        partition_t m_k  = 0;
        weight_t m_l_max = 0;

        const HeavyEdgeMatcherConfiguration* config = nullptr;
        RandomEngine* random_engine                 = nullptr;
        StatisticCollector* m_stat_collector        = nullptr;

        u32 mark = 0;
        std::vector<u32> used;

    public:
        HeavyEdgeMatcher() = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const weight_t t_l_max,
                        RandomEngine& t_random_engine,
                        const ISerialMatcherConfiguration& i_config,
                        StatisticCollector& t_stat_collect) override {
            m_n     = t_n;
            m_m     = t_m;
            m_k     = t_k;
            m_l_max = t_l_max;

            config           = dynamic_cast<const HeavyEdgeMatcherConfiguration*>(&i_config);
            random_engine    = &t_random_engine;
            m_stat_collector = &t_stat_collect;

            mark = 0;
            used.resize(m_n, mark);
        }

        void match(const size_t level,
                   const graph_t& g,
                   Matching& matching) override {
            mark += 1;

            if (config->match_pendant_vertices_first) {
                // first check vertices with degree 1
                forall_gu(g, u)
                    {
                        if (used[u] == mark) { continue; }
                        if (g.size(u) != 1) { continue; }

                        vertex_t v = g.neighbor(u, 0);

                        if (used[v] == mark) { continue; }

                        weight_t u_w = g.weight(u);
                        weight_t v_w = g.weight(v);

                        if (u_w + v_w > m_l_max) { continue; }

                        matching.add(u, v);
                    }
                endfor
            }

            // check all other vertices
            forall_gu(g, u)
                {
                    if (used[u] == mark) { continue; }

                    weight_t u_w        = g.weight(u);
                    vertex_t best_v     = u;
                    weight_t max_weight = 0;

                    forall_guivw(g, u, j, v, w)
                        {
                            if (used[v] == mark) { continue; }

                            weight_t v_w = g.weight(v);

                            if (u_w + v_w > m_l_max) { continue; }

                            if (w > max_weight) {
                                best_v     = v;
                                max_weight = w;
                            }
                        }
                    endfor

                    if (best_v != u) {
                        used[u]      = mark;
                        used[best_v] = mark;

                        matching.add(u, best_v);
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

#endif //HEIPROMAP_HEAVY_EDGE_MATCHER_H
