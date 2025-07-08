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

#ifndef HEIPROMAP_RANDOM_EDGE_MATCHER_H
#define HEIPROMAP_RANDOM_EDGE_MATCHER_H

#include <vector>

#include "../../commons/definitions.h"
#include "../../commons/random_engine.h"
#include "../../commons/statistic_collector.h"

namespace HeiProMap {
    class RandomEdgeMatcherConfiguration {
    public:
    };

    class RandomEdgeMatcher {
        vertex_t    m_n     = 0;
        vertex_t    m_m     = 0;
        partition_t m_k     = 0;
        weight_t    m_l_max = 0;

        const RandomEdgeMatcherConfiguration *config        = nullptr;
        RandomEngine                         *random_engine = nullptr;

        u32               mark = 0;
        AlignedArray<u32> used;

    public:
        RandomEdgeMatcher() = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const weight_t t_l_max,
                        RandomEngine &t_random_engine,
                        const RandomEdgeMatcherConfiguration &i_config) {
            m_n     = t_n;
            m_m     = t_m;
            m_k     = t_k;
            m_l_max = t_l_max;

            config        = dynamic_cast<const RandomEdgeMatcherConfiguration *>(&i_config);
            random_engine = &t_random_engine;

            mark = 0;
            used.initialize(m_n, 0);
        }

        template<typename PartitionManagerT>
        void match(const size_t level,
                   const graph_t &g,
                   PartitionManagerT &p_manager,
                   Matching &matching) {
            mark += 1;

            forall_gu(g, u)
                {
                    if (used[u] == mark) { continue; }

                    weight_t u_w = g.weight(u);

                    f32      counter = 0;
                    vertex_t chosen_v;

                    forall_guiv(g, u, j, v)
                        {
                            if (used[v] == mark) { continue; }
                            if (p_manager[u] != p_manager[v]) { continue; }
                            weight_t v_w = g.weight(v);

                            if (u_w + v_w > m_l_max) { continue; }

                            counter += 1.0;
                            // choose with probability 1/counter as it ensures uniform distribution
                            if (random_engine->get_f32() <= 1.0f / counter) {
                                chosen_v = v;
                            }
                        }
                    endfor
                    if (counter > 0) {
                        used[u]        = mark;
                        used[chosen_v] = mark;

                        matching.add(u, chosen_v);
                    }
                }
            endfor

            /*
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
             */
        }

        JSONString get_stats() {
            return {};
        }
    };
}

#endif //HEIPROMAP_RANDOM_EDGE_MATCHER_H
