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
#include "../interfaces/ISerialMatcher.h"

namespace HeiProMap {
    class HeavyEdgeMatcher final : public ISerialMatcher {
        u32 mark = 0;
        std::vector<u32> used;

        weight_t l_max = 0;

    public:
        HeavyEdgeMatcher() = default;

        void initialize(const vertex_t n, const vertex_t m, const weight_t t_l_max) override {
            mark = 0;
            used.resize(n, mark);
            this->l_max = l_max;
        }

        template <typename TSerialGraph, typename TSerialActiveVertexManager>
        void match(TSerialGraph& g,
                   TSerialActiveVertexManager& av_manager,
                   std::vector<EdgeUV>& matches) {
            mark += 1;
            matches.clear();

            // first check vertices with degree 1
            for (vertex_t u : av_manager) {
                ASSERT(av_manager.is_active(u));

                if (used[u] != mark && g.size(u) == 1) {
                    vertex_t v = g.neighbor(u, 0);

                    if (used[v] != mark && g.get_weight(u) + g.get_weight(v) <= l_max) {
                        used[u] = mark;
                        used[v] = mark;

                        matches.emplace_back(v, u); // pull u into v
                    }
                }
            }

            // check all other vertices
            for (vertex_t u : av_manager) {
                ASSERT(av_manager.is_active(u));

                if (used[u] != mark) {
                    size_t best_idx;
                    weight_t max_weight = 0;

                    for (size_t i = 0; i < g.size(u); ++i) {
                        vertex_t v  = g.neighbor(u, i);
                        weight_t ew = g.get_weight(u, i);
                        ASSERT(u != v);
                        ASSERT(av_manager.is_active(v));
                        if (used[v] != mark && g.get_weight(u) + g.get_weight(v) <= l_max) {
                            if (ew > max_weight) {
                                best_idx   = i;
                                max_weight = ew;
                            }
                        }
                    }

                    if (max_weight != 0) {
                        vertex_t v = g.neighbor(u, best_idx);
                        used[u]    = mark;
                        used[v]    = mark;

                        if (g.size(u) > g.size(v)) {
                            matches.emplace_back(u, v);
                        } else {
                            matches.emplace_back(v, u);
                        }
                    }
                }
            }

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

                if (hit[e.u] == 2) {
                    ASSERT(false);
                }
                if (hit[e.v] == 2) {
                    ASSERT(false);
                }
            }
#endif
        }
    };
}

#endif //HEIPROMAP_HEAVY_EDGE_MATCHER_H
