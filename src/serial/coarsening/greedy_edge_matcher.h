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
    class GreedyEdgeMatcher final : public ISerialMatcher {
        u32              mark = 0;
        std::vector<u32> used;

        weight_t l_max = 0;

    public:
        GreedyEdgeMatcher() = default;

        void initialize(const size_t n, const weight_t l_max) override {
            mark = 0;
            used.resize(n, 0);

            this->l_max = l_max;
        }

        template<typename TSerialGraph, typename TSerialActiveVertexManager>
        void match(TSerialGraph &g,
                   TSerialActiveVertexManager &av_manager,
                   std::vector<EdgeUV> &matches) {
            std::vector<EdgeUVW> edges;
            edges.reserve(g.get_m());
            for (vertex_t u: av_manager) {
                ASSERT(av_manager.is_active(u));

                for (size_t i = 0; i < g.size(u); ++i) {
                    const vertex_t v      = g.neighbor(u, i);
                    const weight_t ew     = g.get_weight(u, i);
                    const f64      rating = (f64) ew / (f64) (g.size(u) * g.size(v));
                    edges.emplace_back(u, v, rating);
                }
            }
            std::sort(edges.begin(), edges.end(), std::greater<>());

            mark += 1;
            matches.clear();
            for (const auto &e: edges) {
                if (used[e.u] != mark && used[e.v] != mark) {
                    used[e.u] = mark;
                    used[e.v] = mark;
                    if (g.size(e.u) > g.size(e.v)) {
                        matches.emplace_back(e.u, e.v);
                    } else {
                        matches.emplace_back(e.v, e.u);
                    }
                }
            }

#if ASSERT_ENABLED
            for (const EdgeUV &e: matches) {
                ASSERT(e.u != e.v);
                ASSERT(av_manager.is_active(e.u));
                ASSERT(av_manager.is_active(e.v));
            }
#endif

#if ASSERT_ENABLED
            std::vector<u8> hit(g.get_n(), 0);
            for (auto &e: matches) {
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

#endif //HEIPROMAP_GREEDY_EDGE_MATCHER_H
