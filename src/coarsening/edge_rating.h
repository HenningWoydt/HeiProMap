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

#ifndef HEIPROMAP_EDGE_RATING_H
#define HEIPROMAP_EDGE_RATING_H

#include "../definitions.h"
#include "../datastructures/csr_graph.h"

namespace HeiProMap {

    template<bool t_uniform_v_weights, bool t_uniform_e_weights, EdgeRatingFunction t_rating_function>
    [[nodiscard]] inline __attribute__((always_inline)) f32 compute_edge_rating(const graph_t &g,
                                                                                const vertex_t u,
                                                                                const vertex_t v,
                                                                                const weight_t u_w,
                                                                                const weight_t v_w,
                                                                                const weight_t ew,
                                                                                const f32 alpha = 0.f,
                                                                                const f32 beta = 0.f) noexcept {
        if constexpr (t_rating_function == EdgeRatingFunction::WEIGHT) {
            return static_cast<f32>(ew);
        } else if constexpr (t_rating_function == EdgeRatingFunction::EXPANSION) {
            return static_cast<f32>(ew) / static_cast<f32>(u_w + v_w);
        } else if constexpr (t_rating_function == EdgeRatingFunction::EXPANSIONSTAR) {
            return static_cast<f32>(ew) / static_cast<f32>(u_w * v_w);
        } else if constexpr (t_rating_function == EdgeRatingFunction::EXPANSIONSTARSTAR) {
            return static_cast<f32>(ew * ew) / static_cast<f32>(u_w * v_w);
        } else if constexpr (t_rating_function == EdgeRatingFunction::INNEROUTER) {
            weight_t out_v = 0;
            for (u64 i_edge = g.neighborhoods[v]; i_edge < g.neighborhoods[v + 1]; ++i_edge) {
                out_v += t_uniform_e_weights ? 1 : g.edges_w[i_edge];
            }
            weight_t out_u = 0;
            for (u64 i_edge = g.neighborhoods[u]; i_edge < g.neighborhoods[u + 1]; ++i_edge) {
                out_u += t_uniform_e_weights ? 1 : g.edges_w[i_edge];
            }
            return static_cast<f32>(ew) / static_cast<f32>(out_v + out_u - (2 * ew));
        } else if constexpr (t_rating_function == EdgeRatingFunction::MESSAGE ||
                             t_rating_function == EdgeRatingFunction::MESSAGESUM) {
            // Contracting (u, v) ends the edge between them and merges the edges to every common neighbor, so the
            // coarse vertex has that many neighbors (messages) fewer; the volume of the edge is not cut any more.
            u64 common = 0;
            for (u64 i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t x = g.edges_v[i];
                if (x == v) { continue; }
                for (u64 j = g.neighborhoods[v]; j < g.neighborhoods[v + 1]; ++j) {
                    if (g.edges_v[j] == x) {
                        ++common;
                        break;
                    }
                }
            }
            const f32 saved = beta * static_cast<f32>(ew) + alpha * static_cast<f32>(1 + common);
            if constexpr (t_rating_function == EdgeRatingFunction::MESSAGE) {
                return saved / static_cast<f32>(u_w * v_w);
            } else {
                return saved / static_cast<f32>(u_w + v_w);
            }
        }
    }

}

#endif //HEIPROMAP_EDGE_RATING_H
