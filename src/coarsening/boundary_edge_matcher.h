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

#ifndef HEIPROMAP_BOUNDARY_EDGE_MATCHER_H
#define HEIPROMAP_BOUNDARY_EDGE_MATCHER_H

#include <vector>

#include "../definitions.h"
#include "../utility/random_engine.h"

namespace HeiProMap {
    class BoundaryEdgeMatcherConfiguration {
    public:
    };

    class BoundaryEdgeMatcher {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;

        const BoundaryEdgeMatcherConfiguration *config = nullptr;
        RandomEngine *random_engine = nullptr;

        u32 mark = 0;
        AlignedArray<u32> used;

    public:
        BoundaryEdgeMatcher() = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        RandomEngine &t_random_engine,
                        const BoundaryEdgeMatcherConfiguration &i_config) {
            ScopedTimer _t("coarsening", "RandomEdgeMatcher", "initialize");

            m_n = t_n;
            m_m = t_m;
            m_k = t_k;

            config = dynamic_cast<const BoundaryEdgeMatcherConfiguration *>(&i_config);
            random_engine = &t_random_engine;

            mark = 0;
            used.initialize(m_n, 0);
        }

        void match([[maybe_unused]] const size_t level,
                   const graph_t &g,
                   PartitionManager &p_manager,
                   BoundaryVertexManager &bv_manager,
                   Mapping &mapping,
                   f64 imbalance) {
            ScopedTimer _t("coarsening", "BoundaryEdgeMatcher", "match");

            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            Matching matching;
            matching.initialize(g.n);

            mark += 1;
            u64 n_matched = 0;

            for (partition_t id = 0; id < bv_manager.get_k(); ++id) {
                for (size_t i = 0; i < bv_manager.size(id); ++i) {
                    vertex_t u = bv_manager.get(id, i);
                    weight_t u_w = g.v_weights[u];

                    if (used[u] == mark) { continue; }

                    f32 counter = 0;
                    vertex_t chosen_v = 0;

                    forall_guiv(g, u, j, v)
                        {
                            if (used[v] == mark) { continue; }
                            if (p_manager[u] != p_manager[v]) { continue; }
                            weight_t v_w = g.v_weights[v];

                            if (u_w + v_w > lmax) { continue; }

                            counter += 1.0;
                            // choose with probability 1/counter as it ensures uniform distribution
                            if (random_engine->get_f32() <= 1.0f / counter) {
                                chosen_v = v;
                            }
                        }
                    endfor
                    if (counter > 0) {
                        used[u] = mark;
                        used[chosen_v] = mark;

                        matching.add(u, chosen_v);
                        n_matched += 2;
                    }
                }
            }

            /*
            // 3) Third: all remaining vertices
            forall_gu(g, u)
                {
                    weight_t u_w = g.v_weights[u];

                    if (used[u] == mark) { continue; }

                    f32 counter = 0;
                    vertex_t chosen_v = 0;

                    forall_guiv(g, u, j, v)
                        {
                            if (used[v] == mark) { continue; }
                            if (p_manager[u] != p_manager[v]) { continue; }
                            weight_t v_w = g.v_weights[v];

                            if (u_w + v_w > lmax) { continue; }

                            counter += 1.0;
                            // choose with probability 1/counter as it ensures uniform distribution
                            if (random_engine->get_f32() <= 1.0f / counter) {
                                chosen_v = v;
                            }
                        }
                    endfor
                    if (counter > 0) {
                        used[u] = mark;
                        used[chosen_v] = mark;

                        matching.add(u, chosen_v);
                        n_matched += 2;
                    }
                }
            endfor
            */

            matching.set_translation();
            mapping.set_coarse_n(matching.get_n_coarse_nodes());
            for (vertex_t u = 0; u < matching.get_n(); ++u) {
                mapping.set(u, matching.get_n(u));
            }
        }
    };
}

#endif //HEIPROMAP_BOUNDARY_EDGE_MATCHER_H
