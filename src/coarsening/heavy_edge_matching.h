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

#ifndef HEIPROMAP_HEAVY_EDGE_MATCHING_H
#define HEIPROMAP_HEAVY_EDGE_MATCHING_H

#include <algorithm>
#include <vector>
#include <numeric>
#include <cmath>

#include "../definitions.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../utility/matching.h"
#include "../utility/random_engine.h"
#include "../utility/mapping.h"
#include "../utility/profiler.h"

namespace HeiProMap {

    class HeavyEdgeMatchingConfiguration {
    public:
        EdgeRatingFunction rating_function = EdgeRatingFunction::EXPANSIONSTAR;
    };

    class HeavyEdgeMatching {
    public:
        void match(const graph_t &g,
                   const p_manager_t &p_manager,
                   Mapping &mapping,
                   f64 imbalance,
                   u64 seed,
                   const HeavyEdgeMatchingConfiguration &config) {
            auto dispatch_with_rating = [&](auto rating_func_const) {
                constexpr EdgeRatingFunction rating_func = rating_func_const;
                if (g.uniform_v_weights && g.uniform_e_weights) {
                    match_templated<true, true, rating_func>(g, p_manager, mapping, imbalance, seed);
                } else if (g.uniform_v_weights) {
                    match_templated<true, false, rating_func>(g, p_manager, mapping, imbalance, seed);
                } else if (g.uniform_e_weights) {
                    match_templated<false, true, rating_func>(g, p_manager, mapping, imbalance, seed);
                } else {
                    match_templated<false, false, rating_func>(g, p_manager, mapping, imbalance, seed);
                }
            };

            switch (config.rating_function) {
                case EdgeRatingFunction::WEIGHT:
                    dispatch_with_rating(std::integral_constant<EdgeRatingFunction, EdgeRatingFunction::WEIGHT>{});
                    break;
                case EdgeRatingFunction::EXPANSION:
                    dispatch_with_rating(std::integral_constant<EdgeRatingFunction, EdgeRatingFunction::EXPANSION>{});
                    break;
                case EdgeRatingFunction::EXPANSIONSTAR:
                    dispatch_with_rating(std::integral_constant<EdgeRatingFunction, EdgeRatingFunction::EXPANSIONSTAR>{});
                    break;
                case EdgeRatingFunction::EXPANSIONSTARSTAR:
                    dispatch_with_rating(std::integral_constant<EdgeRatingFunction, EdgeRatingFunction::EXPANSIONSTARSTAR>{});
                    break;
                case EdgeRatingFunction::INNEROUTER:
                    dispatch_with_rating(std::integral_constant<EdgeRatingFunction, EdgeRatingFunction::INNEROUTER>{});
                    break;
            }
        }

    private:
        template<bool t_uniform_v_weights, bool t_uniform_e_weights, EdgeRatingFunction t_rating_function>
        void match_templated(const graph_t &g,
                            const p_manager_t &p_manager,
                            Mapping &mapping,
                            f64 imbalance,
                            u64 seed) {
            HEIPROMAP_PROFILE_SCOPE("coarsening", "HeavyEdgeMatching", "match");

            mapping.initialize(g.n);

            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            Matching matching;
            matching.initialize(g.n);

            std::vector<vertex_t> permutation(g.n);
            std::iota(permutation.begin(), permutation.end(), 0);
            
            // Randomly shuffle the vertices to avoid bias
            RandomEngine random_engine(seed);
            std::shuffle(permutation.begin(), permutation.end(), random_engine.generator);

            for (vertex_t u : permutation) {
                if (matching.is_matched(u)) continue;

                partition_t u_id = p_manager[u];
                weight_t u_w = t_uniform_v_weights ? 1 : g.v_weights[u];

                vertex_t best_v = u;
                f32 best_rating = -1.0f;

                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    vertex_t v = g.edges_v[j];
                    if (matching.is_matched(v)) continue;
                    
                    if (u_id != p_manager[v]) continue;

                    weight_t v_w = t_uniform_v_weights ? 1 : g.v_weights[v];
                    if (u_w + v_w > lmax) continue;

                    weight_t ew = t_uniform_e_weights ? 1 : g.edges_w[j];
                    
                    f32 edge_rating;
                    if constexpr (t_uniform_v_weights && t_uniform_e_weights) {
                        edge_rating = 1.0f;
                    } else {
                        if constexpr (t_rating_function == EdgeRatingFunction::WEIGHT) {
                            edge_rating = (f32) ew;
                        } else if constexpr (t_rating_function == EdgeRatingFunction::EXPANSION) {
                            edge_rating = (f32) ew / (f32) (u_w + v_w);
                        } else if constexpr (t_rating_function == EdgeRatingFunction::EXPANSIONSTAR) {
                            edge_rating = (f32) ew / (f32) (u_w * v_w);
                        } else if constexpr (t_rating_function == EdgeRatingFunction::EXPANSIONSTARSTAR) {
                            edge_rating = (f32) (ew * ew) / (f32) (u_w * v_w);
                        } else if constexpr (t_rating_function == EdgeRatingFunction::INNEROUTER) {
                            weight_t out_v = 0;
                            for (u64 i = g.neighborhoods[v]; i < g.neighborhoods[v + 1]; ++i) {
                                out_v += g.edges_w[i];
                            }
                            weight_t out_u = 0;
                            for (u64 i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                                out_u += g.edges_w[i];
                            }
                            edge_rating = (f32) ew / ((f32) (out_v + out_u - (2 * ew)));
                        }
                    }

                    if (edge_rating > best_rating) {
                        best_rating = edge_rating;
                        best_v = v;
                    }
                }

                if (best_v != u) {
                    matching.add(u, best_v);
                }
            }

            // Finalize matching and create mapping
            {
                HEIPROMAP_PROFILE_SCOPE("coarsening", "HeavyEdgeMatching", "finalize");
                matching.set_translation();
                mapping.set_coarse_n(matching.get_n_coarse_nodes());
                for (vertex_t u = 0; u < matching.get_n(); ++u) {
                    mapping.set(u, matching.get_n(u));
                }
            }
        }
    };
}

#endif // HEIPROMAP_HEAVY_EDGE_MATCHING_H
