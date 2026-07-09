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

#ifndef HEIPROMAP_SWAP_REFINEMENT_H
#define HEIPROMAP_SWAP_REFINEMENT_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <numeric>

#include "../definitions.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/distance_oracle.h"
#include "../datastructures/partition_manager.h"
#include "../utility/aligned_array.h"
#include "../utility/profiler.h"
#include "../utility/random_engine.h"
#include "../utility/utils.h"
#include "../utility/qap.h"
#include "../utility/functions.h"
#include "../utility/small_map.h"

namespace HeiProMap {
    class SwapRefinementConfiguration final {
    public:
        explicit SwapRefinementConfiguration(const std::string &t_name) {
            name = t_name;
        }

        std::string name;
        bool enabled = false;
        u64 max_iteration = 5;
    };

    class SwapRefinement final {
        std::vector<RandomEngine> m_rnd_engines;
        const SwapRefinementConfiguration *m_config = nullptr;

    public:
        SwapRefinement() = default;

        ~SwapRefinement() = default;

        void initialize(const u64 seed,
                        const SwapRefinementConfiguration &i_config) {
            m_config = &i_config;

            m_rnd_engines.clear();
            m_rnd_engines.emplace_back(seed);
        }

        void refine(graph_t &g,
                    d_oracle_t &d_oracle,
                    p_manager_t &p_manager,
                    const AlignedArray<weight_t> &lmax_constraints) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "SwapRefinement", "refine");

            if (p_manager.k < 2) return;
            if (g.uniform_v_weights && g.uniform_e_weights) refine_impl<true, true>(g, d_oracle, p_manager, lmax_constraints);
            else if (g.uniform_v_weights) refine_impl<true, false>(g, d_oracle, p_manager, lmax_constraints);
            else if (g.uniform_e_weights) refine_impl<false, true>(g, d_oracle, p_manager, lmax_constraints);
            else refine_impl<false, false>(g, d_oracle, p_manager, lmax_constraints);
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void refine_impl(graph_t &g,
                         d_oracle_t &d_oracle,
                         p_manager_t &p_manager,
                         const AlignedArray<weight_t> &lmax_constraints) {
            if (m_rnd_engines.empty()) return;
            RandomEngine &random_engine = m_rnd_engines[0];

            std::vector<vertex_t> vertices(g.n);
            std::iota(vertices.begin(), vertices.end(), 0);

            for (u64 iteration = 0; iteration < m_config->max_iteration; ++iteration) {
                bool improved = false;
                std::shuffle(vertices.begin(), vertices.end(), random_engine.generator);

                // Phase 1: Simple Label Propagation (Single Moves)
                for (vertex_t u : vertices) {
                    partition_t u_id = p_manager[u];
                    const weight_t u_w = t_uniform_v_weights ? 1 : g.v_weights[u];

                    partition_t best_target = u_id;
                    weight_t best_gain = 0;

                    // Collect target blocks from neighbors
                    FlatMap<partition_t, bool> target_blocks;
                    for (size_t e = g.neighborhoods[u]; e < g.neighborhoods[u+1]; ++e) {
                        partition_t v_id = p_manager[g.edges_v[e]];
                        if (v_id != u_id) target_blocks[v_id] = true;
                    }

                    for (auto const& kv : target_blocks) {
                        partition_t target_id = kv.first;
                        if (p_manager.get_bweight(target_id) + u_w > lmax_constraints[target_id]) continue;
                        
                        weight_t gain = get_u_qap_delta(g, u, u_id, target_id, p_manager, d_oracle);
                        if (gain > best_gain) {
                            best_gain = gain;
                            best_target = target_id;
                        }
                    }

                    if (best_target != u_id) {
                        p_manager.move_serial(u, u_w, u_id, best_target);
                        improved = true;
                    }
                }

                // Phase 2: Simple Swap (Vertex Exchanges)
                std::shuffle(vertices.begin(), vertices.end(), random_engine.generator);
                for (vertex_t u : vertices) {
                    partition_t u_id = p_manager[u];
                    const weight_t u_w = t_uniform_v_weights ? 1 : g.v_weights[u];

                    for (size_t e = g.neighborhoods[u]; e < g.neighborhoods[u+1]; ++e) {
                        vertex_t v = g.edges_v[e];
                        partition_t v_id = p_manager[v];
                        if (u_id == v_id) continue;

                        const weight_t v_w = t_uniform_v_weights ? 1 : g.v_weights[v];

                        if (p_manager.get_bweight(u_id) - u_w + v_w > lmax_constraints[u_id]) continue;
                        if (p_manager.get_bweight(v_id) - v_w + u_w > lmax_constraints[v_id]) continue;

                        const weight_t gain_u = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                        const weight_t gain_v = get_u_qap_delta(g, v, v_id, u_id, p_manager, d_oracle);

                        weight_t w_uv = t_uniform_e_weights ? 1 : g.edges_w[e];
                        const weight_t dist_ij = d_oracle.get(u_id, v_id);

                        const weight_t swap_gain = gain_u + gain_v - 2 * w_uv * dist_ij;

                        if (swap_gain > 0) {
                            p_manager.move_serial(u, u_w, u_id, v_id);
                            p_manager.move_serial(v, v_w, v_id, u_id);
                            improved = true;
                            break; // u has swapped, move to next u
                        }
                    }
                }

                if (!improved) break;
            }
        }
    };
}

#endif // HEIPROMAP_SWAP_REFINEMENT_H
