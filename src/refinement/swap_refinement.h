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

#include "../definitions.h"
#include "../datastructures/block_conn.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/distance_oracle.h"
#include "../datastructures/partition_manager.h"
#include "../datastructures/quotient_graph.h"
#include "../utility/aligned_array.h"
#include "../utility/profiler.h"
#include "../utility/random_engine.h"
#include "../utility/utils.h"
#include "../utility/qap.h"
#include "../utility/functions.h"

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
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
        u64 m_threads = 1;

        std::vector<RandomEngine> m_rnd_engines;
        const SwapRefinementConfiguration *m_config = nullptr;

    public:
        SwapRefinement() = default;

        ~SwapRefinement() = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const u64 t_threads,
                        const u64 seed,
                        const SwapRefinementConfiguration &i_config) {
            m_n = t_n;
            m_m = t_m;
            m_k = t_k;
            m_threads = t_threads;
            m_config = &i_config;

            m_rnd_engines.clear();
            m_rnd_engines.reserve(m_threads);
            for (u64 t = 0; t < m_threads; ++t) {
                m_rnd_engines.emplace_back(seed + t);
            }
        }

        void refine(graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph,
                    block_conn_t &block_conn,
                    const AlignedArray<weight_t> &lmax_constraints,
                    bool uniform_v_weights,
                    bool uniform_e_weights) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "SwapRefinement", "refine");

            if (m_k < 2) return;
            if (uniform_v_weights && uniform_e_weights) refine_impl<true, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
            else if (uniform_v_weights) refine_impl<true, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
            else if (uniform_e_weights) refine_impl<false, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
            else refine_impl<false, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void refine_impl(graph_t &g,
                         d_oracle_t &d_oracle,
                         bv_manager_t &bv_manager,
                         p_manager_t &p_manager,
                         q_graph_t &q_graph,
                         block_conn_t &block_conn,
                         const AlignedArray<weight_t> &lmax_constraints) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "SwapRefinement", "refine");
            if (m_rnd_engines.empty()) return;
            RandomEngine &random_engine = m_rnd_engines[0];

            for (u64 iteration = 0; iteration < m_config->max_iteration; ++iteration) {
                bool moved = false;

                for (partition_t i = 0; i < m_k; ++i) {
                    for (partition_t j = i + 1; j < m_k; ++j) {
                        if (!q_graph.has_edge(i, j)) continue;

                        std::vector<vertex_t> b_i;
                        for (size_t k = 0; k < bv_manager.size(i); ++k) {
                            vertex_t u = bv_manager.get(i, k);
                            if (is_connected_to(g, p_manager, u, j)) { b_i.push_back(u); }
                        }
                        if (b_i.empty()) continue;

                        std::vector<vertex_t> b_j;
                        for (size_t k = 0; k < bv_manager.size(j); ++k) {
                            vertex_t v = bv_manager.get(j, k);
                            if (is_connected_to(g, p_manager, v, i)) { b_j.push_back(v); }
                        }
                        if (b_j.empty()) continue;

                        if (b_i.size() > 1) fast_shuffle_unchecked(b_i.data(), b_i.data() + b_i.size(), random_engine.generator);
                        if (b_j.size() > 1) fast_shuffle_unchecked(b_j.data(), b_j.data() + b_j.size(), random_engine.generator);

                        const weight_t dist_ij = d_oracle.get(i, j);

                        for (vertex_t u : b_i) {
                            if (p_manager[u] != i) continue;
                            const weight_t u_w = t_uniform_v_weights ? 1 : g.v_weights[u];
                            const weight_t gain_u = get_u_qap_delta_t<t_uniform_e_weights>(g, u, i, j, p_manager, d_oracle, block_conn);

                            for (vertex_t v : b_j) {
                                if (p_manager[v] != j) continue;
                                const weight_t v_w = t_uniform_v_weights ? 1 : g.v_weights[v];

                                if (p_manager.get_bweight(i) - u_w + v_w > lmax_constraints[i]) continue;
                                if (p_manager.get_bweight(j) - v_w + u_w > lmax_constraints[j]) continue;

                                const weight_t gain_v = get_u_qap_delta_t<t_uniform_e_weights>(g, v, j, i, p_manager, d_oracle, block_conn);

                                weight_t w_uv = 0;
                                for (size_t k = g.neighborhoods[u]; k < g.neighborhoods[u+1]; ++k) {
                                    if (g.edges_v[k] == v) {
                                        w_uv = t_uniform_e_weights ? 1 : g.edges_w[k];
                                        break;
                                    }
                                }

                                const weight_t swap_gain = gain_u + gain_v - 2 * w_uv * dist_ij;

                                if (swap_gain > 0) {
                                    bv_manager.move(g, p_manager, u, i, j);
                                    q_graph.move(g, p_manager, u, i, j);
                                    block_conn.move(g, u, i, j);
                                    p_manager.move_serial(u, u_w, i, j);

                                    bv_manager.move(g, p_manager, v, j, i);
                                    q_graph.move(g, p_manager, v, j, i);
                                    block_conn.move(g, v, j, i);
                                    p_manager.move_serial(v, v_w, j, i);

                                    moved = true;
                                    goto next_u;
                                }
                            }
                            next_u:;
                        }
                    }
                }
                if (!moved) break;
            }
        }
    };
}

#endif // HEIPROMAP_SWAP_REFINEMENT_H
