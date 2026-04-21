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

#ifndef HEIPROMAP_TWO_VERTEX_LABEL_PROPAGATION_REFINEMENT_H
#define HEIPROMAP_TWO_VERTEX_LABEL_PROPAGATION_REFINEMENT_H

#include <iostream>

#include "../definitions.h"
#include "../utility/utils.h"
#include "ISerialRefiner.h"

/*
namespace HeiProMap {
    class TwoVertexLabelPropagationConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit TwoVertexLabelPropagationConfiguration(const std::string &t_name) : ISerialRefinerConfiguration(t_name) {
        }

        u64 max_iteration = 25; // how many iterations to run the algorithm at most
    };

    class TwoVertexLabelPropagationRefinement final : public ISerialRefiner {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
        u64 m_seed = 0;

        AlignedArray<u32> vertex_used;
        u32 vertex_marker = 0;

        AlignedArray<u32> block_used;
        u32 block_marker = 0;

        AlignedArray<vertex_t> curr_boundary;
        size_t curr_boundary_size = 0;

        AlignedArray<partition_t> u_move_ids;
        size_t u_move_ids_size = 0;

        AlignedArray<partition_t> v_move_ids;
        size_t v_move_ids_size = 0;

        RandomEngine random_engine = RandomEngine(0);
        const TwoVertexLabelPropagationConfiguration *config = nullptr;

    public:
        TwoVertexLabelPropagationRefinement() = default;

        ~TwoVertexLabelPropagationRefinement() override = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const ISerialRefinerConfiguration &i_config) override {
            m_n = t_n;
            m_m = t_m;
            m_k = t_k;

            config = dynamic_cast<const TwoVertexLabelPropagationConfiguration *>(&i_config);

            vertex_used.initialize(m_n, 0);
            block_used.initialize(m_k, 0);

            curr_boundary.initialize(m_n);
            curr_boundary_size = 0;

            u_move_ids.initialize(m_k);
            u_move_ids_size = 0;

            v_move_ids.initialize(m_k);
            v_move_ids_size = 0;
        }

        void refine(graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph,
                    block_conn_t &block_conn,
                    f64 imbalance,
                    bool uniform_v_weights,
                    bool uniform_e_weights) override {
            if (uniform_v_weights && uniform_e_weights)      refine_impl<true, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, imbalance);
            else if (uniform_v_weights)                      refine_impl<true, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, imbalance);
            else if (uniform_e_weights)                      refine_impl<false, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, imbalance);
            else                                             refine_impl<false, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, imbalance);
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void refine_impl(graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph,
                    block_conn_t &block_conn,
                    f64 imbalance) {
            ScopedTimer _t("refinement", "TwoVertexLabelPropagationRefinement", "refine");

            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            bool move_occurred = true;
            for (u64 iteration = 0; iteration < config->max_iteration && move_occurred; ++iteration) {
                move_occurred = false;

                curr_boundary_size = 0;
                for (partition_t id = 0; id < bv_manager.get_k(); ++id) {
                    for (size_t i = 0; i < bv_manager.size(id); ++i) { const vertex_t u = bv_manager.get(id, i); {
                            curr_boundary[curr_boundary_size++] = u;
                        }
                    }
                }
                std::shuffle(curr_boundary.get_ptr(), curr_boundary.get_ptr() + curr_boundary_size, random_engine.generator);

                vertex_marker += 1;
                for (size_t i = 0; i < curr_boundary_size; ++i) {
                    vertex_t u = curr_boundary[i];
                    if (vertex_used[u] == vertex_marker) { continue; } // vertex was used
                    if (!bv_manager.is_boundary(u)) { continue; }      // vertex is not boundary

                    weight_t u_weight = g.v_weights[u];
                    partition_t u_id = p_manager[u];

                    // get all connected partitions to u
                    block_marker += 1;
                    block_used[u_id] = block_marker;
                    u_move_ids_size = 0;
                    for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) { const vertex_t neighbor = g.edges_v[i]; {
                            partition_t neighbor_id = p_manager[neighbor];
                            if (block_used[neighbor_id] != block_marker) {
                                u_move_ids[u_move_ids_size++] = neighbor_id;
                                block_used[neighbor_id] = block_marker;
                            }
                        }
                    }

                    partition_t best_u_move_id = 0;
                    vertex_t best_v = 0;
                    partition_t best_v_id = 0;
                    weight_t best_v_weight = 0;
                    partition_t best_v_move_id = 0;
                    weight_t best_qap_delta = -1;

                    // check all neighbors v
                    for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) { const vertex_t v = g.edges_v[i]; {
                            if (vertex_used[v] == vertex_marker) { continue; } // vertex was used
                            if (!bv_manager.is_boundary(v)) { continue; }      // vertex is not boundary

                            weight_t v_weight = g.v_weights[v];
                            partition_t v_id = p_manager[v];

                            // get all connected partitions to v
                            block_marker += 1;
                            block_used[v_id] = block_marker;
                            v_move_ids_size = 0;
                            for (size_t i = g.neighborhoods[v]; i < g.neighborhoods[v + 1]; ++i) { const vertex_t neighbor = g.edges_v[i]; {
                                    partition_t neighbor_id = p_manager[neighbor];

                                    if (block_used[neighbor_id] != block_marker) {
                                        v_move_ids[v_move_ids_size++] = neighbor_id;
                                        block_used[neighbor_id] = block_marker;
                                    }
                                }
                            }

                            // check if moving u to u_ids and v to v_ids simultaneously would improve the score
                            for (size_t j = 0; j < u_move_ids_size; ++j) {
                                for (size_t l = 0; l < v_move_ids_size; ++l) {
                                    partition_t u_move_id = u_move_ids[j];
                                    partition_t v_move_id = v_move_ids[l];

                                    weight_t u_move_id_weight = p_manager.get_bweight(u_move_id);
                                    weight_t v_move_id_weight = p_manager.get_bweight(v_move_id);

                                    if (u_move_id == v_id && u_move_id_weight + u_weight - v_weight > lmax) { continue; }
                                    if (u_move_id != v_id && u_move_id_weight + u_weight > lmax) { continue; }
                                    if (v_move_id == u_id && v_move_id_weight + v_weight - u_weight > lmax) { continue; }
                                    if (v_move_id != u_id && v_move_id_weight + v_weight > lmax) { continue; }
                                    if (u_move_id == v_move_id && u_move_id_weight + u_weight + v_weight > lmax) { continue; }

                                    // no overloading is happening, now compute the qap_delta
                                    weight_t qap_delta = get_qap_delta(g, u, u_id, u_move_id, v, v_id, v_move_id, p_manager, d_oracle);

                                    if (qap_delta > best_qap_delta) {
                                        best_u_move_id = u_move_id;
                                        best_v = v;
                                        best_v_id = v_id;
                                        best_v_weight = v_weight;
                                        best_v_move_id = v_move_id;
                                        best_qap_delta = qap_delta;
                                    }
                                }
                            }
                        }
                    }

                    if (best_qap_delta > 0) {
                        bv_manager.move(g, p_manager, u, u_id, best_u_move_id);
                        q_graph.move(g, p_manager, u, u_id, best_u_move_id);
                        p_manager.move(u, u_weight, u_id, best_u_move_id);

                        bv_manager.move(g, p_manager, best_v, best_v_id, best_v_move_id);
                        q_graph.move(g, p_manager, best_v, best_v_id, best_v_move_id);
                        p_manager.move(best_v, best_v_weight, best_v_id, best_v_move_id);

                        vertex_used[best_v] = vertex_marker;
                        move_occurred = true;
                    }
                    vertex_used[u] = vertex_marker;
                }
            }
        }
    };
}

*/

#endif //HEIPROMAP_TWO_VERTEX_LABEL_PROPAGATION_REFINEMENT_H
