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

#ifndef HEIPROMAP_THREE_VERTEX_LABEL_PROPAGATION_REFINEMENT_H
#define HEIPROMAP_THREE_VERTEX_LABEL_PROPAGATION_REFINEMENT_H

#include <iostream>

#include "../definitions.h"
#include "../utility/utils.h"
#include "ISerialRefiner.h"

namespace HeiProMap {
    class ThreeVertexLabelPropagationConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit ThreeVertexLabelPropagationConfiguration(const std::string &t_name) : ISerialRefinerConfiguration(t_name) {}

        u64 last_n_levels = 2;
        u64 max_iteration = 25; // how many iterations to run the algorithm at most
    };

    class ThreeVertexLabelPropagationRefinement final : public ISerialRefiner {
        vertex_t                 m_n         = 0;
        vertex_t                 m_m         = 0;
        partition_t              m_k         = 0;
        f64                      m_imbalance = 0.0;
        weight_t                 m_lmax      = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t>    m_distance;
        u64                      m_seed      = 0;

        AlignedArray<u32> vertex_used;
        u32               vertex_marker = 0;

        AlignedArray<u32> block_used;
        u32               block_marker = 0;

        AlignedArray<vertex_t> curr_boundary;
        size_t                 curr_boundary_size = 0;

        AlignedArray<partition_t> u_move_ids;
        size_t                    u_move_ids_size = 0;

        AlignedArray<partition_t> v_move_ids;
        size_t                    v_move_ids_size = 0;

        RandomEngine                                   *random_engine = nullptr;
        const ThreeVertexLabelPropagationConfiguration *config        = nullptr;

    public:
        ThreeVertexLabelPropagationRefinement() = default;

        ~ThreeVertexLabelPropagationRefinement() override = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const f64 t_imbalance,
                        const weight_t t_lmax,
                        const std::vector<partition_t> &t_hierarchy,
                        const std::vector<weight_t> &t_distance,
                        RandomEngine &t_random_engine,
                        const ISerialRefinerConfiguration &i_config) override {
            m_n         = t_n;
            m_m         = t_m;
            m_k         = t_k;
            m_lmax      = t_lmax;
            m_imbalance = t_imbalance;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;

            random_engine = &t_random_engine;
            config        = dynamic_cast<const ThreeVertexLabelPropagationConfiguration *>(&i_config);

            vertex_used.initialize(m_n, 0);
            block_used.initialize(m_k, 0);

            curr_boundary.initialize(m_n);
            curr_boundary_size = 0;

            u_move_ids.initialize(m_k);
            u_move_ids_size = 0;

            v_move_ids.initialize(m_k);
            v_move_ids_size = 0;
        }

        void refine(const u64 level,
                    const u64 max_level,
                    graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph) override {
            ScopedTimer _t("refinement", "ThreeVertexLabelPropagationRefinement", "refine");

            if (level + config->last_n_levels < max_level) {
                return;
            }

            bool     move_occurred = true;
            for (u64 iteration     = 0; iteration < config->max_iteration && move_occurred; ++iteration) {
                move_occurred = false;

                curr_boundary_size = 0;
                for (partition_t id = 0; id < bv_manager.get_k(); ++id) {
                    forall_bv_id_iu(bv_manager, id, i, u)
                    {
                        curr_boundary[curr_boundary_size++] = u;
                    }
                    endfor
                }
                std::shuffle(curr_boundary.get_ptr(), curr_boundary.get_ptr() + curr_boundary_size, random_engine->generator);

                vertex_marker += 1;
                for (size_t i = 0; i < curr_boundary_size; ++i) {
                    vertex_t u = curr_boundary[i];
                    if (vertex_used[u] == vertex_marker) { continue; }
                    if (!bv_manager.is_boundary(u)) { continue; }

                    // collect all vertices u, uu, uuu
                    std::vector<vertex_t> vertices;
                    vertices.push_back(u);

                    std::vector<u8> inserted(g.n, 0);
                    inserted[u] = 1;

                    forall_guiv(g, u, j, uu)
                        {
                            if (uu >= u) { continue; }
                            if (inserted[uu] == 0) {
                                vertices.push_back(uu);
                                inserted[uu] = 1;
                            }
                            forall_guiv(g, uu, k, uuu)
                                {
                                    if (uuu >= uu || uuu >= u) { continue; }
                                    if (inserted[uuu] == 0) {
                                        vertices.push_back(uuu);
                                        inserted[uuu] = 1;
                                    }
                                }
                            endfor
                        }
                    endfor

                    // collect for all vertices the connected partitions
                    std::vector<u8>                       id_inserted(m_k);
                    std::vector<std::vector<partition_t>> partitions(vertices.size());

                    for (size_t j = 0; j < vertices.size(); ++j) {
                        vertex_t    v    = vertices[j];
                        partition_t v_id = p_manager[v];
                        std::fill(id_inserted.begin(), id_inserted.end(), 0);

                        forall_guiv(g, v, k, vv)
                            {
                                partition_t vv_id = p_manager[vv];
                                if (v_id == vv_id || id_inserted[vv_id] == 1) { continue; }

                                partitions[j].push_back(vv_id);
                                id_inserted[vv_id] = 1;
                            }
                        endfor
                    }

                    vertex_t    best_v = 0, best_vv = 0, best_vvv = 0;
                    weight_t    best_v_weight = 0, best_vv_weight = 0, best_vvv_weight = 0;
                    partition_t best_v_id = 0, best_vv_id = 0, best_vvv_id = 0;
                    partition_t best_new_v_id = 0, best_new_vv_id = 0, best_new_vvv_id = 0;
                    s64         best_qap_delta = -1;

                    // search for the best triple combination that does not overload
                    for (size_t j = 0; j < vertices.size(); ++j) {
                        for (size_t k = j + 1; k < vertices.size(); ++k) {
                            for (size_t l = k + 1; l < vertices.size(); ++l) {
                                // get the three vertices
                                vertex_t v   = vertices[j];
                                vertex_t vv  = vertices[k];
                                vertex_t vvv = vertices[l];

                                weight_t v_weight   = g.v_weights[v];
                                weight_t vv_weight  = g.v_weights[vv];
                                weight_t vvv_weight = g.v_weights[vvv];

                                partition_t v_id   = p_manager[v];
                                partition_t vv_id  = p_manager[vv];
                                partition_t vvv_id = p_manager[vvv];

                                for (partition_t new_v_id: partitions[j]) {
                                    for (partition_t new_vv_id: partitions[k]) {
                                        for (partition_t new_vvv_id: partitions[l]) {
                                            // make temporary moves
                                            p_manager.move(v, v_weight, v_id, new_v_id);
                                            p_manager.move(vv, vv_weight, vv_id, new_vv_id);
                                            p_manager.move(vvv, vvv_weight, vvv_id, new_vvv_id);

                                            bool overloaded = p_manager.get_bweight(new_v_id) > m_lmax || p_manager.get_bweight(new_vv_id) > m_lmax || p_manager.get_bweight(new_vvv_id) > m_lmax;

                                            // revert the moves
                                            p_manager.move(v, v_weight, new_v_id, v_id);
                                            p_manager.move(vv, vv_weight, new_vv_id, vv_id);
                                            p_manager.move(vvv, vvv_weight, new_vvv_id, vvv_id);

                                            if (overloaded) { continue; }

                                            // get the qap delta
                                            s64 qap_delta = get_qap_delta(g, v, vv, vvv, v_id, vv_id, vvv_id, new_v_id, new_vv_id, new_vvv_id, p_manager, d_oracle);

                                            if (qap_delta > best_qap_delta) {
                                                best_qap_delta = qap_delta;

                                                best_v   = v;
                                                best_vv  = vv;
                                                best_vvv = vvv;

                                                best_v_weight   = v_weight;
                                                best_vv_weight  = vv_weight;
                                                best_vvv_weight = vvv_weight;

                                                best_v_id   = v_id;
                                                best_vv_id  = vv_id;
                                                best_vvv_id = vvv_id;

                                                best_new_v_id   = new_v_id;
                                                best_new_vv_id  = new_vv_id;
                                                best_new_vvv_id = new_vvv_id;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    if (best_qap_delta > 0) {
                        bv_manager.move(g, p_manager, best_v, best_v_id, best_new_v_id);
                        q_graph.move(g, p_manager, best_v, best_v_id, best_new_v_id);
                        p_manager.move(best_v, best_v_weight, best_v_id, best_new_v_id);

                        bv_manager.move(g, p_manager, best_vv, best_vv_id, best_new_vv_id);
                        q_graph.move(g, p_manager, best_vv, best_vv_id, best_new_vv_id);
                        p_manager.move(best_vv, best_vv_weight, best_vv_id, best_new_vv_id);

                        bv_manager.move(g, p_manager, best_vvv, best_vvv_id, best_new_vvv_id);
                        q_graph.move(g, p_manager, best_vvv, best_vvv_id, best_new_vvv_id);
                        p_manager.move(best_vvv, best_vvv_weight, best_vvv_id, best_new_vvv_id);

                        vertex_used[best_v]   = vertex_marker;
                        vertex_used[best_vv]  = vertex_marker;
                        vertex_used[best_vvv] = vertex_marker;
                        move_occurred = true;
                    }
                }
            }
        }

        void refine_layer([[maybe_unused]] const u64 level,
                          [[maybe_unused]] const u64 max_level,
                          [[maybe_unused]] graph_t &g,
                          [[maybe_unused]] d_oracle_t &d_oracle,
                          [[maybe_unused]] bv_manager_t &bv_manager,
                          [[maybe_unused]] p_manager_t &p_manager,
                          [[maybe_unused]] q_graph_t &q_graph,
                          [[maybe_unused]] size_t layer) override {}
    };
}

#endif //HEIPROMAP_THREE_VERTEX_LABEL_PROPAGATION_REFINEMENT_H
