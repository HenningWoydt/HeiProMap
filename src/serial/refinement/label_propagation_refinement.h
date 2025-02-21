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

#ifndef HEIPROMAP_LABEL_PROPAGATION_REFINEMENT_H
#define HEIPROMAP_LABEL_PROPAGATION_REFINEMENT_H

#include <random>

#include "../../definitions.h"
#include "../interfaces/ISerialActiveVertexManager.h"
#include "../interfaces/ISerialBoundaryVertexManager.h"
#include "../interfaces/ISerialDistanceOracle.h"
#include "../interfaces/ISerialGraph.h"
#include "../interfaces/ISerialPartitionManager.h"
#include "../interfaces/ISerialQuotientGraph.h"
#include "../interfaces/ISerialRefiner.h"

namespace HeiProMap {
    struct LabelPropagationConfiguration {
        u64 max_iteration = 25; // how many iterations to run the algorithm at most
    };

    class LabelPropagationRefinement final : public ISerialRefiner {
        vertex_t m_n    = 0;
        vertex_t m_m    = 0;
        partition_t m_k = 0;
        weight_t m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;
        u64 m_seed = 0;

        std::vector<u32> block_used;
        u32 block_marker = 0;

        std::vector<vertex_t> vertices_this_round;
        std::vector<vertex_t> vertices_next_round;
        std::vector<u32> vertex_added;
        u32 vertex_marker = 0;

        std::mt19937 gen;
        std::uniform_real_distribution<float> dis;

    public:
        LabelPropagationRefinement() = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const weight_t t_lmax,
                        const std::vector<partition_t>& t_hierarchy,
                        const std::vector<weight_t>& t_distance,
                        const u64 t_seed) override {
            m_n         = t_n;
            m_m         = t_m;
            m_k         = t_k;
            m_lmax      = t_lmax;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;
            m_seed      = t_seed;

            block_used.resize(t_k, block_marker);
            vertex_added.resize(t_n, vertex_marker);

            gen.seed(m_seed);
            dis = std::uniform_real_distribution<float>(0.0f, 1.0f);
        }

        template <typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle, typename TSerialQuotientGraph>
        void refine(LabelPropagationConfiguration& config,
                    TSerialGraph& g,
                    [[maybe_unused]] TSerialActiveVertexManager& av_manager,
                    TSerialBoundaryVertexManager& bv_manager,
                    TSerialPartitionManager& p_manager,
                    TSerialDistanceOracle& d_oracle,
                    TSerialQuotientGraph& q_graph) {
            static_assert(std::is_base_of_v<ISerialGraph, TSerialGraph>, "TSerialGraph must inherit from ISerialGraph");
            static_assert(std::is_base_of_v<ISerialActiveVertexManager, TSerialActiveVertexManager>, "TSerialActiveVertexManager must inherit from ISerialActiveVertexManager");
            static_assert(std::is_base_of_v<ISerialBoundaryVertexManager, TSerialBoundaryVertexManager>, "TSerialBoundaryVertexManager must inherit from ISerialBoundaryVertexManager");
            static_assert(std::is_base_of_v<ISerialPartitionManager, TSerialPartitionManager>, "TSerialPartitionManager must inherit from ISerialPartitionManager");
            static_assert(std::is_base_of_v<ISerialDistanceOracle, TSerialDistanceOracle>, "TSerialDistanceOracle must inherit from ISerialDistanceOracle");
            static_assert(std::is_base_of_v<ISerialQuotientGraph, TSerialQuotientGraph>, "TSerialQuotientGraph must inherit from ISerialQuotientGraph");

            vertices_this_round.clear();
            for (vertex_t u : bv_manager) { vertices_this_round.push_back(u); }

            for (u64 iteration = 0; iteration < config.max_iteration && !vertices_this_round.empty(); ++iteration) {
                vertex_marker += 1;

                for (vertex_t u : vertices_this_round) {
                    if (!is_boundary(g, p_manager, u)) { continue; }

                    weight_t u_weight = g.get_weight(u);
                    partition_t u_id  = p_manager[u];

                    // make the move that reduces qap the most
                    partition_t best_u_id = u_id;
                    s64 best_qap_delta    = -1;
                    u32 counter           = 0;
                    bool all_overloaded   = true;

                    block_marker += 1;
                    for (const auto& [v, w] : g[u]) {
                        partition_t v_id = p_manager[v];
                        if (v_id == u_id || v_id == best_u_id) { continue; }

                        if (block_used[v_id] != block_marker) {
                            if (p_manager.get_bweight(v_id) + u_weight <= m_lmax) {
                                all_overloaded = false;
                                s64 qap_delta  = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                                if (qap_delta < 0) {
                                    block_used[v_id] = block_marker;
                                    continue;
                                }

                                if (qap_delta > best_qap_delta) {
                                    // straight better so chose this
                                    best_qap_delta = qap_delta;
                                    best_u_id      = v_id;
                                    counter        = 1;
                                } else if (qap_delta == best_qap_delta) {
                                    if (p_manager.get_bweight(v_id) < p_manager.get_bweight(best_u_id)) {
                                        // same, but another partition has smaller weight
                                        best_u_id = v_id;
                                        counter   = 1;
                                    } else {
                                        // same delta and same weight, so uniform choosing
                                        counter += 1;
                                        if (dis(gen) < 1.0f / (f32)counter) {
                                            best_u_id = v_id;
                                        }
                                    }
                                }
                                block_used[v_id] = block_marker;
                            }
                        }
                    }

                    if (best_u_id != u_id) {
                        // choose if positive, if 0-gain choose 50% of the time
                        if (best_qap_delta > 0 || dis(gen) < 0.5) {
                            bv_manager.move(g, p_manager, u, u_id, best_u_id);
                            q_graph.move(g, p_manager, u, u_id, best_u_id);
                            p_manager.move(u, u_weight, u_id, best_u_id);

                            // u gets checked again
                            if (vertex_added[u] != vertex_marker) {
                                vertices_next_round.push_back(u);
                                vertex_added[u] = vertex_marker;
                            }

                            // neighborhood of u also gets checked
                            for (const auto& [v, w] : g[u]) {
                                if (vertex_added[v] != vertex_marker) {
                                    vertices_next_round.push_back(v);
                                    vertex_added[v] = vertex_marker;
                                }
                            }
                        }
                    }

                    if (all_overloaded) {
                        // u gets checked again
                        if (vertex_added[u] != vertex_marker) {
                            vertices_next_round.push_back(u);
                            vertex_added[u] = vertex_marker;
                        }
                    }
                }
                // swap in vertices of next round
                std::swap(vertices_this_round, vertices_next_round);
                vertices_next_round.clear();
            }
        }
    };
}

#endif //HEIPROMAP_LABEL_PROPAGATION_REFINEMENT_H
