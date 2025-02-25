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

#ifndef HEIPROMAP_LABEL_PROPAGATION_REFINEMENT_FARAJ20_H
#define HEIPROMAP_LABEL_PROPAGATION_REFINEMENT_FARAJ20_H

#include <random>

#include "../../definitions.h"
#include "../interfaces/ISerialRefiner.h"

namespace HeiProMap {
    struct LabelPropagationFaraj20Configuration {
        u64 max_iteration = 25; // how many iterations to run the algorithm at most
    };

    /**
     * Executes label propagation refinement as described in
     * > Marcelo Fonseca Faraj, Alexander van der Grinten, Henning Meyerhenke, Jesper Larsson Träff, and Christian Schulz.
     * > High-quality Hierarchical Process Mapping.
     * > In 18th International Symposium on Experimental Algorithms, SEA 2020, June 16-18, 2020, Catania, Italy, volume 160 of LIPIcs, pages 4:1–4:15.
     */
    class LabelPropagationRefinementFaraj20 final : public ISerialRefiner {
        vertex_t m_n    = 0;
        vertex_t m_m    = 0;
        partition_t m_k = 0;
        weight_t m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;
        u64 m_seed = 0;

        u32* vertex_used = nullptr;
        u32 vertex_marker = 0;

        u32* block_used = nullptr;
        u32 block_marker = 0;

        std::mt19937 gen;
        std::uniform_real_distribution<float> dis;

    public:
        LabelPropagationRefinementFaraj20() = default;

        ~LabelPropagationRefinementFaraj20() override {
            free(vertex_used);
            free(block_used);
        }

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

            vertex_t m_n_64 = round_up_64(m_n);
            vertex_used = (u32*) aligned_alloc(64, m_n_64 * sizeof(u32));
            std::fill_n(vertex_used, m_n_64, vertex_marker);

            block_used = (u32*) aligned_alloc(64, m_n_64 * sizeof(u32));
            std::fill_n(vertex_used, m_n_64, vertex_marker);

            gen.seed(m_seed);
            dis = std::uniform_real_distribution<float>(0.0f, 1.0f);
        }

        template <typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle, typename TSerialQuotientGraph>
        void refine(LabelPropagationFaraj20Configuration& config,
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

            bool move_occurred = true;
            for (u64 iteration = 0; iteration < config.max_iteration && move_occurred; ++iteration) {
                move_occurred = false;

                vertex_marker += 1;
                for (vertex_t u : bv_manager) {
                    if (!is_boundary(g, p_manager, u)) { continue; }
                    if (vertex_used[u] == vertex_marker) { continue; } // we already used u in this iteration

                    weight_t u_weight = g.get_weight(u);
                    partition_t u_id  = p_manager[u];

                    // make the move that reduces qap the most
                    partition_t best_u_id = u_id;
                    s64 best_qap_delta    = -1;
                    u32 counter           = 0;

                    block_marker += 1;
                    for (const auto& [v, w] : g[u]) {
                        partition_t v_id = p_manager[v];
                        if (v_id == u_id || v_id == best_u_id) { continue; }

                        if (block_used[v_id] != block_marker) {
                            if (p_manager.get_bweight(v_id) + u_weight <= m_lmax) {
                                s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);

                                if (qap_delta > best_qap_delta) {
                                    best_qap_delta = qap_delta;
                                    best_u_id      = v_id;
                                    counter        = 1;
                                } else if (qap_delta == best_qap_delta && qap_delta != -1) {
                                    counter += 1;
                                    // choose with probability 1/counter as it ensures uniform distribution
                                    if (dis(gen) < 1.0f / (f32)counter) {
                                        best_u_id = v_id;
                                    }
                                }
                            }
                            block_used[v_id] = block_marker;
                        }
                    }

                    if (best_u_id != u_id) {
                        // choose if positive, if 0-gain choose 50% of the time
                        if (best_qap_delta > 0 || dis(gen) < 0.5) {
                            bv_manager.move(g, p_manager, u, u_id, best_u_id);
                            q_graph.move(g, p_manager, u, u_id, best_u_id);
                            p_manager.move(u, u_weight, u_id, best_u_id);
                            vertex_used[u] = vertex_marker;
                            move_occurred  = true;
                        }
                    }
                }
            }
        }
    };
}

#endif //HEIPROMAP_LABEL_PROPAGATION_REFINEMENT_FARAJ20_H
