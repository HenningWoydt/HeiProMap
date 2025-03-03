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
#include <random>

#include "../../definitions.h"
#include "../../macros.h"
#include "../interfaces/ISerialActiveVertexManager.h"
#include "../interfaces/ISerialBoundaryVertexManager.h"
#include "../interfaces/ISerialDistanceOracle.h"
#include "../interfaces/ISerialGraph.h"
#include "../interfaces/ISerialPartitionManager.h"
#include "../interfaces/ISerialQuotientGraph.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/utils.h"

namespace HeiProMap {
    struct TwoVertexLabelPropagationConfiguration {
        u64 max_iteration = 25; // how many iterations to run the algorithm at most
    };

    class TwoVertexLabelPropagationRefinement final : public ISerialRefiner {
        vertex_t m_n    = 0;
        vertex_t m_m    = 0;
        partition_t m_k = 0;
        weight_t m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;
        u64 m_seed = 0;

        u32* vertex_used  = nullptr;
        u32 vertex_marker = 0;

        u32* block_used  = nullptr;
        u32 block_marker = 0;

        vertex_t* curr_boundary   = nullptr;
        size_t curr_boundary_size = 0;

        partition_t* u_move_ids = nullptr;
        size_t u_move_ids_size  = 0;

        partition_t* v_move_ids = nullptr;
        size_t v_move_ids_size  = 0;

        std::mt19937 gen;
        std::uniform_real_distribution<float> dis;

    public:
        TwoVertexLabelPropagationRefinement() = default;

        ~TwoVertexLabelPropagationRefinement() override {
            free(vertex_used);
            free(block_used);
            free(curr_boundary);
            free(u_move_ids);
            free(v_move_ids);
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

            vertex_t m_n_64    = round_up_64(m_n);
            partition_t m_k_64 = round_up_64(m_k);

            vertex_used = (u32*)aligned_alloc(64, m_n_64 * sizeof(u32));
            std::fill_n(vertex_used, m_n_64, vertex_marker);

            block_used = (u32*)aligned_alloc(64, m_k_64 * sizeof(u32));
            std::fill_n(block_used, m_k_64, block_marker);

            curr_boundary      = (vertex_t*)aligned_alloc(64, m_n_64 * sizeof(vertex_t));
            curr_boundary_size = 0;

            u_move_ids      = (vertex_t*)aligned_alloc(64, m_k_64 * sizeof(vertex_t));
            u_move_ids_size = 0;

            v_move_ids      = (vertex_t*)aligned_alloc(64, m_k_64 * sizeof(vertex_t));
            v_move_ids_size = 0;


            gen.seed(m_seed);
            dis = std::uniform_real_distribution<float>(0.0f, 1.0f);
        }

        template <typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle, typename TSerialQuotientGraph>
        void refine(TwoVertexLabelPropagationConfiguration& config,
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
                auto sp                = std::chrono::high_resolution_clock::now();
                u64 max_possible_moves = 0;
                s64 max_possible_gain  = 0;

                move_occurred = false;

                curr_boundary_size = 0;
                for (vertex_t u : bv_manager) { curr_boundary[curr_boundary_size++] = u; }
                std::shuffle(curr_boundary, curr_boundary + curr_boundary_size, gen);

                vertex_marker += 1;
                for (size_t i = 0; i < curr_boundary_size; ++i) {
                    vertex_t u = curr_boundary[i];
                    if (vertex_used[u] == vertex_marker) { continue; }
                    if (!bv_manager.is_boundary(u)) { continue; }

                    weight_t u_weight = g.get_weight(u);
                    partition_t u_id  = p_manager[u];

                    // get all connected partitions to u
                    block_marker += 1;
                    block_used[u_id] = block_marker;
                    u_move_ids_size = 0;
                    for (const auto [neighbor, neighbor_w] : g[u]) {
                        partition_t neighbor_id = p_manager[neighbor];
                        if (block_used[neighbor_id] != block_marker) {
                            u_move_ids[u_move_ids_size++] = neighbor_id;
                            block_used[neighbor_id]       = block_marker;
                        }
                    }

                    partition_t best_u_move_id;
                    vertex_t best_v;
                    partition_t best_v_id;
                    weight_t best_v_weight;
                    partition_t best_v_move_id;
                    s64 best_qap_delta = -1;

                    // get all connected partitions to v
                    for (const auto [v, w] : g[u]) {
                        if (vertex_used[v] == vertex_marker) { continue; }
                        if (!bv_manager.is_boundary(v)) { continue; }

                        weight_t v_weight = g.get_weight(v);
                        partition_t v_id  = p_manager[v];

                        block_marker += 1;
                        block_used[v_id] = block_marker;
                        v_move_ids_size = 0;
                        for (const auto [neighbor, neighbor_w] : g[v]) {
                            partition_t neighbor_id = p_manager[neighbor];

                            if (block_used[neighbor_id] != block_marker) {
                                v_move_ids[v_move_ids_size++] = neighbor_id;
                                block_used[neighbor_id]       = block_marker;
                            }
                        }

                        // check if moving u to u_ids and v to v_ids simultaneously would improve the score
                        for (size_t j = 0; j < u_move_ids_size; ++j) {
                            for (size_t l = 0; l < v_move_ids_size; ++l) {
                                partition_t u_move_id = u_move_ids[j];
                                partition_t v_move_id = v_move_ids[l];

                                weight_t u_move_id_weight = p_manager.get_bweight(u_move_id);
                                weight_t v_move_id_weight = p_manager.get_bweight(v_move_id);

                                if (u_move_id == v_id && u_move_id_weight + u_weight - v_weight > m_lmax) { continue; }
                                if (u_move_id != v_id && u_move_id_weight + u_weight > m_lmax) { continue; }
                                if (v_move_id == u_id && v_move_id_weight + v_weight - u_weight > m_lmax) { continue; }
                                if (v_move_id != u_id && v_move_id_weight + v_weight > m_lmax) { continue; }

                                // no overloading is happening, now compute the qap_delta
                                s64 qap_delta = get_qap_delta(g, u, u_id, u_move_id, v, v_id, v_move_id, p_manager, d_oracle);

                                if (qap_delta > best_qap_delta) {
                                    best_u_move_id = u_move_id;
                                    best_v         = v;
                                    best_v_id      = v_id;
                                    best_v_weight  = v_weight;
                                    best_v_move_id = v_move_id;
                                    best_qap_delta = qap_delta;
                                }
                            }
                        }
                    }

                    if (best_qap_delta != -1) {
                        max_possible_moves += 1;
                        max_possible_gain += best_qap_delta;

                        bv_manager.move(g, p_manager, u, u_id, best_u_move_id);
                        q_graph.move(g, p_manager, u, u_id, best_u_move_id);
                        p_manager.move(u, u_weight, u_id, best_u_move_id);

                        bv_manager.move(g, p_manager, best_v, best_v_id, best_v_move_id);
                        q_graph.move(g, p_manager, best_v, best_v_id, best_v_move_id);
                        p_manager.move(best_v, best_v_weight, best_v_id, best_v_move_id);

                        vertex_used[u]      = 1;
                        vertex_used[best_v] = 1;
                        move_occurred       = true;
                    }
                }
                auto ep     = std::chrono::high_resolution_clock::now();
                f64 seconds = get_seconds(sp, ep);
                // std::cout << "Two Vertex Label Propagation - Iteration " << iteration << ": " << max_possible_moves << " possible moves with " << max_possible_gain << " max gain in " << seconds << " seconds!" << std::endl;
            }
        }
    };
}

#endif //HEIPROMAP_TWO_VERTEX_LABEL_PROPAGATION_REFINEMENT_H
