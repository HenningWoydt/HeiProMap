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

#include <queue>
#include <random>

#include "../../definitions.h"
#include "../../macros.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/utils.h"

namespace HeiProMap {
    /**
     * Executes label propagation refinement as described in
     * > Marcelo Fonseca Faraj, Alexander van der Grinten, Henning Meyerhenke, Jesper Larsson Träff, and Christian Schulz.
     * > High-quality Hierarchical Process Mapping.
     * > In 18th International Symposium on Experimental Algorithms, SEA 2020, June 16-18, 2020, Catania,Italy, volume 160 of LIPIcs, pages 4:1–4:15.
     */
    class LabelPropagationRefinementFaraj20 final : public ISerialRefiner {
        std::vector<partition_t> hierarchy;
        std::vector<weight_t>    distance;
        partition_t              k    = 0;
        weight_t                 lmax = 0;

        std::vector<s32> used;
        s32              mark = -1;

        std::random_device                    rd;
        std::mt19937                          gen;
        std::uniform_real_distribution<float> dis;

    public:
        LabelPropagationRefinementFaraj20() : gen(rd()), dis(0.0f, 1.0f) {}

        void initialize(const vertex_t n,
                        std::vector<partition_t> &t_hierarchy,
                        std::vector<weight_t> &t_distance,
                        const weight_t t_lmax) override {
            hierarchy = t_hierarchy;
            distance  = t_distance;
            k         = prod<partition_t>(hierarchy);
            lmax      = t_lmax;

            used.resize(n, -1);
        }

        /*
        template<typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle, typename TSerialQuotientGraph>
        void refine(TSerialGraph &g,
                    [[maybe_unused]] TSerialActiveVertexManager &av_manager,
                    TSerialBoundaryVertexManager &bv_manager,
                    TSerialPartitionManager &p_manager,
                    TSerialDistanceOracle &d_oracle,
                    TSerialQuotientGraph &q_graph) {
            static_assert(std::is_base_of_v<ISerialGraph, TSerialGraph>, "TSerialGraph must inherit from ISerialGraph");
            static_assert(std::is_base_of_v<ISerialActiveVertexManager, TSerialActiveVertexManager>, "TSerialActiveVertexManager must inherit from ISerialActiveVertexManager");
            static_assert(std::is_base_of_v<ISerialBoundaryVertexManager, TSerialBoundaryVertexManager>, "TSerialBoundaryVertexManager must inherit from ISerialBoundaryVertexManager");
            static_assert(std::is_base_of_v<ISerialPartitionManager, TSerialPartitionManager>, "TSerialPartitionManager must inherit from ISerialPartitionManager");
            static_assert(std::is_base_of_v<ISerialDistanceOracle, TSerialDistanceOracle>, "TSerialDistanceOracle must inherit from ISerialDistanceOracle");
            static_assert(std::is_base_of_v<ISerialQuotientGraph, TSerialQuotientGraph>, "TSerialQuotientGraph must inherit from ISerialQuotientGraph");

            bool     global_move_occurred = true;
            u64      global_max_iteration = 3;
            for (u64 global_iteration     = 0; global_iteration < global_max_iteration && global_move_occurred; ++global_iteration) {
                global_move_occurred = false;

                for (size_t distance_i = 0; distance_i < distance.size(); ++distance_i) {
                    weight_t dist = distance[distance.size() - 1 - distance_i];

                    bool     local_move_occurred = true;
                    u64      local_max_iteration = local_max_iterations[distance.size() - 1 - distance_i];
                    for (u64 local_iteration     = 0; local_iteration < local_max_iteration && local_move_occurred; ++local_iteration) {
                        mark += 1;
                        local_move_occurred = false;

                        for (vertex_t u: bv_manager) {
                            if (used[u] == mark) { continue; } // we already used u in this iteration

                            weight_t    u_weight = g.get_weight(u);
                            partition_t u_id     = p_manager[u];

                            // make the move that reduces qap the most
                            partition_t best_u_id      = u_id;
                            partition_t gain_0_id      = u_id;
                            weight_t    best_qap_delta = 0;

                            for (const auto &[v, w]: g[u]) {
                                partition_t v_id = p_manager[v];
                                if (v_id == u_id || v_id == best_u_id || v_id == gain_0_id) {
                                    // v_id was already considered, so quik check
                                    continue;
                                }

                                if (p_manager.get_bweight(v_id) + u_weight <= lmax && d_oracle.get(u_id, v_id) == dist) {
                                    weight_t qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);

                                    if (qap_delta == 0) { gain_0_id = v_id; }

                                    if (qap_delta > best_qap_delta) {
                                        best_qap_delta = qap_delta;
                                        best_u_id      = v_id;
                                    }
                                }
                            }

                            if (best_u_id != u_id) {
                                bv_manager.move(g, p_manager, u, u_id, best_u_id);
                                q_graph.move(g, p_manager, u, u_id, best_u_id);
                                p_manager.move(u, u_weight, u_id, best_u_id);
                                used[u] = mark;
                                local_move_occurred = true;
                            } else if (gain_0_id != u_id && dis(gen) < 0.5) {
                                // if no positive gain, then random neutral swaps
                                bv_manager.move(g, p_manager, u, u_id, gain_0_id);
                                q_graph.move(g, p_manager, u, u_id, gain_0_id);
                                p_manager.move(u, u_weight, u_id, gain_0_id);
                                used[u] = mark;
                                local_move_occurred = true;
                            }
                        }
                    }
                    global_move_occurred |= local_move_occurred;
                }
            }
        }
         */

        template<typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle, typename TSerialQuotientGraph>
        void refine(TSerialGraph &g,
                    [[maybe_unused]] TSerialActiveVertexManager &av_manager,
                    TSerialBoundaryVertexManager &bv_manager,
                    TSerialPartitionManager &p_manager,
                    TSerialDistanceOracle &d_oracle,
                    TSerialQuotientGraph &q_graph) {
            static_assert(std::is_base_of_v<ISerialGraph, TSerialGraph>, "TSerialGraph must inherit from ISerialGraph");
            static_assert(std::is_base_of_v<ISerialActiveVertexManager, TSerialActiveVertexManager>, "TSerialActiveVertexManager must inherit from ISerialActiveVertexManager");
            static_assert(std::is_base_of_v<ISerialBoundaryVertexManager, TSerialBoundaryVertexManager>, "TSerialBoundaryVertexManager must inherit from ISerialBoundaryVertexManager");
            static_assert(std::is_base_of_v<ISerialPartitionManager, TSerialPartitionManager>, "TSerialPartitionManager must inherit from ISerialPartitionManager");
            static_assert(std::is_base_of_v<ISerialDistanceOracle, TSerialDistanceOracle>, "TSerialDistanceOracle must inherit from ISerialDistanceOracle");
            static_assert(std::is_base_of_v<ISerialQuotientGraph, TSerialQuotientGraph>, "TSerialQuotientGraph must inherit from ISerialQuotientGraph");

            bool     global_move_occurred = true;
            u64      global_max_i         = 3;
            for (u64 global_i             = 0; global_i < global_max_i && global_move_occurred; ++global_i) {
                global_move_occurred = false;

                mark += 1;
                for (vertex_t u: bv_manager) {
                    if (used[u] == mark) { continue; } // we already used u in this iteration

                    weight_t    u_weight = g.get_weight(u);
                    partition_t u_id     = p_manager[u];

                    // make the move that reduces qap the most
                    partition_t best_u_id      = u_id;
                    partition_t gain_0_id      = u_id;
                    weight_t    best_qap_delta = 0;

                    for (const auto &[v, w]: g[u]) {
                        partition_t v_id = p_manager[v];
                        if (v_id == u_id || v_id == best_u_id || v_id == gain_0_id) {
                            // v_id was already considered, so quik check
                            continue;
                        }

                        if (p_manager.get_bweight(v_id) + u_weight <= lmax) {
                            weight_t qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);

                            if (qap_delta == 0) { gain_0_id = v_id; }
                            if (qap_delta > best_qap_delta) {
                                best_qap_delta = qap_delta;
                                best_u_id      = v_id;
                            }
                        }
                    }

                    if (best_u_id != u_id) {
                        bv_manager.move(g, p_manager, u, u_id, best_u_id);
                        q_graph.move(g, p_manager, u, u_id, best_u_id);
                        p_manager.move(u, u_weight, u_id, best_u_id);
                        used[u] = mark;
                        global_move_occurred = true;
                    } else if (gain_0_id != u_id && dis(gen) < 0.5) {
                        // if no positive gain, then random neutral swaps
                        bv_manager.move(g, p_manager, u, u_id, gain_0_id);
                        q_graph.move(g, p_manager, u, u_id, gain_0_id);
                        p_manager.move(u, u_weight, u_id, gain_0_id);
                        used[u] = mark;
                        global_move_occurred = true;
                    }
                }
            }
        }
    };
}

#endif //HEIPROMAP_LABEL_PROPAGATION_REFINEMENT_FARAJ20_H
