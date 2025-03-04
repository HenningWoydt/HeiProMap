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

#ifndef HEIPROMAP_HIERARCHY_AWARE_CYCLE_REFINEMENT_H
#define HEIPROMAP_HIERARCHY_AWARE_CYCLE_REFINEMENT_H

#include <algorithm>
#include <random>

#include "../datastructures/distance_oracle.h"
#include "../datastructures/functions.h"
#include "../interfaces/ISerialActiveVertexManager.h"
#include "../interfaces/ISerialBoundaryVertexManager.h"
#include "../interfaces/ISerialQuotientGraph.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/qap.h"
#include "../utility/utils.h"

namespace HeiProMap {
    struct HierarchyAwareCyclesConfiguration {
        u64 max_iteration = 1; // how many iterations to run the algorithm at most
    };

    inline partition_t get_island_id(partition_t u_id, const std::vector<std::vector<partition_t>> &island_ids) {
        for (partition_t i = 0; i < island_ids.size(); ++i) {
            if (std::find(island_ids[i].begin(), island_ids[i].end(), u_id) != island_ids[i].end()) {
                return i;
            }
        }
        std::cout << "Error u_id " << u_id << " not found!" << std::endl;
        exit(EXIT_FAILURE);
    }

    /**
     * Since the top level of the hierarchy is the most important, try to optimize it the most.
     * Aggregate all partitions of the islands and then try to find moves between the islands instead of individual partitions.
     * If moves between the islands have been found, then try to distribute it onto the individual partitions.
     */
    class HierarchyAwareCycleRefinement final : public ISerialRefiner {
        vertex_t                 m_n    = 0;
        vertex_t                 m_m    = 0;
        partition_t              m_k    = 0;
        weight_t                 m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t>    m_distance;
        u64                      m_seed = 0;

        std::vector<s32> used;
        s32              mark = -1;

        std::random_device                    rd;
        std::mt19937                          gen;
        std::uniform_real_distribution<float> dis;

    public:
        HierarchyAwareCycleRefinement() : gen(rd()), dis(0.0f, 1.0f) {}

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const weight_t t_lmax,
                        const std::vector<partition_t> &t_hierarchy,
                        const std::vector<weight_t> &t_distance,
                        const u64 t_seed) override {
            m_n         = t_n;
            m_m         = t_m;
            m_k         = t_k;
            m_lmax      = t_lmax;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;
            m_seed      = t_seed;

            used.resize(t_n, -1);
        }

        template<typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle, typename TSerialQuotientGraph>
        void refine([[maybe_unused]] HierarchyAwareCyclesConfiguration &config,
                    [[maybe_unused]] TSerialGraph &g,
                    [[maybe_unused]] TSerialActiveVertexManager &av_manager,
                    [[maybe_unused]] TSerialBoundaryVertexManager &bv_manager,
                    [[maybe_unused]] TSerialPartitionManager &p_manager,
                    [[maybe_unused]] TSerialDistanceOracle &d_oracle,
                    [[maybe_unused]] TSerialQuotientGraph &q_graph) {
            static_assert(std::is_base_of_v<ISerialGraph, TSerialGraph>, "TSerialGraph must inherit from ISerialGraph");
            static_assert(std::is_base_of_v<ISerialActiveVertexManager, TSerialActiveVertexManager>, "TSerialActiveVertexManager must inherit from ISerialActiveVertexManager");
            static_assert(std::is_base_of_v<ISerialBoundaryVertexManager, TSerialBoundaryVertexManager>, "TSerialBoundaryVertexManager must inherit from ISerialBoundaryVertexManager");
            static_assert(std::is_base_of_v<ISerialPartitionManager, TSerialPartitionManager>, "TSerialPartitionManager must inherit from ISerialPartitionManager");
            static_assert(std::is_base_of_v<ISerialDistanceOracle, TSerialDistanceOracle>, "TSerialDistanceOracle must inherit from ISerialDistanceOracle");
            static_assert(std::is_base_of_v<ISerialQuotientGraph, TSerialQuotientGraph>, "TSerialQuotientGraph must inherit from ISerialQuotientGraph");

            auto sp = std::chrono::high_resolution_clock::now();

            size_t                                n_islands      = m_hierarchy.back();
            size_t                                ids_per_island = m_k / n_islands;
            std::vector<std::vector<partition_t>> island_ids(n_islands, std::vector<partition_t>(ids_per_island));
            std::vector<weight_t>                 islands_weight(n_islands, 0.0);
            weight_t                              island_lmax    = (weight_t) ids_per_island * m_lmax;

            size_t      id = 0;
            for (size_t i  = 0; i < n_islands; i++) {
                for (size_t j = 0; j < ids_per_island; j++) {
                    island_ids[i][j] = id;
                    islands_weight[i] += p_manager.get_bweight(id);
                    id += 1;
                }
            }

            for (u64 iteration = 0; iteration < config.max_iteration; ++iteration) {
                std::vector<Move> unavailable_moves;

                // insert all boundary vertices
                for (vertex_t u: bv_manager) {
                    partition_t u_id = p_manager[u];

                    // make the move that reduces qap the most
                    partition_t best_id        = u_id;
                    weight_t    best_id_weight = 0;
                    s64         best_qap_delta = 0;
                    f32         counter        = 0;

                    // find all connected partitions to u
                    for (const auto [v, w]: g[u]) {
                        partition_t v_id        = p_manager[v];
                        weight_t    v_id_weight = p_manager.get_bweight(v_id);

                        if (u_id == v_id) { continue; }

                        s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);

                        if (qap_delta > best_qap_delta || (qap_delta == best_qap_delta && v_id_weight < best_id_weight)) {
                            best_id        = v_id;
                            best_id_weight = v_id_weight;
                            best_qap_delta = qap_delta;
                            counter        = 1.0;
                        } else if (qap_delta == best_qap_delta && qap_delta != 0) {
                            counter += 1.0;
                            // choose with probability 1/counter as it ensures uniform distribution
                            if (dis(gen) < 1.0f / counter) {
                                best_id = v_id;
                            }
                        }
                    }

                    if (best_id != u_id) {
                        partition_t u_island_id    = get_island_id(u_id, island_ids);
                        partition_t best_island_id = get_island_id(best_id, island_ids);
                        if (u_island_id == best_island_id) { continue; }
                        weight_t u_weight = g.get_weight(u);

                        bool different_island       = u_island_id != best_island_id;
                        bool would_overload_best_id = best_id_weight + u_weight > m_lmax;
                        bool would_overload_island  = islands_weight[best_island_id] + u_weight > island_lmax;

                        // std::cout << u << " " << u_id << " " << u_island_id << " " << u_weight << " " << best_island_id << " " << best_id_weight << " " << islands_weight[best_island_id] << " " << island_lmax << std::endl;
                        // std::cout << different_island << " " << would_overload_v_id << " " << would_overload_island << std::endl;

                        if (different_island && would_overload_best_id && !would_overload_island) {
                            unavailable_moves.emplace_back(u, best_id, best_qap_delta);
                        }
                    }
                }

                s64 total_qap_delta = 0;
                for (const auto &move: unavailable_moves) {
                    total_qap_delta += move.qap_delta;
                }
                auto ep      = std::chrono::high_resolution_clock::now();
                f64  seconds = get_seconds(sp, ep);
                std::cout << "Hierarchy Aware Cycles: " << unavailable_moves.size() << " moves unavailable, but would improve by " << total_qap_delta << " in " << seconds << " seconds!" << std::endl;
            }
        }
    };
}

#endif //HEIPROMAP_HIERARCHY_AWARE_CYCLE_REFINEMENT_H
