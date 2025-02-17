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
        u64 max_iteration = 2; // how many iterations to run the algorithm at most
    };

    partition_t get_island_id(partition_t u_id, const std::vector<std::vector<partition_t>> &island_ids) {
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
     * If moves between the islands have been found, then try to distribute it onto the individual partzitions.
     */
    class HierarchyAwareCycleRefinement final : public ISerialRefiner {
        std::vector<partition_t> hierarchy;
        std::vector<weight_t> distance;
        partition_t k = 0;
        weight_t lmax = 0;

        std::vector<s32> used;
        s32 mark = -1;

        std::random_device rd;
        std::mt19937 gen;
        std::uniform_real_distribution<float> dis;

    public:
        HierarchyAwareCycleRefinement() : gen(rd()), dis(0.0f, 1.0f) {}

        void initialize(const vertex_t n,
                        std::vector<partition_t>& t_hierarchy,
                        std::vector<weight_t>& t_distance,
                        const weight_t t_lmax) override {
            hierarchy = t_hierarchy;
            distance  = t_distance;
            k         = prod<partition_t>(hierarchy);
            lmax      = t_lmax;

            used.resize(n, -1);
        }

        template <typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle, typename TSerialQuotientGraph>
        void refine([[maybe_unused]] HierarchyAwareCyclesConfiguration& config,
                    [[maybe_unused]] TSerialGraph& g,
                    [[maybe_unused]] TSerialActiveVertexManager& av_manager,
                    [[maybe_unused]] TSerialBoundaryVertexManager& bv_manager,
                    [[maybe_unused]] TSerialPartitionManager& p_manager,
                    [[maybe_unused]] TSerialDistanceOracle& d_oracle,
                    [[maybe_unused]] TSerialQuotientGraph& q_graph) {
            static_assert(std::is_base_of_v<ISerialGraph, TSerialGraph>, "TSerialGraph must inherit from ISerialGraph");
            static_assert(std::is_base_of_v<ISerialActiveVertexManager, TSerialActiveVertexManager>, "TSerialActiveVertexManager must inherit from ISerialActiveVertexManager");
            static_assert(std::is_base_of_v<ISerialBoundaryVertexManager, TSerialBoundaryVertexManager>, "TSerialBoundaryVertexManager must inherit from ISerialBoundaryVertexManager");
            static_assert(std::is_base_of_v<ISerialPartitionManager, TSerialPartitionManager>, "TSerialPartitionManager must inherit from ISerialPartitionManager");
            static_assert(std::is_base_of_v<ISerialDistanceOracle, TSerialDistanceOracle>, "TSerialDistanceOracle must inherit from ISerialDistanceOracle");
            static_assert(std::is_base_of_v<ISerialQuotientGraph, TSerialQuotientGraph>, "TSerialQuotientGraph must inherit from ISerialQuotientGraph");

            size_t n_islands = hierarchy.back();
            size_t ids_per_island = k / n_islands;
            std::vector<std::vector<partition_t>> island_ids(n_islands, std::vector<partition_t>(ids_per_island));
            std::vector<weight_t> islands_weight(n_islands, 0.0);
            weight_t island_lmax = (weight_t) ids_per_island * lmax;

            size_t id = 0;
            for (size_t i = 0; i < n_islands; i++) {
                for (size_t j = 0; j < ids_per_island; j++) {
                    island_ids[i][j] = id;
                    islands_weight[i] += p_manager.get_bweight(id);
                    id += 1;
                }
            }

            for (u64 iteration = 0; iteration < config.max_iteration; ++iteration) {

                // insert all boundary vertices
                for (vertex_t u : bv_manager) {
                    partition_t u_id  = p_manager[u];
                    partition_t u_island_id = get_island_id(u_id, island_ids);
                    weight_t u_weight = g.get_weight(u);

                    // find all connected partitions to u
                    for (const auto& [v, w] : g[u]) {
                        // for (size_t idx = 0; idx < g.size(u); ++idx) {
                        // const vertex_t v = g.neighbor(u, idx);
                        // const weight_t w = g.get_weight(u, idx);
                        partition_t v_id = p_manager[v];
                        partition_t v_island_id = get_island_id(v_id, island_ids);

                        if (u_island_id != v_island_id) {

                        }
                    }
                }

            }

        }
    };
}

#endif //HEIPROMAP_HIERARCHY_AWARE_CYCLE_REFINEMENT_H
