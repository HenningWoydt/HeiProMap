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

#ifndef HEIPROMAP_SIMPLE_REBALANCER_H
#define HEIPROMAP_SIMPLE_REBALANCER_H

#include "../datastructures/distance_oracle.h"
#include "../datastructures/functions.h"
#include "../interfaces/ISerialActiveVertexManager.h"
#include "../interfaces/ISerialBoundaryVertexManager.h"
#include "../interfaces/ISerialQuotientGraph.h"

namespace HeiProMap {
    struct SimpleRebalancerConfiguration {};

    /**
     * Aims to rebalance the partition, such that it adheres to the balancing constraint.
     */
    class SimpleRebalancer {
        vertex_t m_n;
        partition_t m_k;
        weight_t m_lmax;

    public:
        SimpleRebalancer(vertex_t t_n, partition_t t_k, weight_t t_lmax) {
            m_n    = t_n;
            m_k    = t_k;
            m_lmax = t_lmax;
        }

        template <typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle, typename TSerialQuotientGraph>
        void rebalance([[maybe_unused]] SimpleRebalancerConfiguration& config,
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

            for (partition_t id = 0; id < m_k; ++id) {
                if (p_manager.get_bweight(id) > m_lmax) {
                    std::cout << "overloaded: " << id << " " << p_manager.get_bweight(id) << std::endl;
                }
            }
        }
    };
}

#endif //HEIPROMAP_SIMPLE_REBALANCER_H
