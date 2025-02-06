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

#ifndef HEIPROMAP_QAP_H
#define HEIPROMAP_QAP_H

#include "../../definitions.h"
#include "../../macros.h"
#include "../datastructures/distance_oracle.h"

namespace HeiProMap {
    template <typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle>
    weight_t get_qap(TSerialGraph& g,
                     TSerialActiveVertexManager& av_manager,
                     TSerialPartitionManager& p_manager,
                     TSerialDistanceOracle& d_oracle) {
        static_assert(std::is_base_of_v<ISerialGraph, TSerialGraph>, "TGraph must inherit from IGraph");
        static_assert(std::is_base_of_v<ISerialActiveVertexManager, TSerialActiveVertexManager>, "TActiveVertexManager must inherit from IActiveVertexManager");
        static_assert(std::is_base_of_v<ISerialPartitionManager, TSerialPartitionManager>, "TPartitionManager must inherit from IPartitionManager");
        static_assert(std::is_base_of_v<ISerialDistanceOracle, TSerialDistanceOracle>, "TDistanceOracle must inherit from IDistanceOracle");

        weight_t qap = 0;

        for (vertex_t u : av_manager) {
            ASSERT(av_manager.is_active(u));

            partition_t u_id = p_manager[u];

            for (size_t i = 0; i < g.size(u); ++i) {
                vertex_t v       = g.neighbor(u, i);
                weight_t ew      = g.get_weight(u, i);
                partition_t v_id = p_manager[v];
                weight_t d       = d_oracle.get(u_id, v_id);
                qap += (d * ew);
            }
        }

        return qap;
    }

    template <typename TSerialGraph, typename TSerialPartitionManager, typename TSerialDistanceOracle>
    weight_t get_u_qap(TSerialGraph& g,
                       vertex_t u,
                       TSerialPartitionManager& p_manager,
                       TSerialDistanceOracle& d_oracle) {
        static_assert(std::is_base_of_v<ISerialGraph, TSerialGraph>, "TGraph must inherit from IGraph");
        static_assert(std::is_base_of_v<ISerialPartitionManager, TSerialPartitionManager>, "TPartitionManager must inherit from IPartitionManager");
        static_assert(std::is_base_of_v<ISerialDistanceOracle, TSerialDistanceOracle>, "TDistanceOracle must inherit from IDistanceOracle");

        weight_t qap     = 0;
        partition_t u_id = p_manager[u];

        for (size_t i = 0; i < g.size(u); ++i) {
            vertex_t v       = g.neighbor(u, i);
            weight_t ew      = g.get_weight(u, i);
            partition_t v_id = p_manager[v];
            weight_t d       = d_oracle.get(u_id, v_id);
            qap += (d * ew);
        }

        return qap;
    }

    template <typename TSerialGraph, typename TSerialPartitionManager, typename TSerialDistanceOracle>
    weight_t get_u_qap(TSerialGraph& g,
                       vertex_t u,
                       partition_t id,
                       TSerialPartitionManager& p_manager,
                       TSerialDistanceOracle& d_oracle) {
        static_assert(std::is_base_of_v<ISerialGraph, TSerialGraph>, "TGraph must inherit from IGraph");
        static_assert(std::is_base_of_v<ISerialPartitionManager, TSerialPartitionManager>, "TPartitionManager must inherit from IPartitionManager");
        static_assert(std::is_base_of_v<ISerialDistanceOracle, TSerialDistanceOracle>, "TDistanceOracle must inherit from IDistanceOracle");

        weight_t qap = 0;
        for (size_t i = 0; i < g.size(u); ++i) {
            vertex_t v       = g.neighbor(u, i);
            weight_t ew      = g.get_weight(u, i);
            partition_t v_id = p_manager[v];
            weight_t d       = d_oracle.get(id, v_id);
            qap += (d * ew);
        }

        return qap;
    }

    template <typename TSerialGraph, typename TSerialPartitionManager, typename TSerialDistanceOracle>
    s64 get_u_qap_delta_by_idx(TSerialGraph& g,
                            vertex_t u,
                            partition_t old_id,
                            partition_t new_id,
                            TSerialPartitionManager& p_manager,
                            TSerialDistanceOracle& d_oracle) {
        static_assert(std::is_base_of_v<ISerialGraph, TSerialGraph>, "TGraph must inherit from IGraph");
        static_assert(std::is_base_of_v<ISerialPartitionManager, TSerialPartitionManager>, "TPartitionManager must inherit from IPartitionManager");
        static_assert(std::is_base_of_v<ISerialDistanceOracle, TSerialDistanceOracle>, "TDistanceOracle must inherit from IDistanceOracle");

        s64 qap_delta = 0;
        for (size_t i = 0; i < g.size(u); ++i) {
            vertex_t v       = g.neighbor(u, i);
            weight_t w       = g.get_weight(u, i);
            partition_t v_id = p_manager[v];

            weight_t old_d, new_d;
            d_oracle.get(v_id, old_id, new_id, old_d, new_d);
            qap_delta += (old_d - new_d) * w;
        }

        return qap_delta;
    }

    template <typename TSerialGraph, typename TSerialPartitionManager, typename TSerialDistanceOracle>
    s64 get_u_qap_delta(TSerialGraph& g,
                        vertex_t u,
                        partition_t old_id,
                        partition_t new_id,
                        TSerialPartitionManager& p_manager,
                        TSerialDistanceOracle& d_oracle) {
        static_assert(std::is_base_of_v<ISerialGraph, TSerialGraph>, "TGraph must inherit from IGraph");
        static_assert(std::is_base_of_v<ISerialPartitionManager, TSerialPartitionManager>, "TPartitionManager must inherit from IPartitionManager");
        static_assert(std::is_base_of_v<ISerialDistanceOracle, TSerialDistanceOracle>, "TDistanceOracle must inherit from IDistanceOracle");

        s64 qap_delta = 0;
        for (const auto& [v, w] : g[u]) {
            partition_t v_id = p_manager[v];

            weight_t old_d, new_d;
            d_oracle.get(v_id, old_id, new_id, old_d, new_d);
            qap_delta += (old_d - new_d) * w;
        }

        return qap_delta;
    }
}

#endif //HEIPROMAP_QAP_H
