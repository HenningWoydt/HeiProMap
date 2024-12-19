#ifndef HEIDELBERGPROCESSMAPPING_QAP_H
#define HEIDELBERGPROCESSMAPPING_QAP_H

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
        for (size_t i = 0; i < g.size(u); ++i) {
            vertex_t v       = g.neighbor(u, i);
            weight_t ew      = g.get_weight(u, i);
            partition_t v_id = p_manager[v];

            weight_t old_d = d_oracle.get(old_id, v_id);
            weight_t new_d = d_oracle.get(new_id, v_id);
            qap_delta += (old_d * ew) - (new_d * ew);
        }

        return qap_delta;
    }
}

#endif //HEIDELBERGPROCESSMAPPING_QAP_H
