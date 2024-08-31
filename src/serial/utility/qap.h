#ifndef SERIALPROCESSMAPPING_QAP_H
#define SERIALPROCESSMAPPING_QAP_H

#include "../../definitions.h"
#include "../../macros.h"
#include "utils.h"
#include "../datastructures/graph.h"
#include "../datastructures/translation_table.h"
#include "../datastructures/distance_oracle.h"
#include "../datastructures/partition_manager.h"

namespace HeiProMap {
    template<typename TGraph, typename TActiveVertexManager, typename TPartitionManager, typename TDistanceOracle>
    weight_t get_qap(TGraph &g,
                     TActiveVertexManager &av_manager,
                     TPartitionManager &p_manager,
                     TDistanceOracle &d_oracle){
        static_assert(std::is_base_of<IGraph, TGraph>::value, "TGraph must inherit from IGraph");
        static_assert(std::is_base_of<IActiveVertexManager, TActiveVertexManager>::value, "TActiveVertexManager must inherit from IActiveVertexManager");
        static_assert(std::is_base_of<IPartitionManager, TPartitionManager>::value, "TPartitionManager must inherit from IPartitionManager");
        static_assert(std::is_base_of<IDistanceOracle, TDistanceOracle>::value, "TDistanceOracle must inherit from IDistanceOracle");

        weight_t qap = 0;

        for (av_manager.reset_iterator(); av_manager.available(); av_manager.next()) {
            vertex_t u = av_manager.get();
            ASSERT(av_manager.is_active(u));

            partition_t u_id = p_manager[u];

            for (size_t i = 0; i < g.size(u); ++i) {
                vertex_t v = g.neighbor(u, i);
                weight_t ew = g.get_weight(u, i);
                partition_t v_id = p_manager[v];
                weight_t d = d_oracle.get(u_id, v_id);
                qap += (d * ew);
            }
        }

        return qap;
    }

    template<typename TGraph, typename TPartitionManager, typename TDistanceOracle>
    weight_t get_u_qap(TGraph &g,
                       vertex_t u,
                       TPartitionManager &p_manager,
                       TDistanceOracle &d_oracle){
        static_assert(std::is_base_of<IGraph, TGraph>::value, "TGraph must inherit from IGraph");
        static_assert(std::is_base_of<IPartitionManager, TPartitionManager>::value, "TPartitionManager must inherit from IPartitionManager");
        static_assert(std::is_base_of<IDistanceOracle, TDistanceOracle>::value, "TDistanceOracle must inherit from IDistanceOracle");

        weight_t qap = 0;
        partition_t u_id = p_manager[u];

        for (size_t i = 0; i < g.size(u); ++i) {
            vertex_t v = g.neighbor(u, i);
            weight_t ew = g.get_weight(u, i);
            partition_t v_id = p_manager[v];
            weight_t d = d_oracle.get(u_id, v_id);
            qap += (d * ew);
        }

        return qap;
    }

    template<typename TGraph, typename TPartitionManager, typename TDistanceOracle>
    weight_t get_u_qap(TGraph &g,
                       vertex_t u,
                       partition_t id,
                       TPartitionManager &p_manager,
                       TDistanceOracle &d_oracle){
        static_assert(std::is_base_of<IGraph, TGraph>::value, "TGraph must inherit from IGraph");
        static_assert(std::is_base_of<IPartitionManager, TPartitionManager>::value, "TPartitionManager must inherit from IPartitionManager");
        static_assert(std::is_base_of<IDistanceOracle, TDistanceOracle>::value, "TDistanceOracle must inherit from IDistanceOracle");

        weight_t qap = 0;
        for (size_t i = 0; i < g.size(u); ++i) {
            vertex_t v = g.neighbor(u, i);
            weight_t ew = g.get_weight(u, i);
            partition_t v_id = p_manager[v];
            weight_t d = d_oracle.get(id, v_id);
            qap += (d * ew);
        }

        return qap;
    }

    template<typename TGraph, typename TPartitionManager, typename TDistanceOracle>
    s64 get_u_qap_delta(TGraph &g,
                        vertex_t u,
                        partition_t old_id,
                        partition_t new_id,
                        TPartitionManager &p_manager,
                        TDistanceOracle &d_oracle) {
        static_assert(std::is_base_of<IGraph, TGraph>::value, "TGraph must inherit from IGraph");
        static_assert(std::is_base_of<IPartitionManager, TPartitionManager>::value, "TPartitionManager must inherit from IPartitionManager");
        static_assert(std::is_base_of<IDistanceOracle, TDistanceOracle>::value, "TDistanceOracle must inherit from IDistanceOracle");

        s64 qap_delta = 0;
        for (size_t i = 0; i < g.size(u); ++i) {
            vertex_t v = g.neighbor(u, i);
            weight_t ew = g.get_weight(u, i);
            partition_t v_id = p_manager[v];

            weight_t old_d = d_oracle.get(old_id, v_id);
            weight_t new_d = d_oracle.get(new_id, v_id);
            qap_delta += (old_d * ew) - (new_d * ew);
        }

        return qap_delta;
    }
}

#endif //SERIALPROCESSMAPPING_QAP_H
