#include "qap.h"

namespace HeiProMap {

    weight_t get_qap(IGraph &g,
                     IActiveVertexManager &av_manager,
                     IPartitionManager &p_manager,
                     IDistanceOracle &d_oracle) {
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

    weight_t get_u_qap(IGraph &g,
                       vertex_t u,
                       IPartitionManager &p_manager,
                       IDistanceOracle &d_oracle) {
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

    weight_t get_u_qap(IGraph &g,
                       vertex_t u,
                       partition_t id,
                       IPartitionManager &p_manager,
                       IDistanceOracle &d_oracle) {

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

    s64 get_u_qap_delta(IGraph &g,
                        vertex_t u,
                        partition_t old_id,
                        partition_t new_id,
                        IPartitionManager &p_manager,
                        IDistanceOracle &d_oracle) {
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
