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

            for (size_t i = 0; i < g.n_neighbors(u); ++i) {
                const EdgeW &e = g.neighbor(u, i);

                partition_t v_id = p_manager[e.v];
                weight_t d = d_oracle.get(u_id, v_id);
                qap += (d * e.w);
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

        for (size_t i = 0; i < g.n_neighbors(u); ++i) {
            const EdgeW &e = g.neighbor(u, i);
            partition_t v_id = p_manager[e.v];
            weight_t d = d_oracle.get(u_id, v_id);
            qap += (d * e.w);
        }

        return qap;
    }

    weight_t get_u_qap(IGraph &g,
                       vertex_t u,
                       partition_t id,
                       IPartitionManager &p_manager,
                       IDistanceOracle &d_oracle) {

        weight_t qap = 0;
        for (size_t i = 0; i < g.n_neighbors(u); ++i) {
            const EdgeW &e = g.neighbor(u, i);
            partition_t v_id = p_manager[e.v];
            weight_t d = d_oracle.get(id, v_id);
            qap += (d * e.w);
        }

        return qap;
    }
}
