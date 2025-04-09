#ifndef HEIPROMAP_ZEROGAINPERTURBATOR_H
#define HEIPROMAP_ZEROGAINPERTURBATOR_H

#include "../serial_definitions_3.h"
#include "../utility/qap.h"

namespace HeiProMap {

    class ZeroGainPerturbator {
    public:
        ZeroGainPerturbator() = default;

        void perturbate(const u64 level,
                        const u64 max_level,
                        const graph_t &g,
                        const d_oracle_t &d_oracle,
                        bv_manager_t &bv_manager,
                        p_manager_t &p_manager,
                        q_graph_t &q_graph,
                        weight_t lmax,
                        RandomEngine &randomEngine) {
            u64      max_iteration = 1;
            for (u64 iteration     = 0; iteration < max_iteration; ++iteration) {
                forall_gu(g, u)
                    {
                        partition_t u_id = p_manager[u];
                        weight_t    u_w  = g.weight(u);

                        forall_guiv(g, u, i, v)
                            {
                                partition_t v_id = p_manager[v];
                                if (u_id == v_id) { continue; }
                                if (p_manager.get_bweight(v_id) + u_w > lmax) { continue; }

                                s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);

                                if (qap_delta >= 0 || ((-1.0 / (f64) (qap_delta)) < randomEngine.get_f64())) {
                                    bv_manager.move(g, p_manager, u, u_id, v_id);
                                    q_graph.move(g, p_manager, u, u_id, v_id);
                                    p_manager.move(u, u_w, u_id, v_id);
                                    break;
                                }
                            }
                        endfor
                    }
                endfor
            }
        }
    };

}

#endif //HEIPROMAP_ZEROGAINPERTURBATOR_H
