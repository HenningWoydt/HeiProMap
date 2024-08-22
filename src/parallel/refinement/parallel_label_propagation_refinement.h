#ifndef SERIALPROCESSMAPPING_PARALLEL_LABEL_PROPAGATION_REFINEMENT_H
#define SERIALPROCESSMAPPING_PARALLEL_LABEL_PROPAGATION_REFINEMENT_H

/*

#include <queue>
#include <omp.h>
#include <atomic>

#include "../../utility/definitions.h"
#include "../../utility/macros.h"
#include "../../utility/utils.h"
#include "../../datastructures/graph.h"
#include "../../utility/qap.h"
#include "../../datastructures/distance_oracle.h"
#include "../../datastructures/iterators/boundary_vertex_iterator.h"

namespace HeiProMap {

    class ParallelLabelPropagationRefinement {
    private:
        Graph *p_g = nullptr;
        std::vector<u64> hierarchy;
        std::vector<u64> distance;
        u64 k = 0;
        f64 imbalance = 0;
        u64 lmax = 0;
        DistanceOracle *dist_o = nullptr;

        std::vector<s32> used;
        std::vector<u8> neighbor_changed;
        s32 mark = -1;

        u64 n_threads = 1;

    public:
        ParallelLabelPropagationRefinement() = default;

        void initialize(Graph *t_g,
                        std::vector<u64> &t_hierarchy,
                        std::vector<u64> &t_distance,
                        u64 t_k,
                        f64 t_imbalance,
                        u64 t_lmax,
                        DistanceOracle *t_dist_o,
                        u64 t_n_threads) {
            p_g = t_g;
            hierarchy = t_hierarchy;
            distance = t_distance;
            k = t_k;
            imbalance = t_imbalance;
            lmax = t_lmax;
            dist_o = t_dist_o;

            n_threads = t_n_threads;

            used.resize(t_g->get_n(), -1);
        }

        void refine(PartitionManager &pm) {
            ASSERT(p_g != nullptr);
            Graph &g = *p_g;

            std::vector<vertex_t> best_us(n_threads);
            std::vector<partition_t> best_ids(n_threads);
            std::vector<s64> best_qap_deltas(n_threads);
            std::vector<vertex_t> boundary_vertices;

            u64 max_iterations = 1000;
            bool move_occurred = true;

            for(u64 iteration = 0; iteration < max_iterations && move_occurred; ++iteration) {
                move_occurred = false;
                std::fill(best_qap_deltas.begin(), best_qap_deltas.end(), -std::numeric_limits<s64>::max());

                boundary_vertices.clear();
                for (BoundaryVertexIterator bvi(pm.get_bvm()); bvi.not_end(); bvi.next()) {
                    vertex_t u = bvi.get();
                    boundary_vertices.emplace_back(u);
                }

#pragma omp parallel default(none) shared(g, pm, boundary_vertices, best_us, best_ids, best_qap_deltas) num_threads(n_threads)
                {
                    u64 thread_id = omp_get_thread_num();
                    vertex_t local_best_u;
                    partition_t local_best_id;
                    s64 local_best_qap_delta = -std::numeric_limits<s64>::max();
#pragma omp for
                    for (size_t i = 0; i < boundary_vertices.size(); ++i) {
                        vertex_t u = boundary_vertices[i];
                        if (g.get_vertex_state(u) == 1 && pm.is_boundary_vertex(u)) {
                            // this vertex can be moved
                            partition_t u_id = pm[u];
                            weight_t u_weight = g.get_vertex_weight(u);
                            s64 u_qap = get_u_qap(g, u, pm, *dist_o);

                            for (EdgeW &e: g[u]) {
                                vertex_t v = e.v;
                                partition_t v_id = pm[v];
                                weight_t v_id_weight = pm.get_pweight(v_id);

                                if (u_id != v_id && v_id_weight + u_weight <= lmax) {
                                    // check if v_id is a better option
                                    s64 qap_delta = u_qap - get_u_qap(g, u, v_id, pm, *dist_o);
                                    if (qap_delta > local_best_qap_delta && qap_delta > 0) {
                                        local_best_u = u;
                                        local_best_id = v_id;
                                        local_best_qap_delta = qap_delta;
                                    }
                                }
                            }
                        }

                        best_us[thread_id] = local_best_u;
                        best_ids[thread_id] = local_best_id;
                        best_qap_deltas[thread_id] = local_best_qap_delta;
                    }
                }

                size_t max_idx = argmax(best_qap_deltas);

                if(best_qap_deltas[max_idx] > 0){
                    // change the partition of the best vertex
                    pm.move(best_us[max_idx], best_ids[max_idx]);
                    move_occurred = true;
                }
            }
        }
    };
}

 */

#endif //SERIALPROCESSMAPPING_PARALLEL_LABEL_PROPAGATION_REFINEMENT_H
