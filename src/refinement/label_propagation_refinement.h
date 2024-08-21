#ifndef SERIALPROCESSMAPPING_LABEL_PROPAGATION_REFINEMENT_H
#define SERIALPROCESSMAPPING_LABEL_PROPAGATION_REFINEMENT_H

#include <queue>

#include "../utility/definitions.h"
#include "../utility/macros.h"
#include "../utility/utils.h"
#include "../datastructures/graph.h"
#include "../utility/qap.h"
#include "../datastructures/iterators/active_vertex_iterator.h"
#include "../datastructures/distance_oracle.h"
#include "../datastructures/iterators/boundary_vertex_iterator.h"

namespace SPM {

    class LabelPropagationRefinement {
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

    public:
        LabelPropagationRefinement() = default;

        void initialize(Graph *t_g,
                        std::vector<u64> &t_hierarchy,
                        std::vector<u64> &t_distance,
                        u64 t_k,
                        f64 t_imbalance,
                        u64 t_lmax,
                        DistanceOracle *t_dist_o) {
            p_g = t_g;
            hierarchy = t_hierarchy;
            distance = t_distance;
            k = t_k;
            imbalance = t_imbalance;
            lmax = t_lmax;
            dist_o = t_dist_o;

            used.resize(t_g->get_n(), -1);
        }

        void refine_basic(PartitionManager &pm) {
            ASSERT(p_g != nullptr);
            Graph &g = *p_g;

            bool move_occurred = true;
            u64 max_iteration = 10;
            for (u64 iteration = 0; iteration < max_iteration && move_occurred; ++iteration) {
                mark += 1;
                move_occurred = false;

                for (BoundaryVertexIterator bvi(pm.get_bvm()); bvi.not_end(); bvi.next()) {
                    vertex_t u = bvi.get();
                    if (used[u] == mark) {
                        // we already used u in this iteration
                        continue;
                    }

                    weight_t u_weight = g.get_vertex_weight(u);
                    partition_t old_u_id = pm[u];

                    u64 best_qap = std::numeric_limits<u64>::max();
                    partition_t best_u_id = old_u_id;

                    // make the move that reduces qap the most
                    for (EdgeW &e: g[u]) {
                        vertex_t v = e.v;
                        partition_t v_id = pm[v];

                        // only different partitions and do not overload
                        if (old_u_id != v_id && pm.get_pweight(v_id) + u_weight <= lmax) {

                            if (best_qap == std::numeric_limits<u64>::max()) {
                                // initialize with current qap
                                best_qap = get_u_qap(g, u, pm, *dist_o);
                            }

                            // determine new qap
                            u64 new_qap = get_u_qap(g, u, v_id, pm, *dist_o);
                            if (new_qap < best_qap) {
                                best_qap = new_qap;
                                best_u_id = v_id;
                            }
                        }
                    }

                    if (old_u_id != best_u_id) {
                        pm.move(u, best_u_id);
                        used[u] = mark;
                        move_occurred = true;
                    }
                }
            }
        }

        void refine_pq(PartitionManager &pm) {
            ASSERT(p_g != nullptr);
            Graph &g = *p_g;

            bool move_occurred = true;
            u64 max_iteration = 10;
            for (u64 iteration = 0; iteration < max_iteration && move_occurred; ++iteration) {
                move_occurred = false;

                std::priority_queue<Move> pq;
                // Step 1: Calculate all possible moves for boundary vertices and put them into the priority queue
                for (BoundaryVertexIterator bvi(pm.get_bvm()); bvi.not_end(); bvi.next()) {
                    vertex_t u = bvi.get();
                    weight_t u_weight = g.get_vertex_weight(u);
                    vertex_t u_id = pm[u];
                    s64 qap = get_u_qap(g, u, pm, *dist_o);

                    ASSERT(g.get_vertex_state(u) == 1);

                    bool found = false;
                    Move m(u, 0, -std::numeric_limits<s64>::max());
                    for (EdgeW &e: g[u]) {
                        vertex_t v = e.v;
                        vertex_t v_id = pm[v];

                        if (u_id != v_id && pm.get_pweight(v_id) + u_weight <= lmax) {
                            pm[u] = v_id;
                            s64 qap_delta = qap - get_u_qap(g, u, pm, *dist_o);
                            if (qap_delta >= 0 && qap_delta > m.qap_delta) {
                                m.p_id = v_id;
                                m.qap_delta = qap_delta;
                                found = true;
                                break;
                            }
                        }
                    }
                    pm[u] = u_id;

                    if (found) {
                        pq.push(m);
                    }
                }

                // Step 2: Process all items in the queue, and always execute the move with the best gain
                while (!pq.empty()) {
                    Move m = pq.top();
                    pq.pop();

                    vertex_t u = m.u;
                    weight_t u_weight = g.get_vertex_weight(u);
                    vertex_t u_id = pm[u];
                    s64 qap = get_u_qap(g, u, pm, *dist_o);

                    m.qap_delta = -std::numeric_limits<s64>::max();
                    bool found = false;
                    for (EdgeW &e: g[u]) {
                        vertex_t v = e.v;
                        vertex_t v_id = pm[v];

                        if (u_id != v_id && pm.get_pweight(v_id) + u_weight <= lmax) {
                            pm[u] = v_id;
                            s64 qap_delta = qap - get_u_qap(g, u, pm, *dist_o);
                            if (qap_delta >= 0 && qap_delta > m.qap_delta) {
                                m.p_id = v_id;
                                m.qap_delta = qap_delta;
                                found = true;
                            }
                        }
                    }
                    pm[u] = u_id;

                    if (!found) {
                        // no valid moves, go to next move
                        continue;
                    }

                    // compare with top in the queue and insert if smaller
                    if (!pq.empty() && m.qap_delta < pq.top().qap_delta) {
                        pq.push(m);
                        continue;
                    }

                    // execute move
                    move_occurred = true;
                    pm.move(u, m.p_id);
                }
            }

        }

        void refine_exhaustive(PartitionManager &pm) {
            ASSERT(p_g != nullptr);
            Graph &g = *p_g;

            // vector to mark which entries need to be updated
            std::vector<u8> to_update(g.get_n(), 1);
            std::vector<MovePQ> moves(g.get_n());

            bool move_occurred = true;
            while (move_occurred) {
                move_occurred = false;

                s64 best_qap_delta = -std::numeric_limits<s64>::max();
                vertex_t best_u;
                vertex_t old_u_id;
                vertex_t new_u_id;

                // check all boundary vertices, to get the global best move
                for (BoundaryVertexIterator bvi(pm.get_bvm()); bvi.not_end(); bvi.next()) {
                    vertex_t u = bvi.get();
                    weight_t u_weight = g.get_vertex_weight(u);
                    vertex_t u_id = pm[u];

                    if (to_update[u] == 1 || true) {
                        s64 curr_qap = get_u_qap(g, u, pm, *dist_o);
                        moves[u].qap_delta = -std::numeric_limits<s64>::max();
                        for (EdgeW &e: g[u]) {
                            vertex_t v = e.v;
                            vertex_t v_id = pm[v];

                            bool different_partition = u_id != v_id;
                            bool no_overload = pm.get_pweight(v_id) + u_weight <= lmax;

                            if (different_partition && no_overload) {
                                pm[u] = v_id;
                                s64 delta = curr_qap - get_u_qap(g, u, pm, *dist_o);
                                if (delta >= 0 && delta > moves[u].qap_delta) {
                                    moves[u].p_id = v_id;
                                    moves[u].qap_delta = delta;
                                }
                            }
                        }
                        pm[u] = u_id; // reset partition of u
                        to_update[u] = 0; // u was updated
                    }

                    if (moves[u].qap_delta > best_qap_delta) {
                        best_qap_delta = moves[u].qap_delta;
                        best_u = u;
                        old_u_id = pm[u];
                        new_u_id = moves[u].p_id;
                        move_occurred = true;
                    }
                }

                if (move_occurred) {
                    // change the partition of the best vertex
                    pm.move(best_u, new_u_id);

                    // all neighbors need to be updated in the next iteration
                    for (EdgeW &e: g[best_u]) {
                        to_update[e.v] = 1;
                    }
                }
            }
        }

        void refine(PartitionManager &pm) {
            ASSERT(p_g != nullptr);
            Graph &g = *p_g;

            std::vector<partition_t> moves_to_check;
            std::vector<u64> qaps;

            bool global_move_occurred = true;
            u64 global_max_iteration = 2;
            for (u64 global_iteration = 0; global_iteration < global_max_iteration && global_move_occurred; ++global_iteration) {
                global_move_occurred = false;

                std::vector<u64> local_max_iterations = {1, 3, 5};
                for (size_t i = 0; i < distance.size(); ++i) {
                    u64 dist = distance[distance.size() - 1 - i];

                    bool local_move_occurred = true;
                    for (u64 local_iteration = 0; local_iteration < local_max_iterations[distance.size() - 1 - i] && local_move_occurred; ++local_iteration) {
                        mark += 1;
                        local_move_occurred = false;

                        for (BoundaryVertexIterator bvi(pm.get_bvm()); bvi.not_end(); bvi.next()) {
                            vertex_t u = bvi.get();
                            if (used[u] == mark) {
                                // we already used u in this iteration
                                continue;
                            }

                            weight_t u_weight = g.get_vertex_weight(u);
                            vertex_t u_id = pm[u];

                            // make the move that reduces qap the most
                            moves_to_check.clear();
                            for (EdgeW &e: g[u]) {
                                vertex_t v = e.v;
                                vertex_t v_id = pm[v];

                                if (u_id != v_id && pm.get_pweight(v_id) + u_weight <= lmax && (*dist_o).get(u_id, v_id) == dist) {
                                    // initialize with current qap
                                    moves_to_check.push_back(v_id);
                                }
                            }

                            vertex_t best_u_id = u_id;
                            if(!moves_to_check.empty()){
                                u64 curr_qap = get_u_qap(g, u, pm, *dist_o);
                                get_u_qap(g, u, moves_to_check, qaps, pm, *dist_o);
                                size_t best_idx = argmin(qaps);
                                if(qaps[best_idx] < curr_qap){
                                    best_u_id = moves_to_check[best_idx];
                                }
                            }

                            if (best_u_id != u_id) {
                                pm.move(u, best_u_id);
                                used[u] = mark;
                                local_move_occurred = true;
                            }
                        }
                    }
                    global_move_occurred |= local_move_occurred;
                }
            }
        }
    };
}

#endif //SERIALPROCESSMAPPING_LABEL_PROPAGATION_REFINEMENT_H
