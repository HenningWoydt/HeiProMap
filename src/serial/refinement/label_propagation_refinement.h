#ifndef HEIDELBERGPROCESSMAPPING_LABEL_PROPAGATION_REFINEMENT_H
#define HEIDELBERGPROCESSMAPPING_LABEL_PROPAGATION_REFINEMENT_H

#include <queue>
#include <random>

#include "../../definitions.h"
#include "../../macros.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/utils.h"

namespace HeiProMap {
    class LabelPropagationRefinement final : public ISerialRefiner {
        std::vector<partition_t> hierarchy;
        std::vector<weight_t> distance;
        partition_t k = 0;
        weight_t lmax = 0;

        std::vector<s32> used;
        std::vector<u8> neighbor_changed;
        s32 mark = -1;

        std::vector<u64> local_max_iterations = {1, 3, 9};

        std::random_device rd;
        std::mt19937 gen;
        std::uniform_real_distribution<float> dis;

    public:
        LabelPropagationRefinement() : gen(rd()), dis(0.0f, 1.0f) {
        }

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

        /*
        void refine_basic(PartitionManager &pm) {
            ASSERT(m_p_g != nullptr);
            Graph &g = *m_p_g;

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
            ASSERT(m_p_g != nullptr);
            Graph &g = *m_p_g;

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
            ASSERT(m_p_g != nullptr);
            Graph &g = *m_p_g;

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
        */

        template <typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle>
        void refine(TSerialGraph& g,
                    TSerialActiveVertexManager& av_manager,
                    TSerialBoundaryVertexManager& bv_manager,
                    TSerialPartitionManager& p_manager,
                    TSerialDistanceOracle& d_oracle) {
            bool global_move_occurred = true;
            u64 global_max_iteration = 3;
            for (u64 global_iteration = 0; global_iteration < global_max_iteration && global_move_occurred; ++global_iteration) {
                global_move_occurred = false;

                for (size_t distance_i = 0; distance_i < distance.size(); ++distance_i) {
                    weight_t dist = distance[distance.size() - 1 - distance_i];

                    bool local_move_occurred = true;
                    for (u64 local_iteration = 0; local_iteration < local_max_iterations[distance.size() - 1 - distance_i] && local_move_occurred; ++local_iteration) {
                        mark += 1;
                        local_move_occurred = false;

                        std::vector<partition_t> gain_0_ids;
                        for (vertex_t u : bv_manager) {
                            if (used[u] == mark) { continue; } // we already used u in this iteration

                            gain_0_ids.clear(); // clear old ids

                            weight_t u_weight = g.get_weight(u);
                            partition_t u_id = p_manager[u];
                            weight_t u_qap = std::numeric_limits<weight_t>::max();

                            // make the move that reduces qap the most
                            partition_t best_u_id = u_id;
                            weight_t best_qap_delta = 0;
                            for (size_t i = 0; i < g.size(u); ++i) {
                                vertex_t v = g.neighbor(u, i);
                                partition_t v_id = p_manager[v];

                                if (u_id != v_id && p_manager.get_bweight(v_id) + u_weight <= lmax && d_oracle.get(u_id, v_id) == dist) {

                                    if (u_qap == std::numeric_limits<weight_t>::max()) { u_qap = get_u_qap(g, u, p_manager, d_oracle); }

                                    weight_t qap_delta = u_qap - get_u_qap(g, u, v_id, p_manager, d_oracle);
                                    if (qap_delta == 0) {
                                        gain_0_ids.emplace_back(v_id);
                                    }
                                    if (qap_delta > best_qap_delta) {
                                        best_qap_delta = qap_delta;
                                        best_u_id = v_id;
                                    }
                                }
                            }

                            if (best_u_id != u_id) {
                                bv_manager.move(g, p_manager, u, u_id, best_u_id);
                                p_manager.move(u, u_weight, u_id, best_u_id);
                                used[u] = mark;
                                local_move_occurred = true;
                            } else if (!gain_0_ids.empty()) {
                                if (dis(gen) < 0.5) {
                                    // if no positive gain, then random neutral swaps
                                    best_u_id = gain_0_ids.back();
                                    bv_manager.move(g, p_manager, u, u_id, best_u_id);
                                    p_manager.move(u, u_weight, u_id, best_u_id);
                                    used[u] = mark;
                                    local_move_occurred = true;
                                }
                            }
                        }
                    }
                    global_move_occurred |= local_move_occurred;
                }
            }
        }
    };
}

#endif //HEIDELBERGPROCESSMAPPING_LABEL_PROPAGATION_REFINEMENT_H
