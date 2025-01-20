#ifndef HEIDELBERGPROCESSMAPPING_QUOTIENT_GRAPH_REFINEMENT_H
#define HEIDELBERGPROCESSMAPPING_QUOTIENT_GRAPH_REFINEMENT_H

#include <queue>

#include "../datastructures/distance_oracle.h"
#include "../datastructures/indexed_max_heap.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/qap.h"
#include "../utility/utils.h"

namespace HeiProMap {

    class QuotientGraphRefinement final : public ISerialRefiner {
    private:
        std::vector<partition_t> hierarchy;
        std::vector<weight_t> distance;
        partition_t k = 0;
        weight_t lmax = 0;

        // indexed max heaps
        IndexedMaxHeap<s64> imh_u;
        IndexedMaxHeap<s64> imh_v;
        std::vector<vertex_t> moves;
        std::vector<u8> is_good;

        std::vector<s32> used;
        s32 mark = -1;

    public:
        QuotientGraphRefinement() = default;

        void initialize(const vertex_t n,
                        std::vector<partition_t>& t_hierarchy,
                        std::vector<weight_t>& t_distance,
                        const weight_t t_lmax) override {
            hierarchy = t_hierarchy;
            distance  = t_distance;
            k         = prod<partition_t>(hierarchy);
            lmax      = t_lmax;

            used.resize(n, -1);

            // indexed max heaps
            imh_u = IndexedMaxHeap<s64>(n);
            imh_v = IndexedMaxHeap<s64>(n);
        }

        template <typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle, typename TSerialQuotientGraph>
        void refine(TSerialGraph& g,
                    TSerialActiveVertexManager& av_manager,
                    TSerialBoundaryVertexManager& bv_manager,
                    TSerialPartitionManager& p_manager,
                    TSerialDistanceOracle& d_oracle,
                    TSerialQuotientGraph& q_graph) {
            ASSERT(p_g != nullptr);
            Graph &g = *p_g;
            QuotientGraph &qg = pm.get_qg();

            bool global_move_occurred = true;
            u64 global_max_iteration = 1;
            for (u64 global_iteration = 0; global_iteration < global_max_iteration && global_move_occurred; ++global_iteration) {
                global_move_occurred = false;

                for (partition_t u_id = 0; u_id < k; ++u_id) {
                    for (partition_t v_id = u_id + 1; v_id < k; ++v_id) {
                        if (!qg.has_edge(u_id, v_id)) {
                            // no boundary between u_id and v_id
                            continue;
                        }

                        bool local_move_occurred = true;
                        u64 local_max_iteration = 1;
                        for (u64 local_iteration = 0; local_iteration < local_max_iteration && local_move_occurred; ++local_iteration) {
                            local_move_occurred = false;
                            mark += 1;

                            moves.clear();
                            is_good.clear();

                            // add all u boundary vertices
                            imh_u.clear();
                            for (vertex_t u : bv_manager[u_id]) {
                                for (size_t i = 0; i < g.size(u); ++i) {
                                    vertex_t v = g.neighbor(u, i);
                                    if (p_manager[v] == v_id) {
                                        // u is connected to block v_id
                                        s64 qap_delta = get_u_qap(g, u, p_manager, d_oracle) - get_u_qap(g, u, v_id, p_manager, d_oracle);
                                        imh_u.push(u, qap_delta);
                                        break;
                                    }
                                }
                            }

                            // add all v boundary vertices
                            imh_v.clear();
                            for (vertex_t v : bv_manager[v_id]) {
                                for (size_t i = 0; i < g.size(v); ++i) {
                                    vertex_t u = g.neighbor(v, i);
                                    if (p_manager[u] == u_id) {
                                        // v is connected to block u_id
                                        s64 qap_delta = get_u_qap(g, v, p_manager, d_oracle) - get_u_qap(g, v, u_id, p_manager, d_oracle);
                                        imh_v.push(v, qap_delta);
                                        break;
                                    }
                                }
                            }

                            // start moving vertices
                            while (!(imh_u.empty() && imh_v.empty())) {

                                // get stats for u
                                bool u_overloads_v = true;
                                bool u_is_overloaded = p_manager.get_pweight(u_id) > lmax;
                                bool u_qap_improve = false;
                                vertex_t u;
                                s64 u_qap_delta;
                                if (!imh_u.empty()) {
                                    u = imh_u.top_key();
                                    u_qap_delta = imh_u.top();
                                    u_qap_improve = u_qap_delta > 0;
                                    u_overloads_v = p_manager.get_pweight(v_id) + g.get_vertex_weight(u) > lmax;
                                }

                                // get stats for v
                                bool v_overloads_u = true;
                                bool v_is_overloaded = p_manager.get_pweight(v_id) > lmax;
                                bool v_qap_improve = false;
                                vertex_t v;
                                s64 v_qap_delta;
                                if (!imh_v.empty()) {
                                    v = imh_v.top_key();
                                    v_qap_delta = imh_v.top();
                                    v_qap_improve = v_qap_delta > 0;
                                    v_overloads_u = p_manager.get_pweight(u_id) + g.get_vertex_weight(v) > lmax;
                                }

                                bool move = false;
                                vertex_t v_to_move;
                                partition_t id_to_move;
                                u8 is_good_end;
                                // determine which vertex to move
                                if (u_overloads_v && v_overloads_u) {
                                    // move no vertex, since we don't want to mess up balancing even more
                                    // remove both vertices, since they are bad
                                } else if (u_overloads_v && !v_overloads_u) {
                                    // we can move v to u without overloading
                                    if (v_qap_improve) {
                                        // and additionally have a qap improvement
                                        move = true;
                                        is_good_end = 1;
                                        v_to_move = v;
                                        id_to_move = u_id;
                                    } else if (u_is_overloaded) {
                                        // qap gets worse, but u is currently overloaded, so execute the move
                                        move = true;
                                        is_good_end = 0; // but this is not a good end
                                        v_to_move = v;
                                        id_to_move = u_id;
                                    }
                                    // no qap improvement and u is not overloaded, no reason to move
                                } else if (!u_overloads_v && v_overloads_u) {
                                    // we can move u to v without overloading
                                    if (u_qap_improve) {
                                        // and additionally have a qap improvement
                                        move = true;
                                        is_good_end = 1;
                                        v_to_move = u;
                                        id_to_move = v_id;
                                    } else if (v_is_overloaded) {
                                        // qap gets worse, but v is currently overloaded, so execute the move
                                        move = true;
                                        is_good_end = 0; // but this is not a good move
                                        v_to_move = u;
                                        id_to_move = v_id;
                                    }
                                    // no qap improvement and v is not overloaded, no reason to move
                                } else {
                                    // we could move either way
                                    if (u_qap_improve && v_qap_improve) {
                                        if (u_qap_delta >= v_qap_delta) {
                                            move = true;
                                            is_good_end = 1;
                                            v_to_move = u;
                                            id_to_move = v_id;
                                        } else {
                                            move = true;
                                            is_good_end = 1;
                                            v_to_move = v;
                                            id_to_move = u_id;
                                        }
                                    } else if (u_qap_improve && !v_qap_improve) {
                                        move = true;
                                        is_good_end = 1;
                                        v_to_move = u;
                                        id_to_move = v_id;
                                    } else if (!u_qap_improve && v_qap_improve) {
                                        move = true;
                                        is_good_end = 1;
                                        v_to_move = v;
                                        id_to_move = u_id;
                                    } else {
                                        if (u_is_overloaded) {
                                            // no improvement, but u is overloaded so move
                                            move = true;
                                            is_good_end = 0;
                                            v_to_move = u;
                                            id_to_move = v_id;
                                        } else if (v_is_overloaded) {
                                            // no improvement, but v is overloaded so move
                                            move = true;
                                            is_good_end = 0;
                                            v_to_move = v;
                                            id_to_move = u_id;
                                        }
                                        move = false; // move non
                                    }
                                }

                                if (move) {
                                    moves.push_back(v_to_move);
                                    is_good.push_back(is_good_end);
                                    p_manager.move(v_to_move, id_to_move);
                                    used[v_to_move] = mark;

                                    if (id_to_move == u_id) {
                                        // moving v to v_id
                                        imh_v.pop();
                                    } else {
                                        // moving u to v_id
                                        imh_u.pop();
                                    }

                                    // we have to push or update the neighbors that were not moved already
                                    for (size_t i = 0; i < g.size(v_to_move); i++) {
                                        vertex_t ev = g.neighbor(v_to_move, i);
                                        if (p_manager[ev] == u_id && used[ev] != mark) {
                                            s64 qap_delta = get_u_qap(g, u, p_manager, d_oracle) - get_u_qap(g, u, v_id, p_manager, d_oracle);
                                            imh_u.push_update(ev, qap_delta);
                                        }
                                        if (p_manager[ev] == v_id && used[ev] != mark) {
                                            s64 qap_delta = get_u_qap(g, v, p_manager, d_oracle) - get_u_qap(g, v, u_id, p_manager, d_oracle);
                                            imh_v.push_update(ev, qap_delta);
                                        }
                                    }
                                } else {
                                    // remove both vertices, since they are both bad
                                    if (!imh_u.empty()) { imh_u.pop(); }
                                    if (!imh_v.empty()) { imh_v.pop(); }
                                }
                            }
                        }

                        std::cout << moves.size() << " ";

                        // revert to last good end
                        for (size_t i = 0; i < moves.size(); ++i) {
                            if (is_good.back() == 0) {
                                vertex_t u = moves.back();
                                if (p_manager[u] == u_id) {
                                    p_manager.move(u, v_id);
                                } else {
                                    p_manager.move(u, u_id);
                                }

                                moves.pop_back();
                                is_good.pop_back();
                            } else {
                                break;
                            }
                        }

                        std::cout << moves.size() << std::endl;
                    }
                }
            }
        }
    };


}

#endif //HEIDELBERGPROCESSMAPPING_QUOTIENT_GRAPH_REFINEMENT_H
