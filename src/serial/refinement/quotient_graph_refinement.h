#ifndef SERIALPROCESSMAPPING_QUOTIENT_GRAPH_REFINEMENT_H
#define SERIALPROCESSMAPPING_QUOTIENT_GRAPH_REFINEMENT_H

/*

#include <queue>

#include "../utility/definitions.h"
#include "../utility/macros.h"
#include "../utility/utils.h"
#include "../datastructures/graph.h"
#include "../utility/qap.h"
#include "../datastructures/distance_oracle.h"
#include "../datastructures/iterators/boundary_vertex_iterator.h"
#include "../datastructures/iterators/block_boundary_vertex_iterator.h"
#include "../datastructures/indexed_max_heap.h"

namespace HeiProMap {

    class QuotientGraphRefinement {

    private:
        Graph *p_g = nullptr;
        std::vector<u64> hierarchy;
        std::vector<u64> distance;
        u64 k = 0;
        f64 imbalance = 0;
        u64 lmax = 0;
        DistanceOracle *dist_o = nullptr;

        // indexed max heaps
        IndexedMaxHeap<s64> imh_u;
        IndexedMaxHeap<s64> imh_v;
        std::vector<vertex_t> moves;
        std::vector<u8> is_good;

        std::vector<s32> used;
        s32 mark = -1;

    public:
        QuotientGraphRefinement() = default;

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

            // indexed max heaps
            imh_u = IndexedMaxHeap<s64>(t_g->get_n());
            imh_v = IndexedMaxHeap<s64>(t_g->get_n());

            used.resize(t_g->get_n(), -1);
        }

        void refine(PartitionManager &pm) {
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

                        // std::cout << u_id << " " << v_id << std::endl;

                        bool local_move_occurred = true;
                        u64 local_max_iteration = 1;
                        for (u64 local_iteration = 0; local_iteration < local_max_iteration && local_move_occurred; ++local_iteration) {
                            local_move_occurred = false;
                            mark += 1;

                            moves.clear();
                            is_good.clear();

                            // add all u boundary vertices
                            imh_u.clear();
                            for (BlockBoundaryVertexIterator bbvi(pm.get_bvm(), u_id); bbvi.not_end(); bbvi.next()) {
                                vertex_t u = bbvi.get();
                                for (EdgeW &e: g[u]) {
                                    if (pm[e.v] == v_id) {
                                        // u is connected to block v_id
                                        s64 qap_delta = get_u_qap(g, u, pm, *dist_o) - get_u_qap(g, u, v_id, pm, *dist_o);
                                        imh_u.push(u, qap_delta);
                                        break;
                                    }
                                }
                            }

                            // add all v boundary vertices
                            imh_v.clear();
                            for (BlockBoundaryVertexIterator bbvi(pm.get_bvm(), v_id); bbvi.not_end(); bbvi.next()) {
                                vertex_t v = bbvi.get();
                                for (EdgeW &e: g[v]) {
                                    if (pm[e.v] == u_id) {
                                        // v is connected to block u_id
                                        s64 qap_delta = get_u_qap(g, v, pm, *dist_o) - get_u_qap(g, v, u_id, pm, *dist_o);
                                        imh_v.push(v, qap_delta);
                                        break;
                                    }
                                }
                            }

                            std::cout << "a " << imh_u.size() << " " << imh_v.size() << std::endl;

                            // start moving vertices
                            while (!(imh_u.empty() && imh_v.empty())) {
                                // std::cout << imh_u.size() << " " << imh_v.size() << std::endl;

                                // get stats for u
                                bool u_overloads_v = true;
                                bool u_is_overloaded = pm.get_pweight(u_id) > lmax;
                                bool u_qap_improve = false;
                                vertex_t u;
                                s64 u_qap_delta;
                                if (!imh_u.empty()) {
                                    u = imh_u.top_key();
                                    u_qap_delta = imh_u.top();
                                    u_qap_improve = u_qap_delta > 0;
                                    u_overloads_v = pm.get_pweight(v_id) + g.get_vertex_weight(u) > lmax;
                                }

                                // get stats for v
                                bool v_overloads_u = true;
                                bool v_is_overloaded = pm.get_pweight(v_id) > lmax;
                                bool v_qap_improve = false;
                                vertex_t v;
                                s64 v_qap_delta;
                                if (!imh_v.empty()) {
                                    v = imh_v.top_key();
                                    v_qap_delta = imh_v.top();
                                    v_qap_improve = v_qap_delta > 0;
                                    v_overloads_u = pm.get_pweight(u_id) + g.get_vertex_weight(v) > lmax;
                                }

                                // std::cout << u_overloads_v << " " << u_is_overloaded << " " << u_qap_improve << " " << u << " " << u_qap_delta << std::endl;
                                // std::cout << v_overloads_u << " " << v_is_overloaded << " " << v_qap_improve << " " << v << " " << v_qap_delta << std::endl;

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
                                    pm.move(v_to_move, id_to_move);
                                    used[v_to_move] = mark;

                                    if (id_to_move == u_id) {
                                        // moving v to v_id
                                        imh_v.pop();
                                    } else {
                                        // moving u to v_id
                                        imh_u.pop();
                                    }

                                    // we have to push or update the neighbors that were not moved already
                                    for (EdgeW &e: g[v_to_move]) {
                                        if (pm[e.v] == u_id && used[e.v] != mark) {
                                            s64 qap_delta = get_u_qap(g, u, pm, *dist_o) - get_u_qap(g, u, v_id, pm, *dist_o);
                                            imh_u.push_update(e.v, qap_delta);
                                        }
                                        if (pm[e.v] == v_id && used[e.v] != mark) {
                                            s64 qap_delta = get_u_qap(g, v, pm, *dist_o) - get_u_qap(g, v, u_id, pm, *dist_o);
                                            imh_v.push_update(e.v, qap_delta);
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
                                if (pm[u] == u_id) {
                                    pm.move(u, v_id);
                                } else {
                                    pm.move(u, u_id);
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

 */

#endif //SERIALPROCESSMAPPING_QUOTIENT_GRAPH_REFINEMENT_H
