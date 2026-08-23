/*******************************************************************************
 * MIT License
 *
 * This file is part of HeiProMap.
 *
 * Copyright (C) 2025 Henning Woydt <henning.woydt@informatik.uni-heidelberg.de>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#ifndef HEIPROMAP_REBALANCER_H
#define HEIPROMAP_REBALANCER_H

#include <limits>
#include <vector>
#include <queue>

#include "../definitions.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/distance_oracle.h"
#include "../datastructures/partition_manager.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/quotient_graph.h"
#include "../datastructures/block_conn.h"
#include "../utility/random_engine.h"
#include "../utility/qap.h"
#include "../utility/indexed_max_heap.h"

namespace HeiProMap {
    struct RebalanceMovePayload {
        partition_t id;
        weight_t qap_delta;

        bool operator>(const RebalanceMovePayload &other) const {
            return qap_delta > other.qap_delta;
        }
        bool operator<=(const RebalanceMovePayload &other) const {
            return qap_delta <= other.qap_delta;
        }
    };

    struct RebalancerMove {
        vertex_t u;
        partition_t best_id;
        weight_t best_qap;
        u64 state_id;

        RebalancerMove(const vertex_t t_u, const partition_t t_id, const weight_t t_qap, const u64 t_state_id) {
            u = t_u;
            best_id = t_id;
            best_qap = t_qap;
            state_id = t_state_id;
        }

        bool operator>(const RebalancerMove &m) const { return best_qap > m.best_qap; }
        bool operator<(const RebalancerMove &m) const { return best_qap < m.best_qap; }
        bool operator<=(const RebalancerMove &m) const { return best_qap <= m.best_qap; }
    };

    class Rebalancer {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
    public:
        void update_k(partition_t k) { m_k = k; }
    private:

        RandomEngine random_engine;

    public:
        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const u64 seed) {
            HEIPROMAP_PROFILE_SCOPE("rebalance", "Rebalancer", "initialize");

            m_n = t_n;
            m_m = t_m;
            m_k = t_k;

            random_engine = RandomEngine(seed);
        }

        RebalancerMove get_best_move(vertex_t u,
                                     const graph_t &g,
                                     const p_manager_t &p_manager,
                                     [[maybe_unused]] const q_graph_t &q_graph,
                                     const d_oracle_t &d_oracle,
                                     const u64 state_id,
                                     weight_t lmax) {
            RebalancerMove move(u, m_k, -std::numeric_limits<weight_t>::max(), state_id);

            partition_t u_id = p_manager[u];
            weight_t u_weight = g.v_weights[u];

            for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                const vertex_t v = g.edges_v[j];

                partition_t v_id = p_manager[v];
                if (p_manager.get_bweight(v_id) + u_weight > lmax) { continue; }

                weight_t qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                if (qap_delta > move.best_qap || (qap_delta == move.best_qap && p_manager.get_bweight(v_id) < p_manager.get_bweight(move.best_id))) {
                    move.best_qap = qap_delta;
                    move.best_id = v_id;
                }
            }

            if (move.best_id != m_k) { return move; }

            for (size_t i = 0; i < 50; ++i) {
                partition_t v_id = random_engine.get_u64() % m_k;
                if (v_id == u_id) { continue; }
                if (p_manager.get_bweight(v_id) + u_weight > lmax) { continue; }

                weight_t qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                if (qap_delta > move.best_qap || (qap_delta == move.best_qap && p_manager.get_bweight(v_id) < p_manager.get_bweight(move.best_id))) {
                    move.best_qap = qap_delta;
                    move.best_id = v_id;
                }
            }

            return move;
        }

        RebalancerMove get_local_best_move(vertex_t u,
                                           const graph_t &g,
                                           const p_manager_t &p_manager,
                                           const d_oracle_t &d_oracle,
                                           const u64 state_id,
                                           weight_t lmax) const {
            RebalancerMove move(u, m_k, -std::numeric_limits<weight_t>::max(), state_id);

            partition_t u_id = p_manager[u];
            weight_t u_weight = g.v_weights[u];

            for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                const vertex_t v = g.edges_v[j];

                partition_t v_id = p_manager[v];
                if (p_manager.get_bweight(v_id) + u_weight > lmax) { continue; }

                weight_t qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                if (qap_delta > move.best_qap || (qap_delta == move.best_qap && p_manager.get_bweight(v_id) < p_manager.get_bweight(move.best_id))) {
                    move.best_qap = qap_delta;
                    move.best_id = v_id;
                }
            }

            return move;
        }

        void rebalance(const graph_t &g,
                       p_manager_t &p_manager,
                       bv_manager_t &bv_manager,
                       q_graph_t &q_graph,
                       d_oracle_t &d_oracle,
                       block_conn_t &block_conn,
                       f64 imbalance) {
            fill_empty_blocks(g, p_manager, bv_manager, q_graph, d_oracle, block_conn, imbalance);

            HEIPROMAP_PROFILE_SCOPE("rebalance", "Rebalancer", "allocate");

            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            AlignedArray<vertex_t> boundary;
            boundary.initialize(g.n);
            size_t boundary_size = 0;

            AlignedArray<RebalancerMove> moves;
            std::priority_queue<RebalancerMove> global_queue;

            AlignedArray<u64> offsets;
            offsets.initialize(m_k + 1);
            size_t offsets_size = 0;
            AlignedArray<u64> cursor;
            cursor.initialize(m_k + 1);

            AlignedArray<u64> state_ids;
            state_ids.initialize(g.n, 0);

            bool move_made = true;
            u64 iter = 0;
            while (move_made && iter < 100) {
                iter++;
                HEIPROMAP_PROFILE_SCOPE("rebalance", "Rebalancer", "get_boundary");
                move_made = false;

                offsets_size = 1;
                offsets[0] = 0;
                for (partition_t id = 0; id < m_k; ++id) {
                    if (p_manager.get_bweight(id) > lmax) {
                        offsets[offsets_size] = offsets[offsets_size - 1] + bv_manager.size(id);
                    } else {
                        offsets[offsets_size] = offsets[offsets_size - 1];
                    }
                    offsets_size += 1;
                }

                for (size_t i = 0; i < offsets_size; ++i) {
                    cursor[i] = offsets[i];
                }
                boundary_size = offsets[offsets_size - 1];

                for (partition_t id = 0; id < m_k; ++id) {
                    if (p_manager.get_bweight(id) > lmax) {
                        for (size_t i = 0; i < bv_manager.size(id); ++i) {
                            const vertex_t u = bv_manager.get(id, i);

                            boundary[cursor[id]] = u;
                            cursor[id] += 1;
                        }
                    }
                }

                HEIPROMAP_PROFILE_SCOPE("rebalance", "Rebalancer", "fill_heaps");

                moves.initialize(boundary_size);
                for (u64 i = 0; i < boundary_size; ++i) {
                    vertex_t u = boundary[i];
                    state_ids[u] += 1;
                    RebalancerMove move = get_local_best_move(u, g, p_manager, d_oracle, state_ids[u], lmax);
                    moves[i] = move;
                }

                HEIPROMAP_PROFILE_SCOPE("rebalance", "Rebalancer", "global_heap");

                for (u64 i = 0; i < boundary_size; ++i) {
                    if (moves[i].best_id != m_k) {
                        global_queue.push(moves[i]);
                    }
                }

                HEIPROMAP_PROFILE_SCOPE("rebalance", "Rebalancer", "process_heap");

                u64 inner_iter = 0;
                while (!global_queue.empty() && inner_iter < (u64) g.n * 10) {
                    inner_iter++;
                    RebalancerMove move = global_queue.top();
                    global_queue.pop();
                    vertex_t u = move.u;
                    weight_t u_weight = g.v_weights[u];
                    partition_t u_id = p_manager[u];
                    partition_t best_id = move.best_id;

                    if (move.best_id == m_k) { continue; }                 // not a valid destination
                    if (bv_manager.is_boundary(u) == false) { continue; }  // not a boundary vertex
                    if (p_manager.get_bweight(u_id) <= lmax) { continue; } // dont need to move anymore
                    if (state_ids[u] != move.state_id) { continue; }

                    if (p_manager.get_bweight(best_id) + u_weight > lmax) {
                        // best_id is overloaded, recompute
                        state_ids[u] += 1;
                        RebalancerMove new_move = get_local_best_move(u, g, p_manager, d_oracle, state_ids[u], lmax);
                        if (new_move.best_id != m_k) {
                            global_queue.push(new_move);
                        }
                        continue;
                    }

                    // move
                    bv_manager.move(g, p_manager, u, u_id, best_id);
                    q_graph.move(g, p_manager, u, u_id, best_id);
                    block_conn.move(g, u, u_id, best_id);
                    p_manager.move(u, u_weight, u_id, best_id);
                    move_made = true;

                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const vertex_t v = g.edges_v[j];

                        if (bv_manager.is_boundary(v) == false) { continue; }
                        if (p_manager.get_bweight(p_manager[v]) <= lmax) { continue; }

                        state_ids[v] += 1;
                        RebalancerMove new_move = get_local_best_move(v, g, p_manager, d_oracle, state_ids[v], lmax);
                        if (new_move.best_id != m_k) {
                            global_queue.push(new_move);
                        }
                    }
                }
            }
        }

        void rebalance_last_layer(const graph_t &g,
                                  p_manager_t &p_manager,
                                  bv_manager_t &bv_manager,
                                  q_graph_t &q_graph,
                                  d_oracle_t &d_oracle,
                                  block_conn_t &block_conn,
                                  f64 imbalance) {
            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            rebalance(g, p_manager, bv_manager, q_graph, d_oracle, block_conn, imbalance);
            HEIPROMAP_PROFILE_SCOPE("rebalance", "LL-Rebalancer", "allocate");

            AlignedArray<vertex_t> boundary;
            boundary.initialize(g.n);
            size_t boundary_size = 0;

            AlignedArray<RebalancerMove> moves;
            std::priority_queue<RebalancerMove> global_queue;

            AlignedArray<u64> offsets;
            offsets.initialize(m_k + 1);
            size_t offsets_size = 0;

            AlignedArray<u64> cursor;
            cursor.initialize(m_k + 1);

            AlignedArray<u64> state_ids;
            state_ids.initialize(g.n, 0);

            bool move_made = true;
            u64 iter = 0;
            while (move_made && iter < 100) {
                iter++;
                HEIPROMAP_PROFILE_SCOPE("rebalance", "LL-Rebalancer", "get_boundary");
                move_made = false;

                offsets_size = 1;
                offsets[0] = 0;
                for (partition_t id = 0; id < m_k; ++id) {
                    if (p_manager.get_bweight(id) > lmax) {
                        offsets[offsets_size] = offsets[offsets_size - 1] + bv_manager.size(id);
                    } else {
                        offsets[offsets_size] = offsets[offsets_size - 1];
                    }
                    offsets_size += 1;
                }

                for (size_t i = 0; i < offsets_size; ++i) {
                    cursor[i] = offsets[i];
                }
                boundary_size = offsets[offsets_size - 1];


                for (partition_t id = 0; id < m_k; ++id) {
                    if (p_manager.get_bweight(id) > lmax) {
                        for (size_t i = 0; i < bv_manager.size(id); ++i) {
                            const vertex_t u = bv_manager.get(id, i);

                            boundary[cursor[id]] = u;
                            cursor[id] += 1;
                        }
                    }
                }

                HEIPROMAP_PROFILE_SCOPE("rebalance", "LL-Rebalancer", "fill_heaps");

                moves.initialize(boundary_size);
                for (u64 i = 0; i < boundary_size; ++i) {
                    vertex_t u = boundary[i];
                    state_ids[u] += 1;
                    RebalancerMove move = get_best_move(u, g, p_manager, q_graph, d_oracle, state_ids[u], lmax);
                    moves[i] = move;
                }

                HEIPROMAP_PROFILE_SCOPE("rebalance", "LL-Rebalancer", "global_heap");

                global_queue = std::priority_queue<RebalancerMove>();
                for (u64 i = 0; i < boundary_size; ++i) {
                    if (moves[i].best_id != m_k) {
                        global_queue.push(moves[i]);
                    }
                }

                HEIPROMAP_PROFILE_SCOPE("rebalance", "LL-Rebalancer", "process_heap");

                u64 inner_iter = 0;
                while (!global_queue.empty() && inner_iter < (u64) g.n * 10) {
                    inner_iter++;
                    RebalancerMove move = global_queue.top();
                    global_queue.pop();
                    vertex_t u = move.u;
                    weight_t u_weight = g.v_weights[u];
                    partition_t u_id = p_manager[u];
                    partition_t best_id = move.best_id;

                    if (move.best_id == m_k) { continue; }                 // not a valid destination
                    if (bv_manager.is_boundary(u) == false) { continue; }  // not a boundary vertex
                    if (p_manager.get_bweight(u_id) <= lmax) { continue; } // dont need to move anymore
                    if (state_ids[u] != move.state_id) { continue; }

                    if (p_manager.get_bweight(best_id) + u_weight > lmax) {
                        // best_id is overloaded, recompute
                        state_ids[u] += 1;
                        RebalancerMove new_move = get_best_move(u, g, p_manager, q_graph, d_oracle, state_ids[u], lmax);
                        if (new_move.best_id != m_k) {
                            global_queue.push(new_move);
                        }
                        continue;
                    }

                    // move
                    bv_manager.move(g, p_manager, u, u_id, best_id);
                    q_graph.move(g, p_manager, u, u_id, best_id);
                    block_conn.move(g, u, u_id, best_id);
                    p_manager.move(u, u_weight, u_id, best_id);
                    move_made = true;

                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const vertex_t v = g.edges_v[j];

                        if (bv_manager.is_boundary(v) == false) { continue; }
                        if (p_manager.get_bweight(p_manager[v]) <= lmax) { continue; }

                        state_ids[v] += 1;
                        RebalancerMove new_move = get_best_move(v, g, p_manager, q_graph, d_oracle, state_ids[v], lmax);
                        if (new_move.best_id != m_k) {
                            global_queue.push(new_move);
                        }
                    }
                }
            }
        }

        void fill_empty_blocks(const graph_t &g,
                               p_manager_t &p_manager,
                               bv_manager_t &bv_manager,
                               q_graph_t &q_graph,
                               d_oracle_t &d_oracle,
                               block_conn_t &block_conn,
                               f64 imbalance) {
            // return;
            HEIPROMAP_PROFILE_SCOPE("rebalance", "Rebalancer", "fill_empty_blocks");

            std::vector<bool> blocks_to_fill_lookup(m_k, false);
            std::vector<partition_t> blocks_to_fill;
            for (partition_t id = 0; id < m_k; ++id) {
                if (p_manager.get_bweight(id) == 0) {
                    blocks_to_fill_lookup[id] = true;
                    blocks_to_fill.push_back(id);
                }
            }

            if (blocks_to_fill.empty()) { return; }

            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));
            weight_t fill_value = lmax / 2;

            IndexedMaxHeap<RebalanceMovePayload> heap;
            heap.initialize(g.n);

            for (partition_t id = 0; id < m_k; ++id) {
                for (size_t i = 0; i < bv_manager.size(id); ++i) {
                    const vertex_t u = bv_manager.get(id, i);

                    weight_t u_w = g.v_weights[u];

                    partition_t best_id = NO_ID;
                    weight_t best_qap = -std::numeric_limits<weight_t>::max();

                    for (partition_t move_id: blocks_to_fill) {
                        if (p_manager.get_bweight(move_id) + u_w <= lmax) {
                            weight_t qap = get_u_qap_delta(g, u, id, move_id, p_manager, d_oracle, block_conn);

                            if (qap > best_qap) {
                                best_qap = qap;
                                best_id = move_id;
                            }
                        }
                    }

                    if (best_id != NO_ID) {
                        heap.push(u, {best_id, best_qap});
                    }
                }
            }

            u64 empty_iter = 0;
            while (!heap.empty() && empty_iter < (u64) g.n * 10) {
                empty_iter++;
                vertex_t u = heap.top_key();
                partition_t u_id = p_manager[u];
                weight_t u_w = g.v_weights[u];
                partition_t new_id = heap.top().id;
                heap.pop();

                if (p_manager.get_bweight(new_id) >= fill_value) { continue; }
                if (p_manager.get_bweight(new_id) + u_w > lmax) { continue; }

                bv_manager.move(g, p_manager, u, u_id, new_id);
                q_graph.move(g, p_manager, u, u_id, new_id);
                block_conn.move(g, u, u_id, new_id);
                p_manager.move(u, u_w, u_id, new_id);

                if (p_manager.get_bweight(new_id) >= fill_value) {
                    blocks_to_fill_lookup[new_id] = false;
                    // blocks_to_fill.erase(); // todo
                }

                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    const vertex_t v = g.edges_v[i];
                    if (!bv_manager.is_boundary(v)) { continue; }

                    partition_t v_id = p_manager[v];
                    if (blocks_to_fill_lookup[v_id] == true) { continue; }

                    weight_t v_w = g.v_weights[v];
                    partition_t best_id = NO_ID;
                    weight_t best_qap = -std::numeric_limits<weight_t>::max();

                    for (partition_t move_id: blocks_to_fill) {
                        if (blocks_to_fill_lookup[move_id] == false) { continue; }
                        if (p_manager.get_bweight(move_id) + v_w <= lmax) {
                            weight_t qap = get_u_qap_delta(g, v, v_id, move_id, p_manager, d_oracle, block_conn);

                            if (qap > best_qap) {
                                best_qap = qap;
                                best_id = move_id;
                            }
                        }
                    }

                    if (best_id != NO_ID) {
                        heap.push_update(v, {best_id, best_qap});
                    }
                }
            }
        }
    };
}

#endif //HEIPROMAP_REBALANCER_H
