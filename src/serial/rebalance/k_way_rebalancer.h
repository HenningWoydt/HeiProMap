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

#ifndef HEIPROMAP_K_WAY_REBALANCER_H
#define HEIPROMAP_K_WAY_REBALANCER_H

#include <queue>

#include "../serial_definitions_1.h"
#include "../serial_definitions_3.h"
#include "../../commons/definitions.h"
#include "../../commons/random_engine.h"
#include "../../commons/statistic_collector.h"

namespace HeiProMap {
    class KWayRebalancer {
        vertex_t                 m_n    = 0;
        vertex_t                 m_m    = 0;
        partition_t              m_k    = 0;
        weight_t                 m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t>    m_distance;

        AlignedArray<u32> vertex_used;
        u32               vertex_mark = 0;

        AlignedArray<u32> block_used;
        u32               block_marker = 0;

        RandomEngine *random_engine = nullptr;

    public:
        KWayRebalancer() = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const weight_t t_lmax,
                        const std::vector<partition_t> &t_hierarchy,
                        const std::vector<weight_t> &t_distance,
                        RandomEngine &t_random_engine) {
            m_n         = t_n;
            m_m         = t_m;
            m_k         = t_k;
            m_lmax      = t_lmax;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;

            vertex_used.initialize(m_n, 0);
            vertex_mark = 0;

            block_used.initialize(m_k, 0);
            block_marker = 0;

            random_engine = &t_random_engine;
        }

        void rebalance(const u64 level,
                       const u64 max_level,
                       const graph_t &g,
                       d_oracle_t &d_oracle,
                       bv_manager_t &bv_manager,
                       p_manager_t &p_manager,
                       q_graph_t &q_graph) {
            while (p_manager.is_overloaded()) {
                std::vector<std::priority_queue<KWayFMMove>> heaps(m_k);

                // for each block collect the boundary vertices
                for (partition_t u_id = 0; u_id < m_k; ++u_id) {
                    forall_bv_id_iu(bv_manager, u_id, i, u)
                        {
                            weight_t u_weight = g.weight(u);

                            block_marker += 1;
                            forall_guiv(g, u, j, v)
                                {
                                    partition_t v_id = p_manager[v];
                                    if (v_id == u_id) { continue; }
                                    if (block_used[v_id] == block_marker) { continue; }
                                    if (p_manager.get_bweight(v_id) + u_weight > m_lmax) { continue; }

                                    s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                                    heaps[u_id].emplace(u, u_id, v_id, qap_delta);
                                    block_used[v_id] = block_marker;
                                }
                            endfor
                        }
                    endfor
                }

                vertex_mark += 1;
                while (p_manager.is_overloaded()) {
                    // determine the most overloaded block
                    partition_t o_id     = std::numeric_limits<partition_t>::max();
                    weight_t    o_weight = -std::numeric_limits<weight_t>::max();

                    for (partition_t id = 0; id < m_k; ++id) {
                        if (p_manager.get_bweight(id) > m_lmax && p_manager.get_bweight(id) > o_weight && !heaps[id].empty()) {
                            o_id     = id;
                            o_weight = p_manager.get_bweight(id);
                        }
                    }

                    if (o_weight == -std::numeric_limits<weight_t>::max()) {
                        break;
                    }

                    while (!heaps[o_id].empty()) {
                        const KWayFMMove move = heaps[o_id].top();
                        heaps[o_id].pop();

                        vertex_t    vertex        = move.u;
                        partition_t vertex_id     = p_manager[vertex];
                        weight_t    vertex_weight = g.weight(vertex);
                        partition_t move_id       = move.to_move_id;

                        if (vertex_used[vertex] == vertex_mark) { continue; }
                        if (vertex_id != move.u_id) { continue; }
                        if (!bv_manager.is_boundary(vertex)) { continue; }

                        bool is_connected_move_id;
                        s64  temp_qap_delta = get_u_qap_delta_and_is_connected_to(g, vertex, vertex_id, move_id, is_connected_move_id, p_manager, d_oracle);
                        if (!is_connected_move_id) { continue; }
                        if (temp_qap_delta != move.qap_delta) { continue; }

                        // make move in structures
                        vertex_used[vertex] = vertex_mark;
                        bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                        q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                        p_manager.move(vertex, vertex_weight, vertex_id, move_id);

                        // we have to push or update the neighbors
                        forall_guiv(g, vertex, i, neighbor)
                            {
                                if (vertex_used[neighbor] == vertex_mark) { continue; }
                                if (!bv_manager.is_boundary(neighbor)) { continue; }

                                partition_t neighbor_id     = p_manager[neighbor];
                                weight_t    neighbor_weight = g.weight(neighbor);

                                block_marker += 1;
                                forall_guiv(g, neighbor, j, v)
                                    {
                                        partition_t v_id = p_manager[v];
                                        if (v_id == neighbor_id) { continue; }
                                        if (block_used[v_id] == block_marker) { continue; }
                                        if (p_manager.get_bweight(v_id) + neighbor_weight > m_lmax) { continue; }

                                        s64 qap_delta = get_u_qap_delta(g, neighbor, neighbor_id, v_id, p_manager, d_oracle);
                                        heaps[neighbor_id].emplace(neighbor, neighbor_id, v_id, qap_delta);
                                        block_used[v_id] = block_marker;
                                    }
                                endfor
                            }
                        endfor
                    }
                }
            }
        }
    };
}

#endif //HEIPROMAP_K_WAY_REBALANCER_H
