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

#ifndef HEIPROMAP_QAP_H
#define HEIPROMAP_QAP_H

#include "../definitions.h"
#include "macros.h"
#include "../datastructures/distance_oracle.h"

namespace HeiProMap {
    template<typename GraphT, typename PartitionManagerT, typename DistanceOracleT>
    inline weight_t get_qap(const GraphT &g,
                            const PartitionManagerT &p_manager,
                            const DistanceOracleT &d_oracle) {
        ScopedTimer _t("misc", "misc", "get_qap");

        weight_t qap = 0;

        forall_gu(g, u)
            {
                partition_t u_id = p_manager[u];

                forall_guivw(g, u, i, v, w)
                    {
                        partition_t v_id = p_manager[v];
                        weight_t d = d_oracle.get(u_id, v_id);
                        qap += (d * w);
                    }
                endfor
            }
        endfor

        return qap;
    }

    template<typename GraphT, typename DistanceOracleT>
    inline weight_t get_qap(const GraphT &g,
                            const AlignedArray<partition_t> &partition,
                            const DistanceOracleT &d_oracle) {
        ScopedTimer _t("misc", "misc", "get_qap");

        weight_t qap = 0;

        forall_gu(g, u)
            {
                partition_t u_id = partition[u];

                forall_guivw(g, u, i, v, w)
                    {
                        partition_t v_id = partition[v];
                        weight_t d = d_oracle.get(u_id, v_id);
                        qap += (d * w);
                    }
                endfor
            }
        endfor

        return qap;
    }

    template<typename GraphT, typename PartitionManagerT>
    inline weight_t get_edge_cut(GraphT &g,
                                 PartitionManagerT &p_manager) {
        ScopedTimer _t("misc", "misc", "get_qap");

        weight_t edge_cut = 0;

        forall_gu(g, u)
            {
                forall_guivw(g, u, i, v, w)
                    {
                        if (p_manager[u] == p_manager[v]) { continue; }
                        edge_cut += w;
                    }
                endfor
            }
        endfor

        return edge_cut;
    }

    template<typename GraphT, typename PartitionManagerT, typename DistanceOracleT>
    inline weight_t get_qap(GraphT &g,
                            std::vector<vertex_t> &vertices,
                            PartitionManagerT &p_manager,
                            DistanceOracleT &d_oracle) {
        weight_t qap = 0;

        for (vertex_t u: vertices) {
            partition_t u_id = p_manager[u];

            forall_guivw(g, u, i, v, w)
                {
                    partition_t v_id = p_manager[v];
                    weight_t d = d_oracle.get(u_id, v_id);
                    qap += (d * w);
                }
            endfor
        }

        return qap;
    }

    template<typename GraphT, typename PartitionManagerT, typename DistanceOracleT, typename BlockConnT>
    inline weight_t get_qap(GraphT &g,
                            PartitionManagerT &p_manager,
                            DistanceOracleT &d_oracle,
                            BlockConnT &block_conn) {
        ScopedTimer _t("misc", "misc", "get_qap");

        weight_t qap = 0;
        weight_t local_qap = 0;

        forall_gu(g, u)
            {
                partition_t u_id = p_manager[u];

                forall_bc_ui_id_idw(block_conn, u, i, id, id_w)
                    {
                        weight_t d = d_oracle.get(u_id, id);
                        local_qap += (d * id_w);
                    }
                endfor
            }
        endfor
        qap += local_qap;

        return qap;
    }

    template<typename GraphT, typename PartitionManagerT, typename DistanceOracleT>
    inline std::vector<weight_t> get_qap_per_layer(const GraphT &g,
                                                   const PartitionManagerT &p_manager,
                                                   DistanceOracleT &d_oracle,
                                                   const partition_t l) {
        ScopedTimer _t("misc", "misc", "get_qap_per_layer");

        std::vector<weight_t> final_qap(l, 0);
        forall_gu(g, u)
            {
                partition_t u_id = p_manager[u];

                forall_guivw(g, u, i, v, w)
                    {
                        partition_t v_id = p_manager[v];
                        weight_t d = d_oracle.get(u_id, v_id);
                        partition_t layer_id = d_oracle.get_h(u_id, v_id);
                        final_qap[layer_id] += (d * w);
                    }
                endfor
            }
        endfor

        return final_qap;
    }


    template<typename GraphT, typename PartitionManagerT, typename DistanceOracleT>
    inline weight_t get_u_qap_delta(const GraphT &g,
                               const vertex_t u,
                               const partition_t old_id,
                               const partition_t new_id,
                               const PartitionManagerT &p_manager,
                               DistanceOracleT &d_oracle) {
        weight_t qap_delta = 0;

        forall_guivw(g, u, i, v, w)
            {
                partition_t v_id = p_manager[v];
                weight_t old_d = d_oracle.get(v_id, old_id);
                weight_t new_d = d_oracle.get(v_id, new_id);

                qap_delta += (old_d - new_d) * w;
            }
        endfor

        return qap_delta;
    }

    template<typename GraphT, typename PartitionManagerT, typename DistanceOracleT, typename BlockConnT>
    inline weight_t get_u_qap_delta(const GraphT &g,
                               const vertex_t u,
                               const partition_t old_id,
                               const partition_t new_id,
                               const PartitionManagerT &p_manager,
                               DistanceOracleT &d_oracle,
                               BlockConnT &block_conn) {
        weight_t qap_delta = 0;

        forall_bc_ui_id_idw(block_conn, u, i, id, idw)
            {
                weight_t old_d = d_oracle.get(id, old_id);
                weight_t new_d = d_oracle.get(id, new_id);

                qap_delta += (old_d - new_d) * idw;
            }
        endfor

        return qap_delta;
    }

    template<typename GraphT, typename PartitionManagerT, typename DistanceOracleT, typename BlockConnT>
    inline std::pair<partition_t, weight_t> get_u_qap_delta(const GraphT &g,
                                                            const vertex_t u,
                                                            const partition_t old_id,
                                                            const std::vector<partition_t> &new_ids,
                                                            const PartitionManagerT &p_manager,
                                                            DistanceOracleT &d_oracle,
                                                            BlockConnT &block_conn,
                                                            std::vector<weight_t> &qap_deltas) {
        forall_bc_ui_id_idw(block_conn, u, i, id, idw)
            {
                weight_t old_d = d_oracle.get(id, old_id);

                for (size_t j = 0; j < new_ids.size(); j++) {
                    qap_deltas[j] += (old_d - d_oracle.get(id, new_ids[j])) * idw;
                }
            }
        endfor

        std::pair<partition_t, weight_t> best = {new_ids[0], qap_deltas[0]};
        for (size_t j = 1; j < new_ids.size(); j++) {
            if (qap_deltas[j] > best.second) {
                best.first = new_ids[j];
                best.second = qap_deltas[j];
            }
        }
        return best;
    }

    template<typename GraphT, typename PartitionManagerT>
    inline weight_t get_u_edge_cut_delta(const GraphT &g,
                                    const vertex_t u,
                                    const partition_t old_id,
                                    const partition_t new_id,
                                    const PartitionManagerT &p_manager) {
        weight_t edge_cut_delta = 0;

        forall_guivw(g, u, i, v, w)
            {
                partition_t v_id = p_manager[v];

                edge_cut_delta -= (v_id != new_id) * w;
                edge_cut_delta += (v_id != old_id) * w;
            }
        endfor

        return edge_cut_delta;
    }

    template<typename GraphT, typename PartitionManagerT, typename BlockConnT>
    inline weight_t get_u_edge_cut_delta(const GraphT &g,
                                    const vertex_t u,
                                    const partition_t old_id,
                                    const partition_t new_id,
                                    const PartitionManagerT &p_manager,
                                    const BlockConnT &block_conn) {
        weight_t edge_cut_delta = 0;

        forall_bc_ui_id_idw(block_conn, u, i, id, idw)
            {
                edge_cut_delta -= (id != new_id) * idw;
                edge_cut_delta += (id != old_id) * idw;
            }
        endfor

        return edge_cut_delta;
    }

    /*
    inline void get_u_qap_delta(const graph_t &g,
                                const vertex_t u,
                                const partition_t old_id,
                                const partition_t *blocks,
                                weight_t *blocks_qap_delta,
                                const size_t blocks_size,
                                const p_manager_t &p_manager,
                                d_oracle_t &d_oracle) {
        blocks = ASSUME_ALIGNED(partition_t*, blocks, 64);
        blocks_qap_delta = ASSUME_ALIGNED(weight_t *, blocks_qap_delta, 64);

        // reset all to 0
        std::fill_n(blocks_qap_delta, blocks_size, 0);

#pragma GCC unroll 8
        forall_guivw(g, u, j, v, w)
            {
                partition_t v_id = p_manager[v];
                weight_t old_d = d_oracle.get(v_id, old_id);

                for (size_t i = 0; i < blocks_size; ++i) {
                    weight_t new_d = d_oracle.get(v_id, blocks[i]);
                    blocks_qap_delta[i] += (old_d - new_d) * w;
                }
            }
        endfor
    }
    */

    /*
    inline weight_t get_qap_delta(const graph_t &g,
                             const vertex_t u,
                             const partition_t u_old_id,
                             const partition_t u_new_id,
                             const vertex_t v,
                             const partition_t v_old_id,
                             const partition_t v_new_id,
                             const p_manager_t &p_manager,
                             d_oracle_t &d_oracle) {
        weight_t qap_delta = 0;

        // process u
#pragma GCC unroll 8
        forall_guivw(g, u, i, neighbor, w)
            {
                if (neighbor != v) {
                    partition_t neighbor_id = p_manager[neighbor];

                    weight_t old_d = d_oracle.get(neighbor_id, u_old_id);
                    weight_t new_d = d_oracle.get(neighbor_id, u_new_id);
                    qap_delta += (old_d - new_d) * w;
                } else {
                    weight_t old_d = d_oracle.get(v_old_id, u_old_id);
                    weight_t new_d = d_oracle.get(v_new_id, u_new_id);
                    qap_delta += (old_d - new_d) * w;
                }
            }
        endfor

        // process v
#pragma GCC unroll 8
        forall_guivw(g, v, i, neighbor, w)
            {
                if (neighbor == u) { continue; }

                partition_t neighbor_id = p_manager[neighbor];

                weight_t old_d = d_oracle.get(neighbor_id, v_old_id);
                weight_t new_d = d_oracle.get(neighbor_id, v_new_id);
                qap_delta += (old_d - new_d) * w;
            }
        endfor

        return qap_delta;
    }
    */

    /*
    inline weight_t get_qap_delta(const graph_t &g,
                             const vertex_t v,
                             const vertex_t vv,
                             const vertex_t vvv,
                             const partition_t v_id,
                             const partition_t vv_id,
                             const partition_t vvv_id,
                             const partition_t new_v_id,
                             const partition_t new_vv_id,
                             const partition_t new_vvv_id,
                             const p_manager_t &p_manager,
                             d_oracle_t &d_oracle) {
        weight_t qap_delta = 0;

        // process v
#pragma GCC unroll 8
        forall_guivw(g, v, i, neighbor, w)
            {
                if (neighbor == vv) {
                    weight_t old_d = d_oracle.get(v_id, vv_id);
                    weight_t new_d = d_oracle.get(new_v_id, new_vv_id);
                    qap_delta += (old_d - new_d) * w;
                } else if (neighbor == vvv) {
                    weight_t old_d = d_oracle.get(v_id, vvv_id);
                    weight_t new_d = d_oracle.get(new_v_id, new_vvv_id);
                    qap_delta += (old_d - new_d) * w;
                } else {
                    partition_t neighbor_id = p_manager[neighbor];

                    weight_t old_d = d_oracle.get(neighbor_id, v_id);
                    weight_t new_d = d_oracle.get(neighbor_id, new_v_id);
                    qap_delta += (old_d - new_d) * w;
                }
            }
        endfor

        // process vv
#pragma GCC unroll 8
        forall_guivw(g, vv, i, neighbor, w)
            {
                if (neighbor == v) { continue; }
                if (neighbor == vvv) {
                    weight_t old_d = d_oracle.get(vv_id, vvv_id);
                    weight_t new_d = d_oracle.get(new_vv_id, new_vvv_id);
                    qap_delta += (old_d - new_d) * w;
                } else {
                    partition_t neighbor_id = p_manager[neighbor];

                    weight_t old_d = d_oracle.get(neighbor_id, vv_id);
                    weight_t new_d = d_oracle.get(neighbor_id, new_vv_id);
                    qap_delta += (old_d - new_d) * w;
                }
            }
        endfor

        // process vvv
#pragma GCC unroll 8
        forall_guivw(g, vvv, i, neighbor, w)
            {
                if (neighbor == v) { continue; }
                if (neighbor == vv) { continue; }

                partition_t neighbor_id = p_manager[neighbor];

                weight_t old_d = d_oracle.get(neighbor_id, vvv_id);
                weight_t new_d = d_oracle.get(neighbor_id, new_vvv_id);
                qap_delta += (old_d - new_d) * w;
            }
        endfor

        return qap_delta;
    }
    */


    inline weight_t get_u_qap_delta_and_is_boundary(const graph_t &g,
                                               const vertex_t u,
                                               const partition_t old_id,
                                               const partition_t new_id,
                                               bool &is_boundary_old_id,
                                               bool &is_boundary_new_id,
                                               const p_manager_t &p_manager,
                                               d_oracle_t &d_oracle) {
        is_boundary_old_id = false;
        is_boundary_new_id = false;

        weight_t qap_delta = 0;

#pragma GCC unroll 8
        forall_guivw(g, u, i, v, w)
            {
                partition_t v_id = p_manager[v];

                is_boundary_old_id |= (v_id != old_id);
                is_boundary_new_id |= (v_id != new_id);

                weight_t old_d = d_oracle.get(v_id, old_id);
                weight_t new_d = d_oracle.get(v_id, new_id);
                qap_delta += (old_d - new_d) * w;
            }
        endfor

        return qap_delta;
    }

    template<typename GraphT, typename PartitionManagerT, typename DistanceOracleT, typename BlockConnT>
    inline weight_t get_u_qap_delta_and_is_boundary(const GraphT &g,
                                               const vertex_t u,
                                               const partition_t old_id,
                                               const partition_t new_id,
                                               bool &is_boundary_old_id,
                                               bool &is_boundary_new_id,
                                               const PartitionManagerT &p_manager,
                                               DistanceOracleT &d_oracle,
                                               BlockConnT &block_conn) {
        is_boundary_old_id = false;
        is_boundary_new_id = false;

        weight_t qap_delta = 0;

        forall_bc_ui_id_idw(block_conn, u, i, id, idw)
            {
                is_boundary_old_id |= (id != old_id);
                is_boundary_new_id |= (id != new_id);

                weight_t old_d = d_oracle.get(id, old_id);
                weight_t new_d = d_oracle.get(id, new_id);
                qap_delta += (old_d - new_d) * idw;
            }
        endfor


        return qap_delta;
    }

    template<typename GraphT, typename PartitionManagerT, typename DistanceOracleT>
    inline weight_t get_u_qap_delta_and_is_connected_to(const GraphT &g,
                                                   const vertex_t u,
                                                   const partition_t old_id,
                                                   const partition_t new_id,
                                                   bool &is_connected_to_new_id,
                                                   const PartitionManagerT &p_manager,
                                                   DistanceOracleT &d_oracle) {
        is_connected_to_new_id = false;

        weight_t qap_delta = 0;

        forall_guivw(g, u, i, v, w)
            {
                partition_t v_id = p_manager[v];

                is_connected_to_new_id |= (v_id == new_id);

                weight_t old_d = d_oracle.get(v_id, old_id);
                weight_t new_d = d_oracle.get(v_id, new_id);

                qap_delta += (old_d - new_d) * w;
            }
        endfor

        return qap_delta;
    }

    template<typename GraphT, typename PartitionManagerT, typename DistanceOracleT, typename BlockConnT>
    inline weight_t get_u_qap_delta_and_is_connected_to(const GraphT &g,
                                                   const vertex_t u,
                                                   const partition_t old_id,
                                                   const partition_t new_id,
                                                   bool &is_connected_to_new_id,
                                                   const PartitionManagerT &p_manager,
                                                   DistanceOracleT &d_oracle,
                                                   BlockConnT &block_conn) {
        is_connected_to_new_id = false;

        weight_t qap_delta = 0;

        forall_bc_ui_id_idw(block_conn, u, i, id, idw)
            {
                is_connected_to_new_id |= (id == new_id);

                weight_t old_d = d_oracle.get(id, old_id);
                weight_t new_d = d_oracle.get(id, new_id);

                qap_delta += (old_d - new_d) * idw;
            }
        endfor

        return qap_delta;
    }

    template<typename GraphT, typename PartitionManagerT>
    inline weight_t get_u_edge_cut_delta_and_is_connected_to(const GraphT &g,
                                                        const vertex_t u,
                                                        const partition_t old_id,
                                                        const partition_t new_id,
                                                        bool &is_connected_to_new_id,
                                                        const PartitionManagerT &p_manager) {
        is_connected_to_new_id = false;

        weight_t edge_cut_delta = 0;

        forall_guivw(g, u, i, v, w)
            {
                partition_t v_id = p_manager[v];

                is_connected_to_new_id |= (v_id == new_id);

                edge_cut_delta -= (v_id != new_id) * w;
                edge_cut_delta += (v_id != old_id) * w;
            }
        endfor

        return edge_cut_delta;
    }

    template<typename GraphT, typename PartitionManagerT, typename BlockConnT>
    inline weight_t get_u_edge_cut_delta_and_is_connected_to(const GraphT &g,
                                                        const vertex_t u,
                                                        const partition_t old_id,
                                                        const partition_t new_id,
                                                        bool &is_connected_to_new_id,
                                                        const PartitionManagerT &p_manager,
                                                        BlockConnT &block_conn) {
        is_connected_to_new_id = false;

        weight_t edge_cut_delta = 0;

        forall_bc_ui_id_idw(block_conn, u, i, id, idw)
            {
                is_connected_to_new_id |= (id == new_id);

                edge_cut_delta -= (id != new_id) * idw;
                edge_cut_delta += (id != old_id) * idw;
            }
        endfor

        return edge_cut_delta;
    }

    /*
    template<typename PartitionManagerT, typename DistanceOracleT>
    std::vector<weight_t> qap_per_layer(graph_t &g,
                                        PartitionManagerT &p_manager,
                                        DistanceOracleT &d_oracle,
                                        size_t l) {
        std::vector<weight_t> qap(l, 0);

        forall_gu(g, u)
            {
                partition_t u_id = p_manager[u];

                forall_guivw(g, u, i, v, w)
                    {
                        partition_t v_id = p_manager[v];
                        weight_t distance = d_oracle.get(u_id, v_id);
                        weight_t level = d_oracle.get_h(u_id, v_id);
                        qap[level] += w * distance;
                    }
                endfor
            }
        endfor

        return qap;
    }
    */
}

#endif //HEIPROMAP_QAP_H
