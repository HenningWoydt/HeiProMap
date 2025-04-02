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

#include "../../commons/definitions.h"
#include "../../commons/macros.h"
#include "../datastructures/distance_oracle.h"

namespace HeiProMap {
    inline weight_t get_qap(const graph_t& g,
                            const p_manager_t& p_manager,
                            const d_oracle_t& d_oracle) {
        weight_t qap = 0;

        forall_gu(g, u)
            {
                partition_t u_id = p_manager[u];

#pragma GCC unroll 8
                for (size_t i = 0; i < g.size(u); ++i) {
                    vertex_t v       = g.neighbor(u, i);
                    weight_t ew      = g.weight(u, i);
                    partition_t v_id = p_manager[v];
                    weight_t d       = d_oracle.get(u_id, v_id);
                    qap += (d * ew);
                }
            }
        endfor

        return qap;
    }

    inline s64 get_u_qap_delta_size_1(const graph_t& g,
                                      const vertex_t u,
                                      const partition_t old_id,
                                      const partition_t new_id,
                                      const p_manager_t& p_manager,
                                      const d_oracle_t& d_oracle) {
        vertex_t v       = g.neighbor(u, 0);
        weight_t w       = g.weight(u, 0);
        partition_t v_id = p_manager[v];
        weight_t old_d   = d_oracle.get(v_id, old_id);
        weight_t new_d   = d_oracle.get(v_id, new_id);
        return (old_d - new_d) * w;
    }

    inline s64 get_u_qap_delta_size_2(const graph_t& g,
                                      const vertex_t u,
                                      const partition_t old_id,
                                      const partition_t new_id,
                                      const p_manager_t& p_manager,
                                      const d_oracle_t& d_oracle) {
        vertex_t v0 = g.neighbor(u, 0);
        vertex_t v1 = g.neighbor(u, 1);

        weight_t w0 = g.weight(u, 0);
        weight_t w1 = g.weight(u, 1);

        partition_t v0_id = p_manager[v0];
        partition_t v1_id = p_manager[v1];

        weight_t old_d0 = d_oracle.get(old_id, v0_id);
        weight_t old_d1 = d_oracle.get(old_id, v1_id);

        weight_t new_d0 = d_oracle.get(new_id, v0_id);
        weight_t new_d1 = d_oracle.get(new_id, v1_id);

        s64 qap_delta = 0;
        qap_delta += (old_d0 - new_d0) * w0;
        qap_delta += (old_d1 - new_d1) * w1;

        return qap_delta;
    }

    inline s64 get_u_qap_delta_size_3(const graph_t& g,
                                      const vertex_t u,
                                      const partition_t old_id,
                                      const partition_t new_id,
                                      const p_manager_t& p_manager,
                                      const d_oracle_t& d_oracle) {
        vertex_t v0 = g.neighbor(u, 0);
        vertex_t v1 = g.neighbor(u, 1);
        vertex_t v2 = g.neighbor(u, 2);

        weight_t w0 = g.weight(u, 0);
        weight_t w1 = g.weight(u, 1);
        weight_t w2 = g.weight(u, 2);

        partition_t v0_id = p_manager[v0];
        partition_t v1_id = p_manager[v1];
        partition_t v2_id = p_manager[v2];

        weight_t old_d0 = d_oracle.get(old_id, v0_id);
        weight_t old_d1 = d_oracle.get(old_id, v1_id);
        weight_t old_d2 = d_oracle.get(old_id, v2_id);

        weight_t new_d0 = d_oracle.get(new_id, v0_id);
        weight_t new_d1 = d_oracle.get(new_id, v1_id);
        weight_t new_d2 = d_oracle.get(new_id, v2_id);

        s64 qap_delta = 0;
        qap_delta += (old_d0 - new_d0) * w0;
        qap_delta += (old_d1 - new_d1) * w1;
        qap_delta += (old_d2 - new_d2) * w2;

        return qap_delta;
    }

    inline s64 get_u_qap_delta_size_4(const graph_t& g,
                                      const vertex_t u,
                                      const partition_t old_id,
                                      const partition_t new_id,
                                      const p_manager_t& p_manager,
                                      const d_oracle_t& d_oracle) {
        vertex_t v0 = g.neighbor(u, 0);
        vertex_t v1 = g.neighbor(u, 1);
        vertex_t v2 = g.neighbor(u, 2);
        vertex_t v3 = g.neighbor(u, 3);

        weight_t w0 = g.weight(u, 0);
        weight_t w1 = g.weight(u, 1);
        weight_t w2 = g.weight(u, 2);
        weight_t w3 = g.weight(u, 3);

        partition_t v0_id = p_manager[v0];
        partition_t v1_id = p_manager[v1];
        partition_t v2_id = p_manager[v2];
        partition_t v3_id = p_manager[v3];

        weight_t old_d0 = d_oracle.get(old_id, v0_id);
        weight_t old_d1 = d_oracle.get(old_id, v1_id);
        weight_t old_d2 = d_oracle.get(old_id, v2_id);
        weight_t old_d3 = d_oracle.get(old_id, v3_id);

        weight_t new_d0 = d_oracle.get(new_id, v0_id);
        weight_t new_d1 = d_oracle.get(new_id, v1_id);
        weight_t new_d2 = d_oracle.get(new_id, v2_id);
        weight_t new_d3 = d_oracle.get(new_id, v3_id);

        s64 qap_delta = 0;
        qap_delta += (old_d0 - new_d0) * w0;
        qap_delta += (old_d1 - new_d1) * w1;
        qap_delta += (old_d2 - new_d2) * w2;
        qap_delta += (old_d3 - new_d3) * w3;

        return qap_delta;
    }

    inline s64 get_u_qap_delta_size_5(const graph_t& g,
                                      const vertex_t u,
                                      const partition_t old_id,
                                      const partition_t new_id,
                                      const p_manager_t& p_manager,
                                      const d_oracle_t& d_oracle) {
        vertex_t v0 = g.neighbor(u, 0);
        vertex_t v1 = g.neighbor(u, 1);
        vertex_t v2 = g.neighbor(u, 2);
        vertex_t v3 = g.neighbor(u, 3);
        vertex_t v4 = g.neighbor(u, 4);

        weight_t w0 = g.weight(u, 0);
        weight_t w1 = g.weight(u, 1);
        weight_t w2 = g.weight(u, 2);
        weight_t w3 = g.weight(u, 3);
        weight_t w4 = g.weight(u, 4);

        partition_t v0_id = p_manager[v0];
        partition_t v1_id = p_manager[v1];
        partition_t v2_id = p_manager[v2];
        partition_t v3_id = p_manager[v3];
        partition_t v4_id = p_manager[v4];

        weight_t old_d0 = d_oracle.get(old_id, v0_id);
        weight_t old_d1 = d_oracle.get(old_id, v1_id);
        weight_t old_d2 = d_oracle.get(old_id, v2_id);
        weight_t old_d3 = d_oracle.get(old_id, v3_id);
        weight_t old_d4 = d_oracle.get(old_id, v4_id);

        weight_t new_d0 = d_oracle.get(new_id, v0_id);
        weight_t new_d1 = d_oracle.get(new_id, v1_id);
        weight_t new_d2 = d_oracle.get(new_id, v2_id);
        weight_t new_d3 = d_oracle.get(new_id, v3_id);
        weight_t new_d4 = d_oracle.get(new_id, v4_id);

        s64 qap_delta = 0;
        qap_delta += (old_d0 - new_d0) * w0;
        qap_delta += (old_d1 - new_d1) * w1;
        qap_delta += (old_d2 - new_d2) * w2;
        qap_delta += (old_d3 - new_d3) * w3;
        qap_delta += (old_d4 - new_d4) * w4;

        return qap_delta;
    }

    inline s64 get_u_qap_delta_size_6(const graph_t& g,
                                      const vertex_t u,
                                      const partition_t old_id,
                                      const partition_t new_id,
                                      const p_manager_t& p_manager,
                                      const d_oracle_t& d_oracle) {
        vertex_t v0 = g.neighbor(u, 0);
        vertex_t v1 = g.neighbor(u, 1);
        vertex_t v2 = g.neighbor(u, 2);
        vertex_t v3 = g.neighbor(u, 3);
        vertex_t v4 = g.neighbor(u, 4);
        vertex_t v5 = g.neighbor(u, 5);

        weight_t w0 = g.weight(u, 0);
        weight_t w1 = g.weight(u, 1);
        weight_t w2 = g.weight(u, 2);
        weight_t w3 = g.weight(u, 3);
        weight_t w4 = g.weight(u, 4);
        weight_t w5 = g.weight(u, 5);

        partition_t v0_id = p_manager[v0];
        partition_t v1_id = p_manager[v1];
        partition_t v2_id = p_manager[v2];
        partition_t v3_id = p_manager[v3];
        partition_t v4_id = p_manager[v4];
        partition_t v5_id = p_manager[v5];

        weight_t old_d0 = d_oracle.get(old_id, v0_id);
        weight_t old_d1 = d_oracle.get(old_id, v1_id);
        weight_t old_d2 = d_oracle.get(old_id, v2_id);
        weight_t old_d3 = d_oracle.get(old_id, v3_id);
        weight_t old_d4 = d_oracle.get(old_id, v4_id);
        weight_t old_d5 = d_oracle.get(old_id, v5_id);

        weight_t new_d0 = d_oracle.get(new_id, v0_id);
        weight_t new_d1 = d_oracle.get(new_id, v1_id);
        weight_t new_d2 = d_oracle.get(new_id, v2_id);
        weight_t new_d3 = d_oracle.get(new_id, v3_id);
        weight_t new_d4 = d_oracle.get(new_id, v4_id);
        weight_t new_d5 = d_oracle.get(new_id, v5_id);

        s64 qap_delta = 0;
        qap_delta += (old_d0 - new_d0) * w0;
        qap_delta += (old_d1 - new_d1) * w1;
        qap_delta += (old_d2 - new_d2) * w2;
        qap_delta += (old_d3 - new_d3) * w3;
        qap_delta += (old_d4 - new_d4) * w4;
        qap_delta += (old_d5 - new_d5) * w5;

        return qap_delta;
    }

    inline s64 get_u_qap_delta_size_7(const graph_t& g,
                                      const vertex_t u,
                                      const partition_t old_id,
                                      const partition_t new_id,
                                      const p_manager_t& p_manager,
                                      const d_oracle_t& d_oracle) {
        vertex_t v0 = g.neighbor(u, 0);
        vertex_t v1 = g.neighbor(u, 1);
        vertex_t v2 = g.neighbor(u, 2);
        vertex_t v3 = g.neighbor(u, 3);
        vertex_t v4 = g.neighbor(u, 4);
        vertex_t v5 = g.neighbor(u, 5);
        vertex_t v6 = g.neighbor(u, 6);

        weight_t w0 = g.weight(u, 0);
        weight_t w1 = g.weight(u, 1);
        weight_t w2 = g.weight(u, 2);
        weight_t w3 = g.weight(u, 3);
        weight_t w4 = g.weight(u, 4);
        weight_t w5 = g.weight(u, 5);
        weight_t w6 = g.weight(u, 6);

        partition_t v0_id = p_manager[v0];
        partition_t v1_id = p_manager[v1];
        partition_t v2_id = p_manager[v2];
        partition_t v3_id = p_manager[v3];
        partition_t v4_id = p_manager[v4];
        partition_t v5_id = p_manager[v5];
        partition_t v6_id = p_manager[v6];

        weight_t old_d0 = d_oracle.get(old_id, v0_id);
        weight_t old_d1 = d_oracle.get(old_id, v1_id);
        weight_t old_d2 = d_oracle.get(old_id, v2_id);
        weight_t old_d3 = d_oracle.get(old_id, v3_id);
        weight_t old_d4 = d_oracle.get(old_id, v4_id);
        weight_t old_d5 = d_oracle.get(old_id, v5_id);
        weight_t old_d6 = d_oracle.get(old_id, v6_id);

        weight_t new_d0 = d_oracle.get(new_id, v0_id);
        weight_t new_d1 = d_oracle.get(new_id, v1_id);
        weight_t new_d2 = d_oracle.get(new_id, v2_id);
        weight_t new_d3 = d_oracle.get(new_id, v3_id);
        weight_t new_d4 = d_oracle.get(new_id, v4_id);
        weight_t new_d5 = d_oracle.get(new_id, v5_id);
        weight_t new_d6 = d_oracle.get(new_id, v6_id);

        s64 qap_delta = 0;
        qap_delta += (old_d0 - new_d0) * w0;
        qap_delta += (old_d1 - new_d1) * w1;
        qap_delta += (old_d2 - new_d2) * w2;
        qap_delta += (old_d3 - new_d3) * w3;
        qap_delta += (old_d4 - new_d4) * w4;
        qap_delta += (old_d5 - new_d5) * w5;
        qap_delta += (old_d6 - new_d6) * w6;

        return qap_delta;
    }

    inline s64 get_u_qap_delta_size_8(const graph_t& g,
                                      const vertex_t u,
                                      const partition_t old_id,
                                      const partition_t new_id,
                                      const p_manager_t& p_manager,
                                      const d_oracle_t& d_oracle) {
        vertex_t v0 = g.neighbor(u, 0);
        vertex_t v1 = g.neighbor(u, 1);
        vertex_t v2 = g.neighbor(u, 2);
        vertex_t v3 = g.neighbor(u, 3);
        vertex_t v4 = g.neighbor(u, 4);
        vertex_t v5 = g.neighbor(u, 5);
        vertex_t v6 = g.neighbor(u, 6);
        vertex_t v7 = g.neighbor(u, 7);

        weight_t w0 = g.weight(u, 0);
        weight_t w1 = g.weight(u, 1);
        weight_t w2 = g.weight(u, 2);
        weight_t w3 = g.weight(u, 3);
        weight_t w4 = g.weight(u, 4);
        weight_t w5 = g.weight(u, 5);
        weight_t w6 = g.weight(u, 6);
        weight_t w7 = g.weight(u, 7);

        partition_t v0_id = p_manager[v0];
        partition_t v1_id = p_manager[v1];
        partition_t v2_id = p_manager[v2];
        partition_t v3_id = p_manager[v3];
        partition_t v4_id = p_manager[v4];
        partition_t v5_id = p_manager[v5];
        partition_t v6_id = p_manager[v6];
        partition_t v7_id = p_manager[v7];

        weight_t old_d0 = d_oracle.get(old_id, v0_id);
        weight_t old_d1 = d_oracle.get(old_id, v1_id);
        weight_t old_d2 = d_oracle.get(old_id, v2_id);
        weight_t old_d3 = d_oracle.get(old_id, v3_id);
        weight_t old_d4 = d_oracle.get(old_id, v4_id);
        weight_t old_d5 = d_oracle.get(old_id, v5_id);
        weight_t old_d6 = d_oracle.get(old_id, v6_id);
        weight_t old_d7 = d_oracle.get(old_id, v7_id);

        weight_t new_d0 = d_oracle.get(new_id, v0_id);
        weight_t new_d1 = d_oracle.get(new_id, v1_id);
        weight_t new_d2 = d_oracle.get(new_id, v2_id);
        weight_t new_d3 = d_oracle.get(new_id, v3_id);
        weight_t new_d4 = d_oracle.get(new_id, v4_id);
        weight_t new_d5 = d_oracle.get(new_id, v5_id);
        weight_t new_d6 = d_oracle.get(new_id, v6_id);
        weight_t new_d7 = d_oracle.get(new_id, v7_id);

        s64 qap_delta = 0;
        qap_delta += (old_d0 - new_d0) * w0;
        qap_delta += (old_d1 - new_d1) * w1;
        qap_delta += (old_d2 - new_d2) * w2;
        qap_delta += (old_d3 - new_d3) * w3;
        qap_delta += (old_d4 - new_d4) * w4;
        qap_delta += (old_d5 - new_d5) * w5;
        qap_delta += (old_d6 - new_d6) * w6;
        qap_delta += (old_d7 - new_d7) * w7;

        return qap_delta;
    }

    inline s64 get_u_qap_delta(const graph_t& g,
                               const vertex_t u,
                               const partition_t old_id,
                               const partition_t new_id,
                               const p_manager_t& p_manager,
                               const d_oracle_t& d_oracle) {
        switch (g.size(u)) {
        case 1:
            return get_u_qap_delta_size_1(g, u, old_id, new_id, p_manager, d_oracle);
        case 2:
            return get_u_qap_delta_size_2(g, u, old_id, new_id, p_manager, d_oracle);
        case 3:
            return get_u_qap_delta_size_3(g, u, old_id, new_id, p_manager, d_oracle);
        case 4:
            return get_u_qap_delta_size_4(g, u, old_id, new_id, p_manager, d_oracle);
        case 5:
            return get_u_qap_delta_size_5(g, u, old_id, new_id, p_manager, d_oracle);
        case 6:
            return get_u_qap_delta_size_6(g, u, old_id, new_id, p_manager, d_oracle);
        case 7:
            return get_u_qap_delta_size_7(g, u, old_id, new_id, p_manager, d_oracle);
        case 8:
            return get_u_qap_delta_size_8(g, u, old_id, new_id, p_manager, d_oracle);
        default: ;
        }

        s64 qap_delta = 0;

#pragma GCC unroll 8
        forall_guivw(g, u, i, v, w)
            {
                partition_t v_id = p_manager[v];
                weight_t old_d   = d_oracle.get(v_id, old_id);
                weight_t new_d   = d_oracle.get(v_id, new_id);

                qap_delta += (old_d - new_d) * w;
            }
        endfor

        return qap_delta;
    }

    inline s64 get_u_edge_cut_delta(const graph_t& g,
                                    const vertex_t u,
                                    const partition_t old_id,
                                    const partition_t new_id,
                                    const p_manager_t& p_manager) {
        s64 edge_cut_delta = 0;

#pragma GCC unroll 8
        forall_guivw(g, u, i, v, w)
            {
                partition_t v_id = p_manager[v];

                edge_cut_delta -= (v_id != new_id) * w;
                edge_cut_delta += (v_id != old_id) * w;
            }
        endfor

        return edge_cut_delta;
    }

    inline void get_u_qap_delta(const graph_t& g,
                                const vertex_t u,
                                const partition_t old_id,
                                const partition_t* blocks,
                                s64* blocks_qap_delta,
                                const size_t blocks_size,
                                const p_manager_t& p_manager,
                                const d_oracle_t& d_oracle) {
        blocks           = ASSUME_ALIGNED(partition_t*, blocks, 64);
        blocks_qap_delta = ASSUME_ALIGNED(s64 *, blocks_qap_delta, 64);

        // reset all to 0
        std::fill_n(blocks_qap_delta, blocks_size, 0);

#pragma GCC unroll 8
        forall_guivw(g, u, j, v, w)
            {
                partition_t v_id = p_manager[v];
                weight_t old_d   = d_oracle.get(v_id, old_id);

                for (size_t i = 0; i < blocks_size; ++i) {
                    weight_t new_d = d_oracle.get(v_id, blocks[i]);
                    blocks_qap_delta[i] += (old_d - new_d) * w;
                }
            }
        endfor
    }

    inline s64 get_qap_delta(const graph_t& g,
                             const vertex_t u,
                             const partition_t u_old_id,
                             const partition_t u_new_id,
                             const vertex_t v,
                             const partition_t v_old_id,
                             const partition_t v_new_id,
                             const p_manager_t& p_manager,
                             const d_oracle_t& d_oracle) {
        s64 qap_delta = 0;

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

    inline s64 get_qap_delta(const graph_t& g,
                             const vertex_t v,
                             const vertex_t vv,
                             const vertex_t vvv,
                             const partition_t v_id,
                             const partition_t vv_id,
                             const partition_t vvv_id,
                             const partition_t new_v_id,
                             const partition_t new_vv_id,
                             const partition_t new_vvv_id,
                             const p_manager_t& p_manager,
                             const d_oracle_t& d_oracle) {
        s64 qap_delta = 0;

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

    inline s64 get_u_qap_delta_and_is_boundary(const graph_t& g,
                                               const vertex_t u,
                                               const partition_t old_id,
                                               const partition_t new_id,
                                               bool& is_boundary_old_id,
                                               bool& is_boundary_new_id,
                                               const p_manager_t& p_manager,
                                               const d_oracle_t& d_oracle) {
        is_boundary_old_id = false;
        is_boundary_new_id = false;

        s64 qap_delta = 0;

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

    inline s64 get_u_qap_delta_and_is_connected_to(const graph_t& g,
                                                   const vertex_t u,
                                                   const partition_t old_id,
                                                   const partition_t new_id,
                                                   bool& is_connected_to_new_id,
                                                   const p_manager_t& p_manager,
                                                   const d_oracle_t& d_oracle) {
        is_connected_to_new_id = false;

        s64 qap_delta = 0;

#pragma GCC unroll 8
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

    inline s64 get_u_edge_cut_delta_and_is_connected_to(const graph_t& g,
                                                        const vertex_t u,
                                                        const partition_t old_id,
                                                        const partition_t new_id,
                                                        bool& is_connected_to_new_id,
                                                        const p_manager_t& p_manager) {
        is_connected_to_new_id = false;

        s64 edge_cut_delta = 0;

#pragma GCC unroll 8
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
}

#endif //HEIPROMAP_QAP_H
