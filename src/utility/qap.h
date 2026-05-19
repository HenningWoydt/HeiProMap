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

        for (vertex_t u = 0; u < g.n; ++u) {
            partition_t u_id = p_manager[u];

            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t v = g.edges_v[i];
                const weight_t w = g.edges_w[i];
                partition_t v_id = p_manager[v];
                weight_t d = d_oracle.get(u_id, v_id);
                qap += (d * w);
            }
        }

        return qap;
    }

    template<typename GraphT, typename DistanceOracleT>
    inline weight_t get_qap(const GraphT &g,
                            const AlignedArray<partition_t> &partition,
                            const DistanceOracleT &d_oracle) {
        ScopedTimer _t("misc", "misc", "get_qap");

        weight_t qap = 0;

        for (vertex_t u = 0; u < g.n; ++u) {
            partition_t u_id = partition[u];

            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t v = g.edges_v[i];
                const weight_t w = g.edges_w[i];
                partition_t v_id = partition[v];
                weight_t d = d_oracle.get(u_id, v_id);
                qap += (d * w);
            }
        }

        return qap;
    }

    template<typename GraphT, typename PartitionManagerT>
    inline weight_t get_edge_cut(GraphT &g,
                                 PartitionManagerT &p_manager) {
        ScopedTimer _t("misc", "misc", "get_qap");

        weight_t edge_cut = 0;

        for (vertex_t u = 0; u < g.n; ++u) {
            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t v = g.edges_v[i];
                const weight_t w = g.edges_w[i];
                if (p_manager[u] == p_manager[v]) { continue; }
                edge_cut += w;
            }
        }

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

            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t v = g.edges_v[i];
                const weight_t w = g.edges_w[i];
                partition_t v_id = p_manager[v];
                weight_t d = d_oracle.get(u_id, v_id);
                qap += (d * w);
            }
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

        for (vertex_t u = 0; u < g.n; ++u) {
            partition_t u_id = p_manager[u];

            for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) {
                const partition_t id = block_conn.get_id(i);
                const weight_t id_w = block_conn.get_w(i);
                weight_t d = d_oracle.get(u_id, id);
                local_qap += (d * id_w);
            }
        }
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
        for (vertex_t u = 0; u < g.n; ++u) {
            partition_t u_id = p_manager[u];

            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t v = g.edges_v[i];
                const weight_t w = g.edges_w[i];
                partition_t v_id = p_manager[v];
                weight_t d = d_oracle.get(u_id, v_id);
                partition_t layer_id = d_oracle.get_h(u_id, v_id);
                final_qap[layer_id] += (d * w);
            }
        }

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

        for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
            const vertex_t v = g.edges_v[i];
            const weight_t w = g.edges_w[i];
            partition_t v_id = p_manager[v];
            weight_t old_d = d_oracle.get(v_id, old_id);
            weight_t new_d = d_oracle.get(v_id, new_id);

            qap_delta += (old_d - new_d) * w;
        }

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

        for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) {
            const partition_t id = block_conn.get_id(i);
            const weight_t idw = block_conn.get_w(i);
            weight_t old_d = d_oracle.get(id, old_id);
            weight_t new_d = d_oracle.get(id, new_id);

            qap_delta += (old_d - new_d) * idw;
        }

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
        for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) {
            const partition_t id = block_conn.get_id(i);
            const weight_t idw = block_conn.get_w(i);
            weight_t old_d = d_oracle.get(id, old_id);

            for (size_t j = 0; j < new_ids.size(); j++) {
                qap_deltas[j] += (old_d - d_oracle.get(id, new_ids[j])) * idw;
            }
        }

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

        for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
            const vertex_t v = g.edges_v[i];
            const weight_t w = g.edges_w[i];
            partition_t v_id = p_manager[v];

            edge_cut_delta -= (v_id != new_id) * w;
            edge_cut_delta += (v_id != old_id) * w;
        }

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

        for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) {
            const partition_t id = block_conn.get_id(i);
            const weight_t idw = block_conn.get_w(i);
            edge_cut_delta -= (id != new_id) * idw;
            edge_cut_delta += (id != old_id) * idw;
        }

        return edge_cut_delta;
    }


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
        for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
            const vertex_t v = g.edges_v[i];
            const weight_t w = g.edges_w[i];
            partition_t v_id = p_manager[v];

            is_boundary_old_id |= (v_id != old_id);
            is_boundary_new_id |= (v_id != new_id);

            weight_t old_d = d_oracle.get(v_id, old_id);
            weight_t new_d = d_oracle.get(v_id, new_id);
            qap_delta += (old_d - new_d) * w;
        }

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

        for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) {
            const partition_t id = block_conn.get_id(i);
            const weight_t idw = block_conn.get_w(i);
            is_boundary_old_id |= (id != old_id);
            is_boundary_new_id |= (id != new_id);

            weight_t old_d = d_oracle.get(id, old_id);
            weight_t new_d = d_oracle.get(id, new_id);
            qap_delta += (old_d - new_d) * idw;
        }


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

        for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
            const vertex_t v = g.edges_v[i];
            const weight_t w = g.edges_w[i];
            partition_t v_id = p_manager[v];

            is_connected_to_new_id |= (v_id == new_id);

            weight_t old_d = d_oracle.get(v_id, old_id);
            weight_t new_d = d_oracle.get(v_id, new_id);

            qap_delta += (old_d - new_d) * w;
        }

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

        for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) {
            const partition_t id = block_conn.get_id(i);
            const weight_t idw = block_conn.get_w(i);
            is_connected_to_new_id |= (id == new_id);

            weight_t old_d = d_oracle.get(id, old_id);
            weight_t new_d = d_oracle.get(id, new_id);

            qap_delta += (old_d - new_d) * idw;
        }

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

        for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
            const vertex_t v = g.edges_v[i];
            const weight_t w = g.edges_w[i];
            partition_t v_id = p_manager[v];

            is_connected_to_new_id |= (v_id == new_id);

            edge_cut_delta -= (v_id != new_id) * w;
            edge_cut_delta += (v_id != old_id) * w;
        }

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

        for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) {
            const partition_t id = block_conn.get_id(i);
            const weight_t idw = block_conn.get_w(i);
            is_connected_to_new_id |= (id == new_id);

            edge_cut_delta -= (id != new_id) * idw;
            edge_cut_delta += (id != old_id) * idw;
        }

        return edge_cut_delta;
    }

    template<bool t_uniform_e_weights, typename GraphT, typename PartitionManagerT, typename DistanceOracleT, typename BlockConnT>
    inline weight_t get_u_qap_delta_t(const GraphT &g,
                                      const vertex_t u,
                                      const partition_t old_id,
                                      const partition_t new_id,
                                      const PartitionManagerT &p_manager,
                                      DistanceOracleT &d_oracle,
                                      BlockConnT &block_conn) {
        weight_t qap_delta = 0;

        for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) {
            const partition_t id = block_conn.get_id(i);
            const weight_t idw = block_conn.get_w(i);

            qap_delta += (d_oracle.get(id, old_id) - d_oracle.get(id, new_id)) * idw;
        }

        return qap_delta;
    }

    template<bool t_uniform_e_weights, typename GraphT, typename PartitionManagerT, typename DistanceOracleT>
    inline weight_t get_u_qap_delta_and_is_connected_to_t(const GraphT &g,
                                                          const vertex_t u,
                                                          const partition_t old_id,
                                                          const partition_t new_id,
                                                          bool &is_connected_to_new_id,
                                                          const PartitionManagerT &p_manager,
                                                          DistanceOracleT &d_oracle) {
        is_connected_to_new_id = false;
        weight_t qap_delta = 0;

        if constexpr (t_uniform_e_weights) {
            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t v = g.edges_v[i];
                partition_t v_id = p_manager[v];
                is_connected_to_new_id |= (v_id == new_id);
                qap_delta += d_oracle.get(v_id, old_id) - d_oracle.get(v_id, new_id);
            }
        } else {
            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t v = g.edges_v[i];
                const weight_t w = g.edges_w[i];
                partition_t v_id = p_manager[v];
                is_connected_to_new_id |= (v_id == new_id);
                qap_delta += (d_oracle.get(v_id, old_id) - d_oracle.get(v_id, new_id)) * w;
            }
        }

        return qap_delta;
    }

    template<bool t_uniform_e_weights, typename GraphT, typename PartitionManagerT, typename BlockConnT>
    inline weight_t get_u_edge_cut_delta_t(const GraphT &g,
                                           const vertex_t u,
                                           const partition_t old_id,
                                           const partition_t new_id,
                                           const PartitionManagerT &p_manager,
                                           const BlockConnT &block_conn) {
        weight_t edge_cut_delta = 0;

        for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) {
            const partition_t id = block_conn.get_id(i);
            const weight_t idw = block_conn.get_w(i); {
                edge_cut_delta -= (id != new_id) * idw;
                edge_cut_delta += (id != old_id) * idw;
            }
        }

        return edge_cut_delta;
    }

    template<bool t_uniform_e_weights, typename GraphT, typename PartitionManagerT>
    inline weight_t get_u_edge_cut_delta_and_is_connected_to_t(const GraphT &g,
                                                               const vertex_t u,
                                                               const partition_t old_id,
                                                               const partition_t new_id,
                                                               bool &is_connected_to_new_id,
                                                               const PartitionManagerT &p_manager) {
        is_connected_to_new_id = false;
        weight_t edge_cut_delta = 0;

        if constexpr (t_uniform_e_weights) {
            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t v = g.edges_v[i]; {
                    partition_t v_id = p_manager[v];
                    is_connected_to_new_id |= (v_id == new_id);
                    edge_cut_delta -= (v_id != new_id);
                    edge_cut_delta += (v_id != old_id);
                }
            }
        } else {
            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t v = g.edges_v[i];
                const weight_t w = g.edges_w[i]; {
                    partition_t v_id = p_manager[v];
                    is_connected_to_new_id |= (v_id == new_id);
                    edge_cut_delta -= (v_id != new_id) * w;
                    edge_cut_delta += (v_id != old_id) * w;
                }
            }
        }

        return edge_cut_delta;
    }
}

#endif //HEIPROMAP_QAP_H
