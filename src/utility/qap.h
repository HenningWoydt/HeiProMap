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

#include <algorithm>
#include <iomanip>
#include <vector>

#include "../definitions.h"
#include "macros.h"
#include "aligned_array.h"
#include "../datastructures/csr_graph.h"
#include "../distance_oracles/matrix_distance_oracle.h"
#include "../datastructures/dyn_graph.h"
#include "../datastructures/partition_manager.h"

namespace HeiProMap {
    template<typename DistanceOracleT>
    inline weight_t get_qap(const DynGraph &g,
                            const std::vector<partition_t> &partition,
                            const DistanceOracleT &d_oracle) {
        HEIPROMAP_PROFILE_SCOPE("misc", "misc", "get_qap");

        weight_t qap = 0;

        for (vertex_t u = 0; u < g.n; ++u) {
            if (!g.vertex_exists(u)) continue;
            partition_t u_id = partition[u];

            for (const auto &edge: g.neighbors[u]) {
                const vertex_t v = edge.u;
                const weight_t w = edge.w;
                partition_t v_id = partition[v];
                weight_t d = d_oracle.get(u_id, v_id);
                qap += (d * w);
            }
        }

        return qap;
    }

    template<typename GraphT, typename PartitionManagerT, typename DistanceOracleT>
    inline weight_t get_qap(const GraphT &g,
                            const PartitionManagerT &p_manager,
                            const DistanceOracleT &d_oracle) {
        HEIPROMAP_PROFILE_SCOPE("misc", "misc", "get_qap");

        weight_t qap = 0;

        for (vertex_t u = 0; u < g.n; ++u) {
            partition_t u_id = p_manager[u];

            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t v = g.edges_v[i];
                const weight_t w = g.edge_w(i);
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
        HEIPROMAP_PROFILE_SCOPE("misc", "misc", "get_qap");

        weight_t qap = 0;

        for (vertex_t u = 0; u < g.n; ++u) {
            partition_t u_id = partition[u];

            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t v = g.edges_v[i];
                const weight_t w = g.edge_w(i);
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
        HEIPROMAP_PROFILE_SCOPE("misc", "misc", "get_qap");

        weight_t edge_cut = 0;

        for (vertex_t u = 0; u < g.n; ++u) {
            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t v = g.edges_v[i];
                const weight_t w = g.edge_w(i);
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
                const weight_t w = g.edge_w(i);
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
        HEIPROMAP_PROFILE_SCOPE("misc", "misc", "get_qap");

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
        HEIPROMAP_PROFILE_SCOPE("misc", "misc", "get_qap_per_layer");

        std::vector<weight_t> final_qap(l, 0);
        for (vertex_t u = 0; u < g.n; ++u) {
            partition_t u_id = p_manager[u];

            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t v = g.edges_v[i];
                const weight_t w = g.edge_w(i);
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
            const weight_t w = g.edge_w(i);
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
            const weight_t w = g.edge_w(i);
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


    template<typename GraphT, typename DistanceOracleT>
    inline weight_t get_u_qap_delta_and_is_boundary(const GraphT &g,
                                                    const vertex_t u,
                                                    const partition_t old_id,
                                                    const partition_t new_id,
                                                    bool &is_boundary_old_id,
                                                    bool &is_boundary_new_id,
                                                    const p_manager_t &p_manager,
                                                    DistanceOracleT &d_oracle) {
        is_boundary_old_id = false;
        is_boundary_new_id = false;

        weight_t qap_delta = 0;

        #pragma GCC unroll 8
        for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
            const vertex_t v = g.edges_v[i];
            const weight_t w = g.edge_w(i);
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
            const weight_t w = g.edge_w(i);
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
            const weight_t w = g.edge_w(i);
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

    template<typename GraphT, typename PartitionManagerT, typename DistanceOracleT, typename BlockConnT>
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
                const weight_t w = g.edge_w(i);
                partition_t v_id = p_manager[v];
                is_connected_to_new_id |= (v_id == new_id);
                qap_delta += (d_oracle.get(v_id, old_id) - d_oracle.get(v_id, new_id)) * w;
            }
        }

        return qap_delta;
    }

    template<typename GraphT, typename PartitionManagerT, typename BlockConnT>
    inline weight_t get_u_edge_cut_delta_t(const GraphT &g,
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
                const vertex_t v = g.edges_v[i];

                partition_t v_id = p_manager[v];
                is_connected_to_new_id |= (v_id == new_id);
                edge_cut_delta -= (v_id != new_id);
                edge_cut_delta += (v_id != old_id);
            }
        } else {
            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t v = g.edges_v[i];
                const weight_t w = g.edge_w(i);
                partition_t v_id = p_manager[v];
                is_connected_to_new_id |= (v_id == new_id);
                edge_cut_delta -= (v_id != new_id) * w;
                edge_cut_delta += (v_id != old_id) * w;
            }
        }

        return edge_cut_delta;
    }

    template<typename QuotientGraphT, typename DistanceOracleT>
    inline weight_t get_qap_from_quotient_graph(const QuotientGraphT &q_graph,
                                                const partition_t k,
                                                const DistanceOracleT &d_oracle) {
        HEIPROMAP_PROFILE_SCOPE("misc", "misc", "get_qap_from_quotient_graph");
        weight_t qap = 0;
        for (partition_t u = 0; u < k; ++u) {
            q_graph.for_each_neighbor(u, [&](const partition_t v, const weight_t w) {
                if (u < v) {
                    qap += (d_oracle.get(u, v) * w * 2);
                }
            });
        }
        return qap;
    }

    struct MessageCostBreakdown {
        // raw sums (independent of alpha/beta)
        f64 message_count    = 0.0;  // |E_pi|: number of directed messages
        f64 halo_transfers   = 0.0;  // sum S(m): total words crossing partition boundary
        f64 hop_count        = 0.0;  // sum h(m): total hops across all messages
        f64 hop_volume       = 0.0;  // sum h(m)*S(m): hop-weighted volume

        // weighted terms
        f64 term_cardinality = 0.0;  // alpha_end  * message_count
        f64 term_cut_volume  = 0.0;  // beta_end   * halo_transfers
        f64 term_hop_count   = 0.0;  // alpha_link * hop_count
        f64 term_hop_bytes   = 0.0;  // beta_link  * hop_volume
        f64 total            = 0.0;

        void print(std::ostream &os) const {
            auto row = [&](const char *label, f64 val) {
                os << std::left << std::setw(36) << label
                   << ": " << std::right << std::setw(16) << static_cast<u64>(val) << "\n";
            };
            row("Message count (|E_pi|)",           message_count);
            row("Halo transfers (sum S)",            halo_transfers);
            row("Hop count (sum h)",                 hop_count);
            row("Hop volume (sum h*S)",              hop_volume);
            row("Cardinality (alpha_end * |E_pi|)",  term_cardinality);
            row("Cut volume  (beta_end  * sum S)",   term_cut_volume);
            row("Hop count   (alpha_link * sum h)",  term_hop_count);
            row("Hop-bytes   (beta_link * sum h*S)", term_hop_bytes);
            row("Total message cost",                total);
        }
    };

    // Computes the full four-term communication cost model over the quotient graph.
    // Each undirected PE pair (p,q) with p<q represents two directed messages (p->q and q->p),
    // so cardinality and volume-dependent terms are doubled to match the directed formulation.
    template<typename QuotientGraphT, typename DistanceOracleT>
    inline MessageCostBreakdown get_message_cost(const QuotientGraphT &q_graph,
                                                 const partition_t k,
                                                 const DistanceOracleT &d_oracle,
                                                 f64 alpha_end,
                                                 f64 beta_end,
                                                 f64 alpha_link,
                                                 f64 beta_link) {
        HEIPROMAP_PROFILE_SCOPE("misc", "misc", "get_message_cost");

        f64 sum_cardinality = 0.0;
        f64 sum_S           = 0.0;
        f64 sum_h           = 0.0;
        f64 sum_hS          = 0.0;

        for (partition_t u = 0; u < k; ++u) {
            q_graph.for_each_neighbor(u, [&](const partition_t v, const weight_t w) {
                if (u < v) {
                    f64 S = static_cast<f64>(w);
                    f64 h = static_cast<f64>(d_oracle.get(u, v));
                    // multiply by 2: each undirected pair = two directed messages
                    sum_cardinality += 2.0;
                    sum_S           += 2.0 * S;
                    sum_h           += 2.0 * h;
                    sum_hS          += 2.0 * h * S;
                }
            });
        }

        MessageCostBreakdown result;
        result.message_count    = sum_cardinality;
        result.halo_transfers   = sum_S;
        result.hop_count        = sum_h;
        result.hop_volume       = sum_hS;
        result.term_cardinality = alpha_end  * sum_cardinality;
        result.term_cut_volume  = beta_end   * sum_S;
        result.term_hop_count   = alpha_link * sum_h;
        result.term_hop_bytes   = beta_link  * sum_hS;
        result.total            = result.term_cardinality + result.term_cut_volume
                                + result.term_hop_count   + result.term_hop_bytes;
        return result;
    }

    // Busiest PE endpoint: a PE sends one message and receives one message per neighbor, so its load is
    // sum over neighbors of (alpha_end + beta_end * S). Reported to compare against the evaluator's resource bound.
    struct EndpointLoad {
        f64 max_load = 0.0;
        partition_t max_pe = 0;
        f64 max_messages = 0.0;  // largest neighbor count of any PE
    };

    template<typename QuotientGraphT>
    inline EndpointLoad get_max_endpoint_load(const QuotientGraphT &q_graph,
                                              const partition_t k,
                                              const f64 alpha_end,
                                              const f64 beta_end) {
        EndpointLoad result;
        for (partition_t u = 0; u < k; ++u) {
            f64 load = 0.0;
            f64 messages = 0.0;
            q_graph.for_each_neighbor(u, [&](const partition_t, const weight_t w) {
                load += alpha_end + beta_end * static_cast<f64>(w);
                messages += 1.0;
            });
            if (load > result.max_load) {
                result.max_load = load;
                result.max_pe = u;
            }
            result.max_messages = std::max(result.max_messages, messages);
        }
        return result;
    }

    // Computes the maximum load on any directed link under strict XY routing.
    // For each directed message (u->v), the load contribution to each link on its
    // route is (alpha_link + beta_link * S(m)). Returns the max over all directed links.
    template<typename QuotientGraphT>
    inline f64 get_max_link_load(const QuotientGraphT &q_graph,
                                 const partition_t k,
                                 const partition_t W,
                                 f64 alpha_link,
                                 f64 beta_link) {
        HEIPROMAP_PROFILE_SCOPE("misc", "misc", "get_max_link_load");

        const partition_t H = (k + W - 1) / W;

        // Directed link id encoding (strict-XY routing):
        //   East  [0,              H*(W-1)):  east_id(r,c)  = r*(W-1) + c,           c in [0, W-2]
        //   West  [H*(W-1),      2*H*(W-1)):  west_id(r,c)  = H*(W-1)+r*(W-1)+(c-1), c in [1, W-1]
        //   South [2*H*(W-1),    2*H*(W-1)+W*(H-1)):  south_id(c,r) = 2*H*(W-1)+c*(H-1)+r,       r in [0, H-2]
        //   North [2*H*(W-1)+W*(H-1), ...):  north_id(c,r) = 2*H*(W-1)+W*(H-1)+c*(H-1)+(r-1), r in [1, H-1]
        const size_t n_links = 2 * H * (W - 1) + 2 * W * (H - 1);
        std::vector<f64> load(n_links, 0.0);

        auto east_id  = [&](partition_t r, partition_t c) -> size_t {
            return static_cast<size_t>(r) * (W - 1) + c;
        };
        auto west_id  = [&](partition_t r, partition_t c) -> size_t {
            return static_cast<size_t>(H) * (W - 1) + static_cast<size_t>(r) * (W - 1) + (c - 1);
        };
        auto south_id = [&](partition_t c, partition_t r) -> size_t {
            return 2 * static_cast<size_t>(H) * (W - 1) + static_cast<size_t>(c) * (H - 1) + r;
        };
        auto north_id = [&](partition_t c, partition_t r) -> size_t {
            return 2 * static_cast<size_t>(H) * (W - 1) + static_cast<size_t>(W) * (H - 1)
                 + static_cast<size_t>(c) * (H - 1) + (r - 1);
        };

        auto route_and_add = [&](partition_t src, partition_t dst, f64 cost) {
            partition_t sx = src % W, sy = src / W;
            partition_t dx = dst % W, dy = dst / W;
            for (partition_t c = sx; c < dx; ++c) load[east_id(sy, c)]  += cost;
            for (partition_t c = sx; c > dx; --c) load[west_id(sy, c)]  += cost;
            for (partition_t r = sy; r < dy; ++r) load[south_id(dx, r)] += cost;
            for (partition_t r = sy; r > dy; --r) load[north_id(dx, r)] += cost;
        };

        for (partition_t u = 0; u < k; ++u) {
            q_graph.for_each_neighbor(u, [&](const partition_t v, const weight_t w) {
                if (u < v) {
                    f64 cost = alpha_link + beta_link * static_cast<f64>(w);
                    route_and_add(u, v, cost);
                    route_and_add(v, u, cost);
                }
            });
        }

        return load.empty() ? 0.0 : *std::max_element(load.begin(), load.end());
    }

    // Full four-term WSE cost: alpha_end*|E_π| + beta_end*ΣS + alpha_link*Σh + beta_link*Σh·S
    // Each edge is counted once (undirected: we iterate all directed half-edges, halve at end).
    template<typename GraphT, typename PartitionManagerT, typename DistanceOracleT>
    inline f64 get_waver_cost(const GraphT &g,
                              const PartitionManagerT &p_manager,
                              const DistanceOracleT &d_oracle,
                              const f64 alpha_end,
                              const f64 beta_end,
                              const f64 alpha_link,
                              const f64 beta_link) {
        // collect unique communicating block pairs for alpha_end
        std::vector<u64> block_pairs;
        block_pairs.reserve(g.m / 4);

        f64 sum_S   = 0.0;
        f64 sum_h   = 0.0;
        f64 sum_hS  = 0.0;

        for (vertex_t u = 0; u < g.n; ++u) {
            const partition_t u_id = p_manager[u];
            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t    v    = g.edges_v[i];
                const partition_t v_id = p_manager[v];
                if (u_id == v_id) { continue; }
                const f64 w = static_cast<f64>(g.edge_w(i));
                const f64 d = static_cast<f64>(d_oracle.get(u_id, v_id));
                sum_S  += w;
                sum_h  += d;
                sum_hS += w * d;
                // pack ordered pair so each undirected pair is unique
                const partition_t lo = u_id < v_id ? u_id : v_id;
                const partition_t hi = u_id < v_id ? v_id : u_id;
                block_pairs.push_back((static_cast<u64>(lo) << 32) | static_cast<u64>(hi));
            }
        }

        // edges were visited twice (both directions) — halve continuous sums
        sum_S  *= 0.5;
        sum_h  *= 0.5;
        sum_hS *= 0.5;

        // count distinct block pairs
        std::sort(block_pairs.begin(), block_pairs.end());
        const f64 n_pairs = static_cast<f64>(std::unique(block_pairs.begin(), block_pairs.end()) - block_pairs.begin());

        return alpha_end * n_pairs + beta_end * sum_S + alpha_link * sum_h + beta_link * sum_hS;
    }

}

#endif //HEIPROMAP_QAP_H
