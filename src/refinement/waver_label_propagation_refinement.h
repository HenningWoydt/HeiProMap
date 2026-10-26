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

#ifndef HEIPROMAP_WAVER_LABEL_PROPAGATION_REFINEMENT_H
#define HEIPROMAP_WAVER_LABEL_PROPAGATION_REFINEMENT_H

#include <cmath>
#include <limits>
#include <vector>

#include "../definitions.h"
#include "../datastructures/block_conn.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../utility/aligned_array.h"
#include "../utility/profiler.h"
#include "../utility/random_engine.h"

namespace HeiProMap {

    class WaverLabelPropagationConfiguration {
    public:
        explicit WaverLabelPropagationConfiguration(std::string t_name) : name(std::move(t_name)) {}

        std::string name;
        bool enabled = false;
        u64 max_iteration = 25;

        // four-term cost model parameters
        f64 alpha_end  = 0.0;   // alpha_send + alpha_recv
        f64 beta_end   = 0.0;   // beta_send  + beta_recv
        f64 alpha_link = 0.0;
        f64 beta_link  = 0.0;
    };

    // Single-threaded label propagation: every boundary vertex moves to the neighboring block with the best
    // four-term cost gain.
    class WaverLabelPropagationRefinement {
        RandomEngine random_engine;
        WaverLabelPropagationConfiguration config = WaverLabelPropagationConfiguration("waver_label_propagation");

        static inline __attribute__((always_inline)) void apply_move(graph_t &g,
                                                                     bv_manager_t &bv_manager,
                                                                     p_manager_t &p_manager,
                                                                     const vertex_t u,
                                                                     const weight_t u_weight,
                                                                     const partition_t from_id,
                                                                     const partition_t to_id) {
            bv_manager.move(g, p_manager, u, from_id, to_id);
            p_manager.move_serial(u, u_weight, from_id, to_id);
        }

        // Four-term cost delta for moving u from old_id to new_id.
        // Returns -Δcost; positive means the move reduces total cost.
        //
        // Δcost = alpha_end  * Δ|E_π|
        //       + beta_end   * ΔΣS
        //       + alpha_link * ΔΣh
        //       + beta_link  * ΔΣh·S
        //
        // Δ|E_π| approximation: alpha_end is paid once per communicating (src,dst)
        // PE pair (due to message coalescing), NOT once per edge. We approximate
        // Δ|E_π| using boolean indicators over u's neighbors:
        //   - any neighbor in old_id → u's A-A edges become B-A cuts → B-A pair appears: -alpha_end
        //   - any neighbor in new_id → u's A-B cuts become B-B internal → A-B pair may disappear: +alpha_end
        // This charges alpha_end at most once per direction across the A-B boundary,
        // replacing the previous per-edge overcounting.
        //
        // ΔΣS, ΔΣh, ΔΣh·S: per-edge (correct — each edge contributes independently).
        template<bool t_uniform_e_weights, typename DistanceOracleT, typename QGraphT>
        inline __attribute__((always_inline)) f64 compute_cost_delta(const graph_t &g,
                                                                     const vertex_t u,
                                                                     const partition_t old_id,
                                                                     const partition_t new_id,
                                                                     const p_manager_t &p_manager,
                                                                     DistanceOracleT &d_oracle,
                                                                     const QGraphT &q_graph,
                                                                     const block_conn_t &block_conn) const {
            const f64 dist_old_new = static_cast<f64>(d_oracle.get(old_id, new_id));
            f64 delta = 0.0;

            // alpha_end: exact Δ|E_π| via block_conn (u's per-block weights) and q_graph (block-pair totals).
            // u moves from old_id (A) to new_id (B). Positive delta = cost reduction.
            for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) {
                const partition_t C       = block_conn.get_id(i);
                const weight_t    u_to_C  = block_conn.get_w(i);

                if (C == old_id) {
                    // A-A internal → B-A cut; new pair if B-A didn't previously exist
                    if (q_graph.get_weight(new_id, old_id) == 0) {
                        delta -= config.alpha_end;
                    }
                } else if (C == new_id) {
                    // A-B cut → B-B internal; pair disappears iff u was the sole A-B connection
                    if (q_graph.get_weight(old_id, new_id) == u_to_C) {
                        delta += config.alpha_end;
                    }
                } else {
                    // A-C cut → B-C cut
                    if (q_graph.get_weight(old_id, C) == u_to_C) {
                        delta += config.alpha_end;  // A-C pair disappears
                    }
                    if (q_graph.get_weight(new_id, C) == 0) {
                        delta -= config.alpha_end;  // B-C pair appears
                    }
                }
            }

            // beta_end, alpha_link, beta_link: per-edge terms
            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t    v    = g.edges_v[i];
                const partition_t v_id = p_manager[v];
                const f64         w    = t_uniform_e_weights ? 1.0 : static_cast<f64>(g.edge_w(i));

                if (v_id == old_id) {
                    delta -= config.beta_end * w
                           + (config.alpha_link + config.beta_link * w) * dist_old_new;
                } else if (v_id == new_id) {
                    delta += config.beta_end * w
                           + (config.alpha_link + config.beta_link * w) * dist_old_new;
                } else {
                    const f64 d_old = static_cast<f64>(d_oracle.get(v_id, old_id));
                    const f64 d_new = static_cast<f64>(d_oracle.get(v_id, new_id));
                    delta += (config.alpha_link + config.beta_link * w) * (d_old - d_new);
                }
            }
            return delta;
        }

    public:
        WaverLabelPropagationRefinement() = default;
        ~WaverLabelPropagationRefinement() = default;

        void initialize(const u64 t_seed, const WaverLabelPropagationConfiguration &t_config) {
            config = t_config;
            random_engine = RandomEngine(t_seed);
        }

        template<typename DistanceOracleT, typename QGraphT>
        void refine(graph_t &g,
                    DistanceOracleT &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    QGraphT &q_graph,
                    block_conn_t &block_conn) {
            if (g.uniform_v_weights && g.uniform_e_weights) {
                refine_impl<true, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn);
            } else if (g.uniform_v_weights) {
                refine_impl<true, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn);
            } else if (g.uniform_e_weights) {
                refine_impl<false, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn);
            } else {
                refine_impl<false, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn);
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights, typename DistanceOracleT, typename QGraphT>
        void refine_impl(graph_t &g,
                         DistanceOracleT &d_oracle,
                         bv_manager_t &bv_manager,
                         p_manager_t &p_manager,
                         QGraphT &q_graph,
                         block_conn_t &block_conn) {
            bool positive_move_occurred = true;
            for (u64 iteration = 0; iteration < config.max_iteration && positive_move_occurred; ++iteration) {
                positive_move_occurred = false;

                HEIPROMAP_PROFILE_SCOPE("refinement", "WaverLPRefinement", "process_vertices");
                for (vertex_t u = 0; u < g.n; ++u) {
                    if (!bv_manager.is_boundary(u)) { continue; }

                    const weight_t u_weight = t_uniform_v_weights ? 1 : g.v_weight(u);
                    const partition_t u_id  = p_manager[u];

                    partition_t best_id    = NO_ID;
                    f64         best_delta = -std::numeric_limits<f64>::max();
                    f32         counter    = 0.0f;

                    partition_t neighbor_blocks[256];
                    size_t n_neighbor_blocks = 0;

                    for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                        const partition_t v_id = p_manager[g.edges_v[i]];
                        if (v_id == u_id) { continue; }

                        bool found = false;
                        for (size_t j = 0; j < n_neighbor_blocks; ++j) {
                            if (neighbor_blocks[j] == v_id) { found = true; break; }
                        }
                        if (!found && n_neighbor_blocks < 256) {
                            neighbor_blocks[n_neighbor_blocks++] = v_id;
                        }
                    }

                    if (p_manager.get_bweight(u_id) - u_weight < p_manager.lmin[u_id]) { continue; }

                    for (size_t bi = 0; bi < n_neighbor_blocks; ++bi) {
                        const partition_t id = neighbor_blocks[bi];
                        if (p_manager.get_bweight(id) + u_weight > p_manager.lmax[id]) { continue; }

                        const f64 delta = compute_cost_delta<t_uniform_e_weights>(g, u, u_id, id, p_manager, d_oracle, q_graph, block_conn);

                        if (delta > best_delta) {
                            best_id    = id;
                            best_delta = delta;
                            counter    = 1.0f;
                        } else if (delta == best_delta) {
                            counter += 1.0f;
                            if (random_engine.get_f32() < 1.0f / counter) {
                                best_id = id;
                            }
                        }
                    }

                    if (best_delta > 0.0 || (best_delta == 0.0 && random_engine.get_f32() < 0.5f)) {
                        q_graph.move(g, p_manager, u, u_id, best_id);
                        block_conn.move(g, u, u_id, best_id);
                        apply_move(g, bv_manager, p_manager, u, u_weight, u_id, best_id);
                        positive_move_occurred |= (best_delta > 0.0);
                    }
                }
            }
        }
    };

}

#endif //HEIPROMAP_WAVER_LABEL_PROPAGATION_REFINEMENT_H
