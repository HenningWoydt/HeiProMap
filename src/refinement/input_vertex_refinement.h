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

#ifndef HEIPROMAP_INPUT_VERTEX_REFINEMENT_H
#define HEIPROMAP_INPUT_VERTEX_REFINEMENT_H

#include <algorithm>
#include <vector>

#include "../datastructures/block_conn.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/large_quotient_graph.h"
#include "../datastructures/partition_manager.h"
#include "../definitions.h"
#include "../utility/profiler.h"
#include "../utility/qap.h"
#include "../utility/utils.h"
#include "waver_message_refinement.h"

namespace HeiProMap {
    // Refines the placement of the individual input vertices after the block -> PE mapping has been propagated back
    // from the quotient graph: runs WaverMessageRefinement on the input graph with its own, input-sized data
    // structures, so the solver's quotient-sized state is left untouched. No PE may exceed the imbalance the user
    // asked for (the same one KaMinPar was given).
    class InputVertexRefinement {
    public:
        struct Params {
            f64 imbalance = 0.03;     // the user's imbalance: lmax = (1 + imbalance) * average block weight
            f64 min_imbalance = 0.5;  // lmin = min_imbalance * average block weight (a block may shrink to half)
            u64 max_iters = 10;
            f64 alpha_end = 0.0;
            f64 beta_end = 0.0;
            f64 alpha_link = 0.0;
            f64 beta_link = 0.0;
        };

        struct Stats {
            bool ran = false;
            f64 cost_before = 0.0;
            f64 cost_after = 0.0;
            f64 max_link_before = 0.0;
            f64 max_link_after = 0.0;
            weight_t max_rows_before = 0;
            weight_t max_rows_after = 0;
            weight_t min_rows_after = 0;
            f64 seconds = 0.0;
        };

        // part[u]: PE of input vertex u; updated in place.
        template<typename DistanceOracleT>
        static Stats refine(graph_t &g,
                            std::vector<partition_t> &part,
                            const DistanceOracleT &d_oracle,
                            const Params &params) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "InputVertexRefinement", "refine");
            const auto sp = get_time_point();
            const partition_t k = d_oracle.current_k();

            PartitionManager pm;
            pm.initialize(g.n, k, g.g_weight);
            pm.activate_all();
            pm.update_limits(d_oracle.current_topology(), g.g_weight, g.n, params.imbalance, params.min_imbalance);
            for (vertex_t u = 0; u < g.n; ++u) { pm.partition[u] = part[u]; }
            pm.recalculate_weights(g);

            BoundaryVertexManager bv;
            bv.initialize(g.n, k);
            for (vertex_t u = 0; u < g.n; ++u) {
                size_t n_different = 0;
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    n_different += (pm[u] != pm[g.edges_v[i]]);
                }
                bv.add_boundary_vertex_from_count(u, pm[u], n_different);
            }
            LargeQuotientGraph q_graph;
            q_graph.compute_from_scratch(g, pm);
            BlockConn block_conn;
            block_conn.compute_from_scratch(g, pm);

            auto max_rows = [&]() {
                weight_t m = 0;
                for (partition_t id = 0; id < k; ++id) { m = std::max(m, pm.get_bweight(id)); }
                return m;
            };

            Stats stats;
            stats.ran = true;
            stats.cost_before = get_message_cost(q_graph, k, d_oracle, params.alpha_end, params.beta_end,
                                                 params.alpha_link, params.beta_link).total;
            stats.max_link_before = get_max_link_load(q_graph, k, d_oracle.current_W(), params.alpha_link, params.beta_link);
            stats.max_rows_before = max_rows();

            WaverMessageRefinementConfiguration config("Input Vertex Refinement");
            config.enabled = true;
            config.max_iters = params.max_iters;
            WaverMessageRefinement::refine(g, pm, bv, q_graph, d_oracle, block_conn, config,
                                           params.alpha_end, params.beta_end, params.alpha_link, params.beta_link);

            stats.cost_after = get_message_cost(q_graph, k, d_oracle, params.alpha_end, params.beta_end,
                                                params.alpha_link, params.beta_link).total;
            stats.max_link_after = get_max_link_load(q_graph, k, d_oracle.current_W(), params.alpha_link, params.beta_link);
            stats.max_rows_after = max_rows();
            stats.min_rows_after = pm.get_bweight(0);
            for (partition_t id = 1; id < k; ++id) { stats.min_rows_after = std::min(stats.min_rows_after, pm.get_bweight(id)); }

            for (vertex_t u = 0; u < g.n; ++u) { part[u] = pm[u]; }
            stats.seconds = get_seconds(sp, get_time_point());
            return stats;
        }
    };
}

#endif //HEIPROMAP_INPUT_VERTEX_REFINEMENT_H
