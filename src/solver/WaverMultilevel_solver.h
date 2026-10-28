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

#ifndef HEIPROMAP_WAVERMULTILEVEL_SOLVER_H
#define HEIPROMAP_WAVERMULTILEVEL_SOLVER_H

#include "WaverBase_solver.h"

#include "../init_partitioning/greedy_partitioner.h"
#include "../init_partitioning/waver_greedy_partitioner.h"
#include "../refinement/simple_label_propagation_refinement.h"
#include "../refinement/waver_label_propagation_refinement.h"

namespace HeiProMap {
    // The coarsest level is mapped by a greedy partition (the best of several attempts), and the coarse levels are also
    // refined by label propagation.
    class WaverMultilevelSolver : public WaverSolverBase {
        WaverLabelPropagationRefinement lp_refine;
        bool lp_active = false;

    public:
        WaverMultilevelSolver(graph_t &&g, graph_t &&t_topology, const WaverHeiProMapConfiguration &t_ac)
            : WaverSolverBase(std::move(g), std::move(t_topology), t_ac) {
            // Note: the non-LARGE_K label propagation allocates a k x k table (10^10 entries at k = 99856),
            // so it is not set up when it is not used.
            lp_active = ac.waver_label_propagation_config.enabled;
            if (lp_active) {
                lp_refine.initialize(random_engine.get_u64(), ac.waver_label_propagation_config);
            }
        }

    protected:
        void initial_partitioning() override {
            greedy_partitioning();
            finish_initial_partition();
        }

        bool uses_label_propagation() const override { return lp_active; }

        void refine_before_message_refinement() override {
            if (!lp_active) { return; }
            auto sp_local = get_time_point();
            lp_refine.refine(graphs.back(), d_oracle, bv_manager, p_manager, q_graph, block_conn);
            timings.lp_refine_ms += get_milli_seconds(sp_local, get_time_point());
        }

    private:
        // Greedy partition of the coarsest level: the best of initial_kappa attempts (the first one in natural vertex
        // order, the others shuffled).
        void greedy_partitioning() {
            auto p = get_time_point();
            HEIPROMAP_PROFILE_SCOPE("partition", "partition", "partition");

            const graph_t &g = graphs.back();

            // p_manager.n: the arrays keep the size of the finest graph. The working copy carries the PE limits the
            // greedy partitioner respects.
            PartitionManager cur_pm;
            cur_pm.initialize(p_manager.n, ac.k, g.g_weight);
            cur_pm.copy_from(p_manager);
            PartitionManager best_pm;
            best_pm.initialize(p_manager.n, ac.k, g.g_weight);
            f64 best_cost = std::numeric_limits<f64>::max();
            AlignedArray<vertex_t> vertices;
            vertices.initialize(g.n);

            for (u64 trial = 0; trial < ac.initial_kappa; ++trial) {
                const u64 trial_seed = trial == 0 ? 0 : 1 + ac.seed + trial;
                waver_greedy_partition(g, d_oracle, kImbalance, trial_seed, cur_pm, vertices, ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);
                const f64 cur_cost = get_waver_cost(g, cur_pm, d_oracle, ac.alpha_end(), ac.beta_end(), ac.alpha_link, ac.beta_link);

                if (trial == 0 || cur_cost < best_cost) {
                    best_cost = cur_cost;
                    best_pm.copy_from(cur_pm);
                }
            }

            p_manager.copy_from(best_pm);
            p_manager.recalculate_weights(g);

            timings.initial_partitioning_ms += get_milli_seconds(p, get_time_point());
        }
    };
}

#endif //HEIPROMAP_WAVERMULTILEVEL_SOLVER_H
