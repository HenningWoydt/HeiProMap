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

#ifndef HEIPROMAP_WAVERSPECTRAL_SOLVER_H
#define HEIPROMAP_WAVERSPECTRAL_SOLVER_H

#include "WaverBase_solver.h"

namespace HeiProMap {
    // The coarsest level is placed by a spectral (Laplacian eigenvector, recursive coordinate bisection) embedding. The
    // embedding follows the graph hierarchy while uncoarsening, and on the finest level it competes with the result.
    class WaverSpectralSolver : public WaverSolverBase {
        SpectralTracker spectral; // the spectral solution, one level per graph of the hierarchy

    public:
        WaverSpectralSolver(graph_t &&g, graph_t &&t_topology, const WaverHeiProMapConfiguration &t_ac)
            : WaverSolverBase(std::move(g), std::move(t_topology), t_ac) {
            SpectralSmoothing smoothing = ac.spectral_smoothing;
            smoothing.threads = static_cast<int>(ac.threads);
            spectral.set_smoothing(smoothing);
        }

    protected:
        void initial_partitioning() override {
            const auto t_start = get_time_point();
            spectral.init(graphs.back(), ac.seed);
            load_partition(spectral_placement());
            timings.spectral_ms += get_milli_seconds(t_start, get_time_point());
            finish_initial_partition();
        }

        void on_graph_uncoarsen(const Mapping &mapping, const graph_t &finer) override {
            const auto t_start = get_time_point();
            spectral.prolong(mapping, finer);
            timings.spectral_ms += get_milli_seconds(t_start, get_time_point());
        }

        void after_uncoarsening() override { compare_with_spectral_placement(); }

    private:
        std::vector<partition_t> spectral_placement() {
            std::vector<partition_t> part;
            spectral.place(graphs.back(), d_oracle.current_topology(), d_oracle.current_W(), part);
            return part;
        }

        // Once on the finest level (the graph with one vertex per PE): the spectral placement
        // of that level is an independent candidate that ignores the partition found while uncoarsening. If it is
        // cheaper than that partition, it is adopted (rebalanced and refined); if the refined placement is still not
        // cheaper, the partition from uncoarsening is kept. The spectral tracker itself is never overwritten.
        void compare_with_spectral_placement() {
            const graph_t &g = graphs.back();
            if (spectral.size() != g.n) { return; }
            HEIPROMAP_PROFILE_SCOPE("misc", "WaverSolver", "compare_with_spectral_placement");

            const auto t_start = get_time_point();
            const std::vector<partition_t> spectral_part = spectral_placement();
            timings.spectral_ms += get_milli_seconds(t_start, get_time_point());
            const f64 multilevel_cost = placement_cost(p_manager.partition);
            if (placement_cost(spectral_part) >= multilevel_cost) { return; }

            const std::vector<partition_t> multilevel_part = snapshot();
            load_partition(spectral_part);
            recompute_datastructures();
            legalize(mappings.size());
            refinement();
            if (current_message_cost() >= multilevel_cost) {
                load_partition(multilevel_part);
                recompute_datastructures();
            }
        }
    };
}

#endif //HEIPROMAP_WAVERSPECTRAL_SOLVER_H
