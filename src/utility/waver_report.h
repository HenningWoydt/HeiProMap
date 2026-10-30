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

#ifndef HEIPROMAP_WAVER_REPORT_H
#define HEIPROMAP_WAVER_REPORT_H

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "../configuration/WaverHeiProMap_configuration.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../definitions.h"
#include "qap.h"
#include "utils.h"

namespace HeiProMap {

    // Accumulated wall time (ms) per solver phase.
    struct WaverTimings {
        f64 prepartition_ms = 0.0;
        f64 spectral_ms = 0.0; // spectral solver: prolongation of the embedding and the placement of the finest level
        f64 coarsening_ms = 0.0;
        f64 contraction_ms = 0.0;
        f64 initial_partitioning_ms = 0.0;
        f64 uncontraction_ms = 0.0;
        f64 rebalance_ms = 0.0;
        f64 refinement_ms = 0.0;
        f64 msg_refine_ms = 0.0;
        f64 swap_refine_ms = 0.0;
        f64 long_msg_refine_ms = 0.0;
        f64 merge_refine_ms = 0.0;
        f64 busy_refine_ms = 0.0;
        f64 crowd_refine_ms = 0.0;
        f64 shift_refine_ms = 0.0;
        f64 vcycle_ms = 0.0;
        f64 smash_ms = 0.0;
        f64 cycle_refine_ms = 0.0;
        f64 sa_refine_ms = 0.0;
        f64 corner_refine_ms = 0.0;
        f64 misc_ms = 0.0;

        f64 total_ms() const {
            return misc_ms + prepartition_ms + spectral_ms + coarsening_ms + contraction_ms + initial_partitioning_ms
                   + uncontraction_ms + rebalance_ms + refinement_ms + corner_refine_ms;
        }
    };

    // Final report of a Waver-HeiProMap run: configuration, cost breakdown, balance and phase times.
    inline void print_waver_summary(std::ostream &os,
                                    const WaverHeiProMapConfiguration &ac,
                                    const CSRGraph &g,
                                    const f64 duration_s,
                                    const weight_t lmax,
                                    const MessageCostBreakdown &cost,
                                    const f64 max_link_load,
                                    const EndpointLoad &endpoint,
                                    const PartitionManager &p_manager,
                                    const WaverTimings &t) {
        size_t n_empty_partitions = 0;
        size_t n_overloaded_partitions = 0;
        weight_t sum_too_much = 0;
        for (partition_t id = 0; id < ac.k; ++id) {
            n_empty_partitions += p_manager.get_bweight(id) == 0;
            n_overloaded_partitions += p_manager.get_bweight(id) > lmax;
            sum_too_much += std::max((weight_t) 0, p_manager.get_bweight(id) - lmax);
        }

        constexpr int label_w = 36;
        constexpr int value_w = 16;
        const auto row = [&](const std::string &label, const auto &value) {
            os << std::left << std::setw(label_w) << label << ": " << std::right << std::setw(value_w) << value << "\n";
        };
        const auto fixed = [](const f64 v, const int precision) {
            std::ostringstream s;
            s << std::fixed << std::setprecision(precision) << v;
            return s.str();
        };
        const auto time_row = [&](const std::string &label, const f64 ms) {
            row(label, fixed(ms, 1));
        };

        row("Graph", ac.graph_in);
        row("Total time (s)", fixed(duration_s, 3));
        row("#Nodes", g.n);
        row("#Edges", g.m);
        row("k", ac.k);
        row("Lmax", lmax);
        row("Imbalance (input)", ac.imbalance);
        row("Threads", ac.threads);
        row("Spectral smoothing", std::to_string(ac.spectral_smoothing.rounds) + " x "
                                      + std::to_string(ac.spectral_smoothing.depth));
        row("  > " + std::to_string(ac.spectral_smoothing.big_n) + " vertices",
            std::to_string(ac.spectral_smoothing.big_rounds) + " x " + std::to_string(ac.spectral_smoothing.big_depth));
        row("alpha_send", ac.alpha_send);
        row("alpha_recv", ac.alpha_recv);
        row("alpha_link", ac.alpha_link);
        row("beta_send", ac.beta_send);
        row("beta_recv", ac.beta_recv);
        row("beta_link", ac.beta_link);
        os << std::string(label_w + 2 + value_w, '-') << "\n";
        cost.print(os);
        row("Max link load", static_cast<u64>(max_link_load));
        row("Max PE endpoint load", static_cast<u64>(endpoint.max_load));
        row("  busiest PE", endpoint.max_pe);
        row("  max neighbors of a PE", static_cast<u64>(endpoint.max_messages));
        // No schedule can finish before its busiest link or its busiest PE endpoint is done with its messages.
        // Model units; the local computation of a PE (rows, nonzeros) is not part of the endpoint load here.
        const bool link_bound = max_link_load >= endpoint.max_load;
        row("Resource LB (max of the two)", static_cast<u64>(std::max(max_link_load, endpoint.max_load)));
        row("  bound by", link_bound ? "link" : "PE endpoint");
        row("  link / endpoint load", fixed(endpoint.max_load > 0.0 ? max_link_load / endpoint.max_load : 0.0, 3));
        row("Max block weight", max(p_manager.get_bweights()));
        row("#empty partitions", n_empty_partitions);
        row("#overloaded partitions", n_overloaded_partitions);
        row("Sum overload weights", sum_too_much);
        os << "------- Time (ms) -------" << "\n";
        time_row("Prepartition", t.prepartition_ms);
        time_row("Spectral (solver)", t.spectral_ms);
        time_row("Coarsening", t.coarsening_ms);
        time_row("Contraction", t.contraction_ms);
        time_row("Init. Part.", t.initial_partitioning_ms);
        time_row("Uncontraction", t.uncontraction_ms);
        time_row("Rebalance", t.rebalance_ms);
        time_row("Refinement", t.refinement_ms);
        time_row("  Msg Refinement", t.msg_refine_ms);
        time_row("  Swap", t.swap_refine_ms);
        time_row("  Long Msg", t.long_msg_refine_ms);
        time_row("  Merge", t.merge_refine_ms);
        time_row("  Cycle", t.cycle_refine_ms);
        time_row("  Busy", t.busy_refine_ms);
        time_row("  Crowd", t.crowd_refine_ms);
        time_row("  Shift", t.shift_refine_ms);
        time_row("  V-cycle (incl.)", t.vcycle_ms);
        time_row("  Smash (incl.)", t.smash_ms);
        time_row("  Sim. Annealing", t.sa_refine_ms);
        time_row("  Corner (final)", t.corner_refine_ms);
        time_row("Misc (setup, rebuilds, other)", t.misc_ms);
        time_row("ALL", t.total_ms());
        os << std::flush;
    }
}

#endif //HEIPROMAP_WAVER_REPORT_H
