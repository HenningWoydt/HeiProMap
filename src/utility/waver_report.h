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

#include "../configuration/WaverHeiProMap_configuration.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../definitions.h"
#include "qap.h"
#include "utils.h"

namespace HeiProMap {

    // Accumulated wall time (ms) per solver phase.
    struct WaverTimings {
        f64 coarsening_ms = 0.0;
        f64 contraction_ms = 0.0;
        f64 initial_partitioning_ms = 0.0;
        f64 uncontraction_ms = 0.0;
        f64 rebalance_ms = 0.0;
        f64 refinement_ms = 0.0;
        f64 lp_refine_ms = 0.0;
        f64 msg_refine_ms = 0.0;
        f64 swap_refine_ms = 0.0;
        f64 long_msg_refine_ms = 0.0;
        f64 sa_refine_ms = 0.0;
        f64 misc_ms = 0.0;

        f64 total_ms() const {
            return misc_ms + coarsening_ms + contraction_ms + initial_partitioning_ms
                   + uncontraction_ms + rebalance_ms + refinement_ms;
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

        os << "Graph                 : " << ac.graph_in << std::endl;
        os << "Total time (s)        : " << duration_s << std::endl;
        os << "#Nodes                : " << g.n << std::endl;
        os << "#Edges                : " << g.m << std::endl;
        os << "k                     : " << ac.k << std::endl;
        os << "Lmax                  : " << lmax << std::endl;
        os << "Imbalance (input)     : " << ac.imbalance << std::endl;
        os << "Threads               : " << ac.threads << std::endl;
        os << "Spectral smoothing    : " << ac.spectral_smoothing.rounds << " x " << ac.spectral_smoothing.depth
           << " (> " << ac.spectral_smoothing.big_n << " vertices: " << ac.spectral_smoothing.big_rounds << " x "
           << ac.spectral_smoothing.big_depth << ")" << std::endl;
        os << "alpha_send            : " << ac.alpha_send << std::endl;
        os << "alpha_recv            : " << ac.alpha_recv << std::endl;
        os << "alpha_link            : " << ac.alpha_link << std::endl;
        os << "beta_send             : " << ac.beta_send << std::endl;
        os << "beta_recv             : " << ac.beta_recv << std::endl;
        os << "beta_link             : " << ac.beta_link << std::endl;
        os << "-----------------------" << std::endl;
        cost.print(os);
        os << std::left << std::setw(36) << "Max link load"
           << ": " << std::right << std::setw(16) << static_cast<u64>(max_link_load) << "\n";
        os << std::left << std::setw(36) << "Max PE endpoint load"
           << ": " << std::right << std::setw(16) << static_cast<u64>(endpoint.max_load) << "\n";
        os << std::left << std::setw(36) << "  busiest PE"
           << ": " << std::right << std::setw(16) << endpoint.max_pe << "\n";
        os << std::left << std::setw(36) << "  max neighbors of a PE"
           << ": " << std::right << std::setw(16) << static_cast<u64>(endpoint.max_messages) << "\n";
        // No schedule can finish before its busiest link or its busiest PE endpoint is done with its messages.
        // Model units; the local computation of a PE (rows, nonzeros) is not part of the endpoint load here.
        const bool link_bound = max_link_load >= endpoint.max_load;
        os << std::left << std::setw(36) << "Resource LB (max of the two)"
           << ": " << std::right << std::setw(16) << static_cast<u64>(std::max(max_link_load, endpoint.max_load)) << "\n";
        os << std::left << std::setw(36) << "  bound by"
           << ": " << std::right << std::setw(16) << (link_bound ? "link" : "PE endpoint") << "\n";
        os << std::left << std::setw(36) << "  link / endpoint load"
           << ": " << std::right << std::setw(16) << std::fixed << std::setprecision(3)
           << (endpoint.max_load > 0.0 ? max_link_load / endpoint.max_load : 0.0) << std::defaultfloat << "\n";
        os << "max block w          : " << max(p_manager.get_bweights()) << std::endl;
        os << "#empty partitions     : " << n_empty_partitions << std::endl;
        os << "#oload partitions     : " << n_overloaded_partitions << std::endl;
        os << "Sum oload weights     : " << sum_too_much << std::endl;
        os << "------- Time -------" << std::endl;
        os << "Coarsening            : " << t.coarsening_ms << std::endl;
        os << "Contraction           : " << t.contraction_ms << std::endl;
        os << "Init. Part.           : " << t.initial_partitioning_ms << std::endl;
        os << "Uncontraction         : " << t.uncontraction_ms << std::endl;
        os << "Rebalance             : " << t.rebalance_ms << std::endl;
        os << "Refinement            : " << t.refinement_ms << std::endl;
        os << "  Label Prop.         : " << t.lp_refine_ms << std::endl;
        os << "  Msg Refinement      : " << t.msg_refine_ms << std::endl;
        os << "  Swap                : " << t.swap_refine_ms << std::endl;
        os << "  Long Msg            : " << t.long_msg_refine_ms << std::endl;
        os << "  Sim. Annealing      : " << t.sa_refine_ms << std::endl;
        os << "Misc                  : " << t.misc_ms << std::endl;
        os << "ALL                   : " << t.total_ms() << std::endl;
    }
}

#endif //HEIPROMAP_WAVER_REPORT_H
