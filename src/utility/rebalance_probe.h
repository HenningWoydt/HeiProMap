#ifndef HEIPROMAP_REBALANCE_PROBE_H
#define HEIPROMAP_REBALANCE_PROBE_H

#include <algorithm>
#include <iostream>
#include <vector>

#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../definitions.h"
#include "spectral_tracker.h"

namespace HeiProMap {

    // Diagnostics of one call of the rebalancer (--rebalance-verbose): the weight above the limits before and after, how
    // many vertices it moved and how far (distance on the PE grid between the old and the new PE), and the message cost.
    class RebalanceProbe {
        std::vector<partition_t> before;
        f64 overload_before = 0.0, cost_before = 0.0;
        bool cost_known = false;

        static f64 overload(const CSRGraph &g, const PartitionManager &pm, const partition_t k) {
            std::vector<f64> w(k, 0.0);
            for (vertex_t u = 0; u < g.n; ++u) { w[pm.partition[u]] += static_cast<f64>(g.v_weight(u)); }
            f64 excess = 0.0;
            for (partition_t p = 0; p < k; ++p) { excess += std::max(0.0, w[p] - static_cast<f64>(pm.get_lmax(p))); }
            return excess;
        }

    public:
        void begin(const CSRGraph &g, const PartitionManager &pm) {
            before.assign(pm.partition.begin(), pm.partition.begin() + g.n);
            overload_before = overload(g, pm, pm.get_k());
        }

        template<typename OracleT>
        void end(const CSRGraph &g, const PartitionManager &pm, const OracleT &oracle, const partition_t k, const u64 level,
                 const f64 alpha_end, const f64 beta_end, const f64 alpha_link, const f64 beta_link) {
            u64 moved = 0;
            f64 distance = 0.0, max_distance = 0.0;
            std::vector<partition_t> after(pm.partition.begin(), pm.partition.begin() + g.n);
            for (vertex_t u = 0; u < g.n; ++u) {
                if (after[u] == before[u]) { continue; }
                ++moved;
                const f64 d = static_cast<f64>(oracle.get(before[u], after[u]));
                distance += d;
                max_distance = std::max(max_distance, d);
            }
            const f64 c_before = SpectralTracker::message_cost(g, before, oracle, alpha_end, beta_end, alpha_link, beta_link);
            const f64 c_after = SpectralTracker::message_cost(g, after, oracle, alpha_end, beta_end, alpha_link, beta_link);
            std::cout << "[rebalance] graph n=" << g.n << " k=" << k << " levels left " << level << ": overload "
                      << overload_before << " -> " << overload(g, pm, k) << ", moved " << moved << " vertices (mean distance "
                      << (moved ? distance / static_cast<f64>(moved) : 0.0) << ", max " << max_distance << "), cost "
                      << c_before << " -> " << c_after << " (" << (c_before > 0 ? 100.0 * (c_after - c_before) / c_before : 0.0)
                      << "%)" << std::endl;
        }
    };

}

#endif //HEIPROMAP_REBALANCE_PROBE_H
