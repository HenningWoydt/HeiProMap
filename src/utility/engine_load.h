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

#ifndef HEIPROMAP_ENGINE_LOAD_H
#define HEIPROMAP_ENGINE_LOAD_H

#include <algorithm>
#include <vector>

#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../definitions.h"
#include "profiler.h"

namespace HeiProMap {

    // The busy cycles of the compute engine of one PE, from the per-message and per-value costs of the evaluator's
    // timed model (wse-spmv costs.json and its cost breakdown), in cycles. The schedule cannot finish before the
    // busiest engine is done, so the maximum over the PEs is the resource lower bound of the makespan.
    //
    // A PE p runs, in this model,
    //   startup:  startup_fixed + row * rows(p) + local_nnz * nonzeros with both endpoints in p
    //   a send    to every PE q that owns a neighbor of one of its rows: send_msg + send_value * |S_pq|, where
    //             S_pq are the distinct rows of p with a neighbor in q (every value is sent once, however many
    //             of the rows of q use it)
    //   a receive from every such PE q: recv_msg + recv_value * |S_qp|
    // The constants are fits of the cost breakdown (total per category over all messages), not exact.
    struct EngineCosts {
        f64 startup_fixed = 179.0; // matvec setup and launch
        f64 row = 20.0;            // row_even / row_odd
        f64 local_nnz = 4.5;       // nnz_pair for two nonzeros
        f64 send_msg = 383.0;      // send_lead + send_assemble + send_setup + dispatch + guard
        f64 send_value = 5.1;      // send_pair / 2 + wire_value
        f64 recv_msg = 272.0;      // recv_prepare + send_term + scatter setup + dispatch + guard
        f64 recv_value = 15.3;     // wire_value + scatter per value
    };

    struct EngineLoadStats {
        f64 max_load = 0.0;
        f64 mean_load = 0.0;
        partition_t busiest_pe = 0;
        // composition of the busiest PE
        u64 rows = 0;
        u64 local_nnz = 0;
        u64 sends = 0;
        u64 receives = 0;
        u64 sent_values = 0;
        u64 received_values = 0;
        f64 startup = 0.0;
    };

    // Engine load of every PE of the partition of g into k blocks.
    inline EngineLoadStats get_engine_load(const CSRGraph &g, const PartitionManager &p_manager, const partition_t k,
                                           const EngineCosts &c = EngineCosts()) {
        HEIPROMAP_PROFILE_SCOPE("misc", "misc", "get_engine_load");

        std::vector<u64> rows(k, 0), local_nnz(k, 0);
        std::vector<u64> sends(k, 0), receives(k, 0), sent_values(k, 0), received_values(k, 0);

        // One key (A << 32 | B) per row u of block A and block B != A holding a neighbor of u.
        std::vector<u64> keys;
        keys.reserve(g.m);
        std::vector<partition_t> blocks;
        for (vertex_t u = 0; u < g.n; ++u) {
            const partition_t a = p_manager[u];
            rows[a] += static_cast<u64>(g.v_weight(u));
            ++local_nnz[a]; // diagonal
            blocks.clear();
            for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                const partition_t b = p_manager[g.edges_v[j]];
                if (b == a) {
                    ++local_nnz[a];
                } else {
                    blocks.push_back(b);
                }
            }
            std::sort(blocks.begin(), blocks.end());
            blocks.erase(std::unique(blocks.begin(), blocks.end()), blocks.end());
            for (const partition_t b : blocks) { keys.push_back((static_cast<u64>(a) << 32) | b); }
        }
        std::sort(keys.begin(), keys.end());
        for (size_t i = 0; i < keys.size();) {
            size_t e = i;
            while (e < keys.size() && keys[e] == keys[i]) { ++e; }
            const auto a = static_cast<partition_t>(keys[i] >> 32);
            const auto b = static_cast<partition_t>(keys[i] & 0xFFFFFFFFull);
            const u64 values = e - i; // |S_ab|
            ++sends[a];
            sent_values[a] += values;
            ++receives[b];
            received_values[b] += values;
            i = e;
        }

        EngineLoadStats s;
        f64 sum = 0.0;
        for (partition_t p = 0; p < k; ++p) {
            const f64 startup = c.startup_fixed + c.row * static_cast<f64>(rows[p])
                                + c.local_nnz * static_cast<f64>(local_nnz[p]);
            const f64 load = startup
                             + c.send_msg * static_cast<f64>(sends[p]) + c.send_value * static_cast<f64>(sent_values[p])
                             + c.recv_msg * static_cast<f64>(receives[p]) + c.recv_value * static_cast<f64>(received_values[p]);
            sum += load;
            if (load > s.max_load) {
                s.max_load = load;
                s.busiest_pe = p;
                s.rows = rows[p];
                s.local_nnz = local_nnz[p];
                s.sends = sends[p];
                s.receives = receives[p];
                s.sent_values = sent_values[p];
                s.received_values = received_values[p];
                s.startup = startup;
            }
        }
        s.mean_load = k > 0 ? sum / static_cast<f64>(k) : 0.0;
        return s;
    }
}

#endif //HEIPROMAP_ENGINE_LOAD_H
