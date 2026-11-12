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

#ifndef HEIPROMAP_GGG_MULTISECTION_H
#define HEIPROMAP_GGG_MULTISECTION_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>
#include <random>
#include <utility>
#include <vector>

#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../definitions.h"
#include "flow_bisection.h"

namespace HeiProMap {
    /**
     * Very simple initial partitioner: multisection along the hierarchy (highest level first), every level is
     * split by recursive bisection and every bisection is greedy graph growing:
     *   - the seed is the vertex farthest (in hops) from a random vertex,
     *   - the region grows by always adding the vertex with the largest gain (edge weight into the region minus
     *     edge weight to the rest) until it reaches its target weight (never above its weight cap),
     *   - vertices that are not reachable (other components) are picked up by a new seed.
     * At the end the overloaded blocks are repaired (cheapest cut increase first).
     */
    class GggMultisectionPartitioner {
    private:
        std::vector<int> loc; // vertex -> index in the current vertex range, -1 = not in it
        std::mt19937_64 rng;

        // Splits verts[0, cnt) into a first part (returned size) of weight in [min1, max1] (as far as possible) and the rest.
        size_t bisect(const CSRGraph &g, vertex_t *verts, const size_t cnt, const weight_t target1, const weight_t max1) {
            for (size_t i = 0; i < cnt; ++i) { loc[verts[i]] = (int) i; }

            std::vector<s64> wdeg(cnt, 0), conn(cnt, 0);
            for (size_t i = 0; i < cnt; ++i) {
                const vertex_t u = verts[i];
                for (size_t e = g.neighborhoods[u]; e < g.neighborhoods[u + 1]; ++e) {
                    if (loc[g.edges_v[e]] >= 0) { wdeg[i] += g.edge_w(e); }
                }
            }

            // seed: the vertex farthest from a random one (last vertex of a BFS)
            size_t seed_idx = (size_t) (rng() % cnt);
            {
                std::vector<char> seen(cnt, 0);
                std::vector<size_t> queue;
                queue.reserve(cnt);
                queue.push_back(seed_idx);
                seen[seed_idx] = 1;
                for (size_t h = 0; h < queue.size(); ++h) {
                    const vertex_t u = verts[queue[h]];
                    for (size_t e = g.neighborhoods[u]; e < g.neighborhoods[u + 1]; ++e) {
                        const int j = loc[g.edges_v[e]];
                        if (j >= 0 && !seen[j]) { seen[j] = 1; queue.push_back((size_t) j); }
                    }
                }
                seed_idx = queue.back();
            }

            std::vector<char> in1(cnt, 0);
            using Entry = std::pair<s64, size_t>; // (gain, index)
            std::priority_queue<Entry> pq;
            auto gain = [&](const size_t i) -> s64 { return 2 * conn[i] - wdeg[i]; };
            pq.push({gain(seed_idx), seed_idx});
            size_t scan = 0, n1 = 0;
            weight_t w1 = 0;
            while (w1 < target1) {
                if (pq.empty()) {
                    // another component: start a new region part
                    while (scan < cnt && in1[scan]) { ++scan; }
                    if (scan >= cnt) { break; }
                    pq.push({gain(scan), scan});
                }
                const auto [key, i] = pq.top();
                pq.pop();
                if (in1[i] || key != gain(i)) { continue; } // stale entry
                const vertex_t u = verts[i];
                const weight_t wu = g.v_weight(u);
                if (w1 + wu > max1) {
                    // does not fit: leave it (when it was the only candidate, the scan above skips to another vertex)
                    if (pq.empty()) {
                        in1[i] = 2; // blocked, never taken
                        while (scan < cnt && in1[scan]) { ++scan; }
                        if (scan >= cnt) { break; }
                        pq.push({gain(scan), scan});
                    }
                    continue;
                }
                in1[i] = 1;
                w1 += wu;
                ++n1;
                for (size_t e = g.neighborhoods[u]; e < g.neighborhoods[u + 1]; ++e) {
                    const int j = loc[g.edges_v[e]];
                    if (j < 0 || in1[j]) { continue; }
                    conn[j] += g.edge_w(e);
                    pq.push({gain((size_t) j), (size_t) j});
                }
            }

            // stable partition: region first
            std::vector<vertex_t> rest;
            rest.reserve(cnt - n1);
            size_t pos = 0;
            for (size_t i = 0; i < cnt; ++i) {
                if (in1[i] == 1) { verts[pos++] = verts[i]; } else { rest.push_back(verts[i]); }
            }
            for (size_t i = 0; i < rest.size(); ++i) { verts[pos + i] = rest[i]; }
            for (size_t i = 0; i < cnt; ++i) { loc[verts[i]] = -1; }
            return n1;
        }

        // recursive bisection of order[b, e) into kk parts of weight <= part_max; part i = order[bounds[i], bounds[i + 1])
        void split(const CSRGraph &g, vertex_t *order, const size_t b, const size_t e, const partition_t kk,
                   const weight_t part_max, std::vector<size_t> &bounds) {
            bounds.assign(kk + 1, e);
            bounds[0] = b;
            struct Task { size_t b, e; partition_t lo, hi; };
            std::vector<Task> stack = {{b, e, 0, kk}};
            while (!stack.empty()) {
                const Task t = stack.back();
                stack.pop_back();
                const partition_t K = t.hi - t.lo;
                if (K == 1 || t.e - t.b < 2) {
                    bounds[t.lo] = t.b; // everything to the first block of the range, the others stay empty
                    for (partition_t i = t.lo + 1; i < t.hi; ++i) { bounds[i] = t.e; }
                    continue;
                }
                const partition_t mid = t.lo + (K + 1) / 2;
                weight_t total = 0;
                for (size_t x = t.b; x < t.e; ++x) { total += g.v_weight(order[x]); }
                const weight_t max1 = (weight_t) (mid - t.lo) * part_max;
                const weight_t max2 = (weight_t) (t.hi - mid) * part_max;
                weight_t target1 = (weight_t) std::llround((f64) total * (f64) (mid - t.lo) / (f64) K);
                target1 = std::max(target1, total - max2);
                target1 = std::min(target1, max1);
                const size_t n1 = bisect(g, order + t.b, t.e - t.b, target1, max1);
                bounds[t.lo] = t.b;
                bounds[mid] = t.b + n1;
                stack.push_back({t.b, t.b + n1, t.lo, mid});
                stack.push_back({t.b + n1, t.e, mid, t.hi});
            }
        }

    public:
        static void partition(const CSRGraph &g, PartitionManager &p_manager, const std::vector<partition_t> &hierarchy,
                              const f64 imbalance, const u64 seed) {
            GggMultisectionPartitioner p;
            p.rng.seed(seed * 0x9E3779B97F4A7C15ULL + 12345);
            p.loc.assign(g.n, -1);

            const size_t l = hierarchy.size();
            std::vector<partition_t> stride(l); // PEs below one part of level j
            partition_t k = 1;
            for (size_t j = 0; j < l; ++j) { stride[j] = k; k *= hierarchy[j]; }
            const weight_t lmax = (weight_t) std::ceil((1.0 + imbalance) * (f64) g.g_weight / (f64) k);

            std::vector<vertex_t> order(g.n);
            for (vertex_t u = 0; u < g.n; ++u) { order[u] = u; }

            struct Item { size_t b, e, level; partition_t offset; };
            std::vector<Item> stack = {{0, (size_t) g.n, l - 1, 0}};
            std::vector<size_t> bounds;
            while (!stack.empty()) {
                const Item item = stack.back();
                stack.pop_back();
                const partition_t kk = hierarchy[item.level];

                // weight cap of a part: its PEs times lmax, tightened so that the lower levels share the slack
                weight_t part_max = (weight_t) stride[item.level] * lmax;
                weight_t item_w = 0;
                for (size_t x = item.b; x < item.e; ++x) { item_w += g.v_weight(order[x]); }
                const f64 local_k = (f64) stride[item.level] * (f64) kk;
                f64 imb = ((1.0 + imbalance) * local_k * (f64) g.g_weight) / ((f64) k * (f64) std::max<weight_t>(1, item_w));
                imb = std::max(0.0, std::pow(imb, 1.0 / (f64) (item.level + 1)) - 1.0);
                part_max = std::min(part_max, (weight_t) std::ceil((1.0 + imb) * (f64) item_w / (f64) kk));

                p.split(g, order.data(), item.b, item.e, kk, part_max, bounds);
                for (partition_t i = 0; i < kk; ++i) {
                    const size_t pb = bounds[i], pe = (i + 1 < kk) ? bounds[i + 1] : item.e;
                    if (item.level == 0) {
                        for (size_t x = pb; x < pe; ++x) { p_manager.set(order[x], g.v_weight(order[x]), item.offset + i); }
                    } else {
                        stack.push_back({pb, pe, item.level - 1, item.offset + i * stride[item.level]});
                    }
                }
            }
            FlowBisectionPartitioner::repair_global(g, p_manager, k, lmax);
        }
    };
} // namespace HeiProMap

#endif // HEIPROMAP_GGG_MULTISECTION_H
