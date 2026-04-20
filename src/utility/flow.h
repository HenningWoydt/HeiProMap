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

#ifndef HEIPROMAP_FLOW_H
#define HEIPROMAP_FLOW_H

#include <iostream>
#include <stack>
#include <vector>

#include "../definitions.h"
#include "macros.h"
#include "random_engine.h"
#include "translation_table.h"
#include "small_translation_table.h"

namespace HeiProMap {
    class ResidualFlowNetwork {
        vertex_t n = 0;
        vertex_t source = 0;
        vertex_t target = 0;

        // CSR representation
        std::vector<vertex_t> offsets;   // size n+3 (n+2 nodes + sentinel)
        std::vector<vertex_t> neighbors; // flat neighbor array

        // temporary edge buffer used during construction
        std::vector<std::pair<vertex_t, vertex_t>> edge_buf;
        std::vector<vertex_t> write_pos;

    public:
        ResidualFlowNetwork() = default;

        void initialize(vertex_t t_n) {
            n = t_n;
            source = n;
            target = n + 1;
            edge_buf.clear();
        }

        void add_directed_edge(vertex_t u, vertex_t v, [[maybe_unused]] weight_t w) { edge_buf.emplace_back(u, v); }
        void add_edge_to_source(vertex_t u, [[maybe_unused]] weight_t w) { edge_buf.emplace_back(u, source); }
        void add_edge_to_target(vertex_t u, [[maybe_unused]] weight_t w) { edge_buf.emplace_back(u, target); }
        void add_edge_from_source(vertex_t u, [[maybe_unused]] weight_t w) { edge_buf.emplace_back(source, u); }
        void add_edge_from_target(vertex_t u, [[maybe_unused]] weight_t w) { edge_buf.emplace_back(target, u); }

        void finalize() {
            vertex_t total_nodes = n + 2;
            offsets.assign(total_nodes + 1, 0);

            for (auto &[u, v]: edge_buf) { offsets[u + 1]++; }
            for (vertex_t i = 1; i <= total_nodes; ++i) { offsets[i] += offsets[i - 1]; }

            neighbors.resize(edge_buf.size());
            // use a copy of offsets as write cursors
            write_pos.assign(offsets.begin(), offsets.begin() + total_nodes);
            for (auto &[u, v]: edge_buf) { neighbors[write_pos[u]++] = v; }
        }

        vertex_t get_n() const { return n; }
        vertex_t get_source() const { return source; }
        vertex_t get_target() const { return target; }

        // CSR access: neighbors of node u are neighbors[offsets[u]..offsets[u+1])
        vertex_t neighbor_count(vertex_t u) const { return offsets[u + 1] - offsets[u]; }
        vertex_t neighbor_at(vertex_t u, size_t i) const { return neighbors[offsets[u] + i]; }
        vertex_t offset_begin(vertex_t u) const { return offsets[u]; }
        vertex_t offset_end(vertex_t u) const { return offsets[u + 1]; }
    };

    class SCCGraph {
        vertex_t n = 0;
        vertex_t source = 0;
        vertex_t target = 0;

        std::vector<std::vector<vertex_t> > edges;
        std::vector<std::vector<vertex_t> > rev_edges;
        vertex_t scc_source = 0;
        vertex_t scc_target = 0;
        std::vector<weight_t> scc_weights;

        std::vector<vertex_t> scc_s_successors;
        std::vector<vertex_t> scc_t_predecessors;

        static constexpr vertex_t UNVISITED = std::numeric_limits<vertex_t>::max();
        std::vector<std::vector<vertex_t> > vertex_per_scc;
        std::vector<vertex_t> index;
        std::vector<vertex_t> scc_id;
        std::vector<vertex_t> S;
        std::vector<vertex_t> P;
        vertex_t idx = 0;
        vertex_t n_scc = 0;

        const ResidualFlowNetwork *temp_g = nullptr;

        std::vector<vertex_t> stack;
        std::vector<u8> best_closure;
        std::vector<u8> is_active;
        std::vector<vertex_t> in_deg;
        std::vector<vertex_t> topo_order;
        std::vector<u8> in_closure;

    public:
        SCCGraph() = default;

        template<typename GraphT>
        void initialize(const ResidualFlowNetwork &residual_flow_network,
                        const GraphT &g,
                        const TranslationTable<vertex_t> &tt) {
            temp_g = &residual_flow_network;

            n = residual_flow_network.get_n() + 2;
            source = residual_flow_network.get_source();
            target = residual_flow_network.get_target();

            index.clear();
            index.resize(n, UNVISITED);
            idx = 0;

            scc_id.clear();
            scc_id.resize(n, UNVISITED);
            n_scc = 0;

            S.clear();
            P.clear();

            vertex_per_scc.clear();

            // === Gabow's Algorithm ===
            for (vertex_t v = 0; v < n; ++v) {
                // if (index[v] == UNVISITED) { dfs(v, residual_flow_network); }
                if (index[v] == UNVISITED) { dfs_non_recursive(v, residual_flow_network); }
            }

            // clear old edges
            edges.clear();
            edges.resize(n_scc);

            rev_edges.clear();
            rev_edges.resize(n_scc);

            // build the graph and the reversed graph
            for (vertex_t u = 0; u < n; ++u) {
                vertex_t scc_u = scc_id[u];
                for (size_t j = 0; j < residual_flow_network.neighbor_count(u); ++j) {
                    vertex_t scc_v = scc_id[residual_flow_network.neighbor_at(u, j)];
                    if (scc_u == scc_v) { continue; }
                    edges[scc_u].emplace_back(scc_v);
                    rev_edges[scc_v].emplace_back(scc_u);
                }
            }

            // make each edge list unique
            for (vertex_t scc_u = 0; scc_u < n_scc; ++scc_u) {
                std::sort(edges[scc_u].begin(), edges[scc_u].end());
                edges[scc_u].erase(unique(edges[scc_u].begin(), edges[scc_u].end()), edges[scc_u].end());
            }

            for (vertex_t scc_u = 0; scc_u < n_scc; ++scc_u) {
                std::sort(rev_edges[scc_u].begin(), rev_edges[scc_u].end());
                rev_edges[scc_u].erase(unique(rev_edges[scc_u].begin(), rev_edges[scc_u].end()), rev_edges[scc_u].end());
            }

            // set special scc
            scc_source = scc_id[source];
            scc_target = scc_id[target];

            // determine the weight of each scc in the graph
            scc_weights.clear();
            scc_weights.resize(n_scc, 0);
            for (vertex_t u = 0; u < n; ++u) {
                if (u != source && u != target) {
                    scc_weights[scc_id[u]] += g.v_weights[tt.get_o(u)];
                }
            }
        }

        template<typename GraphT>
        void initialize(const ResidualFlowNetwork &residual_flow_network,
                        const GraphT &g,
                        const SmallTranslationTable<vertex_t> &tt) {
            temp_g = &residual_flow_network;

            n = residual_flow_network.get_n() + 2;
            source = residual_flow_network.get_source();
            target = residual_flow_network.get_target();

            index.clear();
            index.resize(n, UNVISITED);
            idx = 0;

            scc_id.clear();
            scc_id.resize(n, UNVISITED);
            n_scc = 0;

            S.clear();
            P.clear();

            vertex_per_scc.clear();

            // === Gabow's Algorithm ===
            for (vertex_t v = 0; v < n; ++v) {
                // if (index[v] == UNVISITED) { dfs(v, residual_flow_network); }
                if (index[v] == UNVISITED) { dfs_non_recursive(v, residual_flow_network); }
            }

            // clear old edges
            edges.clear();
            edges.resize(n_scc);

            rev_edges.clear();
            rev_edges.resize(n_scc);

            // build the graph and the reversed graph
            for (vertex_t u = 0; u < n; ++u) {
                vertex_t scc_u = scc_id[u];
                for (size_t j = 0; j < residual_flow_network.neighbor_count(u); ++j) {
                    vertex_t scc_v = scc_id[residual_flow_network.neighbor_at(u, j)];
                    if (scc_u == scc_v) { continue; }
                    edges[scc_u].emplace_back(scc_v);
                    rev_edges[scc_v].emplace_back(scc_u);
                }
            }

            // make each edge list unique
            for (vertex_t scc_u = 0; scc_u < n_scc; ++scc_u) {
                std::sort(edges[scc_u].begin(), edges[scc_u].end());
                edges[scc_u].erase(unique(edges[scc_u].begin(), edges[scc_u].end()), edges[scc_u].end());
            }

            for (vertex_t scc_u = 0; scc_u < n_scc; ++scc_u) {
                std::sort(rev_edges[scc_u].begin(), rev_edges[scc_u].end());
                rev_edges[scc_u].erase(unique(rev_edges[scc_u].begin(), rev_edges[scc_u].end()), rev_edges[scc_u].end());
            }

            // set special scc
            scc_source = scc_id[source];
            scc_target = scc_id[target];

            // determine the weight of each scc in the graph
            scc_weights.clear();
            scc_weights.resize(n_scc, 0);
            for (vertex_t u = 0; u < n; ++u) {
                if (u != source && u != target) {
                    scc_weights[scc_id[u]] += g.weight(tt.get_o(u));
                }
            }
        }

        vertex_t get_n_scc() const { return n_scc; }

        void reduce() {
            std::vector<vertex_t> curr;
            std::vector<vertex_t> next;
            std::vector<u8> seen(n_scc, 0);

            scc_s_successors.clear();
            next.push_back(scc_source);
            seen[scc_source] = 1;
            while (!next.empty()) {
                curr.swap(next);
                next.clear();

                while (!curr.empty()) {
                    vertex_t scc_u = curr.back();
                    curr.pop_back();
                    scc_s_successors.push_back(scc_u);

                    for (vertex_t scc_v: edges[scc_u]) {
                        if (seen[scc_v] == 1 || seen[scc_v] == 2) { continue; }
                        next.push_back(scc_v);
                        seen[scc_v] = 1;
                    }
                    seen[scc_u] = 2;
                }
            }

            curr.clear();
            next.clear();
            seen.clear();
            seen.resize(n_scc, 0);
            scc_t_predecessors.clear();
            next.push_back(scc_target);
            seen[scc_target] = 1;
            while (!next.empty()) {
                curr.swap(next);
                next.clear();

                while (!curr.empty()) {
                    vertex_t scc_u = curr.back();
                    curr.pop_back();
                    scc_t_predecessors.push_back(scc_u);

                    for (vertex_t scc_v: rev_edges[scc_u]) {
                        if (seen[scc_v] == 1 || seen[scc_v] == 2) { continue; }
                        next.push_back(scc_v);
                        seen[scc_v] = 1;
                    }
                    seen[scc_u] = 2;
                }
            }

            // ASSERT(no_duplicates(scc_s_successors));
            // ASSERT(no_duplicates(scc_t_predecessors));
        }

        bool find_best_closure(weight_t left_non_region_weight,
                               weight_t right_non_region_weight,
                               weight_t left_lmax,
                               weight_t right_lmax,
                               f64 avg_weight,
                               size_t repeats,
                               RandomEngine &rnd_engine,
                               std::vector<u8> &is_left) {
            bool closure_found = false;
            f64 best_cost = std::numeric_limits<f64>::max();
            best_closure.resize(n_scc);

            // determine which sccs do not have to be considered
            is_active.resize(n_scc);
            std::fill(is_active.begin(), is_active.end(), 1);
            for (vertex_t scc_u: scc_s_successors) { is_active[scc_u] = 0; }
            for (vertex_t scc_u: scc_t_predecessors) { is_active[scc_u] = 0; }

            // count active SCCs and use exact method for small DAGs
            vertex_t n_active = 0;
            for (vertex_t scc_u = 0; scc_u < n_scc; ++scc_u) {
                if (is_active[scc_u] == 1) { ++n_active; }
            }
            if (n_active <= 24) {
                return find_best_closure_exact(left_non_region_weight, right_non_region_weight, left_lmax, right_lmax, avg_weight, is_active, is_left);
            }

            // determine in degree of each vertex
            in_deg.resize(n_scc);
            std::fill(in_deg.begin(), in_deg.end(), 0);
            for (vertex_t scc_u = 0; scc_u < n_scc; ++scc_u) {
                if (is_active[scc_u] == 0) { continue; }
                for (vertex_t scc_v: edges[scc_u]) {
                    if (is_active[scc_v] == 0) { continue; }
                    in_deg[scc_v] += 1;
                }
            }
            std::vector<vertex_t> temp_in_deg(in_deg);

            for (size_t i = 0; i < repeats; ++i) {
                std::copy(temp_in_deg.begin(), temp_in_deg.end(), in_deg.begin());

                // push roots into stack
                stack.clear();
                for (vertex_t scc_u = 0; scc_u < n_scc; ++scc_u) {
                    if (is_active[scc_u] == 0) { continue; }
                    if (in_deg[scc_u] == 0) { stack.push_back(scc_u); }
                }
                std::shuffle(stack.begin(), stack.end(), rnd_engine.generator);

                // determine the random topological order
                topo_order.clear();
                while (!stack.empty()) {
                    vertex_t scc_u = stack.back();
                    stack.pop_back();

                    topo_order.push_back(scc_u);

                    for (vertex_t scc_v: edges[scc_u]) {
                        if (is_active[scc_v] == 0) { continue; }
                        in_deg[scc_v] -= 1;
                        if (in_deg[scc_v] == 0) {
                            stack.push_back(scc_v);
                            size_t idx = rnd_engine.get_u32() % stack.size();
                            std::swap(stack[idx], stack.back());
                        }
                    }
                }

                // go through the order and determine the best closure
                weight_t closure_weight = 0;
                in_closure.resize(n_scc);
                std::fill(in_closure.begin(), in_closure.end(), 0);
                for (vertex_t scc_u: scc_s_successors) {
                    in_closure[scc_u] = 1;
                    closure_weight += scc_weights[scc_u];
                }

                for (vertex_t scc_u: topo_order) {
                    if (is_active[scc_u] == 0) { continue; }

                    in_closure[scc_u] = 1;
                    closure_weight += scc_weights[scc_u];
                }

                weight_t complement_weight = 0;
                for (vertex_t scc_u: scc_t_predecessors) {
                    in_closure[scc_u] = 0;
                    complement_weight += scc_weights[scc_u];
                }

                if (left_non_region_weight + closure_weight <= left_lmax && right_non_region_weight + complement_weight <= right_lmax) {
                    f64 left = (f64)(left_non_region_weight + closure_weight);
                    f64 right = (f64)(right_non_region_weight + complement_weight);
                    f64 cost = (left - avg_weight) * (left - avg_weight) + (right - avg_weight) * (right - avg_weight);
                    if (cost < best_cost) {
                        closure_found = true;
                        best_cost = cost;
                        best_closure = in_closure;
                    }
                }

                for (vertex_t scc_u: topo_order) {
                    if (is_active[scc_u] == 0) { continue; }

                    in_closure[scc_u] = 0;
                    closure_weight -= scc_weights[scc_u];
                    complement_weight += scc_weights[scc_u];

                    if (left_non_region_weight + closure_weight <= left_lmax && right_non_region_weight + complement_weight <= right_lmax) {
                        f64 left = (f64)(left_non_region_weight + closure_weight);
                        f64 right = (f64)(right_non_region_weight + complement_weight);
                        f64 cost = (left - avg_weight) * (left - avg_weight) + (right - avg_weight) * (right - avg_weight);
                        if (cost < best_cost) {
                            closure_found = true;
                            best_cost = cost;
                            best_closure = in_closure;
                        }
                    }
                }
            }

            is_left.resize(n - 2);
            for (vertex_t u = 0; u < n - 2; ++u) { is_left[u] = best_closure[scc_id[u]]; }

            return closure_found;
        }

        bool find_best_closure_exact(weight_t left_non_region_weight,
                                     weight_t right_non_region_weight,
                                     weight_t left_lmax,
                                     weight_t right_lmax,
                                     f64 avg_weight,
                                     const std::vector<u8> &t_is_active,
                                     std::vector<u8> &is_left) {
            // Collect active SCCs and compute fixed weights
            std::vector<vertex_t> active_nodes;
            for (vertex_t scc_u = 0; scc_u < n_scc; ++scc_u) {
                if (t_is_active[scc_u] == 1) { active_nodes.push_back(scc_u); }
            }
            vertex_t m = active_nodes.size();

            weight_t s_weight = 0;
            for (vertex_t scc_u: scc_s_successors) { s_weight += scc_weights[scc_u]; }
            weight_t t_weight = 0;
            for (vertex_t scc_u: scc_t_predecessors) { t_weight += scc_weights[scc_u]; }
            weight_t total_active_weight = 0;
            for (vertex_t i = 0; i < m; ++i) { total_active_weight += scc_weights[active_nodes[i]]; }

            // Max closure weight that could satisfy left constraint
            weight_t max_w = left_lmax - left_non_region_weight - s_weight;
            if (max_w < 0) { return false; }
            // Min closure weight that could satisfy right constraint
            weight_t min_w = total_active_weight - (right_lmax - right_non_region_weight - t_weight);

            // Build local topo order of active nodes
            std::vector<vertex_t> local_id(n_scc, UNVISITED);
            for (vertex_t i = 0; i < m; ++i) { local_id[active_nodes[i]] = i; }

            std::vector<vertex_t> local_in_deg(m, 0);
            std::vector<std::vector<vertex_t>> local_preds(m);
            for (vertex_t i = 0; i < m; ++i) {
                vertex_t scc_u = active_nodes[i];
                for (vertex_t scc_v: edges[scc_u]) {
                    if (local_id[scc_v] == UNVISITED) { continue; }
                    vertex_t j = local_id[scc_v];
                    local_in_deg[j] += 1;
                    local_preds[j].push_back(i);
                }
            }

            std::vector<vertex_t> local_topo;
            local_topo.reserve(m);
            std::vector<vertex_t> q;
            for (vertex_t i = 0; i < m; ++i) {
                if (local_in_deg[i] == 0) { q.push_back(i); }
            }
            while (!q.empty()) {
                vertex_t i = q.back();
                q.pop_back();
                local_topo.push_back(i);
                vertex_t scc_u = active_nodes[i];
                for (vertex_t scc_v: edges[scc_u]) {
                    if (local_id[scc_v] == UNVISITED) { continue; }
                    vertex_t j = local_id[scc_v];
                    if (--local_in_deg[j] == 0) { q.push_back(j); }
                }
            }

            // DP: dp[w] = 1 means closure weight w is achievable
            // Process nodes in topo order. For each node, we can include it
            // only if all its active predecessors are included.
            // We track a bitmask per weight to know which nodes are "must be included"
            // -- too expensive. Instead: process in topo order, for each node either
            // include (add weight) or exclude (then all successors must be excluded).
            //
            // Better formulation: process in REVERSE topo order.
            // dp[i][w] = can we achieve active-closure-weight w considering nodes
            //            local_topo[i..m-1], where node i is the last in topo order.
            // In reverse topo: if we exclude node i, we can still include or exclude
            // later nodes (earlier in topo). If we include node i, all its predecessors
            // (earlier in reverse topo = later in topo) must also be included -- but
            // they haven't been decided yet.
            //
            // Correct approach: process in topo order with "forced inclusion" tracking.
            // Use a 1D DP array. For each node in topo order:
            //   - If node has no active predecessors excluded, it CAN be included.
            //   - We branch: include (add weight) or exclude (mark successors as blocked).
            //
            // Since tracking blocked status per-node requires exponential states,
            // we use a different formulation:
            //
            // A valid closure = any antichain-bounded downward-closed set.
            // Enumerate by processing topo order and using DP on achievable weights,
            // but only allow including node i if we can prove all preds are included.
            // Since we process in topo order, all preds of node i come before it.
            // We need to know WHICH nodes are included -- that's exponential.
            //
            // Practical exact approach for small m: bitmask DP.
            // For m <= 24, use bitmask enumeration of valid closures.

            if (m == 0) {
                // No active nodes, check if fixed assignment is feasible
                weight_t left = left_non_region_weight + s_weight;
                weight_t right = right_non_region_weight + t_weight + total_active_weight;
                if (left <= left_lmax && right <= right_lmax) {
                    best_closure.resize(n_scc);
                    std::fill(best_closure.begin(), best_closure.end(), 0);
                    for (vertex_t scc_u: scc_s_successors) { best_closure[scc_u] = 1; }
                    is_left.resize(n - 2);
                    for (vertex_t u = 0; u < n - 2; ++u) { is_left[u] = best_closure[scc_id[u]]; }
                    return true;
                }
                return false;
            }

            // Bitmask enumeration for small DAGs
            // A subset S is a valid closure iff for every node i in S,
            // all predecessors of i are also in S.
            // Precompute predecessor masks.
            std::vector<u64> pred_mask(m, 0);
            for (vertex_t i = 0; i < m; ++i) {
                for (vertex_t p: local_preds[i]) {
                    pred_mask[i] |= (1ULL << p);
                }
            }
            // Transitive closure of pred_mask: pred_mask[i] should include
            // all ancestors, not just direct predecessors.
            for (vertex_t idx = 0; idx < local_topo.size(); ++idx) {
                vertex_t i = local_topo[idx];
                for (vertex_t p: local_preds[i]) {
                    pred_mask[i] |= pred_mask[p];
                }
            }

            bool closure_found = false;
            f64 best_cost = std::numeric_limits<f64>::max();
            u64 best_mask = 0;

            u64 total = 1ULL << m;
            for (u64 mask = 0; mask < total; ++mask) {
                // Check closure property: for each included node, all ancestors included
                bool valid = true;
                for (vertex_t i = 0; i < m; ++i) {
                    if ((mask >> i) & 1) {
                        if ((mask & pred_mask[i]) != pred_mask[i]) {
                            valid = false;
                            break;
                        }
                    }
                }
                if (!valid) { continue; }

                weight_t w = 0;
                for (vertex_t i = 0; i < m; ++i) {
                    if ((mask >> i) & 1) { w += scc_weights[active_nodes[i]]; }
                }

                weight_t left_w = left_non_region_weight + s_weight + w;
                weight_t right_w = right_non_region_weight + t_weight + (total_active_weight - w);
                if (left_w <= left_lmax && right_w <= right_lmax) {
                    f64 dl = (f64)left_w - avg_weight;
                    f64 dr = (f64)right_w - avg_weight;
                    f64 cost = dl * dl + dr * dr;
                    if (cost < best_cost) {
                        closure_found = true;
                        best_cost = cost;
                        best_mask = mask;
                    }
                }
            }

            if (!closure_found) { return false; }

            best_closure.resize(n_scc);
            std::fill(best_closure.begin(), best_closure.end(), 0);
            for (vertex_t scc_u: scc_s_successors) { best_closure[scc_u] = 1; }
            for (vertex_t i = 0; i < m; ++i) {
                if ((best_mask >> i) & 1) { best_closure[active_nodes[i]] = 1; }
            }

            is_left.resize(n - 2);
            for (vertex_t u = 0; u < n - 2; ++u) { is_left[u] = best_closure[scc_id[u]]; }
            return closure_found;
        }

    private:
        void dfs(vertex_t v, const ResidualFlowNetwork &residual_flow_network) {
            index[v] = idx++;
            S.push_back(v);
            P.push_back(v);

            for (size_t j = 0; j < residual_flow_network.neighbor_count(v); ++j) {
                vertex_t u = residual_flow_network.neighbor_at(v, j);
                if (index[u] == UNVISITED) {
                    dfs(u, residual_flow_network);
                } else if (scc_id[u] == UNVISITED) {
                    while (!P.empty() && index[P.back()] > index[u]) {
                        P.pop_back();
                    }
                }
            }

            if (!P.empty() && P.back() == v) {
                vertex_per_scc.emplace_back(); // Start a new SCC
                vertex_t w;
                do {
                    w = S.back();
                    S.pop_back();
                    scc_id[w] = n_scc;
                    vertex_per_scc.back().push_back(w); // Add to current SCC
                } while (w != v);
                P.pop_back();
                ++n_scc;
            }
        }

        void dfs_non_recursive(vertex_t start, const ResidualFlowNetwork &residual_flow_network) {
            std::stack<std::pair<vertex_t, size_t> > dfs_stack; // (node, next neighbor index)

            dfs_stack.push({start, 0});
            index[start] = idx++;
            S.push_back(start);
            P.push_back(start);

            while (!dfs_stack.empty()) {
                auto &[v, i] = dfs_stack.top(); // current node and index into its neighbors

                if (i < residual_flow_network.neighbor_count(v)) {
                    vertex_t u = residual_flow_network.neighbor_at(v, i++);
                    if (index[u] == UNVISITED) {
                        index[u] = idx++;
                        S.push_back(u);
                        P.push_back(u);
                        dfs_stack.push({u, 0}); // dive deeper
                    } else if (scc_id[u] == UNVISITED) {
                        while (!P.empty() && index[P.back()] > index[u]) {
                            P.pop_back();
                        }
                    }
                } else {
                    // done processing all neighbors
                    if (!P.empty() && P.back() == v) {
                        vertex_per_scc.emplace_back();
                        vertex_t w;
                        do {
                            w = S.back();
                            S.pop_back();
                            scc_id[w] = n_scc;
                            vertex_per_scc.back().push_back(w);
                        } while (w != v);
                        P.pop_back();
                        ++n_scc;
                    }
                    dfs_stack.pop();
                }
            }
        }
    };
}

#endif //HEIPROMAP_FLOW_H
