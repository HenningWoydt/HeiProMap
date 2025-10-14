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
#include <queue>
#include <stack>
#include <unordered_set>
#include <vector>

#include "definitions.h"
#include "macros.h"
#include "random_engine.h"
#include "translation_table.h"
#include "small_translation_table.h"
#include "../../extern/maxflow/graph.h"
#include "../serial/serial_definitions_1.h"

namespace HeiProMap {

    class ResidualFlowNetwork {
        vertex_t                         n      = 0;
        vertex_t                         source = 0;
        vertex_t                         target = 0;
        std::vector<std::vector<EdgeVW>> edges;

    public:
        ResidualFlowNetwork() = default;

        void initialize(vertex_t t_n) {
            n      = t_n;
            source = n;
            target = n + 1;

            if (edges.size() != n + 2) {
                edges.resize(n + 2);
            }
            for (auto &vec: edges) {
                vec.clear();
            }
        }

        void add_directed_edge(vertex_t u, vertex_t v, weight_t w) {
            // for (auto &e: edges[u]) { ASSERT(e.v != v); }

            edges[u].emplace_back(v, w);
        }

        void add_edge_to_source(vertex_t u, weight_t w) {
            // for (auto &e: edges[u]) { ASSERT(e.v != source); }

            edges[u].emplace_back(source, w);
        }

        void add_edge_to_target(vertex_t u, weight_t w) {
            // for (auto &e: edges[u]) { ASSERT(e.v != target); }

            edges[u].emplace_back(target, w);
        }

        void add_edge_from_source(vertex_t u, weight_t w) {
            // for (auto &e: edges[source]) { ASSERT(e.v != u); }

            edges[source].emplace_back(u, w);
        }

        void add_edge_from_target(vertex_t u, weight_t w) {
            // for (auto &e: edges[target]) { ASSERT(e.v != u); }

            edges[target].emplace_back(u, w);
        }

        vertex_t get_n() const { return n; }

        vertex_t get_source() const { return source; }

        vertex_t get_target() const { return target; }

        const std::vector<EdgeVW> &operator[](vertex_t i) const { return edges[i]; }
    };

    class SCCGraph {
        vertex_t n      = 0;
        vertex_t source = 0;
        vertex_t target = 0;

        std::vector<std::vector<vertex_t>> edges;
        std::vector<std::vector<vertex_t>> rev_edges;
        vertex_t                           scc_source = 0;
        vertex_t                           scc_target = 0;
        std::vector<weight_t>              scc_weights;

        std::vector<vertex_t> scc_s_successors;
        std::vector<vertex_t> scc_t_predecessors;

        static constexpr vertex_t          UNVISITED = std::numeric_limits<vertex_t>::max();
        std::vector<std::vector<vertex_t>> vertex_per_scc;
        std::vector<vertex_t>              index;
        std::vector<vertex_t>              scc_id;
        std::vector<vertex_t>              S;
        std::vector<vertex_t>              P;
        vertex_t                           idx       = 0;
        vertex_t                           n_scc     = 0;

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

            n      = residual_flow_network.get_n() + 2;
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
                for (const auto [v, w]: residual_flow_network[u]) {
                    vertex_t scc_v = scc_id[v];
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

        template<typename GraphT>
        void initialize(const ResidualFlowNetwork &residual_flow_network,
                        const GraphT &g,
                        const SmallTranslationTable<vertex_t> &tt) {
            temp_g = &residual_flow_network;

            n      = residual_flow_network.get_n() + 2;
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
                for (const auto [v, w]: residual_flow_network[u]) {
                    vertex_t scc_v = scc_id[v];
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
            std::vector<u8>       seen(n_scc, 0);

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
                               size_t repeats,
                               RandomEngine &rnd_engine,
                               std::vector<u8> &is_left) {
            bool            closure_found = false;
            weight_t        best_diff     = -1;
            best_closure.resize(n_scc);

            // determine which sccs do not have to be considered
            is_active.resize(n_scc);
            std::fill(is_active.begin(), is_active.end(), 1);
            for (vertex_t   scc_u: scc_s_successors) { is_active[scc_u] = 0; }
            for (vertex_t   scc_u: scc_t_predecessors) { is_active[scc_u] = 0; }

            // determine in degree of each vertex
            in_deg.resize(n_scc);
            std::fill(in_deg.begin(), in_deg.end(), 0);
            for (vertex_t         scc_u = 0; scc_u < n_scc; ++scc_u) {
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
                for (vertex_t         scc_u = 0; scc_u < n_scc; ++scc_u) {
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
                weight_t        closure_weight = 0;
                in_closure.resize(n_scc);
                std::fill(in_closure.begin(), in_closure.end(), 0);
                for (vertex_t   scc_u: scc_s_successors) {
                    in_closure[scc_u] = 1;
                    closure_weight += scc_weights[scc_u];
                }

                for (vertex_t scc_u: topo_order) {
                    if (is_active[scc_u] == 0) { continue; }

                    in_closure[scc_u] = 1;
                    closure_weight += scc_weights[scc_u];
                }

                weight_t      complement_weight = 0;
                for (vertex_t scc_u: scc_t_predecessors) {
                    in_closure[scc_u] = 0;
                    complement_weight += scc_weights[scc_u];
                }

                if (left_non_region_weight + closure_weight <= left_lmax && right_non_region_weight + complement_weight <= right_lmax) {
                    weight_t left  = left_non_region_weight + closure_weight;
                    weight_t right = right_non_region_weight + complement_weight;
                    weight_t diff  = std::min(left_lmax - left, right_lmax - right);
                    if (diff > best_diff) {
                        closure_found = true;
                        best_diff     = diff;
                        best_closure  = in_closure;
                    }
                }

                for (vertex_t scc_u: topo_order) {
                    if (is_active[scc_u] == 0) { continue; }

                    in_closure[scc_u] = 0;
                    closure_weight -= scc_weights[scc_u];
                    complement_weight += scc_weights[scc_u];

                    if (left_non_region_weight + closure_weight <= left_lmax && right_non_region_weight + complement_weight <= right_lmax) {
                        weight_t left  = left_non_region_weight + closure_weight;
                        weight_t right = right_non_region_weight + complement_weight;
                        weight_t diff  = std::min(left_lmax - left, right_lmax - right);
                        if (diff > best_diff) {
                            closure_found = true;
                            best_diff     = diff;
                            best_closure  = in_closure;
                        }
                    }
                }
            }

            is_left.resize(n - 2);
            for (vertex_t u = 0; u < n - 2; ++u) { is_left[u] = best_closure[scc_id[u]]; }

            return closure_found;
        }

        bool find_best_closure_old(weight_t left_non_region_weight,
                               weight_t right_non_region_weight,
                               weight_t left_lmax,
                               weight_t right_lmax,
                               size_t repeats,
                               RandomEngine &rnd_engine,
                               std::vector<u8> &is_left) {
            bool            closure_found = false;
            weight_t        best_diff     = -1;
            best_closure.resize(n_scc);

            // determine which sccs do not have to be considered
            is_active.resize(n_scc);
            std::fill(is_active.begin(), is_active.end(), 1);
            for (vertex_t   scc_u: scc_s_successors) { is_active[scc_u] = 0; }
            for (vertex_t   scc_u: scc_t_predecessors) { is_active[scc_u] = 0; }

            for (size_t i = 0; i < repeats; ++i) {
                // determine in degree of each vertex
                in_deg.resize(n_scc);
                std::fill(in_deg.begin(), in_deg.end(), 0);
                for (vertex_t         scc_u = 0; scc_u < n_scc; ++scc_u) {
                    if (is_active[scc_u] == 0) { continue; }
                    for (vertex_t scc_v: edges[scc_u]) {
                        if (is_active[scc_v] == 0) { continue; }
                        in_deg[scc_v] += 1;
                    }
                }

                // push roots into stack
                stack.clear();
                for (vertex_t         scc_u = 0; scc_u < n_scc; ++scc_u) {
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
                        if (in_deg[scc_v] == 0) { stack.push_back(scc_v); }
                    }
                    std::shuffle(stack.begin(), stack.end(), rnd_engine.generator);
                }

                // go through the order and determine the best closure
                weight_t        closure_weight = 0;
                in_closure.resize(n_scc);
                std::fill(in_closure.begin(), in_closure.end(), 0);
                for (vertex_t   scc_u: scc_s_successors) {
                    in_closure[scc_u] = 1;
                    closure_weight += scc_weights[scc_u];
                }

                for (vertex_t scc_u: topo_order) {
                    if (is_active[scc_u] == 0) { continue; }

                    in_closure[scc_u] = 1;
                    closure_weight += scc_weights[scc_u];
                }

                weight_t      complement_weight = 0;
                for (vertex_t scc_u: scc_t_predecessors) {
                    in_closure[scc_u] = 0;
                    complement_weight += scc_weights[scc_u];
                }

                if (left_non_region_weight + closure_weight <= left_lmax && right_non_region_weight + complement_weight <= right_lmax) {
                    weight_t left  = left_non_region_weight + closure_weight;
                    weight_t right = right_non_region_weight + complement_weight;
                    weight_t diff  = std::min(left_lmax - left, right_lmax - right);
                    if (diff > best_diff) {
                        closure_found = true;
                        best_diff     = diff;
                        best_closure  = in_closure;
                    }
                }

                for (vertex_t scc_u: topo_order) {
                    if (is_active[scc_u] == 0) { continue; }

                    in_closure[scc_u] = 0;
                    closure_weight -= scc_weights[scc_u];
                    complement_weight += scc_weights[scc_u];

                    if (left_non_region_weight + closure_weight <= left_lmax && right_non_region_weight + complement_weight <= right_lmax) {
                        weight_t left  = left_non_region_weight + closure_weight;
                        weight_t right = right_non_region_weight + complement_weight;
                        weight_t diff  = std::min(left_lmax - left, right_lmax - right);
                        if (diff > best_diff) {
                            closure_found = true;
                            best_diff     = diff;
                            best_closure  = in_closure;
                        }
                    }
                }
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

            for (auto [u, w]: residual_flow_network[v]) {
                if (index[u] == UNVISITED) {
                    dfs(u, residual_flow_network);
                } else if (scc_id[u] == UNVISITED) {
                    while (!P.empty() && index[P.back()] > index[u]) {
                        P.pop_back();
                    }
                }
            }

            if (!P.empty() && P.back() == v) {
                vertex_per_scc.emplace_back();  // Start a new SCC
                vertex_t w;
                do {
                    w = S.back();
                    S.pop_back();
                    scc_id[w] = n_scc;
                    vertex_per_scc.back().push_back(w);  // Add to current SCC
                } while (w != v);
                P.pop_back();
                ++n_scc;
            }
        }

        void dfs_non_recursive(vertex_t start, const ResidualFlowNetwork &residual_flow_network) {
            std::stack<std::pair<vertex_t, size_t>> dfs_stack; // (node, next neighbor index)

            dfs_stack.push({start, 0});
            index[start] = idx++;
            S.push_back(start);
            P.push_back(start);

            while (!dfs_stack.empty()) {
                auto       &[v, i]    = dfs_stack.top();  // current node and index into its neighbors
                const auto &neighbors = residual_flow_network[v];

                if (i < neighbors.size()) {
                    vertex_t u = neighbors[i++].v;
                    if (index[u] == UNVISITED) {
                        index[u] = idx++;
                        S.push_back(u);
                        P.push_back(u);
                        dfs_stack.push({u, 0});  // dive deeper
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

    class FlowNetwork {
        vertex_t             n;
        Graph<int, int, int> g;
        vertex_t             source;
        vertex_t             target;

        // std::vector<std::vector<EdgeVW>> adj;

    public:
        FlowNetwork() : g(Graph<int, int, int>(0, 0)) {
            n      = 0;
            source = 0;
            target = 0;
        }

        void initialize(size_t t_n) {
            n = t_n;
            g.reset();
            g.add_node(n + 2);
            source = n;
            target = n + 1;

            /*
            if (adj.size() != n + 2) {
                adj.resize(n + 2);
            }
            for (auto &vec: adj) {
                vec.clear();
            }
             */
        }

        void add(vertex_t u, vertex_t v, weight_t w) {
            ASSERT(u < n);
            ASSERT(v < n);
            ASSERT(w >= 0);
            // for(auto & e : adj[u]){ ASSERT(e.v != v); }
            // for(auto & e : adj[v]){ ASSERT(e.v != u); }

            g.add_edge(u, v, w, w);

            // adj[u].emplace_back(v, w);
            // adj[v].emplace_back(u, w);
        }

        void add_s_edge(vertex_t v, weight_t w) {
            ASSERT(v < n);
            ASSERT(w >= 0);
            // for(auto & e : adj[source]){ ASSERT(e.v != v); }

            g.add_edge(source, v, w, 0);

            // adj[source].emplace_back(v, w);
        }

        void add_t_edge(vertex_t v, weight_t w) {
            ASSERT(v < n);
            ASSERT(w >= 0);
            // for(auto & e : adj[target]){ ASSERT(e.v != v); }

            g.add_edge(v, target, w, 0);

            // adj[v].emplace_back(target, w);
        }

        void build_residual_network(ResidualFlowNetwork &residual_g) {
            residual_g.initialize(n);

            int                          n_edges = g.get_arc_num();
            Graph<int, int, int>::arc_id arc     = g.get_first_arc();

            int u, v;
            for (int     i = 0; i < n_edges; ++i) {
                g.get_arc_ends(arc, u, v);
                weight_t w = g.get_rcap(arc);

                if (w > 0) {
                    if (u == source) {
                        residual_g.add_edge_from_source(v, w);
                    } else if (v == source) {
                        residual_g.add_edge_to_source(u, w);
                    } else if (u == target) {
                        residual_g.add_edge_from_target(v, w);
                    } else if (v == target) {
                        residual_g.add_edge_to_target(u, w);
                    } else {
                        residual_g.add_directed_edge(u, v, w);
                    }
                }

                arc = g.get_next_arc(arc);
            }
        }

        void solve() {
            /*
            for (auto &vec: adj) {

                for (size_t i = 0; i < vec.size(); ++i) {
                    for (size_t j = i + 1; j < vec.size(); ++j) {
                        ASSERT(vec[i].v != vec[j].v);
                    }
                }
            }
             */

            const int INF = std::numeric_limits<int>::max() / 2;
            g.add_tweights(source, INF, 0);
            g.add_tweights(target, 0, INF);

            g.maxflow();
        }

        void get_cut(std::vector<u8> &is_left) {
            is_left.resize(n);
            for (vertex_t u = 0; u < n; ++u) {
                is_left[u] = g.what_segment(u) == Graph<int, int, int>::SOURCE;
            }
        }

        void print() const {
            /*
            std::cout << "Flow Network Adjacency List:\n";
            std::cout << "-----------------------------\n";

            for (vertex_t u = 0; u < adj.size(); ++u) {
                if (u == source) {
                    std::cout << "Source (" << u << "): ";
                } else if (u == target) {
                    std::cout << "Target (" << u << "): ";
                } else {
                    std::cout << "Node " << u << ": ";
                }

                for (const auto &edge: adj[u]) {
                    std::cout << edge.v << "(w=" << edge.w << ") ";
                }

                std::cout << "\n";
            }

            std::cout << "-----------------------------\n";
             */
        }
    };

}

#endif //HEIPROMAP_FLOW_H
