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

#ifndef HEIPROMAP_FLOW_BASED_REFINEMENT_H
#define HEIPROMAP_FLOW_BASED_REFINEMENT_H

#include <algorithm>

#include "../../extern/maxflow-v3.04.src/graph.h"

#include "quotient_graph_refinement.h"
#include "../../commons/indexed_update_heap.h"
#include "../../commons/random_engine.h"
#include "../../commons/statistic_collector.h"
#include "../../commons/utils.h"
#include "../datastructures/functions.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/qap.h"

namespace HeiProMap {
    class ResidualFlowNetwork {
        vertex_t n      = 0;
        vertex_t source = 0;
        vertex_t target = 0;
        std::vector<std::vector<EdgeVW>> edges;

    public:
        ResidualFlowNetwork() = default;

        void initialize(vertex_t t_n) {
            n      = t_n;
            source = n;
            target = n + 1;
            edges.clear();
            edges.resize(n + 2);
        }

        void add_directed_edge(vertex_t u, vertex_t v, weight_t w) {
            edges[u].emplace_back(v, w);
        }

        void add_edge_to_source(vertex_t u, weight_t w) {
            edges[u].emplace_back(source, w);
        }

        void add_edge_to_target(vertex_t u, weight_t w) {
            edges[u].emplace_back(target, w);
        }

        void add_edge_from_source(vertex_t u, weight_t w) {
            edges[source].emplace_back(u, w);
        }

        void add_edge_from_target(vertex_t u, weight_t w) {
            edges[target].emplace_back(u, w);
        }

        vertex_t get_n() const { return n; }
        vertex_t get_source() const { return source; }
        vertex_t get_target() const { return target; }

        const std::vector<EdgeVW>& operator[](vertex_t i) const { return edges[i]; }
    };

    class SCCGraph {
        vertex_t n      = 0;
        vertex_t source = 0;
        vertex_t target = 0;

        std::vector<std::vector<vertex_t>> edges;
        std::vector<std::vector<vertex_t>> rev_edges;
        vertex_t scc_source = 0;
        vertex_t scc_target = 0;
        std::vector<weight_t> scc_weight;
        weight_t total_weight = 0;

        std::vector<vertex_t> scc_s_successors;
        weight_t scc_s_succ_weight = 0;
        std::vector<vertex_t> scc_t_predecessors;
        weight_t scc_t_pred_weight = 0;

        static constexpr vertex_t UNVISITED = std::numeric_limits<vertex_t>::max();
        std::vector<vertex_t> index;
        std::vector<vertex_t> scc_id;
        std::vector<vertex_t> S;
        std::vector<vertex_t> P;
        vertex_t idx         = 0;
        vertex_t current_scc = 0;

    public:
        SCCGraph() = default;

        void initialize(const ResidualFlowNetwork& residual_flow_network,
                        const graph_t& g,
                        const TranslationTable<vertex_t>& tt) {
            n      = residual_flow_network.get_n() + 2;
            source = residual_flow_network.get_source();
            target = residual_flow_network.get_target();

            index.clear();
            scc_id.clear();
            S.clear();
            P.clear();
            index.resize(n, UNVISITED);
            scc_id.resize(n, UNVISITED);
            idx         = 0;
            current_scc = 0;

            // === Gabow's Algorithm ===
            for (vertex_t v = 0; v < n; ++v) {
                if (index[v] == UNVISITED) { dfs(v, residual_flow_network); }
            }

            // clear old edges
            edges.clear();
            edges.resize(current_scc);

            rev_edges.clear();
            rev_edges.resize(current_scc);

            // build the graph and the reversed graph
            for (vertex_t u = 0; u < n; ++u) {
                vertex_t scc_u = scc_id[u];
                for (const auto [v, w] : residual_flow_network[u]) {
                    vertex_t scc_v = scc_id[v];
                    if (scc_u == scc_v) { continue; }
                    edges[scc_u].emplace_back(scc_v);
                    rev_edges[scc_v].emplace_back(scc_u);
                }
            }

            // make each edge list unique
            for (vertex_t scc_u = 0; scc_u < current_scc; ++scc_u) {
                std::sort(edges[scc_u].begin(), edges[scc_u].end());
                edges[scc_u].erase(unique(edges[scc_u].begin(), edges[scc_u].end()), edges[scc_u].end());
            }

            // set special scc
            scc_source = scc_id[source];
            scc_target = scc_id[target];

            // determine the weight of each scc in the graph
            scc_weight.clear();
            scc_weight.resize(current_scc, 0);
            for (vertex_t u = 0; u < n; ++u) {
                if (u != source && u != target) {
                    scc_weight[scc_id[u]] += g.weight(tt.get_o(u));
                }
            }

            // get total weight of all scc
            total_weight = 0;
            for (weight_t w : scc_weight) { total_weight += w; }
        }

        void reduce() {
            std::vector<vertex_t> curr;
            std::vector<vertex_t> next;
            std::vector<u8> seen(current_scc, 0);

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

                    for (vertex_t scc_v : edges[scc_u]) {
                        if (seen[scc_v] == 1 || seen[scc_v] == 2) { continue; }
                        next.push_back(scc_v);
                        seen[scc_v] = 1;
                    }
                    seen[scc_u] = 2;
                }
            }

            scc_s_succ_weight = 0;
            for (vertex_t scc : scc_s_successors) { scc_s_succ_weight += scc_weight[scc]; }

            curr.clear();
            next.clear();
            seen.clear();
            seen.resize(current_scc, 0);
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

                    for (vertex_t scc_v : rev_edges[scc_u]) {
                        if (seen[scc_v] == 1 || seen[scc_v] == 2) { continue; }
                        next.push_back(scc_v);
                        seen[scc_v] = 1;
                    }
                    seen[scc_u] = 2;
                }
            }

            scc_t_pred_weight = 0;
            for (vertex_t scc : scc_t_predecessors) { scc_t_pred_weight += scc_weight[scc]; }

            ASSERT(no_duplicates(scc_s_successors));
            ASSERT(no_duplicates(scc_t_predecessors));
        }

        bool find_best_closure(weight_t left_non_region_weight,
                               weight_t right_non_region_weight,
                               weight_t both_region_weight,
                               weight_t lmax,
                               size_t repeats,
                               RandomEngine& rnd_engine,
                               std::vector<u8>& is_left) {
            bool closure_found = false;
            weight_t best_diff = std::numeric_limits<weight_t>::max();
            std::vector<u8> best_closure(current_scc);

            // determine which sccs do not have to be considered
            std::vector<u8> is_active(current_scc, 1);
            for (vertex_t scc : scc_s_successors) { is_active[scc] = 0; }
            for (vertex_t scc : scc_t_predecessors) { is_active[scc] = 0; }

            for (size_t i = 0; i < repeats; ++i) {
                // determine in degree of each vertex
                std::vector<vertex_t> in_deg(current_scc, 0);
                for (vertex_t scc = 0; scc < current_scc; ++scc) {
                    for (vertex_t scc_n : edges[scc]) { in_deg[scc_n] += 1; }
                }

                // push roots into stack
                std::vector<vertex_t> stack;
                for (vertex_t scc = 0; scc < current_scc; ++scc) {
                    if (in_deg[scc] > 0) { stack.push_back(scc); }
                }
                std::shuffle(stack.begin(), stack.end(), rnd_engine.gen);

                // determine the random topological order
                std::vector<vertex_t> topo_order;
                while (!stack.empty()) {
                    vertex_t u = stack.back();
                    stack.pop_back();

                    topo_order.push_back(u);

                    for (vertex_t scc : edges[u]) {
                        in_deg[scc] -= 1;
                        if (in_deg[scc] == 0) { stack.push_back(scc); }
                    }
                    std::shuffle(stack.begin(), stack.end(), rnd_engine.gen);
                }

                // go through the order and determine the best closure
                weight_t closure_weight = scc_s_succ_weight;
                std::vector<u8> in_closure(current_scc, 0);
                for (vertex_t scc : scc_s_successors) { in_closure[scc] = 1; }

                for (vertex_t u : topo_order) {
                    if (is_active[u] == 0) { continue; }

                    in_closure[u] = 1;
                    closure_weight += scc_weight[u];
                }

                weight_t complement_weight = scc_t_pred_weight;

                if (right_non_region_weight + complement_weight > lmax) { return false; }

                if (left_non_region_weight + closure_weight <= lmax && right_non_region_weight + complement_weight <= lmax) {
                    weight_t diff = std::abs((both_region_weight / 2) - (left_non_region_weight + closure_weight));
                    if (diff < best_diff) {
                        closure_found = true;
                        best_diff     = diff;
                        best_closure  = in_closure;
                    }
                }

                for (vertex_t u : topo_order) {
                    if (is_active[u] == 0) { continue; }

                    in_closure[u] = 0;
                    closure_weight -= scc_weight[u];
                    complement_weight += scc_weight[u];

                    if (left_non_region_weight + closure_weight > lmax) { break; }

                    if (left_non_region_weight + closure_weight <= lmax && right_non_region_weight + complement_weight <= lmax) {
                        weight_t diff = std::abs((both_region_weight / 2) - (left_non_region_weight + closure_weight));
                        if (diff < best_diff) {
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

        void print() {
            std::cout << "SCC Graph" << std::endl;
            std::cout << "Source  : " << scc_source << std::endl;
            std::cout << "Target  : " << scc_target << std::endl;
            std::cout << "n       : " << n << std::endl;
            std::cout << "curr_scc: " << current_scc << std::endl;

            for (vertex_t scc_u = 0; scc_u < current_scc; ++scc_u) {
                std::cout << scc_u << " (" << scc_weight[scc_u] << ") : ";
                for (vertex_t scc_v : edges[scc_u]) {
                    std::cout << scc_v << ", ";
                }
                std::cout << std::endl;
            }
        }

    private:
        void dfs(vertex_t v, const ResidualFlowNetwork& residual_flow_network) {
            index[v] = idx++;
            S.push_back(v);
            P.push_back(v);

            for (auto [u, w] : residual_flow_network[v]) {
                if (index[u] == UNVISITED) {
                    dfs(u, residual_flow_network);
                } else if (scc_id[u] == UNVISITED) {
                    while (!P.empty() && index[P.back()] > index[u]) {
                        P.pop_back();
                    }
                }
            }

            if (!P.empty() && P.back() == v) {
                vertex_t w;
                do {
                    w = S.back();
                    S.pop_back();
                    scc_id[w] = current_scc;
                } while (w != v);
                P.pop_back();
                ++current_scc;
            }
        }
    };

    class FlowNetwork {
        vertex_t n;
        Graph<weight_t, weight_t, weight_t> g;
        vertex_t source;
        vertex_t target;

    public:
        FlowNetwork() : g(Graph<weight_t, weight_t, weight_t>(0, 0)) {
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
        }

        void add(vertex_t u, vertex_t v, weight_t w) {
            ASSERT(u < n);
            ASSERT(v < n);
            ASSERT(w >= 0);
            g.add_edge(u, v, w, w);
        }

        void add_s_edge(vertex_t v, weight_t w) {
            ASSERT(v < n);
            ASSERT(w >= 0);
            g.add_edge(source, v, w, 0);
        }

        void add_t_edge(vertex_t v, weight_t w) {
            ASSERT(v < n);
            ASSERT(w >= 0);
            g.add_edge(v, target, w, 0);
        }

        void build_residual_network(ResidualFlowNetwork& residual_g) {
            residual_g.initialize(n);

            int n_edges                                     = g.get_arc_num();
            Graph<weight_t, weight_t, weight_t>::arc_id arc = g.get_first_arc();
            unsigned int u, v;
            for (int i = 0; i < n_edges; ++i) {
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
            const weight_t INF = std::numeric_limits<weight_t>::max() / 2;
            g.add_tweights(source, INF, 0);
            g.add_tweights(target, 0, INF);

            g.maxflow();
        }

        void get_cut(std::vector<u8>& is_left) {
            is_left.resize(n);
            for (vertex_t u = 0; u < n; ++u) {
                is_left[u] = g.what_segment(u) == Graph<weight_t, weight_t, weight_t>::SOURCE;
            }
        }
    };

    class FlowBasedRefinementConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit FlowBasedRefinementConfiguration(const std::string& t_name) : ISerialRefinerConfiguration(t_name) {}

        u64 max_global_iteration   = 1;
        u64 max_local_iteration    = 3;
        f64 alpha                  = 2.0;
        f64 alpha_upper_bound      = 8.0;
        f64 alpha_modifier         = 2.0;
        bool use_closed_vertex_set = true;
    };

    class FlowBasedRefinement final : public ISerialRefiner {
        vertex_t m_n    = 0;
        vertex_t m_m    = 0;
        partition_t m_k = 0;
        f64 m_imbalance = 0.0;
        weight_t m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;

        // active block scheduling
        u8* active_this_round = nullptr;
        u8* active_next_round = nullptr;
        PairWeight* pairs     = nullptr;
        size_t pairs_size     = 0;

        // array for boundary vertices
        vertex_t* left_boundary   = nullptr;
        size_t left_boundary_size = 0;

        vertex_t* right_boundary   = nullptr;
        size_t right_boundary_size = 0;

        // array for regions
        vertex_t* left_region   = nullptr;
        size_t left_region_size = 0;

        vertex_t* right_region   = nullptr;
        size_t right_region_size = 0;

        u32* is_left_region  = nullptr;
        u32* is_right_region = nullptr;
        u32 is_region_mark   = 0;

        u32* bfs_level = nullptr;

        vertex_t* queue   = nullptr;
        size_t queue_size = 0;

        u32* seen     = nullptr;
        u32 seen_mark = 0;

        // array for penalties
        weight_t* left_penalties  = nullptr;
        weight_t* right_penalties = nullptr;

        //Translation Table for mapping
        TranslationTable<vertex_t> translation_table;

        FlowNetwork flow_network;
        ResidualFlowNetwork residual_flow_network;
        SCCGraph scc_graph;

        RandomEngine* random_engine                    = nullptr;
        const FlowBasedRefinementConfiguration* config = nullptr;
        StatisticCollector* m_stat_collector           = nullptr;

    public:
        FlowBasedRefinement() = default;

        ~FlowBasedRefinement() override {
            free(active_this_round);
            free(active_next_round);

            free(left_boundary);
            free(right_boundary);

            free(left_region);
            free(right_region);
            free(is_left_region);
            free(is_right_region);
            free(queue);
            free(seen);

            free(left_penalties);
            free(right_penalties);
        }

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const f64 t_imbalance,
                        const weight_t t_lmax,
                        const std::vector<partition_t>& t_hierarchy,
                        const std::vector<weight_t>& t_distance,
                        RandomEngine& t_random_engine,
                        const ISerialRefinerConfiguration& i_config,
                        StatisticCollector& t_stat_collect) override {
            m_n         = t_n;
            m_m         = t_m;
            m_k         = t_k;
            m_lmax      = t_lmax;
            m_imbalance = t_imbalance;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;

            random_engine    = &t_random_engine;
            config           = dynamic_cast<const FlowBasedRefinementConfiguration*>(&i_config);
            m_stat_collector = &t_stat_collect;

            vertex_t m_n_64        = round_up_64(m_n);
            partition_t m_k_64     = round_up_64(m_k);
            partition_t m_k_m_k_64 = round_up_64(m_k * m_k);

            // active block scheduling
            active_this_round = (u8*)aligned_alloc(64, m_k_64 * sizeof(u8));
            active_next_round = (u8*)aligned_alloc(64, m_k_64 * sizeof(u8));
            pairs             = (PairWeight*)aligned_alloc(64, m_k_m_k_64 * sizeof(PairWeight));
            pairs_size        = 0;

            left_boundary      = (vertex_t*)aligned_alloc(64, m_n_64 * sizeof(vertex_t));
            left_boundary_size = 0;

            right_boundary      = (vertex_t*)aligned_alloc(64, m_n_64 * sizeof(vertex_t));
            right_boundary_size = 0;

            left_region      = (vertex_t*)aligned_alloc(64, m_n_64 * sizeof(vertex_t));
            left_region_size = 0;

            right_region      = (vertex_t*)aligned_alloc(64, m_n_64 * sizeof(vertex_t));
            right_region_size = 0;

            bfs_level = (u32*)aligned_alloc(64, m_n_64 * sizeof(u32));

            is_left_region = (u32*)aligned_alloc(64, m_n_64 * sizeof(u32));
            std::fill_n(is_left_region, m_n_64, 0);
            is_right_region = (u32*)aligned_alloc(64, m_n_64 * sizeof(u32));
            std::fill_n(is_right_region, m_n_64, 0);
            is_region_mark = 0;

            queue      = (vertex_t*)aligned_alloc(64, m_n_64 * sizeof(vertex_t));
            queue_size = 0;

            seen = (u32*)aligned_alloc(64, m_n_64 * sizeof(u32));
            std::fill_n(seen, m_n_64, 0);
            seen_mark = 0;

            left_penalties  = (weight_t*)aligned_alloc(64, m_n_64 * sizeof(weight_t));
            right_penalties = (weight_t*)aligned_alloc(64, m_n_64 * sizeof(weight_t));

            translation_table.reserve(m_n_64, m_n_64);
        }

        void refine(const u64 level,
                    const u64 max_level,
                    const graph_t& g,
                    const d_oracle_t& d_oracle,
                    bv_manager_t& bv_manager,
                    p_manager_t& p_manager,
                    q_graph_t& q_graph) override {
            std::fill_n(active_this_round, m_k, 1);
            std::fill_n(active_next_round, m_k, 0);

            for (u64 iteration = 0; iteration < config->max_global_iteration; ++iteration) {
                // determine all pairs in the quotient graph
                pairs_size = 0;
                for (partition_t u_id = 0; u_id < m_k; ++u_id) {
                    for (partition_t v_id = u_id + 1; v_id < m_k; ++v_id) {
                        if (q_graph.has_edge(u_id, v_id) && (active_this_round[u_id] || active_this_round[v_id])) {
                            pairs[pairs_size++] = {u_id, v_id, d_oracle.get(u_id, v_id)};
                        }
                    }
                }
                std::sort(pairs, pairs + pairs_size, std::greater<>());

                if (pairs_size == 0) { return; }

                for (size_t i = 0; i < pairs_size; ++i) {
                    partition_t left_id  = pairs[i].id1;
                    partition_t right_id = pairs[i].id2;
                    refine_blocks(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, left_id, right_id);
                }

                std::swap(active_this_round, active_next_round);
                std::fill_n(active_next_round, m_k, 0);
            }
        }

        void refine_blocks(const u64 level,
                           const u64 max_level,
                           const graph_t& g,
                           const d_oracle_t& d_oracle,
                           bv_manager_t& bv_manager,
                           p_manager_t& p_manager,
                           q_graph_t& q_graph,
                           partition_t left_id,
                           partition_t right_id) {
            ASSERT(left_id != right_id);
            ASSERT(max(p_manager.get_bweights()) <= m_lmax);

            f64 alpha             = config->alpha;
            f64 alpha_upper_bound = config->alpha_upper_bound;
            f64 alpha_modifier    = config->alpha_modifier;

            u64 max_local_iteration = config->max_local_iteration;
            u64 iteration           = 0;
            while (iteration < max_local_iteration) {
                iteration += 1;
                // get boundary vertices
                determine_boundary_vertices(g, bv_manager, p_manager, left_id, right_id);

                // calc max weight for each bfs
                weight_t lmax             = std::ceil((1.0 + m_imbalance * alpha) * ((f64)g.weight() / (f64)m_k));
                weight_t left_max_weight  = lmax - p_manager.get_bweight(right_id);
                weight_t right_max_weight = lmax - p_manager.get_bweight(left_id);

                // get both regions
                weight_t left_region_weight;
                weight_t right_region_weight;
                determine_regions(g, p_manager, left_id, left_max_weight, &left_region_weight, right_id, right_max_weight, &right_region_weight);

                if (left_region_size + right_region_size == 0) {
                    // if both regions are empty, increase their sizes
                    if (alpha == alpha_upper_bound) {
                        // the regions will not change next iteration, so just return
                        return;
                    }
                    alpha = std::min(alpha_modifier * alpha, alpha_upper_bound);
                    continue;
                }

                // determine penalties for all vertices
                determine_penalties(g, p_manager, d_oracle, left_id, right_id);

                // build a translation table from graph to flow network
                vertex_t new_u = 0;
                for (size_t i = 0; i < left_region_size; ++i) { translation_table.add(left_region[i], new_u++); }
                for (size_t i = 0; i < right_region_size; ++i) { translation_table.add(right_region[i], new_u++); }

                // build flownetwork
                build_flow_network(g, d_oracle, left_id, right_id);

                // solve the flow network
                flow_network.solve();

                bool qap_normal_calculated = false;
                weight_t qap_normal_change;
                std::vector<u8> is_left;
                if (config->use_closed_vertex_set) {
#if HEAVYASSERT_ENABLED
                    // get the first cut for comparison
                    flow_network.get_cut(is_left);

                    if (cut_is_valid(g, p_manager, left_id, right_id, is_left)) {
                        // make the changes
                        // change_boundary(g, bv_manager, p_manager, q_graph, is_left, left_id, right_id);
                        // HEAVYASSERT(assert_correct_boundary(g, p_manager, bv_manager, m_k));

                        // qap_normal_calculated = true;
                        // qap_normal_change     = get_qap(g, p_manager, d_oracle);

                        // change_boundary(g, bv_manager, p_manager, q_graph, is_left, left_id, right_id);
                        // HEAVYASSERT(assert_correct_boundary(g, p_manager, bv_manager, m_k));
                    }
#endif

                    // build residual network
                    build_residual_flow_network();

                    // build scc graph
                    scc_graph.initialize(residual_flow_network, g, translation_table);

                    // reduce the scc graph
                    scc_graph.reduce();

                    // determine best balanced min cut
                    weight_t left_non_region_weight  = p_manager.get_bweight(left_id) - left_region_weight;
                    weight_t right_non_region_weight = p_manager.get_bweight(right_id) - right_region_weight;
                    weight_t both_region_weight      = p_manager.get_bweight(left_id) + p_manager.get_bweight(right_id);
                    bool closure_found               = scc_graph.find_best_closure(left_non_region_weight, right_non_region_weight, both_region_weight, m_lmax, 10, *random_engine, is_left);

                    if (!closure_found) {
                        alpha = std::max(alpha / alpha_modifier, 1.0);
                        continue;
                    }
                } else {
                    // simply use the first cut found
                    flow_network.get_cut(is_left);

                    if (!cut_is_valid(g, p_manager, left_id, right_id, is_left)) {
                        alpha = std::max(alpha / alpha_modifier, 1.0);
                        continue;
                    }
                }

                // check if the cut actually changes the partition
                if (!cut_changes_partition(is_left)) {
                    // cut is valid, but does not change anything
                    alpha = std::max(alpha / alpha_modifier, 1.0);
                    continue;
                }

                // cut is valid and changes the partition, increase alpha
                alpha = std::min(alpha * alpha_modifier, alpha_upper_bound);

                // make the changes
                change_boundary(g, bv_manager, p_manager, q_graph, is_left, left_id, right_id);
                HEAVYASSERT(assert_correct_boundary(g, p_manager, bv_manager, m_k));

#if HEAVYASSERT_ENABLED
                if (qap_normal_calculated) {
                    weight_t qap_closure_change = get_qap(g, p_manager, d_oracle);
                    if (qap_normal_change != qap_closure_change) {
                        std::cout << "Normal change: " << qap_normal_change << std::endl;
                        std::cout << "Closure change: " << qap_closure_change << std::endl;
                        exit(EXIT_FAILURE);
                    }
                }
#endif
                active_next_round[left_id]  = 1;
                active_next_round[right_id] = 1;
            }
            ASSERT(max(p_manager.get_bweights()) <= m_lmax);
        }

        void determine_boundary_vertices(const graph_t& g,
                                         const bv_manager_t& bv_manager,
                                         const p_manager_t& p_manager,
                                         partition_t left_id,
                                         partition_t right_id) {
            left_boundary_size = 0;
            forall_bv_id_iu(bv_manager, left_id, i, u)
                {
                    forall_guiv(g, u, j, v)
                        {
                            partition_t v_id = p_manager[v];
                            if (v_id == right_id) {
                                left_boundary[left_boundary_size++] = u;
                                break;
                            }
                        }
                    endfor
                }
            endfor

            right_boundary_size = 0;
            forall_bv_id_iu(bv_manager, right_id, i, u)
                {
                    forall_guiv(g, u, j, v)
                        {
                            partition_t v_id = p_manager[v];
                            if (v_id == left_id) {
                                right_boundary[right_boundary_size++] = u;
                                break;
                            }
                        }
                    endfor
                }
            endfor
        }

        void determine_regions(const graph_t& g,
                               const p_manager_t& p_manager,
                               partition_t left_id,
                               weight_t left_max_weight,
                               weight_t* left_region_weight,
                               partition_t right_id,
                               weight_t right_max_weight,
                               weight_t* right_region_weight) {
            is_region_mark += 1;
            seen_mark += 2;
            // seen[u] == seen_mark     means u is processed
            // seen[u] == seen_mark - 1 means u is in the queue

            weight_t left_curr_weight = 0;

            queue_size = 0;
            for (size_t i = 0; i < left_boundary_size; ++i) {
                vertex_t u          = left_boundary[i];
                ASSERT(p_manager[u] == left_id);
                queue[queue_size++] = u;
                seen[u]             = seen_mark - 1;
                bfs_level[u]        = 0;
            }

            size_t queue_idx = 0;
            left_region_size = 0;
            while (queue_idx < queue_size) {
                vertex_t u = queue[queue_idx++];
                if (seen[u] == seen_mark) { continue; }

                if (left_curr_weight + g.weight(u) <= left_max_weight) {
                    left_region[left_region_size++] = u;
                    is_left_region[u]               = is_region_mark;
                    left_curr_weight += g.weight(u);
                    forall_guiv(g, u, i, v)
                        {
                            partition_t v_id = p_manager[v];
                            if (v_id != left_id) { continue; }

                            if (seen[v] != seen_mark && seen[v] != seen_mark - 1) {
                                queue[queue_size++] = v;
                                seen[v]             = seen_mark - 1;
                                bfs_level[v]        = bfs_level[u] + 1;
                            }
                        }
                    endfor
                }
                seen[u] = seen_mark;
            }
            *left_region_weight = left_curr_weight;

            weight_t right_curr_weight = 0;

            queue_size = 0;
            for (size_t i = 0; i < right_boundary_size; ++i) {
                vertex_t u          = right_boundary[i];
                ASSERT(p_manager[u] == right_id);
                queue[queue_size++] = u;
                seen[u]             = seen_mark - 1;
                bfs_level[u]        = 0;
            }

            right_region_size = 0;
            queue_idx         = 0;
            while (queue_idx < queue_size) {
                vertex_t u = queue[queue_idx++];
                if (seen[u] == seen_mark) { continue; }

                if (right_curr_weight + g.weight(u) <= right_max_weight) {
                    right_region[right_region_size++] = u;
                    is_right_region[u]                = is_region_mark;
                    right_curr_weight += g.weight(u);
                    forall_guiv(g, u, i, v)
                        {
                            partition_t v_id = p_manager[v];
                            if (v_id != right_id) { continue; }

                            if (seen[v] != seen_mark && seen[v] != seen_mark - 1) {
                                queue[queue_size++] = v;
                                seen[v]             = seen_mark - 1;
                                bfs_level[v]        = bfs_level[u] + 1;
                            }
                        }
                    endfor
                }
                seen[u] = seen_mark;
            }
            *right_region_weight = right_curr_weight;
        }

        void determine_penalties(const graph_t& g,
                                 const p_manager_t& p_manager,
                                 const d_oracle_t& d_oracle,
                                 partition_t left_id,
                                 partition_t right_id) {
            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u         = left_region[j];
                left_penalties[u]  = 0;
                right_penalties[u] = 0;
                forall_guivw(g, u, i, v, w)
                    {
                        if (is_left_region[v] == is_region_mark || is_right_region[v] == is_region_mark) { continue; } // ignore neighbors that are in the region
                        partition_t v_id = p_manager[v];
                        left_penalties[u] += w * d_oracle.get(left_id, v_id);
                        right_penalties[u] += w * d_oracle.get(right_id, v_id);
                    }
                endfor
                left_penalties[u] *= 2;
                right_penalties[u] *= 2;
            }
            for (size_t j = 0; j < right_region_size; ++j) {
                vertex_t u         = right_region[j];
                left_penalties[u]  = 0;
                right_penalties[u] = 0;
                forall_guivw(g, u, i, v, w)
                    {
                        if (is_right_region[v] == is_region_mark || is_left_region[v] == is_region_mark) { continue; } // ignore neighbors that are in the region
                        partition_t v_id = p_manager[v];
                        left_penalties[u] += w * d_oracle.get(left_id, v_id);
                        right_penalties[u] += w * d_oracle.get(right_id, v_id);
                    }
                endfor
                left_penalties[u] *= 2;
                right_penalties[u] *= 2;
            }
        }

        void build_flow_network(const graph_t& g,
                                const d_oracle_t& d_oracle,
                                partition_t left_id,
                                partition_t right_id) {
            weight_t distance = d_oracle.get(left_id, right_id);

            // build flownetwork
            size_t n = left_region_size + right_region_size;
            flow_network.initialize(n);

            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u = left_region[j];
                forall_guivw(g, u, i, v, w)
                    {
                        weight_t mult = 2;
                        if (is_left_region[v] != is_region_mark && is_right_region[v] != is_region_mark) { continue; } // vertex gets handled by penalties
                        if (is_left_region[v] == is_region_mark && bfs_level[u] < bfs_level[v]) { continue; } // only forward edges allowed
                        if (is_left_region[v] == is_region_mark && bfs_level[u] == bfs_level[v]) { mult = 1; } // the edge from v to u will also be added, therefore only one time the weight

                        vertex_t new_u = translation_table.get_n(u);
                        vertex_t new_v = translation_table.get_n(v);

                        ASSERT(new_u != new_v);

                        flow_network.add(new_u, new_v, mult * w * distance);
                    }
                endfor
            }
            for (size_t j = 0; j < right_region_size; ++j) {
                vertex_t u = right_region[j];
                forall_guivw(g, u, i, v, w)
                    {
                        weight_t mult = 2;
                        if (is_right_region[v] != is_region_mark) { continue; } // vertex gets handled by penalties, or if v belongs to the left region, no edge is made
                        if (bfs_level[u] > bfs_level[v]) { continue; } // only forward edges allowed
                        if (bfs_level[u] == bfs_level[v]) { mult = 1; } // the edge from v to u will also be added, therefore only one time the weight
                        vertex_t new_u = translation_table.get_n(u);
                        vertex_t new_v = translation_table.get_n(v);

                        ASSERT(new_u != new_v);

                        flow_network.add(new_u, new_v, mult * w * distance);
                    }
                endfor
            }

            // add the penalties
            for (size_t i = 0; i < left_region_size; ++i) {
                vertex_t u             = left_region[i];
                vertex_t new_u         = translation_table.get_n(u);
                weight_t left_penalty  = left_penalties[u];
                weight_t right_penalty = right_penalties[u];
                if (left_penalty > 0) { flow_network.add_t_edge(new_u, left_penalty); }
                if (right_penalty > 0) { flow_network.add_s_edge(new_u, right_penalty); }
            }
            for (size_t i = 0; i < right_region_size; ++i) {
                vertex_t u             = right_region[i];
                vertex_t new_u         = translation_table.get_n(u);
                weight_t left_penalty  = left_penalties[u];
                weight_t right_penalty = right_penalties[u];
                if (left_penalty > 0) { flow_network.add_t_edge(new_u, left_penalty); }
                if (right_penalty > 0) { flow_network.add_s_edge(new_u, right_penalty); }
            }
        }

        void build_residual_flow_network() {
            vertex_t n = left_region_size + right_region_size;
            residual_flow_network.initialize(n);
        }

        bool cut_is_valid(const graph_t& g,
                          const p_manager_t& p_manager,
                          partition_t left_id,
                          partition_t right_id,
                          std::vector<u8>& is_left) {
            weight_t left_weight  = p_manager.get_bweight(left_id);
            weight_t right_weight = p_manager.get_bweight(right_id);
            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u        = left_region[j];
                weight_t u_weight = g.weight(u);
                vertex_t new_u    = translation_table.get_n(u);
                if (is_left[new_u] == 0) {
                    left_weight -= u_weight;
                    right_weight += u_weight;
                }
            }

            for (size_t j = 0; j < right_region_size; ++j) {
                vertex_t u        = right_region[j];
                weight_t u_weight = g.weight(u);
                vertex_t new_u    = translation_table.get_n(u);
                if (is_left[new_u] == 1) {
                    right_weight -= u_weight;
                    left_weight += u_weight;
                }
            }

            return left_weight <= m_lmax && right_weight <= m_lmax;
        }

        bool cut_changes_partition(std::vector<u8>& is_left) {
            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u     = left_region[j];
                vertex_t new_u = translation_table.get_n(u);
                if (is_left[new_u] == 0) {
                    return true;
                }
            }
            for (size_t j = 0; j < right_region_size; ++j) {
                vertex_t u     = right_region[j];
                vertex_t new_u = translation_table.get_n(u);
                if (is_left[new_u] == 1) {
                    return true;
                }
            }

            return false;
        }

        void change_boundary(const graph_t& g,
                             bv_manager_t& bv_manager,
                             p_manager_t& p_manager,
                             q_graph_t& q_graph,
                             std::vector<u8>& is_left,
                             partition_t left_id,
                             partition_t right_id) {
            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u     = left_region[j];
                vertex_t new_u = translation_table.get_n(u);

                ASSERT(is_left_region[u] == is_region_mark);
                ASSERT(new_u < left_region_size + right_region_size);

                if (is_left[new_u] == 0) {
                    if (bv_manager.is_boundary(u)) {
                        bv_manager.move(g, p_manager, u, left_id, right_id);
                    } else {
                        bv_manager.add_new(g, p_manager, u, right_id);
                    }

                    q_graph.move(g, p_manager, u, left_id, right_id);
                    p_manager.move(u, g.weight(u), left_id, right_id);
                }
            }
            for (size_t j = 0; j < right_region_size; ++j) {
                vertex_t u     = right_region[j];
                vertex_t new_u = translation_table.get_n(u);

                ASSERT(is_right_region[u] == is_region_mark);
                ASSERT(new_u < left_region_size + right_region_size);

                if (is_left[new_u] == 1) {
                    if (bv_manager.is_boundary(u)) {
                        bv_manager.move(g, p_manager, u, right_id, left_id);
                    } else {
                        bv_manager.add_new(g, p_manager, u, left_id);
                    }

                    q_graph.move(g, p_manager, u, right_id, left_id);
                    p_manager.move(u, g.weight(u), right_id, left_id);
                }
            }
        }

        JSONString get_stats() override {
            std::string stats = "{ \n";
#if COLLECT_METRICS

#endif
            stats.pop_back();
            stats.pop_back();
            stats += "\n}";

            JSONString json_stats;
            json_stats.s = stats;
            return json_stats;
        }
    };
}

#endif //HEIPROMAP_FLOW_BASED_REFINEMENT_H
