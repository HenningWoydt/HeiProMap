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

#ifndef HEIPROMAP_SIMPLE_FLOW_BASED_REFINEMENT_H
#define HEIPROMAP_SIMPLE_FLOW_BASED_REFINEMENT_H

#include <algorithm>
#include <cmath>
#include <vector>
#include <limits>
#include <omp.h>

#include "../definitions.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../datastructures/distance_1_matching.h"
#include "../datastructures/active_block_scheduling.h"
#include "../datastructures/block_adjacency.h"
#include "../utility/aligned_array.h"
#include "../utility/profiler.h"
#include "../utility/random_engine.h"
#include "../utility/translation_table.h"
#include "../utility/push_relabel.h"
#include "../utility/memory_stack.h"
#include "../utility/flow.h"

namespace HeiProMap {

    class SimpleFlowBasedRefinementConfiguration final {
    public:
        explicit SimpleFlowBasedRefinementConfiguration(std::string t_name) : name(std::move(t_name)) {}

        std::string name;
        bool enabled = false;
        u64 max_global_iteration = 1;
        u64 max_local_iteration = 2;
        f64 alpha = 2.0;
        f64 alpha_upper_bound = 8.0;
        f64 alpha_modifier = 2.0;
        bool use_closed_vertex_set = true;
        u64 closed_vertex_sets_repeats = 100;
        bool use_edge_cut = true;
    };

    template<bool LARGE_K>
    class SimpleFlowBasedRefinement {
        vertex_t m_n = 0;
        partition_t m_k = 0;
        u64 m_threads = 1;

        const SimpleFlowBasedRefinementConfiguration *config = nullptr;
        std::vector<RandomEngine> rnd_engines;

        // Block adjacency and independent pair matching (no dense quotient graph)
        BlockAdjacency block_adjacency;
        Distance1Matching<LARGE_K> d1_matcher;
        ActiveBlockScheduling active_block_scheduling;

        // Per-thread state for local flow networks
        std::vector<AlignedArray<u32>> seen_vecs;
        std::vector<AlignedArray<u32>> region_vecs;
        std::vector<u32> seen_marker_vecs;
        std::vector<u32> region_marker_vecs;

        std::vector<std::vector<vertex_t>> left_boundaries;
        std::vector<std::vector<vertex_t>> right_boundaries;
        std::vector<std::vector<vertex_t>> left_regions;
        std::vector<std::vector<vertex_t>> right_regions;
        std::vector<std::vector<vertex_t>> queues;
        std::vector<std::vector<u8>> is_left_vecs;
        std::vector<std::vector<u8>> is_left_2_vecs;
        std::vector<std::vector<u8>> s_connected_vecs;
        std::vector<std::vector<u8>> t_connected_vecs;

        std::vector<PushRelabel<weight_t>> prs;
        std::vector<MemoryStack> pr_mems;
        std::vector<ResidualFlowNetwork> residual_flow_networks;
        std::vector<SCCGraph> scc_graphs;
        std::vector<TranslationTable<vertex_t>> translation_tables;

    public:
        SimpleFlowBasedRefinement() = default;
        ~SimpleFlowBasedRefinement() = default;

        void initialize(vertex_t n, partition_t k, u64 threads, u64 seed,
                        const SimpleFlowBasedRefinementConfiguration &cfg,
                        MemoryPool *pool = nullptr) {
            m_n = n;
            m_k = k;
            m_threads = threads;
            config = &cfg;

            block_adjacency.initialize(k);
            d1_matcher.initialize(k, pool);
            active_block_scheduling.initialize(k, pool);

            rnd_engines.clear();
            rnd_engines.reserve(m_threads);
            for (u64 t = 0; t < m_threads; ++t) {
                rnd_engines.emplace_back(seed + t * 997 + 1);
            }

            seen_vecs.resize(m_threads);
            region_vecs.resize(m_threads);
            seen_marker_vecs.assign(m_threads, 0);
            region_marker_vecs.assign(m_threads, 0);

            left_boundaries.resize(m_threads);
            right_boundaries.resize(m_threads);
            left_regions.resize(m_threads);
            right_regions.resize(m_threads);
            queues.resize(m_threads);
            is_left_vecs.resize(m_threads);
            is_left_2_vecs.resize(m_threads);
            s_connected_vecs.resize(m_threads);
            t_connected_vecs.resize(m_threads);

            prs.resize(m_threads);
            pr_mems.resize(m_threads);
            residual_flow_networks.resize(m_threads);
            scc_graphs.resize(m_threads);
            translation_tables.resize(m_threads);

            for (u64 t = 0; t < m_threads; ++t) {
                seen_vecs[t].set_pool(pool);
                seen_vecs[t].initialize(m_n, 0);
                region_vecs[t].set_pool(pool);
                region_vecs[t].initialize(m_n, 0);
                translation_tables[t].reserve(m_n, m_n);
            }
        }

        template<typename DistanceOracleT>
        void refine(graph_t &g,
                    DistanceOracleT &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager) {
            if (g.uniform_v_weights && g.uniform_e_weights) {
                refine_impl<true, true>(g, d_oracle, bv_manager, p_manager);
            } else if (g.uniform_v_weights) {
                refine_impl<true, false>(g, d_oracle, bv_manager, p_manager);
            } else if (g.uniform_e_weights) {
                refine_impl<false, true>(g, d_oracle, bv_manager, p_manager);
            } else {
                refine_impl<false, false>(g, d_oracle, bv_manager, p_manager);
            }
        }

    private:
        template<bool t_uniform_v_weights, bool t_uniform_e_weights, typename DistanceOracleT>
        void refine_impl(graph_t &g,
                         DistanceOracleT &d_oracle,
                         bv_manager_t &bv_manager,
                         p_manager_t &p_manager) {
            active_block_scheduling.reset(m_k);
            std::vector<std::vector<std::pair<partition_t, partition_t>>> static_matchings;

            for (u64 iteration = 0; iteration < config->max_global_iteration; ++iteration) {
                HEIPROMAP_PROFILE_SCOPE("refinement", "SimpleFlowRefinement", "recompute_adjacency");
                block_adjacency.compute(g, p_manager, bv_manager, m_threads);

                d1_matcher.reset_used_edges();
                bool any_positive_move = false;

                d1_matcher.compute_static_matchings(
                    block_adjacency, p_manager, active_block_scheduling.active_this_round,
                    static_matchings, 512 * 512, 8);

                if (static_matchings.empty()) break;

                for (size_t round = 0; round < static_matchings.size(); ++round) {
                    const auto &matching = static_matchings[round];
                    if (matching.empty()) continue;

                    #pragma omp parallel for num_threads(m_threads) schedule(dynamic)
                    for (size_t j = 0; j < matching.size(); ++j) {
                        const partition_t u_id = matching[j].first;
                        const partition_t v_id = matching[j].second;
                        const u64 tid = omp_get_thread_num();

                        bool improved = refine_pair<t_uniform_v_weights, t_uniform_e_weights>(
                            g, d_oracle, bv_manager, p_manager, u_id, v_id, tid);

                        if (improved) {
                            active_block_scheduling.activate(u_id, v_id);
                            #pragma omp atomic write
                            any_positive_move = true;
                        }
                    }
                }

                if (!any_positive_move) break;
                active_block_scheduling.next_round();
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights, typename DistanceOracleT>
        bool refine_pair(graph_t &g,
                         DistanceOracleT &d_oracle,
                         bv_manager_t &bv_manager,
                         p_manager_t &p_manager,
                         partition_t left_id,
                         partition_t right_id,
                         u64 tid) {
            RandomEngine &random_engine = rnd_engines[tid];
            f64 alpha = config->alpha;
            f64 alpha_upper_bound = config->alpha_upper_bound;
            f64 alpha_modifier = config->alpha_modifier;

            std::vector<vertex_t> &left_boundary = left_boundaries[tid];
            std::vector<vertex_t> &right_boundary = right_boundaries[tid];
            weight_t left_boundary_weight = 0;
            weight_t right_boundary_weight = 0;

            std::vector<vertex_t> &left_region = left_regions[tid];
            std::vector<vertex_t> &right_region = right_regions[tid];

            PushRelabel<weight_t> &pr = prs[tid];
            MemoryStack &pr_mem = pr_mems[tid];
            ResidualFlowNetwork &residual_flow_network = residual_flow_networks[tid];
            SCCGraph &scc_graph = scc_graphs[tid];

            TranslationTable<vertex_t> &translation_table = translation_tables[tid];
            AlignedArray<u32> &seen = seen_vecs[tid];
            AlignedArray<u32> &region_marker = region_vecs[tid];
            u32 &seen_mark = seen_marker_vecs[tid];
            u32 &region_mark = region_marker_vecs[tid];

            std::vector<vertex_t> &queue = queues[tid];
            std::vector<u8> &is_left = is_left_vecs[tid];
            std::vector<u8> &is_left_2 = is_left_2_vecs[tid];
            std::vector<u8> &s_connected = s_connected_vecs[tid];
            std::vector<u8> &t_connected = t_connected_vecs[tid];

            f64 avg_weight = (f64) g.g_weight / (f64) m_k;
            bool total_improved = false;

            for (u64 iter = 0; iter < config->max_local_iteration; ++iter) {
                left_boundary.clear();
                right_boundary.clear();
                left_boundary_weight = 0;
                right_boundary_weight = 0;
                left_region.clear();
                right_region.clear();

                // Find boundary vertices between left_id and right_id
                for (const vertex_t u : bv_manager.boundary(left_id)) {
                    bool adjacent_to_right = false;
                    for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                        if (p_manager[g.edges_v[i]] == right_id) {
                            adjacent_to_right = true;
                            break;
                        }
                    }
                    if (adjacent_to_right) {
                        left_boundary.push_back(u);
                        left_boundary_weight += t_uniform_v_weights ? 1 : g.v_weight(u);
                    }
                }

                for (const vertex_t v : bv_manager.boundary(right_id)) {
                    bool adjacent_to_left = false;
                    for (size_t i = g.neighborhoods[v]; i < g.neighborhoods[v + 1]; ++i) {
                        if (p_manager[g.edges_v[i]] == left_id) {
                            adjacent_to_left = true;
                            break;
                        }
                    }
                    if (adjacent_to_left) {
                        right_boundary.push_back(v);
                        right_boundary_weight += t_uniform_v_weights ? 1 : g.v_weight(v);
                    }
                }

                if (left_boundary.empty() || right_boundary.empty()) break;

                weight_t adapt_lmax_left = (weight_t) std::ceil(avg_weight + alpha * ((f64) p_manager.lmax[left_id] - avg_weight));
                weight_t adapt_lmax_right = (weight_t) std::ceil(avg_weight + alpha * ((f64) p_manager.lmax[right_id] - avg_weight));

                weight_t left_max_weight = adapt_lmax_left - p_manager.get_bweight(right_id);
                weight_t right_max_weight = adapt_lmax_right - p_manager.get_bweight(left_id);

                if (region_mark > std::numeric_limits<u32>::max() - 4) {
                    std::fill_n(region_marker.get_ptr(), m_n, 0);
                    region_mark = 1;
                }
                if (seen_mark > std::numeric_limits<u32>::max() - 4) {
                    std::fill_n(seen.get_ptr(), m_n, 0);
                    seen_mark = 1;
                }
                region_mark += 2;
                u32 left_mark = region_mark - 1;
                u32 right_mark = region_mark;
                seen_mark += 2;

                weight_t left_region_weight = determine_region<t_uniform_v_weights>(
                    g, p_manager, left_id, left_mark, left_max_weight, left_boundary, left_region, seen, seen_mark, region_marker, queue);
                weight_t right_region_weight = determine_region<t_uniform_v_weights>(
                    g, p_manager, right_id, right_mark, right_max_weight, right_boundary, right_region, seen, seen_mark, region_marker, queue);

                if (left_region.empty() && right_region.empty()) {
                    if (alpha >= alpha_upper_bound) break;
                    alpha = std::min(alpha_modifier * alpha, alpha_upper_bound);
                    continue;
                }

                vertex_t new_u = 0;
                for (size_t i = 0; i < left_region.size(); ++i) translation_table.add(left_region[i], new_u++);
                for (size_t i = 0; i < right_region.size(); ++i) translation_table.add(right_region[i], new_u++);

                // Build flow network with or without tree penalties
                if (!d_oracle.last_level_pair(left_id, right_id) || !config->use_edge_cut) {
                    build_flow_network_with_penalties<t_uniform_e_weights>(
                        g, d_oracle, p_manager, left_id, right_id, left_region, right_region,
                        pr, pr_mem, translation_table, region_marker, region_mark, alpha, s_connected, t_connected);
                } else {
                    build_flow_network_no_penalties<t_uniform_e_weights>(
                        g, p_manager, left_id, right_id, left_region, right_region,
                        pr, pr_mem, translation_table, region_marker, region_mark, s_connected, t_connected);
                }

                pr.maxflow();

                size_t flow_n = left_region.size() + right_region.size();
                is_left.resize(flow_n);
                for (size_t i = 0; i < flow_n; ++i) {
                    is_left[i] = (pr.what_label(static_cast<int>(i)) == SOURCE) ? 1 : 0;
                }

                bool is_valid = cut_is_valid<t_uniform_v_weights>(
                    g, p_manager, left_id, right_id, is_left, left_region, right_region, translation_table);

                if (!is_valid && !config->use_closed_vertex_set) {
                    if (alpha <= 1.0) break;
                    alpha = std::max(alpha / alpha_modifier, 1.0);
                    continue;
                }

                if (!is_valid) {
                    build_residual_from_cut(
                        g, pr, left_region, right_region, translation_table, region_marker, region_mark, s_connected, t_connected, residual_flow_network);
                    residual_flow_network.finalize();

                    scc_graph.initialize(residual_flow_network, g, translation_table);
                    scc_graph.reduce();

                    weight_t left_non_region_weight = p_manager.get_bweight(left_id) - left_region_weight;
                    weight_t right_non_region_weight = p_manager.get_bweight(right_id) - right_region_weight;
                    bool closure_found = scc_graph.find_best_closure(
                        left_non_region_weight, right_non_region_weight,
                        p_manager.get_lmax(left_id), p_manager.get_lmax(right_id),
                        avg_weight, config->closed_vertex_sets_repeats, random_engine, is_left_2);

                    if (closure_found) {
                        weight_t mincut_gain = calculate_gain<t_uniform_e_weights>(
                            g, p_manager, d_oracle, is_left, left_id, right_id, left_region, right_region, translation_table, region_marker, region_mark);
                        weight_t closure_gain = calculate_gain<t_uniform_e_weights>(
                            g, p_manager, d_oracle, is_left_2, left_id, right_id, left_region, right_region, translation_table, region_marker, region_mark);

                        if (closure_gain >= mincut_gain) {
                            std::swap(is_left, is_left_2);
                        } else {
                            closure_found = false;
                        }
                    }

                    if (!closure_found) {
                        if (alpha <= 1.0) break;
                        alpha = std::max(alpha / alpha_modifier, 1.0);
                        continue;
                    }
                }

                // Check if cut changes partition
                bool cut_changed = false;
                for (size_t i = 0; i < left_region.size(); ++i) {
                    if (is_left[translation_table.get_n(left_region[i])] == 0) { cut_changed = true; break; }
                }
                if (!cut_changed) {
                    for (size_t i = 0; i < right_region.size(); ++i) {
                        if (is_left[translation_table.get_n(right_region[i])] == 1) { cut_changed = true; break; }
                    }
                }
                if (!cut_changed) {
                    if (alpha <= 1.0) break;
                    alpha = std::max(alpha / alpha_modifier, 1.0);
                    continue;
                }

                weight_t gain = calculate_gain<t_uniform_e_weights>(
                    g, p_manager, d_oracle, is_left, left_id, right_id, left_region, right_region, translation_table, region_marker, region_mark);

                if (gain <= 0) {
                    if (alpha <= 1.0) break;
                    alpha = std::max(alpha / alpha_modifier, 1.0);
                    continue;
                }

                // Apply cut moves to p_manager and bv_manager
                for (size_t i = 0; i < left_region.size(); ++i) {
                    vertex_t u = left_region[i];
                    if (is_left[translation_table.get_n(u)] == 0) {
                        weight_t w = t_uniform_v_weights ? 1 : g.v_weight(u);
                        bv_manager.move(g, p_manager, u, left_id, right_id);
                        p_manager.move_serial(u, w, left_id, right_id);
                    }
                }
                for (size_t i = 0; i < right_region.size(); ++i) {
                    vertex_t u = right_region[i];
                    if (is_left[translation_table.get_n(u)] == 1) {
                        weight_t w = t_uniform_v_weights ? 1 : g.v_weight(u);
                        bv_manager.move(g, p_manager, u, right_id, left_id);
                        p_manager.move_serial(u, w, right_id, left_id);
                    }
                }

                total_improved = true;
                alpha = std::min(alpha * alpha_modifier, alpha_upper_bound);
            }

            return total_improved;
        }

        template<bool t_uniform_v_weights>
        weight_t determine_region(const graph_t &g,
                                  const p_manager_t &p_manager,
                                  partition_t id,
                                  u32 mark,
                                  weight_t max_weight,
                                  const std::vector<vertex_t> &boundary,
                                  std::vector<vertex_t> &region,
                                  AlignedArray<u32> &seen,
                                  u32 &seen_mark,
                                  AlignedArray<u32> &region_marker,
                                  std::vector<vertex_t> &queue) {
            seen_mark += 2;
            weight_t curr_weight = 0;
            queue.clear();

            for (const vertex_t u : boundary) {
                queue.push_back(u);
                seen[u] = seen_mark - 1;
            }

            size_t queue_idx = 0;
            while (queue_idx < queue.size()) {
                vertex_t u = queue[queue_idx++];
                weight_t w = t_uniform_v_weights ? 1 : g.v_weight(u);

                if (curr_weight + w > max_weight) continue;
                curr_weight += w;
                region.push_back(u);
                region_marker[u] = mark;
                seen[u] = seen_mark;

                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    vertex_t v = g.edges_v[i];
                    if (p_manager[v] != id) continue;
                    if (seen[v] != seen_mark && seen[v] != seen_mark - 1) {
                        queue.push_back(v);
                        seen[v] = seen_mark - 1;
                    }
                }
            }
            return curr_weight;
        }

        template<bool t_uniform_e_weights, typename DistanceOracleT>
        void build_flow_network_with_penalties(const graph_t &g,
                                               const DistanceOracleT &d_oracle,
                                               const p_manager_t &p_manager,
                                               partition_t left_id,
                                               partition_t right_id,
                                               std::vector<vertex_t> &left_region,
                                               std::vector<vertex_t> &right_region,
                                               PushRelabel<weight_t> &pr,
                                               MemoryStack &pr_mem,
                                               TranslationTable<vertex_t> &translation_table,
                                               AlignedArray<u32> &region_marker,
                                               u32 &region_mark,
                                               f64,
                                               std::vector<u8> &s_connected,
                                               std::vector<u8> &t_connected) {
            u32 left_mark = region_mark - 1;
            u32 right_mark = region_mark;
            weight_t distance = d_oracle.get(left_id, right_id);

            int n = static_cast<int>(left_region.size() + right_region.size());
            int src = n, snk = n + 1;
            int m = 2 * n;
            for (size_t j = 0; j < left_region.size(); ++j) {
                m += static_cast<int>(g.neighborhoods[left_region[j] + 1] - g.neighborhoods[left_region[j]]);
            }
            for (size_t j = 0; j < right_region.size(); ++j) {
                m += static_cast<int>(g.neighborhoods[right_region[j] + 1] - g.neighborhoods[right_region[j]]);
            }

            pr.init(n + 2, src, snk, m, pr_mem);
            s_connected.assign(n, 0);
            t_connected.assign(n, 0);

            for (size_t j = 0; j < left_region.size(); ++j) {
                vertex_t u = left_region[j];
                int new_u = static_cast<int>(translation_table.get_n(u));
                weight_t left_penalty = 0;
                weight_t right_penalty = 0;

                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    const vertex_t v = g.edges_v[i];
                    const weight_t w = g.edge_w(i);

                    if (region_marker[v] == right_mark) {
                        int new_v = static_cast<int>(translation_table.get_n(v));
                        weight_t cap = t_uniform_e_weights ? distance : w * distance;
                        pr.add_edge(new_u, new_v, cap, cap);
                        continue;
                    }
                    if (region_marker[v] == left_mark) {
                        if (u < v) continue;
                        int new_v = static_cast<int>(translation_table.get_n(v));
                        weight_t cap = t_uniform_e_weights ? distance : w * distance;
                        pr.add_edge(new_u, new_v, cap, cap);
                        continue;
                    }

                    partition_t v_id = p_manager[v];
                    if constexpr (t_uniform_e_weights) {
                        left_penalty += d_oracle.get(left_id, v_id);
                        right_penalty += d_oracle.get(right_id, v_id);
                    } else {
                        left_penalty += w * d_oracle.get(left_id, v_id);
                        right_penalty += w * d_oracle.get(right_id, v_id);
                    }
                }

                if (left_penalty > 0) {
                    pr.add_edge(new_u, snk, left_penalty);
                    t_connected[new_u] = 1;
                }
                if (right_penalty > 0) {
                    pr.add_edge(src, new_u, right_penalty);
                    s_connected[new_u] = 1;
                }
            }

            for (size_t j = 0; j < right_region.size(); ++j) {
                vertex_t u = right_region[j];
                int new_u = static_cast<int>(translation_table.get_n(u));
                weight_t left_penalty = 0;
                weight_t right_penalty = 0;

                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    const vertex_t v = g.edges_v[i];
                    const weight_t w = g.edge_w(i);
                    if (region_marker[v] == left_mark || region_marker[v] == right_mark) {
                        if (region_marker[v] == right_mark && u < v) continue;
                        int new_v = static_cast<int>(translation_table.get_n(v));
                        weight_t cap = t_uniform_e_weights ? distance : w * distance;
                        pr.add_edge(new_u, new_v, cap, cap);
                        continue;
                    }

                    partition_t v_id = p_manager[v];
                    if constexpr (t_uniform_e_weights) {
                        left_penalty += d_oracle.get(left_id, v_id);
                        right_penalty += d_oracle.get(right_id, v_id);
                    } else {
                        left_penalty += w * d_oracle.get(left_id, v_id);
                        right_penalty += w * d_oracle.get(right_id, v_id);
                    }
                }

                if (left_penalty > 0) {
                    pr.add_edge(new_u, snk, left_penalty);
                    t_connected[new_u] = 1;
                }
                if (right_penalty > 0) {
                    pr.add_edge(src, new_u, right_penalty);
                    s_connected[new_u] = 1;
                }
            }
        }

        template<bool t_uniform_e_weights>
        void build_flow_network_no_penalties(const graph_t &g,
                                             const p_manager_t &p_manager,
                                             partition_t left_id,
                                             partition_t right_id,
                                             std::vector<vertex_t> &left_region,
                                             std::vector<vertex_t> &right_region,
                                             PushRelabel<weight_t> &pr,
                                             MemoryStack &pr_mem,
                                             TranslationTable<vertex_t> &translation_table,
                                             AlignedArray<u32> &region_marker,
                                             u32 &region_mark,
                                             std::vector<u8> &s_connected,
                                             std::vector<u8> &t_connected) {
            u32 left_mark = region_mark - 1;
            u32 right_mark = region_mark;

            int n = static_cast<int>(left_region.size() + right_region.size());
            int src = n, snk = n + 1;
            int m = 2 * n;
            for (size_t j = 0; j < left_region.size(); ++j) {
                m += static_cast<int>(g.neighborhoods[left_region[j] + 1] - g.neighborhoods[left_region[j]]);
            }
            for (size_t j = 0; j < right_region.size(); ++j) {
                m += static_cast<int>(g.neighborhoods[right_region[j] + 1] - g.neighborhoods[right_region[j]]);
            }

            pr.init(n + 2, src, snk, m, pr_mem);
            s_connected.assign(n, 0);
            t_connected.assign(n, 0);

            for (size_t j = 0; j < left_region.size(); ++j) {
                vertex_t u = left_region[j];
                int new_u = static_cast<int>(translation_table.get_n(u));
                weight_t w_left = 0, w_right = 0;

                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    const vertex_t v = g.edges_v[i];
                    const weight_t w = g.edge_w(i);
                    if (region_marker[v] == left_mark || region_marker[v] == right_mark) {
                        if (u < v) continue;
                        int new_v = static_cast<int>(translation_table.get_n(v));
                        weight_t cap = t_uniform_e_weights ? 1 : w;
                        pr.add_edge(new_u, new_v, cap, cap);
                        continue;
                    }
                    partition_t v_id = p_manager[v];
                    if (v_id == left_id) w_left += (t_uniform_e_weights ? 1 : w);
                    else if (v_id == right_id) w_right += (t_uniform_e_weights ? 1 : w);
                }

                if (w_left > 0) { pr.add_edge(src, new_u, w_left); s_connected[new_u] = 1; }
                if (w_right > 0) { pr.add_edge(new_u, snk, w_right); t_connected[new_u] = 1; }
            }

            for (size_t j = 0; j < right_region.size(); ++j) {
                vertex_t u = right_region[j];
                int new_u = static_cast<int>(translation_table.get_n(u));
                weight_t w_left = 0, w_right = 0;

                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    const vertex_t v = g.edges_v[i];
                    const weight_t w = g.edge_w(i);
                    if (region_marker[v] == left_mark || region_marker[v] == right_mark) {
                        if (u < v) continue;
                        int new_v = static_cast<int>(translation_table.get_n(v));
                        weight_t cap = t_uniform_e_weights ? 1 : w;
                        pr.add_edge(new_u, new_v, cap, cap);
                        continue;
                    }
                    partition_t v_id = p_manager[v];
                    if (v_id == left_id) w_left += (t_uniform_e_weights ? 1 : w);
                    else if (v_id == right_id) w_right += (t_uniform_e_weights ? 1 : w);
                }

                if (w_left > 0) { pr.add_edge(src, new_u, w_left); s_connected[new_u] = 1; }
                if (w_right > 0) { pr.add_edge(new_u, snk, w_right); t_connected[new_u] = 1; }
            }
        }

        template<bool t_uniform_v_weights>
        bool cut_is_valid(const graph_t &g,
                          const p_manager_t &p_manager,
                          partition_t left_id,
                          partition_t right_id,
                          const std::vector<u8> &is_left,
                          const std::vector<vertex_t> &left_region,
                          const std::vector<vertex_t> &right_region,
                          const TranslationTable<vertex_t> &translation_table) {
            weight_t left_weight = p_manager.get_bweight(left_id);
            weight_t right_weight = p_manager.get_bweight(right_id);

            for (const vertex_t u : left_region) {
                if (is_left[translation_table.get_n(u)] == 0) {
                    weight_t w = t_uniform_v_weights ? 1 : g.v_weight(u);
                    left_weight -= w;
                    right_weight += w;
                }
            }
            for (const vertex_t u : right_region) {
                if (is_left[translation_table.get_n(u)] == 1) {
                    weight_t w = t_uniform_v_weights ? 1 : g.v_weight(u);
                    right_weight -= w;
                    left_weight += w;
                }
            }
            return (left_weight <= p_manager.get_lmax(left_id) && right_weight <= p_manager.get_lmax(right_id));
        }

        template<bool t_uniform_e_weights, typename DistanceOracleT>
        weight_t calculate_gain(const graph_t &g,
                                const p_manager_t &p_manager,
                                const DistanceOracleT &d_oracle,
                                const std::vector<u8> &is_left,
                                partition_t left_id,
                                partition_t right_id,
                                const std::vector<vertex_t> &left_region,
                                const std::vector<vertex_t> &right_region,
                                const TranslationTable<vertex_t> &translation_table,
                                const AlignedArray<u32> &region_marker,
                                u32 region_mark) {
            u32 left_mark = region_mark - 1;
            u32 right_mark = region_mark;
            weight_t gain = 0;

            auto process_vertex = [&](vertex_t u, partition_t old_id, partition_t new_id) {
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    vertex_t v = g.edges_v[i];
                    weight_t w = t_uniform_e_weights ? 1 : g.edge_w(i);
                    partition_t v_id = p_manager[v];
                    partition_t v_new_id = v_id;

                    bool v_moves = false;
                    if (region_marker[v] == left_mark) {
                        if (is_left[translation_table.get_n(v)] == 0) {
                            v_new_id = right_id;
                            v_moves = true;
                        }
                    } else if (region_marker[v] == right_mark) {
                        if (is_left[translation_table.get_n(v)] == 1) {
                            v_new_id = left_id;
                            v_moves = true;
                        }
                    }

                    if (v_moves && u > v) continue;
                    weight_t old_dist = d_oracle.get(old_id, v_id);
                    weight_t new_dist = d_oracle.get(new_id, v_new_id);
                    gain += w * (old_dist - new_dist);
                }
            };

            for (vertex_t u : left_region) {
                if (is_left[translation_table.get_n(u)] == 0) process_vertex(u, left_id, right_id);
            }
            for (vertex_t u : right_region) {
                if (is_left[translation_table.get_n(u)] == 1) process_vertex(u, right_id, left_id);
            }
            return gain;
        }

        void build_residual_from_cut(const graph_t &g,
                                     const PushRelabel<weight_t> &pr,
                                     const std::vector<vertex_t> &left_region,
                                     const std::vector<vertex_t> &right_region,
                                     [[maybe_unused]] const TranslationTable<vertex_t> &translation_table,
                                     [[maybe_unused]] const AlignedArray<u32> &region_marker,
                                     [[maybe_unused]] u32 region_mark,
                                     [[maybe_unused]] const std::vector<u8> &s_connected,
                                     [[maybe_unused]] const std::vector<u8> &t_connected,
                                     ResidualFlowNetwork &residual_flow_network) {
            vertex_t flow_n = left_region.size() + right_region.size();
            residual_flow_network.initialize(flow_n);

            int pr_n = static_cast<int>(flow_n + 2);
            for (int u = 0; u < pr_n; ++u) {
                for (int i = pr.arc_first(u); i < pr.arc_last(u); ++i) {
                    if (pr.arc_cap(i) > 0) {
                        int v = pr.arc_to(i);
                        residual_flow_network.add_directed_edge(static_cast<vertex_t>(u), static_cast<vertex_t>(v), 1);
                    }
                }
            }
        }
    };
}

#endif // HEIPROMAP_SIMPLE_FLOW_BASED_REFINEMENT_H
