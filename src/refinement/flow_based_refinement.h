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
#include <cmath>
#include <queue>
#include <stack>
#include <string>
#include <unordered_set>
#include <vector>

#include <omp.h>

#include "../definitions.h"
#include "../datastructures/block_conn.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/distance_oracle.h"
#include "../datastructures/partition_manager.h"
#include "../datastructures/quotient_graph.h"
#include "../utility/aligned_array.h"
#include "../utility/memory_stack.h"
#include "../utility/profiler.h"
#include "../utility/translation_table.h"
#include "quotient_graph_refinement.h"
#include "../utility/flow.h"
#include "../utility/dinics.h"
#include "../utility/push_relabel.h"
#include "../utility/random_engine.h"
#include "../utility/utils.h"
#include "../utility/assert_state.h"
#include "../utility/qap.h"

namespace HeiProMap {
    enum struct GrowthStrategy {
        BFS,
        HEAVY_FIRST,
        LIGHT_FIRST
    };

    class FlowBasedRefinementConfiguration final {
    public:
        explicit FlowBasedRefinementConfiguration(const std::string &t_name) {
            name = t_name;
        }

        std::string name;
        bool enabled = false;
        bool use_active_block_scheduling = true;
        u64 max_global_iteration = 1;
        u64 max_local_iteration = 3;
        f64 alpha = 2.0;
        f64 alpha_upper_bound = 8.0;
        f64 alpha_modifier = 2.0;
        bool use_closed_vertex_set = true;
        u64 closed_vertex_sets_repeats = 10;
        bool always_include_boundary = false;
        GrowthStrategy growth_strategy = GrowthStrategy::BFS;
        bool use_edge_cut = true;
    };

    class FlowBasedRefinement final {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
        u64 m_threads = 1;

        u64 m_scc_successes = 0;
        u64 m_scc_failures = 0;

        std::vector<AlignedArray<u32> > seen_vecs;
        std::vector<AlignedArray<u32> > region_vecs;
        std::vector<TranslationTable<vertex_t> > translation_tables;
        std::vector<u32> seen_marker_vecs;
        std::vector<u32> region_marker_vecs;

        const FlowBasedRefinementConfiguration *config = nullptr;
        std::vector<RandomEngine> rnd_engines;

        // per-thread persistent storage to avoid repeated allocations
        std::vector<std::vector<vertex_t> > left_boundaries;
        std::vector<std::vector<vertex_t> > right_boundaries;
        std::vector<std::vector<vertex_t> > left_regions;
        std::vector<std::vector<vertex_t> > right_regions;
        std::vector<std::vector<vertex_t> > queues;
        std::vector<std::vector<u8> > is_left_vecs;
        std::vector<std::vector<u8> > is_left_2_vecs;
        std::vector<std::vector<u8> > s_connected_vecs;
        std::vector<std::vector<u8> > t_connected_vecs;
        std::vector<PushRelabel<weight_t> > prs;
        std::vector<MemoryStack> pr_mems;
        std::vector<ResidualFlowNetwork> residual_flow_networks;
        std::vector<SCCGraph> scc_graphs;

    public:
        FlowBasedRefinement() = default;

        ~FlowBasedRefinement() = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const u64 t_threads,
                        const u64 seed,
                        const FlowBasedRefinementConfiguration &i_config) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "FlowBasedRefinement", "initialize");

            m_n = t_n;
            m_m = t_m;
            m_k = t_k;
            m_threads = t_threads;

            seen_vecs.resize(m_threads);
            region_vecs.resize(m_threads);
            translation_tables.resize(m_threads);
            for (u64 t = 0; t < m_threads; ++t) {
                seen_vecs[t].initialize(m_n, 0);
                region_vecs[t].initialize(m_n, 0);
                translation_tables[t].reserve(m_n, m_n);
            }
            seen_marker_vecs.resize(m_threads, 1);
            region_marker_vecs.resize(m_threads, 1);

            config = &i_config;

            rnd_engines.resize(m_threads);
            for (u64 t = 0; t < m_threads; ++t) {
                rnd_engines[t] = RandomEngine(seed + t);
            }

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
        }

        void refine(graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph,
                    block_conn_t &block_conn,
                    const AlignedArray<weight_t> &lmax_constraints) {
            if (g.uniform_v_weights && g.uniform_e_weights) refine_impl<true, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
            else if (g.uniform_v_weights) refine_impl<true, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
            else if (g.uniform_e_weights) refine_impl<false, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
            else refine_impl<false, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void refine_impl(graph_t &g,
                         d_oracle_t &d_oracle,
                         bv_manager_t &bv_manager,
                         p_manager_t &p_manager,
                         q_graph_t &q_graph,
                         block_conn_t &block_conn,
                         const AlignedArray<weight_t> &lmax_constraints) {
            m_scc_successes = 0;
            m_scc_failures = 0;

            // active block scheduling
            AlignedArray<u8> active_this_round;
            AlignedArray<u8> active_next_round;
            AlignedArray<u8> used_this_round;

            HEIPROMAP_PROFILE_SCOPE("refinement", "FlowBasedRefinement", "allocate");
            active_this_round.initialize(m_k, 1);
            active_next_round.initialize(m_k, 0);
            used_this_round.initialize(m_k * m_k);

            std::vector<std::pair<partition_t, partition_t> > matching;

            for (u64 iteration = 0; iteration < config->max_global_iteration; ++iteration) {
                HEIPROMAP_PROFILE_SCOPE("refinement", "FlowBasedRefinement", "reset_used_edges");
                std::fill_n(used_this_round.get_ptr(), m_k * m_k, 0);

                HEIPROMAP_PROFILE_SCOPE("refinement", "FlowBasedRefinement", "matching");
                bool found_matching = q_graph.find_distance_3_matching(active_this_round, used_this_round, matching);

                if (!found_matching) break;

                while (found_matching) {
                    #pragma omp parallel for num_threads(m_threads) schedule(dynamic)
                    for (size_t i = 0; i < matching.size(); ++i) {
                        partition_t u_id = matching[i].first;
                        partition_t v_id = matching[i].second;

                        u64 thread_id = omp_get_thread_num();
                        refine_blocks<t_uniform_v_weights, t_uniform_e_weights>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, u_id, v_id, lmax_constraints, active_next_round, thread_id, seen_marker_vecs[thread_id], region_marker_vecs[thread_id]);
                    }

                    HEIPROMAP_PROFILE_SCOPE("refinement", "FlowBasedRefinement", "matching");
                    found_matching = q_graph.find_distance_3_matching(active_this_round, used_this_round, matching);
                }

                // swap active
                HEIPROMAP_PROFILE_SCOPE("refinement", "FlowBasedRefinement", "swap_active");
                if (config->use_active_block_scheduling) {
                    std::swap(active_this_round, active_next_round);
                    active_next_round.initialize(m_k, 0);
                }
            }

            if (config->use_closed_vertex_set) {
                // std::cout << "[FlowBasedRefinement] SCC closure search summary: "
                //           << "Successes = " << m_scc_successes << ", "
                //           << "Failures = " << m_scc_failures << std::endl;
            }
        }


        template<bool t_uniform_e_weights>
        weight_t compute_flow_cut_capacity(const graph_t &g,
                                           const p_manager_t &p_manager,
                                           const d_oracle_t &d_oracle,
                                           partition_t left_id,
                                           partition_t right_id,
                                           const std::vector<vertex_t> &left_region,
                                           const std::vector<vertex_t> &right_region,
                                           const std::vector<u8> &is_left,
                                           const TranslationTable<vertex_t> &translation_table,
                                           const AlignedArray<u32> &region_marker,
                                           u32 region_mark) {
            u32 left_mark = region_mark - 1;
            u32 right_mark = region_mark;
            weight_t distance = d_oracle.get(left_id, right_id);
            bool use_penalties = !d_oracle.last_level_pair(left_id, right_id) || !config->use_edge_cut;

            weight_t total_cut = 0;

            auto process_region = [&](const std::vector<vertex_t> &region) {
                for (size_t j = 0; j < region.size(); ++j) {
                    vertex_t u = region[j];
                    int new_u = static_cast<int>(translation_table.get_n(u));
                    bool u_src = (is_left[new_u] == 1);

                    weight_t left_penalty = 0;
                    weight_t right_penalty = 0;

                    for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                        const vertex_t v = g.edges_v[i];
                        const weight_t w = g.edges_w[i];

                        if (region_marker[v] == left_mark || region_marker[v] == right_mark) {
                            if (u < v) {
                                int new_v = static_cast<int>(translation_table.get_n(v));
                                bool v_src = (is_left[new_v] == 1);
                                if (u_src != v_src) {
                                    weight_t cap = use_penalties ? (t_uniform_e_weights ? distance : w * distance)
                                                                 : (t_uniform_e_weights ? 1 : w);
                                    total_cut += cap;
                                }
                            }
                        } else {
                            partition_t v_id = p_manager[v];
                            if (use_penalties) {
                                weight_t dist_l = t_uniform_e_weights ? d_oracle.get(left_id, v_id) : w * d_oracle.get(left_id, v_id);
                                weight_t dist_r = t_uniform_e_weights ? d_oracle.get(right_id, v_id) : w * d_oracle.get(right_id, v_id);
                                left_penalty += dist_l;
                                right_penalty += dist_r;
                            } else {
                                weight_t e_w = t_uniform_e_weights ? 1 : w;
                                if (v_id == left_id) right_penalty += e_w;
                                else if (v_id == right_id) left_penalty += e_w;
                            }
                        }
                    }

                    if (u_src) {
                        total_cut += left_penalty;
                    } else {
                        total_cut += right_penalty;
                    }
                }
            };

            process_region(left_region);
            process_region(right_region);

            return total_cut;
        }

        template<bool t_uniform_e_weights>
        weight_t calculate_gain(const graph_t &g,
                                const p_manager_t &p_manager,
                                const d_oracle_t &d_oracle,
                                const std::vector<u8> &is_left,
                                partition_t left_id,
                                partition_t right_id,
                                const std::vector<vertex_t> &left_region,
                                const std::vector<vertex_t> &right_region,
                                const TranslationTable<vertex_t> &translation_table,
                                const AlignedArray<u32> &region_marker,
                                u32 region_mark) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "FlowBasedRefinement", "calculate_gain");

            u32 left_mark = region_mark - 1;
            u32 right_mark = region_mark;

            weight_t gain = 0;

            auto process_vertex = [&](vertex_t u, partition_t old_id, partition_t new_id) {
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    vertex_t v = g.edges_v[i];
                    weight_t w = t_uniform_e_weights ? 1 : g.edges_w[i];
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
                if (is_left[translation_table.get_n(u)] == 0) {
                    process_vertex(u, left_id, right_id);
                }
            }
            for (vertex_t u : right_region) {
                if (is_left[translation_table.get_n(u)] == 1) {
                    process_vertex(u, right_id, left_id);
                }
            }

            return gain;
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void refine_blocks(graph_t &g,
                           d_oracle_t &d_oracle,
                           bv_manager_t &bv_manager,
                           p_manager_t &p_manager,
                           q_graph_t &q_graph,
                           block_conn_t &block_conn,
                           partition_t left_id,
                           partition_t right_id,
                           const AlignedArray<weight_t> &lmax_constraints,
                           AlignedArray<u8> &active_next_round,
                           u64 thread_id,
                           u32 &seen_mark,
                           u32 &region_mark) {
            ASSERT(left_id != right_id);

            RandomEngine &random_engine = rnd_engines[thread_id];

            f64 alpha = config->alpha;
            f64 alpha_upper_bound = config->alpha_upper_bound;
            f64 alpha_modifier = config->alpha_modifier;

            u64 max_local_iteration = config->max_local_iteration;
            u64 iteration = 0;

            std::vector<vertex_t> &left_boundary = left_boundaries[thread_id];
            std::vector<vertex_t> &right_boundary = right_boundaries[thread_id];
            weight_t left_boundary_weight;
            weight_t right_boundary_weight;

            std::vector<vertex_t> &left_region = left_regions[thread_id];
            std::vector<vertex_t> &right_region = right_regions[thread_id];

            PushRelabel<weight_t> &pr = prs[thread_id];
            MemoryStack &pr_mem = pr_mems[thread_id];
            ResidualFlowNetwork &residual_flow_network = residual_flow_networks[thread_id];
            SCCGraph &scc_graph = scc_graphs[thread_id];

            TranslationTable<vertex_t> &translation_table = translation_tables[thread_id];
            AlignedArray<u32> &seen = seen_vecs[thread_id];
            AlignedArray<u32> &region_marker = region_vecs[thread_id];

            std::vector<vertex_t> &queue = queues[thread_id];
            std::vector<u8> &is_left = is_left_vecs[thread_id];
            std::vector<u8> &is_left_2 = is_left_2_vecs[thread_id];
            std::vector<u8> &s_connected = s_connected_vecs[thread_id];
            std::vector<u8> &t_connected = t_connected_vecs[thread_id];

            f64 avg_weight = (f64) g.g_weight / (f64) m_k;

            while (iteration < max_local_iteration) {
                left_boundary.clear();
                right_boundary.clear();
                left_boundary_weight = 0;
                right_boundary_weight = 0;
                left_region.clear();
                right_region.clear();

                iteration += 1;

                // get boundary vertices
                determine_boundary_vertices<t_uniform_v_weights>(g, bv_manager, p_manager, block_conn, left_id, right_id, left_boundary, right_boundary, left_boundary_weight, right_boundary_weight, random_engine);

                // calc max weight for each bfs
                weight_t adapt_lmax_left = (weight_t) std::ceil(avg_weight + alpha * ((f64) lmax_constraints[left_id] - avg_weight));
                weight_t adapt_lmax_right = (weight_t) std::ceil(avg_weight + alpha * ((f64) lmax_constraints[right_id] - avg_weight));

                weight_t left_max_weight = adapt_lmax_left - p_manager.get_bweight(right_id);
                weight_t right_max_weight = adapt_lmax_right - p_manager.get_bweight(left_id);

                // get both regions
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

                weight_t left_region_weight = determine_region<t_uniform_v_weights>(g, p_manager, left_id, left_mark, left_max_weight, left_boundary, left_region, left_boundary_weight, seen, seen_mark, region_marker, region_mark, queue);
                weight_t right_region_weight = determine_region<t_uniform_v_weights>(g, p_manager, right_id, right_mark, right_max_weight, right_boundary, right_region, right_boundary_weight, seen, seen_mark, region_marker, region_mark, queue);

                if (left_region.size() + right_region.size() == 0) {
                    if (alpha >= alpha_upper_bound) { return; }
                    alpha = std::min(alpha_modifier * alpha, alpha_upper_bound);
                    continue;
                }

                // build a translation table from graph to flow network
                vertex_t new_u = 0;
                for (size_t i = 0; i < left_region.size(); ++i) { translation_table.add(left_region[i], new_u++); }
                for (size_t i = 0; i < right_region.size(); ++i) { translation_table.add(right_region[i], new_u++); }

                if (!d_oracle.last_level_pair(left_id, right_id) || !config->use_edge_cut) {
                    build_flow_network_with_penalties<t_uniform_e_weights>(g, d_oracle, p_manager, left_id, right_id, left_region, right_region, pr, pr_mem, translation_table, region_marker, region_mark, alpha, s_connected, t_connected);
                } else {
                    build_flow_network_no_penalties<t_uniform_e_weights>(g, p_manager, left_id, right_id, left_region, right_region, pr, pr_mem, translation_table, region_marker, region_mark, s_connected, t_connected);
                }

                // solve the flow network
                HEIPROMAP_PROFILE_SCOPE("refinement", "FlowBasedRefinement", "solve_maxflow");
                pr.maxflow();

                // get the cut
                HEIPROMAP_PROFILE_SCOPE("refinement", "FlowBasedRefinement", "get_cut");
                size_t flow_n = left_region.size() + right_region.size();
                is_left.resize(flow_n);
                for (size_t i = 0; i < flow_n; ++i) {
                    is_left[i] = (pr.what_label(static_cast<int>(i)) == SOURCE) ? 1 : 0;
                }

                // check if cut is valid
                HEIPROMAP_PROFILE_SCOPE("refinement", "FlowBasedRefinement", "cut_is_valid");
                bool is_valid = cut_is_valid<t_uniform_v_weights>(g, p_manager, left_id, right_id, is_left, lmax_constraints, left_region, right_region, translation_table);

                // not valid and no chance for improvement -> cancel
                if (!is_valid && !config->use_closed_vertex_set) {
                    if (alpha <= 1.0) { return; }
                    alpha = std::max(alpha / alpha_modifier, 1.0);
                    continue;
                }

                // not valid but no chance for improvement
                if (!is_valid) {
                    // build residual network from cut labels

                    HEIPROMAP_PROFILE_SCOPE("refinement", "FlowBasedRefinement", "build_residual_network");
                    build_residual_from_cut(g, pr, left_region, right_region, translation_table, region_marker, region_mark, s_connected, t_connected, residual_flow_network);
                    residual_flow_network.finalize();

                    // build scc graph
                    scc_graph.initialize(residual_flow_network, g, translation_table);

                    // reduce the scc graph
                    HEIPROMAP_PROFILE_SCOPE("refinement", "FlowBasedRefinement", "reduce_scc");
                    scc_graph.reduce();

                    // determine best balanced min cut
                    HEIPROMAP_PROFILE_SCOPE("refinement", "FlowBasedRefinement", "scc_find");
                    weight_t left_non_region_weight = p_manager.get_bweight(left_id) - left_region_weight;
                    weight_t right_non_region_weight = p_manager.get_bweight(right_id) - right_region_weight;
                    bool closure_found = scc_graph.find_best_closure(left_non_region_weight, right_non_region_weight, lmax_constraints[left_id], lmax_constraints[right_id], avg_weight, config->closed_vertex_sets_repeats, random_engine, is_left_2);

                    if (closure_found) {
                        #pragma omp atomic
                        m_scc_successes++;

                        // SANITY CHECK 1: Verify balance constraints
                        bool closure_is_valid = cut_is_valid<t_uniform_v_weights>(g, p_manager, left_id, right_id, is_left_2, lmax_constraints, left_region, right_region, translation_table);
                        ASSERT(closure_is_valid);

                        weight_t mincut_gain = calculate_gain<t_uniform_e_weights>(g, p_manager, d_oracle, is_left, left_id, right_id, left_region, right_region, translation_table, region_marker, region_mark);
                        weight_t closure_gain = calculate_gain<t_uniform_e_weights>(g, p_manager, d_oracle, is_left_2, left_id, right_id, left_region, right_region, translation_table, region_marker, region_mark);

                        if (closure_gain >= mincut_gain) {
                            std::swap(is_left, is_left_2);
                        } else {
                            closure_found = false;
                        }
                    } else {
                        #pragma omp atomic
                        m_scc_failures++;
                    }

                    // still not valid or did not improve gain
                    if (!closure_found) {
                        if (alpha <= 1.0) { return; }
                        alpha = std::max(alpha / alpha_modifier, 1.0);
                        continue;
                    }
                }

                // check if the cut actually changes the partition
                if (!cut_changes_partition(is_left, left_region, right_region, translation_table)) {
                    if (alpha <= 1.0) { return; }
                    alpha = std::max(alpha / alpha_modifier, 1.0);
                    continue;
                }

                // check if cut improves qap
                // weight_t gain = calculate_gain<t_uniform_e_weights>(g, p_manager, d_oracle, is_left, left_id, right_id, left_region, right_region, translation_table, region_marker, region_mark);

                // if (gain == 0) {
                //     // cut changes the partition but has no better qap, so search again but without increasing alpha
                //     continue;
                // }

                // cut is valid, positive and changes the partition, increase alpha
                alpha = std::min(alpha * alpha_modifier, alpha_upper_bound);

                // make the changes
                change_boundary<t_uniform_v_weights>(g, bv_manager, p_manager, q_graph, block_conn, is_left, left_id, right_id, left_region, right_region, translation_table, region_marker, region_mark);

                active_next_round[left_id] = 1;
                active_next_round[right_id] = 1;
            }
        }

        template<bool t_uniform_v_weights>
        void determine_boundary_vertices(const graph_t &g,
                                         const bv_manager_t &bv_manager,
                                         const p_manager_t &p_manager,
                                         const block_conn_t &block_conn,
                                         partition_t left_id,
                                         partition_t right_id,
                                         std::vector<vertex_t> &left_boundary,
                                         std::vector<vertex_t> &right_boundary,
                                         weight_t &left_boundary_weight,
                                         weight_t &right_boundary_weight,
                                         RandomEngine &random_engine) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "FlowBasedRefinement", "determine_boundary_vertices");

            left_boundary_weight = 0;
            for (size_t i = 0; i < bv_manager.size(left_id); ++i) {
                const vertex_t u = bv_manager.get(left_id, i);
                for (size_t j = block_conn.start(u); j < block_conn.end(u); ++j) {
                    const partition_t id = block_conn.get_id(j);
                    if (id == right_id) {
                        left_boundary.push_back(u);
                        left_boundary_weight += t_uniform_v_weights ? 1 : g.v_weights[u];
                        break;
                    }
                }
            }
            fast_shuffle_unchecked(left_boundary.data(), left_boundary.data() + left_boundary.size(), random_engine.generator);

            right_boundary_weight = 0;
            for (size_t i = 0; i < bv_manager.size(right_id); ++i) {
                const vertex_t u = bv_manager.get(right_id, i); {
                    for (size_t j = block_conn.start(u); j < block_conn.end(u); ++j) {
                        const partition_t id = block_conn.get_id(j); {
                            if (id == left_id) {
                                right_boundary.push_back(u);
                                right_boundary_weight += t_uniform_v_weights ? 1 : g.v_weights[u];
                                break;
                            }
                        }
                    }
                }
            }
            fast_shuffle_unchecked(right_boundary.data(), right_boundary.data() + right_boundary.size(), random_engine.generator);
        }

        template<bool t_uniform_v_weights>
        weight_t determine_region(const graph_t &g,
                                  const p_manager_t &p_manager,
                                  partition_t id,
                                  u32 mark,
                                  weight_t max_weight,
                                  const std::vector<vertex_t> &boundary,
                                  std::vector<vertex_t> &region,
                                  [[maybe_unused]] weight_t boundary_weight,
                                  AlignedArray<u32> &seen,
                                  u32 &seen_mark,
                                  AlignedArray<u32> &region_marker,
                                  [[maybe_unused]] u32 &region_mark_ref,
                                  std::vector<vertex_t> &queue) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "FlowBasedRefinement", "determine_regions");

            seen_mark += 2;
            // seen[u] == seen_mark     means u is processed
            // seen[u] == seen_mark - 1 means u is in the queue

            weight_t curr_weight = 0;

            queue.clear();
            for (size_t i = 0; i < boundary.size(); ++i) {
                vertex_t u = boundary[i];
                queue.push_back(u);
                seen[u] = seen_mark - 1;
            }

            size_t boundary_size = queue.size();
            region.clear();

            if (config->growth_strategy == GrowthStrategy::BFS) {
                size_t queue_idx = 0;
                while (queue_idx < queue.size()) {
                    vertex_t u = queue[queue_idx++];
                    if (seen[u] == seen_mark) { continue; }

                    weight_t u_w = t_uniform_v_weights ? 1 : g.v_weights[u];
                    bool include = (curr_weight + u_w <= max_weight);
                    if (config->always_include_boundary && queue_idx <= boundary_size) {
                        include = true;
                    }

                    if (include) {
                        region.push_back(u);
                        region_marker[u] = mark;
                        curr_weight += u_w;
                        for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                            const vertex_t v = g.edges_v[i];
                            partition_t v_id = p_manager[v];
                            if (v_id != id) { continue; }

                            if (seen[v] != seen_mark && seen[v] != seen_mark - 1) {
                                queue.push_back(v);
                                seen[v] = seen_mark - 1;
                            }
                        }
                    }
                    seen[u] = seen_mark;
                }
            } else {
                auto heavy_cmp = [&](vertex_t lhs, vertex_t rhs) {
                    weight_t wl = t_uniform_v_weights ? 1 : g.v_weights[lhs];
                    weight_t wr = t_uniform_v_weights ? 1 : g.v_weights[rhs];
                    return wl < wr;
                };
                auto light_cmp = [&](vertex_t lhs, vertex_t rhs) {
                    weight_t wl = t_uniform_v_weights ? 1 : g.v_weights[lhs];
                    weight_t wr = t_uniform_v_weights ? 1 : g.v_weights[rhs];
                    return wl > wr;
                };

                size_t queue_idx = 0;
                if (config->always_include_boundary) {
                    while (queue_idx < boundary_size) {
                        vertex_t u = queue[queue_idx++];
                        if (seen[u] == seen_mark) continue;

                        region.push_back(u);
                        region_marker[u] = mark;
                        curr_weight += t_uniform_v_weights ? 1 : g.v_weights[u];

                        for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                            vertex_t v = g.edges_v[i];
                            if (p_manager[v] == id && seen[v] != seen_mark && seen[v] != seen_mark - 1) {
                                queue.push_back(v);
                                seen[v] = seen_mark - 1;
                            }
                        }
                        seen[u] = seen_mark;
                    }
                }

                if (queue_idx < queue.size()) {
                    if (config->growth_strategy == GrowthStrategy::HEAVY_FIRST)
                        std::make_heap(queue.begin() + queue_idx, queue.end(), heavy_cmp);
                    else
                        std::make_heap(queue.begin() + queue_idx, queue.end(), light_cmp);
                }

                while (queue_idx < queue.size() && curr_weight < max_weight) {
                    if (config->growth_strategy == GrowthStrategy::HEAVY_FIRST)
                        std::pop_heap(queue.begin() + queue_idx, queue.end(), heavy_cmp);
                    else
                        std::pop_heap(queue.begin() + queue_idx, queue.end(), light_cmp);

                    vertex_t u = queue.back();
                    queue.pop_back();

                    if (seen[u] == seen_mark) continue;

                    weight_t u_w = t_uniform_v_weights ? 1 : g.v_weights[u];
                    if (curr_weight + u_w <= max_weight) {
                        region.push_back(u);
                        region_marker[u] = mark;
                        curr_weight += u_w;
                        for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                            const vertex_t v = g.edges_v[i];
                            if (p_manager[v] == id && seen[v] != seen_mark && seen[v] != seen_mark - 1) {
                                queue.push_back(v);
                                if (config->growth_strategy == GrowthStrategy::HEAVY_FIRST)
                                    std::push_heap(queue.begin() + queue_idx, queue.end(), heavy_cmp);
                                else
                                    std::push_heap(queue.begin() + queue_idx, queue.end(), light_cmp);
                                seen[v] = seen_mark - 1;
                            }
                        }
                    }
                    seen[u] = seen_mark;
                }
            }

            return curr_weight;
        }

        template<bool t_uniform_e_weights>
        void build_flow_network_with_penalties(const graph_t &g,
                                               const d_oracle_t &d_oracle,
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
                                               f64 alpha,
                                               std::vector<u8> &s_connected,
                                               std::vector<u8> &t_connected) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "FlowBasedRefinement", "build_flow_network_with_penalties");

            u32 left_mark = region_mark - 1;
            u32 right_mark = region_mark;

            weight_t distance = d_oracle.get(left_id, right_id);

            int n = static_cast<int>(left_region.size() + right_region.size());
            int src = n, snk = n + 1;

            // loose upper bound: sum of degrees + 2 terminal edges per vertex
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

            // build left region: left-left, left-right, and penalties
            for (size_t j = 0; j < left_region.size(); ++j) {
                vertex_t u = left_region[j];
                int new_u = static_cast<int>(translation_table.get_n(u));
                weight_t left_penalty = 0;
                weight_t right_penalty = 0;

                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    const vertex_t v = g.edges_v[i];
                    const weight_t w = g.edges_w[i]; {
                        if (region_marker[v] == right_mark) {
                            int new_v = static_cast<int>(translation_table.get_n(v));
                            weight_t cap = t_uniform_e_weights ? distance : w * distance;
                            pr.add_edge(new_u, new_v, cap, cap);
                            continue;
                        }

                        if (region_marker[v] == left_mark) {
                            if (u < v) { continue; }
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

            // build right region: only right-right and penalties
            for (size_t j = 0; j < right_region.size(); ++j) {
                vertex_t u = right_region[j];
                int new_u = static_cast<int>(translation_table.get_n(u));
                weight_t left_penalty = 0;
                weight_t right_penalty = 0;

                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    const vertex_t v = g.edges_v[i];
                    const weight_t w = g.edges_w[i]; {
                        if (region_marker[v] == right_mark) {
                            if (u < v) { continue; }
                            int new_v = static_cast<int>(translation_table.get_n(v));
                            weight_t cap = t_uniform_e_weights ? distance : w * distance;
                            pr.add_edge(new_u, new_v, cap, cap);
                            continue;
                        }

                        if (region_marker[v] == left_mark) {
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
            HEIPROMAP_PROFILE_SCOPE("refinement", "FlowBasedRefinement", "build_flow_network_no_penalties");

            u32 left_mark = region_mark - 1;
            u32 right_mark = region_mark;

            int n = static_cast<int>(left_region.size() + right_region.size());
            int src = n, snk = n + 1;

            // loose upper bound: sum of degrees + 2 terminal edges per vertex
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

            // build left region
            for (size_t j = 0; j < left_region.size(); ++j) {
                vertex_t u = left_region[j];
                int new_u = static_cast<int>(translation_table.get_n(u));
                weight_t w_left = 0;
                weight_t w_right = 0;

                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    const vertex_t v = g.edges_v[i];
                    const weight_t w = g.edges_w[i]; {
                        if (region_marker[v] == left_mark || region_marker[v] == right_mark) {
                            if (u < v) { continue; }
                            int new_v = static_cast<int>(translation_table.get_n(v));
                            weight_t cap = t_uniform_e_weights ? 1 : w;
                            pr.add_edge(new_u, new_v, cap, cap);
                            continue;
                        }

                        partition_t v_id = p_manager[v];
                        if (v_id == left_id) { w_left += t_uniform_e_weights ? 1 : w; } else if (v_id == right_id) { w_right += t_uniform_e_weights ? 1 : w; }
                    }
                }

                if (w_left > 0) {
                    pr.add_edge(src, new_u, w_left);
                    s_connected[new_u] = 1;
                }
                if (w_right > 0) {
                    pr.add_edge(new_u, snk, w_right);
                    t_connected[new_u] = 1;
                }
            }

            // build right region
            for (size_t j = 0; j < right_region.size(); ++j) {
                vertex_t u = right_region[j];
                int new_u = static_cast<int>(translation_table.get_n(u));
                weight_t w_left = 0;
                weight_t w_right = 0;

                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    const vertex_t v = g.edges_v[i];
                    const weight_t w = g.edges_w[i]; {
                        if (region_marker[v] == left_mark || region_marker[v] == right_mark) {
                            if (u < v) { continue; }
                            int new_v = static_cast<int>(translation_table.get_n(v));
                            weight_t cap = t_uniform_e_weights ? 1 : w;
                            pr.add_edge(new_u, new_v, cap, cap);
                            continue;
                        }

                        partition_t v_id = p_manager[v];
                        if (v_id == left_id) { w_left += t_uniform_e_weights ? 1 : w; } else if (v_id == right_id) { w_right += t_uniform_e_weights ? 1 : w; }
                    }
                }

                if (w_left > 0) {
                    pr.add_edge(src, new_u, w_left);
                    s_connected[new_u] = 1;
                }
                if (w_right > 0) {
                    pr.add_edge(new_u, snk, w_right);
                    t_connected[new_u] = 1;
                }
            }
        }

        void build_residual_from_cut(const graph_t &g,
                                     PushRelabel<weight_t> &pr,
                                     std::vector<vertex_t> &left_region,
                                     std::vector<vertex_t> &right_region,
                                     TranslationTable<vertex_t> &translation_table,
                                     AlignedArray<u32> &region_marker,
                                     u32 &region_mark,
                                     std::vector<u8> &s_connected,
                                     std::vector<u8> &t_connected,
                                     ResidualFlowNetwork &residual_g) {
            vertex_t flow_n = left_region.size() + right_region.size();
            residual_g.initialize(flow_n);

            int pr_n = static_cast<int>(flow_n + 2);
            for (int u = 0; u < pr_n; ++u) {
                for (int i = pr.arc_first(u); i < pr.arc_last(u); ++i) {
                    if (pr.arc_cap(i) > 0) {
                        int v = pr.arc_to(i);
                        residual_g.add_directed_edge(static_cast<vertex_t>(u), static_cast<vertex_t>(v), 1);
                    }
                }
            }
        }

        template<bool t_uniform_v_weights>
        bool cut_is_valid(const graph_t &g,
                          const p_manager_t &p_manager,
                          partition_t left_id,
                          partition_t right_id,
                          std::vector<u8> &is_left,
                          const AlignedArray<weight_t> &lmax_constraints,
                          std::vector<vertex_t> &left_region,
                          std::vector<vertex_t> &right_region,
                          TranslationTable<vertex_t> &translation_table) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "FlowBasedRefinement", "cut_is_valid");

            weight_t left_weight = p_manager.get_bweight(left_id);
            weight_t right_weight = p_manager.get_bweight(right_id);
            for (size_t j = 0; j < left_region.size(); ++j) {
                vertex_t u = left_region[j];
                weight_t u_weight = t_uniform_v_weights ? 1 : g.v_weights[u];
                vertex_t new_u = translation_table.get_n(u);
                if (is_left[new_u] == 0) {
                    left_weight -= u_weight;
                    right_weight += u_weight;
                }
            }

            for (size_t j = 0; j < right_region.size(); ++j) {
                vertex_t u = right_region[j];
                weight_t u_weight = t_uniform_v_weights ? 1 : g.v_weights[u];
                vertex_t new_u = translation_table.get_n(u);
                if (is_left[new_u] == 1) {
                    right_weight -= u_weight;
                    left_weight += u_weight;
                }
            }

            return left_weight <= lmax_constraints[left_id] && right_weight <= lmax_constraints[right_id];
        }

        bool cut_changes_partition(std::vector<u8> &is_left,
                                   std::vector<vertex_t> &left_region,
                                   std::vector<vertex_t> &right_region,
                                   TranslationTable<vertex_t> &translation_table) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "FlowBasedRefinement", "cut_changes_partition");

            for (size_t j = 0; j < left_region.size(); ++j) {
                vertex_t u = left_region[j];
                vertex_t new_u = translation_table.get_n(u);
                if (is_left[new_u] == 0) {
                    return true;
                }
            }
            for (size_t j = 0; j < right_region.size(); ++j) {
                vertex_t u = right_region[j];
                vertex_t new_u = translation_table.get_n(u);
                if (is_left[new_u] == 1) {
                    return true;
                }
            }

            return false;
        }

        template<bool t_uniform_v_weights>
        std::vector<u8> change_boundary(const graph_t &g,
                                        bv_manager_t &bv_manager,
                                        p_manager_t &p_manager,
                                        q_graph_t &q_graph,
                                        block_conn_t &block_conn,
                                        std::vector<u8> &is_left,
                                        partition_t left_id,
                                        partition_t right_id,
                                        std::vector<vertex_t> &left_region,
                                        std::vector<vertex_t> &right_region,
                                        TranslationTable<vertex_t> &translation_table,
                                        AlignedArray<u32> &region_marker,
                                        u32 &region_mark) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "FlowBasedRefinement", "change_boundary");

            [[maybe_unused]] u32 left_mark = region_mark - 1;
            [[maybe_unused]] u32 right_mark = region_mark;

            std::vector<u8> changed(is_left.size(), 0);

            for (size_t j = 0; j < left_region.size(); ++j) {
                vertex_t u = left_region[j];
                vertex_t new_u = translation_table.get_n(u);

                ASSERT(region_marker[u] == left_mark);
                ASSERT(new_u < left_region.size() + right_region.size());

                if (is_left[new_u] == 0) {
                    changed[new_u] = 1;
                    if (bv_manager.is_boundary(u)) {
                        bv_manager.move(g, p_manager, u, left_id, right_id);
                    } else {
                        bv_manager.add_new(g, p_manager, u, right_id);
                    }

                    q_graph.move(g, p_manager, u, left_id, right_id);
                    block_conn.move(g, u, left_id, right_id);
                    p_manager.move(u, t_uniform_v_weights ? 1 : g.v_weights[u], left_id, right_id);
                }
            }
            for (size_t j = 0; j < right_region.size(); ++j) {
                vertex_t u = right_region[j];
                vertex_t new_u = translation_table.get_n(u);

                ASSERT(region_marker[u] == right_mark);
                ASSERT(new_u < left_region.size() + right_region.size());

                if (is_left[new_u] == 1) {
                    changed[new_u] = 1;
                    if (bv_manager.is_boundary(u)) {
                        bv_manager.move(g, p_manager, u, right_id, left_id);
                    } else {
                        bv_manager.add_new(g, p_manager, u, left_id);
                    }

                    q_graph.move(g, p_manager, u, right_id, left_id);
                    block_conn.move(g, u, right_id, left_id);
                    p_manager.move(u, t_uniform_v_weights ? 1 : g.v_weights[u], right_id, left_id);
                }
            }
            return changed;
        }

        template<bool t_uniform_v_weights>
        void revert_boundary(const graph_t &g,
                             bv_manager_t &bv_manager,
                             p_manager_t &p_manager,
                             q_graph_t &q_graph,
                             block_conn_t &block_conn,
                             std::vector<u8> &changed,
                             partition_t left_id,
                             partition_t right_id,
                             std::vector<vertex_t> &left_region,
                             std::vector<vertex_t> &right_region,
                             TranslationTable<vertex_t> &translation_table) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "FlowBasedRefinement", "revert_boundary");

            for (vertex_t new_u = 0; new_u < left_region.size() + right_region.size(); ++new_u) {
                if (changed[new_u] == 0) { continue; }

                vertex_t old_u = translation_table.get_o(new_u);
                if (p_manager[old_u] == left_id) {
                    if (bv_manager.is_boundary(old_u)) {
                        bv_manager.move(g, p_manager, old_u, left_id, right_id);
                    } else {
                        bv_manager.add_new(g, p_manager, old_u, right_id);
                    }

                    q_graph.move(g, p_manager, old_u, left_id, right_id);
                    block_conn.move(g, old_u, left_id, right_id);
                    p_manager.move(old_u, t_uniform_v_weights ? 1 : g.v_weights[old_u], left_id, right_id);
                } else {
                    if (bv_manager.is_boundary(old_u)) {
                        bv_manager.move(g, p_manager, old_u, right_id, left_id);
                    } else {
                        bv_manager.add_new(g, p_manager, old_u, left_id);
                    }

                    q_graph.move(g, p_manager, old_u, right_id, left_id);
                    block_conn.move(g, old_u, right_id, left_id);
                    p_manager.move(old_u, t_uniform_v_weights ? 1 : g.v_weights[old_u], right_id, left_id);
                }
            }
        }
    };
}

#endif //HEIPROMAP_FLOW_BASED_REFINEMENT_H
