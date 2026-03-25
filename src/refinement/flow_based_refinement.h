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
#include <stack>
#include <unordered_set>

#include <omp.h>


#include "ISerialRefiner.h"
#include "quotient_graph_refinement.h"
#include "../utility/flow.h"
#include "../utility/flow_algorithm_adapters.h"
#include "../utility/random_engine.h"
#include "../utility/utils.h"
#include "../utility/assert_state.h"
#include "../utility/qap.h"

namespace HeiProMap {
    class FlowBasedRefinementConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit FlowBasedRefinementConfiguration(const std::string &t_name) : ISerialRefinerConfiguration(t_name) {
        }

        u64 max_global_iteration = 1;
        u64 max_local_iteration = 3;
        f64 alpha = 2.0;
        f64 alpha_upper_bound = 8.0;
        f64 alpha_modifier = 2.0;
        bool use_closed_vertex_set = true;
        u64 closed_vertex_sets_repeats = 10;
    };

    class FlowBasedRefinement final : public ISerialRefiner {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
        u64 m_threads = 1;

        AlignedArray<weight_t> left_penalties;
        AlignedArray<weight_t> right_penalties;

        const FlowBasedRefinementConfiguration *config = nullptr;

    public:
        FlowBasedRefinement() = default;

        ~FlowBasedRefinement() override = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const u64 t_threads,
                        const ISerialRefinerConfiguration &i_config) override {
            ScopedTimer _t("io", "FlowBasedRefinement", "initialize");

            m_n = t_n;
            m_m = t_m;
            m_k = t_k;
            m_threads = t_threads;

            left_penalties.initialize(m_n);
            right_penalties.initialize(m_n);

            config = dynamic_cast<const FlowBasedRefinementConfiguration *>(&i_config);
        }

        void refine(graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph,
                    block_conn_t &block_conn,
                    f64 imbalance) override {
            RandomEngine random_engine = RandomEngine(0);

            // active block scheduling
            AlignedArray<u8> active_this_round;
            AlignedArray<u8> active_next_round;

            //
            {
                ScopedTimer _t("refinement", "FlowBasedRefinement", "allocate");

                active_this_round.initialize(m_k, 1);
                active_next_round.initialize(m_k, 0);
            }

            AlignedArray<u8> used_this_round;
            //
            {
                ScopedTimer _t("refinement", "FlowBasedRefinement", "allocate");

                used_this_round.initialize(m_k * m_k);
            }

            std::vector<std::pair<partition_t, partition_t> > matching;

            for (u64 iteration = 0; iteration < config->max_global_iteration; ++iteration) {
                //
                {
                    ScopedTimer _t("refinement", "FlowBasedRefinement", "reset_used_edges");

                    std::fill_n(used_this_round.get_ptr(), m_k * m_k, 0);
                }

                bool found_matching = false;
                //
                {
                    ScopedTimer _t("refinement", "FlowBasedRefinement", "matching");
                    found_matching = q_graph.find_distance_3_matching(active_this_round, used_this_round, matching);
                }

                while (found_matching) {
#pragma omp parallel for num_threads(m_threads) schedule(dynamic)
                    for (size_t i = 0; i < matching.size(); ++i) {
                        partition_t u_id = matching[i].first;
                        partition_t v_id = matching[i].second;

                        u64 thread_id = omp_get_thread_num();
                        refine_blocks(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, u_id, v_id, imbalance, active_next_round);
                    }

                    //
                    {
                        ScopedTimer _t("refinement", "FlowBasedRefinement", "matching");
                        found_matching = q_graph.find_distance_3_matching(active_this_round, used_this_round, matching);
                    }
                }

                // swap active
                {
                    ScopedTimer _t("refinement", "FlowBasedRefinement", "swap_active");
                    std::swap(active_this_round, active_next_round);
                    active_next_round.initialize(m_k, 0);
                }
            }
        }

        void refine_blocks(graph_t &g,
                           d_oracle_t &d_oracle,
                           bv_manager_t &bv_manager,
                           p_manager_t &p_manager,
                           q_graph_t &q_graph,
                           block_conn_t &block_conn,
                           partition_t left_id,
                           partition_t right_id,
                           f64 imbalance,
                           AlignedArray<u8> &active_next_round) {
            ASSERT(left_id != right_id);

            RandomEngine random_engine = RandomEngine(left_id * right_id + g.n);

            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) m_k));

            f64 alpha = config->alpha;
            f64 alpha_upper_bound = config->alpha_upper_bound;
            f64 alpha_modifier = config->alpha_modifier;

            u64 max_local_iteration = config->max_local_iteration;
            u64 iteration = 0;

            std::vector<vertex_t> left_boundary;
            std::vector<vertex_t> right_boundary;
            weight_t left_boundary_weight;
            weight_t right_boundary_weight;

            std::vector<vertex_t> left_region;
            std::vector<vertex_t> right_region;

            BKAdapter<int, int, int> flow_network;
            // IBFSAdapter<int, int, int> flow_network;
            // IBFSAdapter<int, int, int> flow_network;
            // HiPrAdapter<int, int, int> flow_network;
            ResidualFlowNetwork residual_flow_network;
            SCCGraph scc_graph;

            TranslationTable<vertex_t> translation_table;
            translation_table.reserve(g.n, g.n);

            std::vector<u32> seen(g.n, 0);
            u32 seen_mark = 1;

            std::vector<u32> region_marker(g.n, 0);
            u32 region_mark = 1;

            while (iteration < max_local_iteration) {
                left_boundary.clear();
                right_boundary.clear();
                left_boundary_weight = 0;
                right_boundary_weight = 0;
                left_region.clear();
                right_region.clear();

                iteration += 1;

                // get boundary vertices
                determine_boundary_vertices(g, bv_manager, p_manager, block_conn, left_id, right_id, left_boundary, right_boundary, left_boundary_weight, right_boundary_weight, random_engine);

                // calc max weight for each bfs
                weight_t adapt_lmax = std::ceil((1.0 + (imbalance * alpha)) * ((f64) g.g_weight / (f64) m_k));
                weight_t left_max_weight = adapt_lmax - p_manager.get_bweight(right_id);
                weight_t right_max_weight = adapt_lmax - p_manager.get_bweight(left_id);

                // get both regions
                region_mark += 2;
                u32 left_mark = region_mark - 1;
                u32 right_mark = region_mark;
                seen_mark += 2;
                vertex_t left_seed_vertex = left_boundary[0];
                vertex_t right_seed_vertex = right_boundary[0];
                forall_guiv(g, left_seed_vertex, i, v)
                    {
                        if (p_manager[v] == right_id) {
                            right_seed_vertex = v;
                            break;
                        }
                    }
                endfor
                weight_t left_region_weight = determine_region(g, p_manager, left_id, left_mark, left_max_weight, left_boundary, left_region, left_seed_vertex, left_boundary_weight, seen, seen_mark, region_marker, region_mark);
                weight_t right_region_weight = determine_region(g, p_manager, right_id, right_mark, right_max_weight, right_boundary, right_region, right_seed_vertex, right_boundary_weight, seen, seen_mark, region_marker, region_mark);

                if (left_region.size() + right_region.size() == 0) {
                    // if both regions are empty, increase their sizes
                    if (alpha == alpha_upper_bound) { return; }
                    alpha = std::min(alpha_modifier * alpha, alpha_upper_bound);
                    continue;
                }

                // determine penalties for all vertices
                determine_penalties(g, p_manager, d_oracle, left_id, right_id, left_region, right_region, region_marker, region_mark);

                // build a translation table from graph to flow network
                vertex_t new_u = 0;
                for (size_t i = 0; i < left_region.size(); ++i) { translation_table.add(left_region[i], new_u++); }
                for (size_t i = 0; i < right_region.size(); ++i) { translation_table.add(right_region[i], new_u++); }

                // build flownetwork
                build_flow_network(g, d_oracle, left_id, right_id, left_region, right_region, flow_network, translation_table, region_marker, region_mark);

                // solve the flow network
                {
                    ScopedTimer _t("refinement", "FlowBasedRefinement", "solve_flow_network");
                    flow_network.solve();
                }

                std::vector<u8> is_left;
                // get the cut
                {
                    ScopedTimer _t("refinement", "FlowBasedRefinement", "get_cut");
                    flow_network.get_cut(is_left);
                }
                bool is_valid;
                // check if cut is valid
                {
                    ScopedTimer _t("refinement", "FlowBasedRefinement", "cut_is_valid");
                    is_valid = cut_is_valid(g, p_manager, left_id, right_id, is_left, lmax, left_region, right_region, translation_table);
                }

                if (!is_valid) {
                    if (!config->use_closed_vertex_set) {
                        if (alpha == 1.0) { return; }
                        if (alpha == alpha_upper_bound) { return; }
                        alpha = std::max(alpha / alpha_modifier, 1.0);
                        continue;
                    }

                    // build residual network
                    {
                        ScopedTimer _t("refinement", "FlowBasedRefinement", "build_residual_network");
                        flow_network.build_residual_network(residual_flow_network);
                    }
                    // build scc graph
                    {
                        ScopedTimer _t("refinement", "FlowBasedRefinement", "build_scc");
                        scc_graph.initialize(residual_flow_network, g, translation_table);
                    }
                    // reduce the scc graph
                    {
                        ScopedTimer _t("refinement", "FlowBasedRefinement", "reduce_scc");
                        scc_graph.reduce();
                    }
                    bool closure_found;
                    // determine best balanced min cut
                    {
                        ScopedTimer _t("refinement", "FlowBasedRefinement", "scc_find");
                        weight_t left_non_region_weight = p_manager.get_bweight(left_id) - left_region_weight;
                        weight_t right_non_region_weight = p_manager.get_bweight(right_id) - right_region_weight;
                        closure_found = scc_graph.find_best_closure(left_non_region_weight, right_non_region_weight, lmax, lmax, config->closed_vertex_sets_repeats, random_engine, is_left);
                    }
                    //
                    if (!closure_found) {
                        if (alpha == 1.0) { return; }
                        alpha = std::max(alpha / alpha_modifier, 1.0);
                        continue;
                    }
                }

                // check if the cut actually changes the partition
                if (!cut_changes_partition(is_left, left_region, right_region, translation_table)) {
                    // cut is valid, but does not change anything
                    if (alpha == 1.0) { return; }
                    if (alpha == alpha_upper_bound) { return; }
                    alpha = std::max(alpha / alpha_modifier, 1.0);
                    continue;
                }
                // cut is valid and changes the partition, increase alpha
                alpha = std::min(alpha * alpha_modifier, alpha_upper_bound);

                // make the changes
                change_boundary(g, bv_manager, p_manager, q_graph, block_conn, is_left, left_id, right_id, left_region, right_region, translation_table, region_marker, region_mark);

                active_next_round[left_id] = 1;
                active_next_round[right_id] = 1;
            }
        }

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
            ScopedTimer _t("refinement", "FlowBasedRefinement", "determine_boundary_vertices");

            left_boundary_weight = 0;
            forall_bv_id_iu(bv_manager, left_id, i, u)
                {
                    forall_bc_ui_id(block_conn, u, j, id)
                        {
                            if (id == right_id) {
                                left_boundary.push_back(u);
                                left_boundary_weight += g.v_weights[u];
                                break;
                            }
                        }
                    endfor
                }
            endfor
            fast_shuffle_unchecked(left_boundary.data(), left_boundary.data() + left_boundary.size(), random_engine.generator);

            right_boundary_weight = 0;
            forall_bv_id_iu(bv_manager, right_id, i, u)
                {
                    forall_bc_ui_id(block_conn, u, j, id)
                        {
                            if (id == left_id) {
                                right_boundary.push_back(u);
                                right_boundary_weight += g.v_weights[u];
                                break;
                            }
                        }
                    endfor
                }
            endfor
            fast_shuffle_unchecked(right_boundary.data(), right_boundary.data() + right_boundary.size(), random_engine.generator);
        }

        weight_t determine_region(const graph_t &g,
                                  const p_manager_t &p_manager,
                                  partition_t id,
                                  u32 mark,
                                  weight_t max_weight,
                                  std::vector<vertex_t> &boundary,
                                  std::vector<vertex_t> &region,
                                  vertex_t seed_vertex,
                                  weight_t boundary_weight,
                                  std::vector<u32> &seen,
                                  u32 &seen_mark,
                                  std::vector<u32> &region_marker,
                                  u32 &region_mark) {
            ScopedTimer _t_sort_pairs("refinement", "FlowBasedRefinement", "determine_regions");

            seen_mark += 2;
            // seen[u] == seen_mark     means u is processed
            // seen[u] == seen_mark - 1 means u is in the queue

            std::vector<vertex_t> queue;
            weight_t curr_weight = 0;

            if (boundary_weight <= max_weight || true) {
                for (size_t i = 0; i < boundary.size(); ++i) {
                    vertex_t u = boundary[i];
                    queue.push_back(u);
                    seen[u] = seen_mark - 1;
                }
            } else {
                queue.push_back(seed_vertex);
                seen[seed_vertex] = seen_mark - 1;
            }

            size_t queue_idx = 0;
            region.clear();
            while (queue_idx < queue.size()) {
                vertex_t u = queue[queue_idx++];
                if (seen[u] == seen_mark) { continue; }

                if (curr_weight + g.v_weights[u] <= max_weight) {
                    region.push_back(u);
                    region_marker[u] = mark;
                    curr_weight += g.v_weights[u];
                    forall_guiv(g, u, i, v)
                        {
                            partition_t v_id = p_manager[v];
                            if (v_id != id) { continue; }

                            if (seen[v] != seen_mark && seen[v] != seen_mark - 1) {
                                queue.push_back(v);
                                seen[v] = seen_mark - 1;
                            }
                        }
                    endfor
                }
                seen[u] = seen_mark;
            }
            return curr_weight;
        }

        void determine_penalties(const graph_t &g,
                                 const p_manager_t &p_manager,
                                 d_oracle_t &d_oracle,
                                 partition_t left_id,
                                 partition_t right_id,
                                 std::vector<vertex_t> &left_region,
                                 std::vector<vertex_t> &right_region,
                                 std::vector<u32> &region_marker,
                                 u32 &region_mark) {
            ScopedTimer _t("refinement", "FlowBasedRefinement", "determine_penalties");

            u32 left_mark = region_mark - 1;
            u32 right_mark = region_mark;

            for (size_t j = 0; j < left_region.size(); ++j) {
                vertex_t u = left_region[j];
                left_penalties[u] = 0;
                right_penalties[u] = 0;
                forall_guivw(g, u, i, v, w)
                    {
                        if (region_marker[v] == left_mark || region_marker[v] == right_mark) { continue; } // ignore neighbors that are in the region, they will be handled later
                        partition_t v_id = p_manager[v];

                        left_penalties[u] += w * d_oracle.get(left_id, v_id);
                        right_penalties[u] += w * d_oracle.get(right_id, v_id);
                    }
                endfor
            }
            for (size_t j = 0; j < right_region.size(); ++j) {
                vertex_t u = right_region[j];
                left_penalties[u] = 0;
                right_penalties[u] = 0;
                forall_guivw(g, u, i, v, w)
                    {
                        if (region_marker[v] == right_mark || region_marker[v] == left_mark) { continue; } // ignore neighbors that are in the region, they will be handled later
                        partition_t v_id = p_manager[v];

                        left_penalties[u] += w * d_oracle.get(left_id, v_id);
                        right_penalties[u] += w * d_oracle.get(right_id, v_id);
                    }
                endfor
            }
        }

        void build_flow_network(const graph_t &g,
                                d_oracle_t &d_oracle,
                                partition_t left_id,
                                partition_t right_id,
                                std::vector<vertex_t> &left_region,
                                std::vector<vertex_t> &right_region,
                                BKAdapter<int, int, int> &flow_network,
                                TranslationTable<vertex_t> &translation_table,
                                std::vector<u32> &region_marker,
                                u32 &region_mark) {
            ScopedTimer _t("refinement", "FlowBasedRefinement", "build_flow_network");

            u32 left_mark = region_mark - 1;
            u32 right_mark = region_mark;

            weight_t distance = d_oracle.get(left_id, right_id);

            // build flownetwork
            size_t n = left_region.size() + right_region.size();
            flow_network.initialize(n);

            // build left region
            for (size_t j = 0; j < left_region.size(); ++j) {
                vertex_t u = left_region[j];
                ASSERT(region_marker[u] == left_mark);

                forall_guivw(g, u, i, v, w)
                    {
                        if (region_marker[v] == right_mark) {
                            vertex_t new_u = translation_table.get_n(u);
                            vertex_t new_v = translation_table.get_n(v);

                            ASSERT(new_u != new_v);

                            flow_network.add(new_u, new_v, w * distance);
                            continue;
                        }

                        if (region_marker[v] != left_mark) { continue; }
                        if (u < v) { continue; } // no double edges

                        vertex_t new_u = translation_table.get_n(u);
                        vertex_t new_v = translation_table.get_n(v);

                        ASSERT(new_u != new_v);

                        flow_network.add(new_u, new_v, w * distance);
                    }
                endfor
            }

            // build right region
            for (size_t j = 0; j < right_region.size(); ++j) {
                vertex_t u = right_region[j];
                ASSERT(region_marker[u] == right_mark);

                forall_guivw(g, u, i, v, w)
                    {
                        if (region_marker[v] != right_mark) { continue; } // vertex gets handled by penalties, or if v belongs to the left region, no edge is made
                        if (u < v) { continue; } // no double edges

                        vertex_t new_u = translation_table.get_n(u);
                        vertex_t new_v = translation_table.get_n(v);

                        ASSERT(new_u != new_v);

                        flow_network.add(new_u, new_v, w * distance);
                    }
                endfor
            }

            // add the penalties
            for (size_t i = 0; i < left_region.size(); ++i) {
                vertex_t u = left_region[i];
                vertex_t new_u = translation_table.get_n(u);
                weight_t left_penalty = left_penalties[u];
                weight_t right_penalty = right_penalties[u];
                if (left_penalty > 0) { flow_network.add_t_edge(new_u, left_penalty); }
                if (right_penalty > 0) { flow_network.add_s_edge(new_u, right_penalty); }
            }
            for (size_t i = 0; i < right_region.size(); ++i) {
                vertex_t u = right_region[i];
                vertex_t new_u = translation_table.get_n(u);
                weight_t left_penalty = left_penalties[u];
                weight_t right_penalty = right_penalties[u];
                if (left_penalty > 0) { flow_network.add_t_edge(new_u, left_penalty); }
                if (right_penalty > 0) { flow_network.add_s_edge(new_u, right_penalty); }
            }
        }

        bool cut_is_valid(const graph_t &g,
                          const p_manager_t &p_manager,
                          partition_t left_id,
                          partition_t right_id,
                          std::vector<u8> &is_left,
                          weight_t lmax,
                          std::vector<vertex_t> &left_region,
                          std::vector<vertex_t> &right_region,
                          TranslationTable<vertex_t> &translation_table) {
            ScopedTimer _t_sort_pairs("refinement", "FlowBasedRefinement", "cut_is_valid");

            weight_t left_weight = p_manager.get_bweight(left_id);
            weight_t right_weight = p_manager.get_bweight(right_id);
            for (size_t j = 0; j < left_region.size(); ++j) {
                vertex_t u = left_region[j];
                weight_t u_weight = g.v_weights[u];
                vertex_t new_u = translation_table.get_n(u);
                if (is_left[new_u] == 0) {
                    left_weight -= u_weight;
                    right_weight += u_weight;
                }
            }

            for (size_t j = 0; j < right_region.size(); ++j) {
                vertex_t u = right_region[j];
                weight_t u_weight = g.v_weights[u];
                vertex_t new_u = translation_table.get_n(u);
                if (is_left[new_u] == 1) {
                    right_weight -= u_weight;
                    left_weight += u_weight;
                }
            }

            return left_weight <= lmax && right_weight <= lmax;
        }

        bool cut_changes_partition(std::vector<u8> &is_left,
                                   std::vector<vertex_t> &left_region,
                                   std::vector<vertex_t> &right_region,
                                   TranslationTable<vertex_t> &translation_table) {
            ScopedTimer _t_sort_pairs("refinement", "FlowBasedRefinement", "cut_changes_partition");

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
                                        std::vector<u32> &region_marker,
                                        u32 &region_mark) {
            ScopedTimer _t("refinement", "FlowBasedRefinement", "change_boundary");

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
                    p_manager.move(u, g.v_weights[u], left_id, right_id);
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
                    p_manager.move(u, g.v_weights[u], right_id, left_id);
                }
            }
            return changed;
        }

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
            ScopedTimer _t_sort_pairs("refinement", "FlowBasedRefinement", "revert_boundary");

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
                    p_manager.move(old_u, g.v_weights[old_u], left_id, right_id);
                } else {
                    if (bv_manager.is_boundary(old_u)) {
                        bv_manager.move(g, p_manager, old_u, right_id, left_id);
                    } else {
                        bv_manager.add_new(g, p_manager, old_u, left_id);
                    }

                    q_graph.move(g, p_manager, old_u, right_id, left_id);
                    block_conn.move(g, old_u, right_id, left_id);
                    p_manager.move(old_u, g.v_weights[old_u], right_id, left_id);
                }
            }
        }
    };
}

#endif //HEIPROMAP_FLOW_BASED_REFINEMENT_H
