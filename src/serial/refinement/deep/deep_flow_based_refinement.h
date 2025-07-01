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

#ifndef HEIPROMAP_DEEP_FLOW_BASED_REFINEMENT_H
#define HEIPROMAP_DEEP_FLOW_BASED_REFINEMENT_H

#include <algorithm>

#include "ISerialDeepRefiner.h"
#include "../flow_based_refinement.h"
#include "../../../commons/random_engine.h"
#include "../../../commons/utils.h"

#define PROFILE(var, statement) \
    do { \
        auto __sp = std::chrono::high_resolution_clock::now(); \
        statement; \
        auto __ep = std::chrono::high_resolution_clock::now(); \
        var += get_seconds(__sp, __ep); \
    } while (0)

namespace HeiProMap {
    class DeepFlowBasedRefinementConfiguration final : public ISerialDeepRefinerConfiguration {
    public:
        explicit DeepFlowBasedRefinementConfiguration(const std::string &t_name) : ISerialDeepRefinerConfiguration(t_name) {}

        u64 max_global_iteration = 1;
        u64 max_local_iteration = 3;
        f64 alpha = 2.0;
        f64 alpha_upper_bound = 8.0;
        f64 alpha_modifier = 2.0;
        bool use_closed_vertex_set = true;
        u64 closed_vertex_sets_repeats = 10;
        u64 max_level = 100;
        u64 min_level = 0;
    };

    class DeepFlowBasedRefinement final : public ISerialDeepRefiner {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
        u64 m_threads = 1;
        f64 m_imbalance = 0.0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;
        std::vector<partition_t> k_rem;

        // active block scheduling
        AlignedArray<u8> active_this_round;
        AlignedArray<u8> active_next_round;

        struct thread_info {
            // array for boundary
            AlignedArray<vertex_t> left_boundary;
            size_t left_boundary_size = 0;

            AlignedArray<vertex_t> right_boundary;
            size_t right_boundary_size = 0;

            // array for regions
            AlignedArray<vertex_t> left_region;
            size_t left_region_size = 0;

            AlignedArray<vertex_t> right_region;
            size_t right_region_size = 0;

            AlignedArray<u32> is_left_region;
            AlignedArray<u32> is_right_region;
            u32 is_region_mark = 0;

            AlignedArray<vertex_t> queue;
            size_t queue_size = 0;

            AlignedArray<u32> seen;
            u32 seen_mark = 0;

            // array for penalties
            AlignedArray<weight_t> left_penalties;
            AlignedArray<weight_t> right_penalties;

            //Translation Table for mapping
            TranslationTable<vertex_t> translation_table;

            FlowNetwork flow_network;
            ResidualFlowNetwork residual_flow_network;
            SCCGraph scc_graph;

            RandomEngine random_engine;

            f64 time_get_boundary = 0;
            f64 time_get_region = 0;
            f64 time_penalties = 0;
            f64 time_build_network = 0;
            f64 time_solve_network = 0;
            f64 time_build_residual_network = 0;
            f64 time_init_scc = 0;
            f64 time_reduce_scc = 0;
            f64 time_search_scc = 0;
            f64 time_change_boundary = 0;
        };
        std::vector<thread_info> thread_infos;

        const DeepFlowBasedRefinementConfiguration *config = nullptr;

    public:
        DeepFlowBasedRefinement() = default;

        ~DeepFlowBasedRefinement() override = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const f64 t_imbalance,
                        const u64 t_threads,
                        const std::vector<partition_t> &t_hierarchy,
                        const std::vector<weight_t> &t_distance,
                        RandomEngine &t_random_engine,
                        const ISerialDeepRefinerConfiguration &i_config) override {
            m_n = t_n;
            m_m = t_m;
            m_k = t_k;
            m_threads = t_threads;

            m_imbalance = t_imbalance;
            m_hierarchy = t_hierarchy;
            m_distance = t_distance;

            partition_t temp_k = 1;
            k_rem.push_back(temp_k);
            for (u64 i = 0; i < m_hierarchy.size(); ++i) {
                temp_k *= m_hierarchy[m_hierarchy.size() - 1 - i];
                k_rem.push_back(k_rem.back() * m_hierarchy[i]);
            }

            config = dynamic_cast<const DeepFlowBasedRefinementConfiguration *>(&i_config);

            // active block scheduling
            active_this_round.initialize(m_k);
            active_next_round.initialize(m_k);

            thread_infos.resize(m_threads);
            for (size_t i = 0; i < m_threads; ++i) {
                thread_infos[i].left_boundary.initialize(m_n);
                thread_infos[i].left_boundary_size = 0;

                thread_infos[i].right_boundary.initialize(m_n);
                thread_infos[i].right_boundary_size = 0;

                thread_infos[i].left_region.initialize(m_n);
                thread_infos[i].left_region_size = 0;

                thread_infos[i].right_region.initialize(m_n);
                thread_infos[i].right_region_size = 0;

                thread_infos[i].is_left_region.initialize(m_n, 0);
                thread_infos[i].is_right_region.initialize(m_n, 0);
                thread_infos[i].is_region_mark = 0;

                thread_infos[i].queue.initialize(m_n);
                thread_infos[i].queue_size = 0;

                thread_infos[i].seen.initialize(m_n, 0);
                thread_infos[i].seen_mark = 0;

                thread_infos[i].left_penalties.initialize(m_n);
                thread_infos[i].right_penalties.initialize(m_n);

                thread_infos[i].translation_table.reserve(m_n, m_n);

                thread_infos[i].random_engine = RandomEngine(t_random_engine.get_u32());
            }
        }

        void refine(u64 level,
                    u64 max_level,
                    const graph_t &g,
                    deep_d_oracle_t &d_oracle,
                    deep_bv_manager_t &bv_manager,
                    deep_p_manager_t &p_manager,
                    deep_q_graph_t &q_graph) override {
            if (!(config->min_level <= level && level < config->max_level)) { return; }

            active_this_round.initialize(m_k, 1);
            active_next_round.initialize(m_k, 0);

            for (u64 iteration = 0; iteration < config->max_global_iteration; ++iteration) {
                std::vector<std::vector<std::pair<partition_t, partition_t>>> matchings = q_graph.get_distance_3_matchings(active_this_round);

                for (auto &matching: matchings) {
#pragma omp parallel for num_threads(m_threads) schedule(dynamic)
                    for (auto [u_id, v_id]: matching) {
                        u64 thread_id = omp_get_thread_num();
                        refine_blocks(g, d_oracle, bv_manager, p_manager, q_graph, u_id, v_id, thread_id);
                    }
                }

                std::swap(active_this_round, active_next_round);
                active_next_round.initialize(m_k, 0);
            }

            // Accumulate times
            f64 total_time_get_boundary = 0;
            f64 total_time_get_region = 0;
            f64 total_time_penalties = 0;
            f64 total_time_build_network = 0;
            f64 total_time_solve_network = 0;
            f64 total_time_build_residual_network = 0;
            f64 total_time_init_scc = 0;
            f64 total_time_reduce_scc = 0;
            f64 total_time_search_scc = 0;
            f64 total_time_change_boundary = 0;

            for (const auto& t : thread_infos) {
                total_time_get_boundary += t.time_get_boundary;
                total_time_get_region += t.time_get_region;
                total_time_penalties += t.time_penalties;
                total_time_build_network += t.time_build_network;
                total_time_solve_network += t.time_solve_network;
                total_time_build_residual_network += t.time_build_residual_network;
                total_time_init_scc += t.time_init_scc;
                total_time_reduce_scc += t.time_reduce_scc;
                total_time_search_scc += t.time_search_scc;
                total_time_change_boundary += t.time_change_boundary;
            }

            // Compute grand total
            f64 total_time =
                    total_time_get_boundary +
                    total_time_get_region +
                    total_time_penalties +
                    total_time_build_network +
                    total_time_solve_network +
                    total_time_build_residual_network +
                    total_time_init_scc +
                    total_time_reduce_scc +
                    total_time_search_scc +
                    total_time_change_boundary;

            // Output
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "\n--- Profiling Summary ---\n";
            std::cout << std::left << std::setw(30) << "Category"
                      << std::right << std::setw(12) << "Time (s)"
                      << std::setw(12) << "Percent\n";
            std::cout << std::string(54, '-') << "\n";

            auto print_line = [&](const std::string& label, f64 time) {
                f64 percent = (total_time > 0.0) ? (100.0 * time / total_time) : 0.0;
                std::cout << std::left << std::setw(30) << label
                          << std::right << std::setw(12) << time
                          << std::setw(11) << std::setprecision(2) << percent << "%\n";
            };

            print_line("get_boundary", total_time_get_boundary);
            print_line("get_region", total_time_get_region);
            print_line("penalties", total_time_penalties);
            print_line("build_network", total_time_build_network);
            print_line("solve_network", total_time_solve_network);
            print_line("build_residual_network", total_time_build_residual_network);
            print_line("init_scc", total_time_init_scc);
            print_line("reduce_scc", total_time_reduce_scc);
            print_line("search_scc", total_time_search_scc);
            print_line("change_boundary", total_time_change_boundary);

            std::cout << std::string(54, '-') << "\n";
            std::cout << std::left << std::setw(30) << "TOTAL"
                      << std::right << std::setw(12) << total_time
                      << std::setw(11) << "100.00%\n";
        }

        void refine_blocks(const graph_t &g,
                           deep_d_oracle_t &d_oracle,
                           deep_bv_manager_t &bv_manager,
                           deep_p_manager_t &p_manager,
                           deep_q_graph_t &q_graph,
                           partition_t left_id,
                           partition_t right_id,
                           u64 thread_id) {
            ASSERT(left_id != right_id);

            size_t &left_boundary_size = thread_infos[thread_id].left_boundary_size;
            size_t &left_region_size = thread_infos[thread_id].left_region_size;
            AlignedArray<vertex_t> &left_boundary = thread_infos[thread_id].left_boundary;
            AlignedArray<vertex_t> &left_region = thread_infos[thread_id].left_region;

            size_t &right_boundary_size = thread_infos[thread_id].right_boundary_size;
            size_t &right_region_size = thread_infos[thread_id].right_region_size;
            AlignedArray<vertex_t> &right_boundary = thread_infos[thread_id].right_boundary;
            AlignedArray<vertex_t> &right_region = thread_infos[thread_id].right_region;

            FlowNetwork &flow_network = thread_infos[thread_id].flow_network;
            ResidualFlowNetwork &residual_flow_network = thread_infos[thread_id].residual_flow_network;
            SCCGraph &scc_graph = thread_infos[thread_id].scc_graph;
            TranslationTable<vertex_t> &translation_table = thread_infos[thread_id].translation_table;

            RandomEngine &random_engine = thread_infos[thread_id].random_engine;

            f64 alpha = config->alpha;
            f64 alpha_upper_bound = config->alpha_upper_bound;
            f64 alpha_modifier = config->alpha_modifier;

            u64 max_local_iteration = config->max_local_iteration;
            u64 iteration = 0;

            while (iteration < max_local_iteration) {
                iteration += 1;

                // get boundary vertices
                PROFILE(thread_infos[thread_id].time_get_boundary, determine_boundary_vertices(g, bv_manager, p_manager, left_id, right_id, thread_id));
                // determine_boundary_vertices(g, bv_manager, p_manager, left_id, right_id, thread_id);

                // calc max weight for each bfs
                weight_t left_lmax = std::ceil((1.0 + (m_imbalance * alpha)) * ((f64) g.weight() / (f64) m_k)) * k_rem[p_manager.get_hierarchy_level(left_id)];
                weight_t right_lmax = std::ceil((1.0 + (m_imbalance * alpha)) * ((f64) g.weight() / (f64) m_k)) * k_rem[p_manager.get_hierarchy_level(right_id)];
                weight_t left_max_add_weight = left_lmax - p_manager.get_bweight(left_id);
                weight_t right_max_add_weight = right_lmax - p_manager.get_bweight(right_id);

                // get both regions
                weight_t left_region_weight;
                weight_t right_region_weight;

                PROFILE(thread_infos[thread_id].time_get_region, determine_regions(g, p_manager, left_id, left_max_add_weight, &left_region_weight, right_id, right_max_add_weight, &right_region_weight, thread_id));
                // determine_regions(g, p_manager, left_id, left_max_add_weight, &left_region_weight, right_id, right_max_add_weight, &right_region_weight, thread_id);

                weight_t left_boundary_max_weight = 0;
                weight_t right_boundary_max_weight = 0;
                for (size_t i = 0; i < left_boundary_size; ++i) { left_boundary_max_weight = std::max(left_boundary_max_weight, g.weight(left_boundary[i])); }
                for (size_t i = 0; i < right_boundary_size; ++i) { right_boundary_max_weight = std::max(right_boundary_max_weight, g.weight(right_boundary[i])); }

                if (left_region_size + right_region_size <= 10) {
                    // if both regions are too small, increase their sizes
                    if (alpha == alpha_upper_bound) { return; }
                    alpha = std::min(alpha_modifier * alpha, alpha_upper_bound);
                    continue;
                }

                // build a translation table from graph to flow network
                vertex_t new_u = 0;
                for (size_t i = 0; i < left_region_size; ++i) { translation_table.add(left_region[i], new_u++); }
                for (size_t i = 0; i < right_region_size; ++i) { translation_table.add(right_region[i], new_u++); }

                if (d_oracle.last_level_pair(left_id, right_id)) {
                    // build flownetwork
                    PROFILE(thread_infos[thread_id].time_build_network, build_flow_network_lowest_level(g, d_oracle, p_manager, left_id, right_id, thread_id));
                    // build_flow_network_lowest_level(g, d_oracle, p_manager, left_id, right_id, thread_id);
                } else {
                    // determine penalties for all vertices
                    PROFILE(thread_infos[thread_id].time_penalties, determine_penalties(g, p_manager, d_oracle, left_id, right_id, thread_id));
                    // determine_penalties(g, p_manager, d_oracle, left_id, right_id, thread_id);

                    // build flownetwork
                    PROFILE(thread_infos[thread_id].time_build_network, build_flow_network(g, d_oracle, left_id, right_id, thread_id));
                    // build_flow_network(g, d_oracle, left_id, right_id, thread_id);
                }

                // solve the flow network
                PROFILE(thread_infos[thread_id].time_solve_network, flow_network.solve());
                // flow_network.solve();

                std::vector<u8> is_left;
                if (config->use_closed_vertex_set) {
                    // build residual network
                    PROFILE(thread_infos[thread_id].time_build_residual_network, flow_network.build_residual_network(residual_flow_network));
                    // flow_network.build_residual_network(residual_flow_network);

                    // build scc graph
                    PROFILE(thread_infos[thread_id].time_init_scc, scc_graph.initialize(residual_flow_network, g, translation_table));
                    // scc_graph.initialize(residual_flow_network, g, translation_table);

                    // reduce the scc graph
                    PROFILE(thread_infos[thread_id].time_reduce_scc, scc_graph.reduce());
                    // scc_graph.reduce();

                    // determine best balanced min cut
                    weight_t left_non_region_weight = p_manager.get_bweight(left_id) - left_region_weight;
                    weight_t right_non_region_weight = p_manager.get_bweight(right_id) - right_region_weight;
                    bool closure_found;
                    PROFILE(thread_infos[thread_id].time_search_scc, closure_found = scc_graph.find_best_closure(left_non_region_weight, right_non_region_weight, p_manager.get_lmax(left_id), p_manager.get_lmax(right_id), config->closed_vertex_sets_repeats, random_engine, is_left));
                    // closure_found = scc_graph.find_best_closure(left_non_region_weight, right_non_region_weight, p_manager.get_lmax(left_id), p_manager.get_lmax(right_id), config->closed_vertex_sets_repeats, random_engine, is_left);

                    if (!closure_found) {
                        if (alpha == 1.0) { return; }
                        alpha = std::max(alpha / alpha_modifier, 1.0);
                        continue;
                    }
                } else {
                    // use the first cut found
                    flow_network.get_cut(is_left);

                    if (!cut_is_valid(g, p_manager, left_id, right_id, is_left, thread_id)) {
                        if (alpha == 1.0) { return; }
                        alpha = std::max(alpha / alpha_modifier, 1.0);
                        continue;
                    }
                }

                // check if the cut actually changes the partition
                if (!cut_changes_partition(is_left, thread_id)) {
                    // the cut is valid, but does not change anything
                    if (alpha == 1.0) { return; }
                    alpha = std::max(alpha / alpha_modifier, 1.0);
                    continue;
                }

                // cut is valid and changes the partition, increase alpha
                alpha = std::min(alpha * alpha_modifier, alpha_upper_bound);

                // make the changes
                PROFILE(thread_infos[thread_id].time_change_boundary, change_boundary(g, bv_manager, p_manager, q_graph, is_left, left_id, right_id, thread_id));
                // change_boundary(g, bv_manager, p_manager, q_graph, is_left, left_id, right_id, thread_id);

                active_next_round[left_id] = 1;
                active_next_round[right_id] = 1;
            }
        }

        void determine_boundary_vertices(const graph_t &g,
                                         const deep_bv_manager_t &bv_manager,
                                         const deep_p_manager_t &p_manager,
                                         partition_t left_id,
                                         partition_t right_id,
                                         u64 thread_id) {
            u32 &seen_mark = thread_infos[thread_id].seen_mark;
            AlignedArray<u32> &seen = thread_infos[thread_id].seen;

            size_t &left_boundary_size = thread_infos[thread_id].left_boundary_size;
            AlignedArray<vertex_t> &left_boundary = thread_infos[thread_id].left_boundary;

            size_t &right_boundary_size = thread_infos[thread_id].right_boundary_size;
            AlignedArray<vertex_t> &right_boundary = thread_infos[thread_id].right_boundary;

            seen_mark += 1;
            left_boundary_size = 0;
            right_boundary_size = 0;
            forall_bv_id_iu(bv_manager, left_id, i, u)
                {
                    bool is_boundary = false;
                    forall_guiv(g, u, j, v)
                        {
                            partition_t v_id = p_manager[v];
                            if (v_id == right_id) {
                                is_boundary = true;
                                if(seen[v] != seen_mark){
                                    right_boundary[right_boundary_size++] = v;
                                    seen[v] = seen_mark;
                                }
                            }
                        }
                    endfor
                    if(is_boundary){
                        left_boundary[left_boundary_size++] = u;
                    }
                }
            endfor
        }

        void determine_regions(const graph_t &g,
                               const deep_p_manager_t &p_manager,
                               partition_t left_id,
                               weight_t left_max_add_weight,
                               weight_t *left_region_weight,
                               partition_t right_id,
                               weight_t right_max_add_weight,
                               weight_t *right_region_weight,
                               u64 thread_id) {
            u32 &is_region_mark = thread_infos[thread_id].is_region_mark;

            u32 &seen_mark = thread_infos[thread_id].seen_mark;
            AlignedArray<u32> &seen = thread_infos[thread_id].seen;

            size_t &queue_size = thread_infos[thread_id].queue_size;
            AlignedArray<vertex_t> &queue = thread_infos[thread_id].queue;

            size_t &left_boundary_size = thread_infos[thread_id].left_boundary_size;
            size_t &left_region_size = thread_infos[thread_id].left_region_size;
            AlignedArray<vertex_t> &left_boundary = thread_infos[thread_id].left_boundary;
            AlignedArray<vertex_t> &left_region = thread_infos[thread_id].left_region;
            AlignedArray<u32> &is_left_region = thread_infos[thread_id].is_left_region;

            size_t &right_boundary_size = thread_infos[thread_id].right_boundary_size;
            size_t &right_region_size = thread_infos[thread_id].right_region_size;
            AlignedArray<vertex_t> &right_boundary = thread_infos[thread_id].right_boundary;
            AlignedArray<vertex_t> &right_region = thread_infos[thread_id].right_region;
            AlignedArray<u32> &is_right_region = thread_infos[thread_id].is_right_region;

            is_region_mark += 1;
            seen_mark += 2;
            // seen[u] == seen_mark means u is processed
            // seen[u] == seen_mark - 1 means u is in the queue

            weight_t left_curr_weight = 0;

            queue_size = 0;
            for (size_t i = 0; i < left_boundary_size; ++i) {
                vertex_t u = left_boundary[i];
                ASSERT(p_manager[u] == left_id);
                queue[queue_size++] = u;
                seen[u] = seen_mark - 1;
            }

            // left region can only gro to right_max_add_weight
            size_t queue_idx = 0;
            left_region_size = 0;
            while (queue_idx < queue_size) {
                vertex_t u = queue[queue_idx++];
                if (seen[u] == seen_mark) { continue; }

                if (left_curr_weight + g.weight(u) <= right_max_add_weight) {
                    left_region[left_region_size++] = u;
                    is_left_region[u] = is_region_mark;
                    left_curr_weight += g.weight(u);
                    forall_guiv(g, u, i, v)
                        {
                            partition_t v_id = p_manager[v];
                            if (v_id != left_id) { continue; }

                            if (seen[v] != seen_mark && seen[v] != seen_mark - 1) {
                                queue[queue_size++] = v;
                                seen[v] = seen_mark - 1;
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
                vertex_t u = right_boundary[i];
                ASSERT(p_manager[u] == right_id);
                queue[queue_size++] = u;
                seen[u] = seen_mark - 1;
            }

            // right region can only gro to right_max_add_weight
            right_region_size = 0;
            queue_idx = 0;
            while (queue_idx < queue_size) {
                vertex_t u = queue[queue_idx++];
                if (seen[u] == seen_mark) { continue; }

                if (right_curr_weight + g.weight(u) <= left_max_add_weight) {
                    right_region[right_region_size++] = u;
                    is_right_region[u] = is_region_mark;
                    right_curr_weight += g.weight(u);
                    forall_guiv(g, u, i, v)
                        {
                            partition_t v_id = p_manager[v];
                            if (v_id != right_id) { continue; }

                            if (seen[v] != seen_mark && seen[v] != seen_mark - 1) {
                                queue[queue_size++] = v;
                                seen[v] = seen_mark - 1;
                            }
                        }
                    endfor
                }
                seen[u] = seen_mark;
            }
            *right_region_weight = right_curr_weight;
        }

        void determine_penalties(const graph_t &g,
                                 const deep_p_manager_t &p_manager,
                                 deep_d_oracle_t &d_oracle,
                                 partition_t left_id,
                                 partition_t right_id,
                                 u64 thread_id) {
            u32 &is_region_mark = thread_infos[thread_id].is_region_mark;

            size_t &left_region_size = thread_infos[thread_id].left_region_size;
            AlignedArray<vertex_t> &left_region = thread_infos[thread_id].left_region;
            AlignedArray<u32> &is_left_region = thread_infos[thread_id].is_left_region;
            AlignedArray<weight_t> &left_penalties = thread_infos[thread_id].left_penalties;

            size_t &right_region_size = thread_infos[thread_id].right_region_size;
            AlignedArray<vertex_t> &right_region = thread_infos[thread_id].right_region;
            AlignedArray<u32> &is_right_region = thread_infos[thread_id].is_right_region;
            AlignedArray<weight_t> &right_penalties = thread_infos[thread_id].right_penalties;

            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u = left_region[j];
                left_penalties[u] = 0;
                right_penalties[u] = 0;
                forall_guivw(g, u, i, v, w)
                    {
                        if (is_left_region[v] == is_region_mark || is_right_region[v] == is_region_mark) { continue; } // ignore neighbors that are in the region, they will be handled later
                        partition_t v_id = p_manager[v];

                        left_penalties[u] += w * d_oracle.get(left_id, v_id);
                        right_penalties[u] += w * d_oracle.get(right_id, v_id);
                    }
                endfor
                // left_penalties[u] *= 2;
                // right_penalties[u] *= 2;
            }
            for (size_t j = 0; j < right_region_size; ++j) {
                vertex_t u = right_region[j];
                left_penalties[u] = 0;
                right_penalties[u] = 0;
                forall_guivw(g, u, i, v, w)
                    {
                        if (is_right_region[v] == is_region_mark || is_left_region[v] == is_region_mark) { continue; } // ignore neighbors that are in the region, they will be handled later
                        partition_t v_id = p_manager[v];

                        left_penalties[u] += w * d_oracle.get(left_id, v_id);
                        right_penalties[u] += w * d_oracle.get(right_id, v_id);
                    }
                endfor
                // left_penalties[u] *= 2;
                // right_penalties[u] *= 2;
            }
        }

        void build_flow_network(const graph_t &g,
                                deep_d_oracle_t &d_oracle,
                                partition_t left_id,
                                partition_t right_id,
                                u64 thread_id) {
            u32 &is_region_mark = thread_infos[thread_id].is_region_mark;

            size_t &left_region_size = thread_infos[thread_id].left_region_size;
            AlignedArray<vertex_t> &left_region = thread_infos[thread_id].left_region;
            AlignedArray<u32> &is_left_region = thread_infos[thread_id].is_left_region;
            AlignedArray<weight_t> &left_penalties = thread_infos[thread_id].left_penalties;

            size_t &right_region_size = thread_infos[thread_id].right_region_size;
            AlignedArray<vertex_t> &right_region = thread_infos[thread_id].right_region;
            AlignedArray<u32> &is_right_region = thread_infos[thread_id].is_right_region;
            AlignedArray<weight_t> &right_penalties = thread_infos[thread_id].right_penalties;

            FlowNetwork &flow_network = thread_infos[thread_id].flow_network;
            TranslationTable<vertex_t> &translation_table = thread_infos[thread_id].translation_table;

            weight_t distance = d_oracle.get(left_id, right_id);

            // build flownetwork
            size_t n = left_region_size + right_region_size;
            flow_network.initialize(n);

            // build the left region
            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u = left_region[j];
                ASSERT(is_left_region[u] == is_region_mark);

                forall_guivw(g, u, i, v, w)
                    {
                        if (is_right_region[v] == is_region_mark) {
                            vertex_t new_u = translation_table.get_n(u);
                            vertex_t new_v = translation_table.get_n(v);

                            ASSERT(new_u != new_v);

                            flow_network.add(new_u, new_v, w * distance);
                            continue;
                        }

                        if (is_left_region[v] != is_region_mark) { continue; }
                        if (u < v) { continue; } // no double edges

                        vertex_t new_u = translation_table.get_n(u);
                        vertex_t new_v = translation_table.get_n(v);

                        ASSERT(new_u != new_v);

                        flow_network.add(new_u, new_v, w * distance);
                    }
                endfor
            }

            // build the right region
            for (size_t j = 0; j < right_region_size; ++j) {
                vertex_t u = right_region[j];
                ASSERT(is_right_region[u] == is_region_mark);

                forall_guivw(g, u, i, v, w)
                    {
                        if (is_right_region[v] != is_region_mark) { continue; } // vertex gets handled by penalties, or if v belongs to the left region, no edge is made
                        if (u < v) { continue; } // no double edges

                        vertex_t new_u = translation_table.get_n(u);
                        vertex_t new_v = translation_table.get_n(v);

                        ASSERT(new_u != new_v);

                        flow_network.add(new_u, new_v, w * distance);
                    }
                endfor
            }

            // add the penalties
            for (size_t i = 0; i < left_region_size; ++i) {
                vertex_t u = left_region[i];
                vertex_t new_u = translation_table.get_n(u);
                weight_t left_penalty = left_penalties[u];
                weight_t right_penalty = right_penalties[u];
                if (left_penalty > 0) { flow_network.add_t_edge(new_u, left_penalty); }
                if (right_penalty > 0) { flow_network.add_s_edge(new_u, right_penalty); }
            }
            for (size_t i = 0; i < right_region_size; ++i) {
                vertex_t u = right_region[i];
                vertex_t new_u = translation_table.get_n(u);
                weight_t left_penalty = left_penalties[u];
                weight_t right_penalty = right_penalties[u];
                if (left_penalty > 0) { flow_network.add_t_edge(new_u, left_penalty); }
                if (right_penalty > 0) { flow_network.add_s_edge(new_u, right_penalty); }
            }
        }

        void build_flow_network_lowest_level(const graph_t &g,
                                             deep_d_oracle_t &d_oracle,
                                             deep_p_manager_t &p_manager,
                                             partition_t left_id,
                                             partition_t right_id,
                                             u64 thread_id) {
            u32 &is_region_mark = thread_infos[thread_id].is_region_mark;

            size_t &left_region_size = thread_infos[thread_id].left_region_size;
            AlignedArray<vertex_t> &left_region = thread_infos[thread_id].left_region;
            AlignedArray<u32> &is_left_region = thread_infos[thread_id].is_left_region;

            size_t &right_region_size = thread_infos[thread_id].right_region_size;
            AlignedArray<vertex_t> &right_region = thread_infos[thread_id].right_region;
            AlignedArray<u32> &is_right_region = thread_infos[thread_id].is_right_region;

            FlowNetwork &flow_network = thread_infos[thread_id].flow_network;
            TranslationTable<vertex_t> &translation_table = thread_infos[thread_id].translation_table;

            // build flownetwork
            size_t n = left_region_size + right_region_size;
            flow_network.initialize(n);

            // build the left region
            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u = left_region[j];
                ASSERT(is_left_region[u] == is_region_mark);

                weight_t w_left = 0;
                weight_t w_right = 0;

                forall_guivw(g, u, i, v, w)
                    {
                        if (is_left_region[v] == is_region_mark) {
                            // v belongs to the left region
                            if (u < v) { continue; } // no double edges
                            vertex_t new_u = translation_table.get_n(u);
                            vertex_t new_v = translation_table.get_n(v);

                            ASSERT(new_u != new_v);

                            flow_network.add(new_u, new_v, w);
                            continue;
                        }
                        if (is_right_region[v] == is_region_mark) {
                            // v belongs to the right region
                            if (u < v) { continue; } // no double edges
                            vertex_t new_u = translation_table.get_n(u);
                            vertex_t new_v = translation_table.get_n(v);

                            ASSERT(new_u != new_v);

                            flow_network.add(new_u, new_v, w);
                            continue;
                        }
                        if (p_manager[v] == left_id) {
                            // v belongs to the left block, but not the left region
                            w_left += w;
                            continue;
                        }
                        if (p_manager[v] == right_id) {
                            // v belongs to the left block, but not the left region
                            w_right += w;
                            continue;
                        }
                    }
                endfor
                vertex_t new_u = translation_table.get_n(u);
                if (w_left > 0) { flow_network.add_s_edge(new_u, w_left); }
                if (w_right > 0) { flow_network.add_t_edge(new_u, w_right); }
            }

            // build the right region
            for (size_t j = 0; j < right_region_size; ++j) {
                vertex_t u = right_region[j];
                ASSERT(is_right_region[u] == is_region_mark);

                weight_t w_left = 0;
                weight_t w_right = 0;

                forall_guivw(g, u, i, v, w)
                    {
                        if (is_left_region[v] == is_region_mark) {
                            // v belongs to the left region
                            if (u < v) { continue; } // no double edges
                            vertex_t new_u = translation_table.get_n(u);
                            vertex_t new_v = translation_table.get_n(v);

                            ASSERT(new_u != new_v);

                            flow_network.add(new_u, new_v, w);
                            continue;
                        }
                        if (is_right_region[v] == is_region_mark) {
                            // v belongs to the right region
                            if (u < v) { continue; } // no double edges
                            vertex_t new_u = translation_table.get_n(u);
                            vertex_t new_v = translation_table.get_n(v);

                            ASSERT(new_u != new_v);

                            flow_network.add(new_u, new_v, w);
                            continue;
                        }
                        if (p_manager[v] == left_id) {
                            // v belongs to the left block, but not the left region
                            w_left += w;
                            continue;
                        }
                        if (p_manager[v] == right_id) {
                            // v belongs to the left block, but not the left region
                            w_right += w;
                            continue;
                        }
                    }
                endfor
                vertex_t new_u = translation_table.get_n(u);
                if (w_left > 0) { flow_network.add_s_edge(new_u, w_left); }
                if (w_right > 0) { flow_network.add_t_edge(new_u, w_right); }
            }
        }

        bool cut_is_valid(const graph_t &g,
                          deep_p_manager_t &p_manager,
                          partition_t left_id,
                          partition_t right_id,
                          std::vector<u8> &is_left,
                          u64 thread_id) {
            size_t &left_region_size = thread_infos[thread_id].left_region_size;
            AlignedArray<vertex_t> &left_region = thread_infos[thread_id].left_region;

            size_t &right_region_size = thread_infos[thread_id].right_region_size;
            AlignedArray<vertex_t> &right_region = thread_infos[thread_id].right_region;

            TranslationTable<vertex_t> &translation_table = thread_infos[thread_id].translation_table;

            weight_t left_weight = p_manager.get_bweight(left_id);
            weight_t right_weight = p_manager.get_bweight(right_id);
            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u = left_region[j];
                weight_t u_weight = g.weight(u);
                vertex_t new_u = translation_table.get_n(u);
                if (is_left[new_u] == 0) {
                    left_weight -= u_weight;
                    right_weight += u_weight;
                }
            }

            for (size_t j = 0; j < right_region_size; ++j) {
                vertex_t u = right_region[j];
                weight_t u_weight = g.weight(u);
                vertex_t new_u = translation_table.get_n(u);
                if (is_left[new_u] == 1) {
                    right_weight -= u_weight;
                    left_weight += u_weight;
                }
            }

            return left_weight <= p_manager.get_lmax(left_id) && right_weight <= p_manager.get_lmax(right_id);
        }

        bool cut_changes_partition(std::vector<u8> &is_left,
                                   u64 thread_id) {
            size_t &left_region_size = thread_infos[thread_id].left_region_size;
            AlignedArray<vertex_t> &left_region = thread_infos[thread_id].left_region;

            size_t &right_region_size = thread_infos[thread_id].right_region_size;
            AlignedArray<vertex_t> &right_region = thread_infos[thread_id].right_region;

            TranslationTable<vertex_t> &translation_table = thread_infos[thread_id].translation_table;

            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u = left_region[j];
                vertex_t new_u = translation_table.get_n(u);
                if (is_left[new_u] == 0) {
                    return true;
                }
            }
            for (size_t j = 0; j < right_region_size; ++j) {
                vertex_t u = right_region[j];
                vertex_t new_u = translation_table.get_n(u);
                if (is_left[new_u] == 1) {
                    return true;
                }
            }

            return false;
        }

        void change_boundary(const graph_t &g,
                             deep_bv_manager_t &bv_manager,
                             deep_p_manager_t &p_manager,
                             deep_q_graph_t &q_graph,
                             std::vector<u8> &is_left,
                             partition_t left_id,
                             partition_t right_id,
                             u64 thread_id) {
            u32 &is_region_mark = thread_infos[thread_id].is_region_mark;
            size_t &left_region_size = thread_infos[thread_id].left_region_size;
            AlignedArray<vertex_t> &left_region = thread_infos[thread_id].left_region;
            AlignedArray<u32> &is_left_region = thread_infos[thread_id].is_left_region;

            size_t &right_region_size = thread_infos[thread_id].right_region_size;
            AlignedArray<vertex_t> &right_region = thread_infos[thread_id].right_region;
            AlignedArray<u32> &is_right_region = thread_infos[thread_id].is_right_region;

            TranslationTable<vertex_t> &translation_table = thread_infos[thread_id].translation_table;

            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u = left_region[j];
                vertex_t new_u = translation_table.get_n(u);

                ASSERT(is_left_region[u] == is_region_mark);
                ASSERT(new_u < left_region_size + right_region_size);

                if (is_left[new_u] == 0) {
                    bv_manager.move(g, p_manager, u, left_id, right_id);
                    q_graph.move(g, p_manager, u, left_id, right_id);
                    p_manager.move(u, g.weight(u), left_id, right_id);
                }
            }
            for (size_t j = 0; j < right_region_size; ++j) {
                vertex_t u = right_region[j];
                vertex_t new_u = translation_table.get_n(u);

                ASSERT(is_right_region[u] == is_region_mark);
                ASSERT(new_u < left_region_size + right_region_size);

                if (is_left[new_u] == 1) {
                    bv_manager.move(g, p_manager, u, right_id, left_id);
                    q_graph.move(g, p_manager, u, right_id, left_id);
                    p_manager.move(u, g.weight(u), right_id, left_id);
                }
            }
        }

        void revert_boundary(const graph_t &g,
                             deep_bv_manager_t &bv_manager,
                             deep_p_manager_t &p_manager,
                             deep_q_graph_t &q_graph,
                             std::vector<u8> &changed,
                             partition_t left_id,
                             partition_t right_id,
                             u64 thread_id) {
            size_t &left_region_size = thread_infos[thread_id].left_region_size;

            size_t &right_region_size = thread_infos[thread_id].right_region_size;

            TranslationTable<vertex_t> &translation_table = thread_infos[thread_id].translation_table;

            for (vertex_t new_u = 0; new_u < left_region_size + right_region_size; ++new_u) {
                if (changed[new_u] == 0) { continue; }

                vertex_t old_u = translation_table.get_o(new_u);
                if (p_manager[old_u] == left_id) {
                    if (bv_manager.is_boundary(old_u)) {
                        bv_manager.move(g, p_manager, old_u, left_id, right_id);
                    } else {
                        bv_manager.add_new(g, p_manager, old_u, right_id);
                    }

                    q_graph.move(g, p_manager, old_u, left_id, right_id);
                    p_manager.move(old_u, g.weight(old_u), left_id, right_id);
                } else {
                    if (bv_manager.is_boundary(old_u)) {
                        bv_manager.move(g, p_manager, old_u, right_id, left_id);
                    } else {
                        bv_manager.add_new(g, p_manager, old_u, left_id);
                    }

                    q_graph.move(g, p_manager, old_u, right_id, left_id);
                    p_manager.move(old_u, g.weight(old_u), right_id, left_id);
                }
            }
        }
    };
}

#endif //HEIPROMAP_DEEP_FLOW_BASED_REFINEMENT_H
