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


#include "ISerialRefiner.h"
#include "quotient_graph_refinement.h"
#include "../utility/flow.h"
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

        // active block scheduling
        AlignedArray<u8> active_this_round;
        AlignedArray<u8> active_next_round;
        AlignedArray<PairWeight> pairs;
        size_t pairs_size = 0;

        // array for boundary
        AlignedArray<vertex_t> left_boundary;
        size_t left_boundary_size = 0;
        weight_t left_boundary_weight = 0;

        AlignedArray<vertex_t> right_boundary;
        size_t right_boundary_size = 0;
        weight_t right_boundary_weight = 0;

        // array for regions
        AlignedArray<vertex_t> left_region;
        size_t left_region_size = 0;

        AlignedArray<vertex_t> right_region;
        size_t right_region_size = 0;

        AlignedArray<u32> region_marker;
        u32 region_mark = 0;

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

        RandomEngine random_engine = RandomEngine(0);
        const FlowBasedRefinementConfiguration *config = nullptr;

    public:
        FlowBasedRefinement() = default;

        ~FlowBasedRefinement() override = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const ISerialRefinerConfiguration &i_config) override {
            ScopedTimer _t("io", "FlowBasedRefinement", "initialize");

            m_n = t_n;
            m_m = t_m;
            m_k = t_k;

            config = dynamic_cast<const FlowBasedRefinementConfiguration *>(&i_config);

            // active block scheduling
            active_this_round.initialize(m_k);
            active_next_round.initialize(m_k);
            size_t size = (size_t) m_k * (size_t) m_k;
            pairs.initialize(size);
            pairs_size = 0;

            left_boundary.initialize(m_n);
            left_boundary_size = 0;

            right_boundary.initialize(m_n);
            right_boundary_size = 0;

            left_region.initialize(m_n);
            left_region_size = 0;

            right_region.initialize(m_n);
            right_region_size = 0;

            region_marker.initialize(m_n, 0);
            region_mark = 0;

            queue.initialize(m_n);
            queue_size = 0;

            seen.initialize(m_n, 0);
            seen_mark = 0;

            left_penalties.initialize(m_n);
            right_penalties.initialize(m_n);

            translation_table.reserve(m_n, m_n);
        }

        void refine(graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph,
                    f64 imbalance) override {
            //
            {
                ScopedTimer _t("refinement", "FlowBasedRefinement", "allocate");
                active_this_round.initialize(m_k, 1);
                active_next_round.initialize(m_k, 0);
            }

            for (u64 iteration = 0; iteration < config->max_global_iteration; ++iteration) {
                // determine all pairs in the quotient graph
                {
                    ScopedTimer _t("refinement", "FlowBasedRefinement", "get_pairs");
                    pairs_size = 0;
                    for (partition_t u_id = 0; u_id < m_k; ++u_id) {
                        for (partition_t v_id = u_id + 1; v_id < m_k; ++v_id) {
                            if (q_graph.has_edge(u_id, v_id) && (active_this_round[u_id] || active_this_round[v_id])) {
                                pairs[pairs_size++] = {u_id, v_id, d_oracle.get(u_id, v_id)};
                            }
                        }
                    }
                    if (pairs_size == 0) { return; }
                }

                // shuffle the pairs
                {
                    ScopedTimer _t_sort_pairs("refinement", "FlowBasedRefinement", "shuffle_pairs");
                    std::shuffle(pairs.get_ptr(), pairs.get_ptr() + pairs_size, random_engine.generator);
                }

                for (size_t i = 0; i < pairs_size; ++i) {
                    partition_t left_id = pairs[i].id1;
                    partition_t right_id = pairs[i].id2;
                    refine_blocks(g, d_oracle, bv_manager, p_manager, q_graph, left_id, right_id, imbalance);
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
                           partition_t left_id,
                           partition_t right_id,
                           f64 imbalance) {
            ASSERT(left_id != right_id);

            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) m_k));

            f64 alpha = config->alpha;
            f64 alpha_upper_bound = config->alpha_upper_bound;
            f64 alpha_modifier = config->alpha_modifier;

            u64 max_local_iteration = config->max_local_iteration;
            u64 iteration = 0;

            while (iteration < max_local_iteration) {
                iteration += 1;

                // get boundary vertices
                determine_boundary_vertices(g, bv_manager, p_manager, left_id, right_id);

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
                forall_guiv(g, left_seed_vertex, i, v) {
                        if (p_manager[v] == right_id) {
                            right_seed_vertex = v;
                            break;
                        }
                    }
                endfor
                weight_t left_region_weight = determine_region(g, p_manager, left_id, left_mark, left_max_weight, left_boundary, left_boundary_size, left_region, left_region_size, left_seed_vertex, left_boundary_weight);
                weight_t right_region_weight = determine_region(g, p_manager, right_id, right_mark, right_max_weight, right_boundary, right_boundary_size, right_region, right_region_size, right_seed_vertex, right_boundary_weight);

                if (left_region_size + right_region_size == 0) {
                    // if both regions are empty, increase their sizes
                    if (alpha == alpha_upper_bound) { return; }
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
                    is_valid = cut_is_valid(g, p_manager, left_id, right_id, is_left, lmax);
                }

                if (!is_valid && config->use_closed_vertex_set) {
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
                if (!cut_changes_partition(is_left)) {
                    // cut is valid, but does not change anything
                    if (alpha == 1.0) { return; }
                    if (alpha == alpha_upper_bound) { return; }
                    alpha = std::max(alpha / alpha_modifier, 1.0);
                    continue;
                }

                // cut is valid and changes the partition, increase alpha
                alpha = std::min(alpha * alpha_modifier, alpha_upper_bound);

                // make the changes
                change_boundary(g, bv_manager, p_manager, q_graph, is_left, left_id, right_id);

                active_next_round[left_id] = 1;
                active_next_round[right_id] = 1;
            }
        }

        void determine_boundary_vertices(const graph_t &g,
                                         const bv_manager_t &bv_manager,
                                         const p_manager_t &p_manager,
                                         partition_t left_id,
                                         partition_t right_id) {
            ScopedTimer _t("refinement", "FlowBasedRefinement", "determine_boundary_vertices");

            left_boundary_size = 0;
            left_boundary_weight = 0;
            forall_bv_id_iu(bv_manager, left_id, i, u) {
                    forall_guiv(g, u, j, v) {
                            partition_t v_id = p_manager[v];
                            if (v_id == right_id) {
                                left_boundary[left_boundary_size++] = u;
                                left_boundary_weight += g.v_weights[u];
                                break;
                            }
                        }
                    endfor
                }
            endfor
            std::shuffle(left_boundary.get_ptr(), left_boundary.get_ptr() + left_boundary_size, random_engine.generator);

            right_boundary_size = 0;
            right_boundary_weight = 0;
            forall_bv_id_iu(bv_manager, right_id, i, u) {
                    forall_guiv(g, u, j, v) {
                            partition_t v_id = p_manager[v];
                            if (v_id == left_id) {
                                right_boundary[right_boundary_size++] = u;
                                right_boundary_weight += g.v_weights[u];
                                break;
                            }
                        }
                    endfor
                }
            endfor
            std::shuffle(right_boundary.get_ptr(), right_boundary.get_ptr() + right_boundary_size, random_engine.generator);
        }

        weight_t determine_region(const graph_t &g,
                                  const p_manager_t &p_manager,
                                  partition_t id,
                                  u32 mark,
                                  weight_t max_weight,
                                  AlignedArray<vertex_t> &boundary,
                                  vertex_t boundary_size,
                                  AlignedArray<vertex_t> &region,
                                  vertex_t &region_size,
                                  vertex_t seed_vertex,
                                  weight_t boundary_weight) {
            ScopedTimer _t_sort_pairs("refinement", "FlowBasedRefinement", "determine_regions");

            seen_mark += 2;
            // seen[u] == seen_mark     means u is processed
            // seen[u] == seen_mark - 1 means u is in the queue

            weight_t curr_weight = 0;

            queue_size = 0;
            if (boundary_weight <= max_weight || true) {
                for (size_t i = 0; i < boundary_size; ++i) {
                    vertex_t u = boundary[i];
                    queue[queue_size++] = u;
                    seen[u] = seen_mark - 1;
                }
            } else {
                queue[queue_size++] = seed_vertex;
                seen[seed_vertex] = seen_mark - 1;
            }

            size_t queue_idx = 0;
            region_size = 0;
            while (queue_idx < queue_size) {
                vertex_t u = queue[queue_idx++];
                if (seen[u] == seen_mark) { continue; }

                if (curr_weight + g.v_weights[u] <= max_weight) {
                    region[region_size++] = u;
                    region_marker[u] = mark;
                    curr_weight += g.v_weights[u];
                    forall_guiv(g, u, i, v) {
                            partition_t v_id = p_manager[v];
                            if (v_id != id) { continue; }

                            if (seen[v] != seen_mark && seen[v] != seen_mark - 1) {
                                queue[queue_size++] = v;
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
                                 partition_t right_id) {
            ScopedTimer _t("refinement", "FlowBasedRefinement", "determine_penalties");

            u32 left_mark = region_mark - 1;
            u32 right_mark = region_mark;

            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u = left_region[j];
                left_penalties[u] = 0;
                right_penalties[u] = 0;
                forall_guivw(g, u, i, v, w) {
                        if (region_marker[v] == left_mark || region_marker[v] == right_mark) { continue; } // ignore neighbors that are in the region, they will be handled later
                        partition_t v_id = p_manager[v];

                        left_penalties[u] += w * d_oracle.get(left_id, v_id);
                        right_penalties[u] += w * d_oracle.get(right_id, v_id);
                    }
                endfor
            }
            for (size_t j = 0; j < right_region_size; ++j) {
                vertex_t u = right_region[j];
                left_penalties[u] = 0;
                right_penalties[u] = 0;
                forall_guivw(g, u, i, v, w) {
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
                                partition_t right_id) {
            ScopedTimer _t("refinement", "FlowBasedRefinement", "build_flow_network");

            u32 left_mark = region_mark - 1;
            u32 right_mark = region_mark;

            weight_t distance = d_oracle.get(left_id, right_id);

            // build flownetwork
            size_t n = left_region_size + right_region_size;
            flow_network.initialize(n);

            // build left region
            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u = left_region[j];
                ASSERT(region_marker[u] == left_mark);

                forall_guivw(g, u, i, v, w) {
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
            for (size_t j = 0; j < right_region_size; ++j) {
                vertex_t u = right_region[j];
                ASSERT(region_marker[u] == right_mark);

                forall_guivw(g, u, i, v, w) {
                        if (region_marker[v] != right_mark) { continue; } // vertex gets handled by penalties, or if v belongs to the left region, no edge is made
                        if (u < v) { continue; }                          // no double edges

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

        bool cut_is_valid(const graph_t &g,
                          const p_manager_t &p_manager,
                          partition_t left_id,
                          partition_t right_id,
                          std::vector<u8> &is_left,
                          weight_t lmax) {
            ScopedTimer _t_sort_pairs("refinement", "FlowBasedRefinement", "cut_is_valid");

            weight_t left_weight = p_manager.get_bweight(left_id);
            weight_t right_weight = p_manager.get_bweight(right_id);
            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u = left_region[j];
                weight_t u_weight = g.v_weights[u];
                vertex_t new_u = translation_table.get_n(u);
                if (is_left[new_u] == 0) {
                    left_weight -= u_weight;
                    right_weight += u_weight;
                }
            }

            for (size_t j = 0; j < right_region_size; ++j) {
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

        bool cut_changes_partition(std::vector<u8> &is_left) {
            ScopedTimer _t_sort_pairs("refinement", "FlowBasedRefinement", "cut_changes_partition");

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

        std::vector<u8> change_boundary(const graph_t &g,
                                        bv_manager_t &bv_manager,
                                        p_manager_t &p_manager,
                                        q_graph_t &q_graph,
                                        std::vector<u8> &is_left,
                                        partition_t left_id,
                                        partition_t right_id) {
            ScopedTimer _t("refinement", "FlowBasedRefinement", "change_boundary");

            [[maybe_unused]] u32 left_mark = region_mark - 1;
            [[maybe_unused]] u32 right_mark = region_mark;

            std::vector<u8> changed(is_left.size(), 0);

            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u = left_region[j];
                vertex_t new_u = translation_table.get_n(u);

                ASSERT(region_marker[u] == left_mark);
                ASSERT(new_u < left_region_size + right_region_size);

                if (is_left[new_u] == 0) {
                    changed[new_u] = 1;
                    if (bv_manager.is_boundary(u)) {
                        bv_manager.move(g, p_manager, u, left_id, right_id);
                    } else {
                        bv_manager.add_new(g, p_manager, u, right_id);
                    }

                    q_graph.move(g, p_manager, u, left_id, right_id);
                    p_manager.move(u, g.v_weights[u], left_id, right_id);
                }
            }
            for (size_t j = 0; j < right_region_size; ++j) {
                vertex_t u = right_region[j];
                vertex_t new_u = translation_table.get_n(u);

                ASSERT(region_marker[u] == right_mark);
                ASSERT(new_u < left_region_size + right_region_size);

                if (is_left[new_u] == 1) {
                    changed[new_u] = 1;
                    if (bv_manager.is_boundary(u)) {
                        bv_manager.move(g, p_manager, u, right_id, left_id);
                    } else {
                        bv_manager.add_new(g, p_manager, u, left_id);
                    }

                    q_graph.move(g, p_manager, u, right_id, left_id);
                    p_manager.move(u, g.v_weights[u], right_id, left_id);
                }
            }
            return changed;
        }

        void revert_boundary(const graph_t &g,
                             bv_manager_t &bv_manager,
                             p_manager_t &p_manager,
                             q_graph_t &q_graph,
                             std::vector<u8> &changed,
                             partition_t left_id,
                             partition_t right_id) {
            ScopedTimer _t_sort_pairs("refinement", "FlowBasedRefinement", "revert_boundary");

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
                    p_manager.move(old_u, g.v_weights[old_u], left_id, right_id);
                } else {
                    if (bv_manager.is_boundary(old_u)) {
                        bv_manager.move(g, p_manager, old_u, right_id, left_id);
                    } else {
                        bv_manager.add_new(g, p_manager, old_u, left_id);
                    }

                    q_graph.move(g, p_manager, old_u, right_id, left_id);
                    p_manager.move(old_u, g.v_weights[old_u], right_id, left_id);
                }
            }
        }
    };
}

#endif //HEIPROMAP_FLOW_BASED_REFINEMENT_H
