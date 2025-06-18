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
#include "../../commons/flow.h"
#include "../../commons/random_engine.h"
#include "../../commons/utils.h"
#include "../utility/assert_state.h"
#include "../utility/qap.h"

namespace HeiProMap {
    class FlowBasedRefinementConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit FlowBasedRefinementConfiguration(const std::string &t_name) : ISerialRefinerConfiguration(t_name) {}

        u64  max_global_iteration       = 1;
        u64  max_local_iteration        = 3;
        f64  alpha                      = 2.0;
        f64  alpha_upper_bound          = 8.0;
        f64  alpha_modifier             = 2.0;
        bool use_closed_vertex_set      = true;
        u64  closed_vertex_sets_repeats = 10;
        u64  max_level                  = 100;
        u64  min_level                  = 0;
    };

    class FlowBasedRefinement final : public ISerialRefiner {
        vertex_t                 m_n         = 0;
        vertex_t                 m_m         = 0;
        partition_t              m_k         = 0;
        f64                      m_imbalance = 0.0;
        weight_t                 m_lmax      = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t>    m_distance;

        // active block scheduling
        AlignedArray<u8>         active_this_round;
        AlignedArray<u8>         active_next_round;
        AlignedArray<PairWeight> pairs;
        size_t                   pairs_size = 0;

        // array for boundary
        AlignedArray<vertex_t> left_boundary;
        size_t                 left_boundary_size = 0;

        AlignedArray<vertex_t> right_boundary;
        size_t                 right_boundary_size = 0;

        // array for regions
        AlignedArray<vertex_t> left_region;
        size_t                 left_region_size = 0;

        AlignedArray<vertex_t> right_region;
        size_t                 right_region_size = 0;

        AlignedArray<u32> is_left_region;
        AlignedArray<u32> is_right_region;
        u32               is_region_mark = 0;

        AlignedArray<vertex_t> queue;
        size_t                 queue_size = 0;

        AlignedArray<u32> seen;
        u32               seen_mark = 0;

        // array for penalties
        AlignedArray<weight_t> left_penalties;
        AlignedArray<weight_t> right_penalties;

        //Translation Table for mapping
        TranslationTable<vertex_t> translation_table;

        FlowNetwork         flow_network;
        ResidualFlowNetwork residual_flow_network;
        SCCGraph            scc_graph;

        RandomEngine                           *random_engine    = nullptr;
        const FlowBasedRefinementConfiguration *config           = nullptr;

    public:
        FlowBasedRefinement() = default;

        ~FlowBasedRefinement() override = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const f64 t_imbalance,
                        const weight_t t_lmax,
                        const std::vector<partition_t> &t_hierarchy,
                        const std::vector<weight_t> &t_distance,
                        RandomEngine &t_random_engine,
                        const ISerialRefinerConfiguration &i_config) override {
            m_n         = t_n;
            m_m         = t_m;
            m_k         = t_k;
            m_lmax      = t_lmax;
            m_imbalance = t_imbalance;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;

            random_engine    = &t_random_engine;
            config           = dynamic_cast<const FlowBasedRefinementConfiguration *>(&i_config);

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

            is_left_region.initialize(m_n, 0);
            is_right_region.initialize(m_n, 0);
            is_region_mark = 0;

            queue.initialize(m_n);
            queue_size = 0;

            seen.initialize(m_n, 0);
            seen_mark = 0;

            left_penalties.initialize(m_n);
            right_penalties.initialize(m_n);

            translation_table.reserve(m_n, m_n);
        }

        void refine(const u64 level,
                    const u64 max_level,
                    graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph) override {

            if (!(config->min_level <= level && level < config->max_level)) { return; }

            active_this_round.initialize(m_k, 1);
            active_next_round.initialize(m_k, 0);

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
                std::sort(pairs.get_ptr(), pairs.get_ptr() + pairs_size, std::greater<>());

                if (pairs_size == 0) { return; }

                for (size_t i = 0; i < pairs_size; ++i) {
                    partition_t left_id  = pairs[i].id1;
                    partition_t right_id = pairs[i].id2;
                    refine_blocks(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, left_id, right_id);
                }

                std::swap(active_this_round, active_next_round);
                active_next_round.initialize(m_k, 0);
            }
        }

        void refine_blocks(const u64 level,
                           const u64 max_level,
                           graph_t &g,
                           d_oracle_t &d_oracle,
                           bv_manager_t &bv_manager,
                           p_manager_t &p_manager,
                           q_graph_t &q_graph,
                           partition_t left_id,
                           partition_t right_id) {
            ASSERT(left_id != right_id);

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
                weight_t lmax             = std::ceil((1.0 + (m_imbalance * alpha)) * ((f64) g.weight() / (f64) m_k));
                weight_t left_max_weight  = lmax - p_manager.get_bweight(right_id);
                weight_t right_max_weight = lmax - p_manager.get_bweight(left_id);

                // get both regions
                weight_t left_region_weight;
                weight_t right_region_weight;
                determine_regions(g, p_manager, left_id, left_max_weight, &left_region_weight, right_id, right_max_weight, &right_region_weight);

                if (left_region_size + right_region_size == 0) {
                    // if both regions are empty, increase their sizes
                    if (alpha == alpha_upper_bound) { return; }
                    alpha = std::min(alpha_modifier * alpha, alpha_upper_bound);
                    continue;
                }

                // determine penalties for all vertices
                determine_penalties(g, p_manager, d_oracle, left_id, right_id);

                // build a translation table from graph to flow network
                vertex_t    new_u = 0;
                for (size_t i     = 0; i < left_region_size; ++i) { translation_table.add(left_region[i], new_u++); }
                for (size_t i     = 0; i < right_region_size; ++i) { translation_table.add(right_region[i], new_u++); }

                // build flownetwork
                build_flow_network(g, d_oracle, left_id, right_id);

                // solve the flow network
                flow_network.solve();

                std::vector<u8> is_left;
                if (config->use_closed_vertex_set) {
#if HEAVYASSERT_ENABLED
                    // get the first cut for comparison
                    flow_network.get_cut(is_left);

                    if (cut_is_valid(g, p_manager, left_id, right_id, is_left)) {
                        // make the changes
                        weight_t qap            = get_qap(g, p_manager, d_oracle);
                        std::vector<u8> changed = change_boundary(g, bv_manager, p_manager, q_graph, is_left, left_id, right_id);
                        HEAVYASSERT(assert_correct_boundary(g, p_manager, bv_manager, m_k));


                        revert_boundary(g, bv_manager, p_manager, q_graph, changed, left_id, right_id);
                        HEAVYASSERT(assert_correct_boundary(g, p_manager, bv_manager, m_k));
                    }
#endif

                    // build residual network
                    flow_network.build_residual_network(residual_flow_network);

                    // build scc graph
                    scc_graph.initialize(residual_flow_network, g, translation_table);

                    // reduce the scc graph
                    scc_graph.reduce();

                    // determine best balanced min cut
                    weight_t left_non_region_weight  = p_manager.get_bweight(left_id) - left_region_weight;
                    weight_t right_non_region_weight = p_manager.get_bweight(right_id) - right_region_weight;
                    bool     closure_found           = scc_graph.find_best_closure(left_non_region_weight, right_non_region_weight, m_lmax, m_lmax, config->closed_vertex_sets_repeats, *random_engine, is_left);
/*
#if HEAVYASSERT_ENABLED
                    std::vector<std::vector<u8>> all_is_left = scc_graph.get_all_closures(left_non_region_weight, right_non_region_weight, m_lmax, 10, *random_engine);
                    std::vector<weight_t> qap_deltas;
                    // make the changes
                    weight_t qap = get_qap(g, p_manager, d_oracle);
                    for (auto& is_left : all_is_left) {
                        std::vector<u8> changed = change_boundary(g, bv_manager, p_manager, q_graph, is_left, left_id, right_id);
                        HEAVYASSERT(assert_correct_boundary(g, p_manager, bv_manager, m_k));

                        qap_deltas.push_back(qap - get_qap(g, p_manager, d_oracle));

                        revert_boundary(g, bv_manager, p_manager, q_graph, changed, left_id, right_id);
                        HEAVYASSERT(assert_correct_boundary(g, p_manager, bv_manager, m_k));
                    }
                    for (size_t i = 0; i < qap_deltas.size(); ++i) {
                        if (qap_deltas[0] != qap_deltas[i]) {
                            print(qap_deltas);
                            std::cout << "original flow cut: " << qap_normal_change << std::endl;
                            print(is_left);
                            for (auto vec : all_is_left) {
                                print(vec);
                            }

                            flow_network.print();
                            residual_flow_network.print();
                            scc_graph.print();

                            exit(1);
                        }
                    }
#endif
 */

                    if (left_region_size == 5 && right_region_size == 5 && scc_graph.get_n_scc() >= 4 && false) {
                        std::cout << "Left-Region" << std::endl;
                        for (size_t i = 0; i < left_region_size; ++i) {
                            vertex_t    u    = left_region[i];
                            weight_t    u_w  = g.weight(u);
                            partition_t u_id = p_manager[u];
                            std::cout << "(" << u << ", " << u_w << ", " << u_id << ") : ";
                            forall_guivw(g, u, j, v, w)
                                {
                                    partition_t v_id = p_manager[v];
                                    std::cout << "(" << v << ", " << w << ", " << d_oracle.get(u_id, v_id) << ", " << v_id << ") ";
                                }
                            endfor
                            std::cout << std::endl;
                        }

                        std::cout << "Right-Region" << std::endl;
                        for (size_t i = 0; i < right_region_size; ++i) {
                            vertex_t    u    = right_region[i];
                            weight_t    u_w  = g.weight(u);
                            partition_t u_id = p_manager[u];
                            std::cout << "(" << u << ", " << u_w << ", " << u_id << ") : ";
                            forall_guivw(g, u, j, v, w)
                                {
                                    partition_t v_id = p_manager[v];
                                    std::cout << "(" << v << ", " << w << ", " << d_oracle.get(u_id, v_id) << ", " << v_id << ") ";
                                }
                            endfor
                            std::cout << std::endl;
                        }

                        std::cout << "Left-Region Penalties" << std::endl;
                        for (size_t i = 0; i < left_region_size; ++i) {
                            vertex_t u         = left_region[i];
                            weight_t l_penalty = left_penalties[u];
                            weight_t r_penalty = right_penalties[u];
                            std::cout << u << " : " << l_penalty << " -- " << r_penalty << std::endl;
                        }

                        std::cout << "Right-Region Penalties" << std::endl;
                        for (size_t i = 0; i < right_region_size; ++i) {
                            vertex_t u         = right_region[i];
                            weight_t l_penalty = left_penalties[u];
                            weight_t r_penalty = right_penalties[u];
                            std::cout << u << " : " << l_penalty << " -- " << r_penalty << std::endl;
                        }

                        flow_network.print();
                        residual_flow_network.print();
                        scc_graph.print();

                        std::cout << "best found closure: ";
                        print(is_left);
                        exit(1);
                    }

                    if (!closure_found) {
                        if (alpha == 1.0) { return; }
                        alpha = std::max(alpha / alpha_modifier, 1.0);
                        continue;
                    }
                } else {
                    // simply use the first cut found
                    flow_network.get_cut(is_left);

                    if (!cut_is_valid(g, p_manager, left_id, right_id, is_left)) {
                        if (alpha == 1.0) { return; }
                        alpha = std::max(alpha / alpha_modifier, 1.0);
                        continue;
                    }
                }

                // check if the cut actually changes the partition
                if (!cut_changes_partition(is_left)) {
                    // cut is valid, but does not change anything
                    if (alpha == 1.0) { return; }
                    alpha = std::max(alpha / alpha_modifier, 1.0);
                    continue;
                }

                // cut is valid and changes the partition, increase alpha
                alpha = std::min(alpha * alpha_modifier, alpha_upper_bound);

                // make the changes
                change_boundary(g, bv_manager, p_manager, q_graph, is_left, left_id, right_id);

                active_next_round[left_id]  = 1;
                active_next_round[right_id] = 1;
            }
        }

        void determine_boundary_vertices(const graph_t &g,
                                         const bv_manager_t &bv_manager,
                                         const p_manager_t &p_manager,
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

        void determine_regions(const graph_t &g,
                               const p_manager_t &p_manager,
                               partition_t left_id,
                               weight_t left_max_weight,
                               weight_t *left_region_weight,
                               partition_t right_id,
                               weight_t right_max_weight,
                               weight_t *right_region_weight) {
            is_region_mark += 1;
            seen_mark += 2;
            // seen[u] == seen_mark     means u is processed
            // seen[u] == seen_mark - 1 means u is in the queue

            weight_t left_curr_weight = 0;

            queue_size = 0;
            for (size_t i = 0; i < left_boundary_size; ++i) {
                vertex_t u = left_boundary[i];
                ASSERT(p_manager[u] == left_id);
                queue[queue_size++] = u;
                seen[u]             = seen_mark - 1;
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
                            }
                        }
                    endfor
                }
                seen[u]    = seen_mark;
            }
            *left_region_weight = left_curr_weight;

            weight_t right_curr_weight = 0;

            queue_size = 0;
            for (size_t i = 0; i < right_boundary_size; ++i) {
                vertex_t u = right_boundary[i];
                ASSERT(p_manager[u] == right_id);
                queue[queue_size++] = u;
                seen[u]             = seen_mark - 1;
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
                            }
                        }
                    endfor
                }
                seen[u]    = seen_mark;
            }
            *right_region_weight = right_curr_weight;
        }

        void determine_penalties(const graph_t &g,
                                 const p_manager_t &p_manager,
                                 d_oracle_t &d_oracle,
                                 partition_t left_id,
                                 partition_t right_id) {
            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u = left_region[j];
                left_penalties[u]  = 0;
                right_penalties[u] = 0;
                forall_guivw(g, u, i, v, w)
                    {
                        if (is_left_region[v] == is_region_mark || is_right_region[v] == is_region_mark) { continue; } // ignore neighbors that are in the region, they will be handled later
                        partition_t v_id = p_manager[v];

                        left_penalties[u] += w * d_oracle.get(left_id, v_id);
                        right_penalties[u] += w * d_oracle.get(right_id, v_id);

                        /*
                        if (v_id != left_id && v_id != right_id) {
                            // peripheral edge
                            left_penalties[u] += w * d_oracle.get(left_id, v_id);
                            right_penalties[u] += w * d_oracle.get(right_id, v_id);
                        } else if (v_id == left_id) {
                            // edge from left region into left block, only right penalty
                            // left_penalties[u] += w * d_oracle.get(left_id, v_id);
                            right_penalties[u] += w * d_oracle.get(right_id, v_id);
                        } else {
                            // edge from left region into right block, only left penalty
                            left_penalties[u] += w * d_oracle.get(left_id, v_id);
                            // right_penalties[u] += w * d_oracle.get(right_id, v_id);
                        }
                         */
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
                        if (is_right_region[v] == is_region_mark || is_left_region[v] == is_region_mark) { continue; } // ignore neighbors that are in the region, they will be handled later
                        partition_t v_id = p_manager[v];

                        left_penalties[u] += w * d_oracle.get(left_id, v_id);
                        right_penalties[u] += w * d_oracle.get(right_id, v_id);

                        /*
                        if (v_id != left_id && v_id != right_id) {
                            // peripheral edge
                            left_penalties[u] += w * d_oracle.get(left_id, v_id);
                            right_penalties[u] += w * d_oracle.get(right_id, v_id);
                        } else if (v_id == left_id) {
                            // edge from right region into left block, only right penalty
                            // left_penalties[u] += w * d_oracle.get(left_id, v_id);
                            right_penalties[u] += w * d_oracle.get(right_id, v_id);
                        } else {
                            // edge from right region into right block, only left penalty
                            left_penalties[u] += w * d_oracle.get(left_id, v_id);
                            // right_penalties[u] += w * d_oracle.get(right_id, v_id);
                        }
                         */
                    }
                endfor
                left_penalties[u] *= 2;
                right_penalties[u] *= 2;
            }
        }

        void build_flow_network(const graph_t &g,
                                d_oracle_t &d_oracle,
                                partition_t left_id,
                                partition_t right_id) {
            weight_t distance = d_oracle.get(left_id, right_id);

            // build flownetwork
            size_t n = left_region_size + right_region_size;
            flow_network.initialize(n);

            // build left region
            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u = left_region[j];
                ASSERT(is_left_region[u] == is_region_mark);

                forall_guivw(g, u, i, v, w)
                    {
                        if (is_right_region[v] == is_region_mark) {
                            vertex_t new_u = translation_table.get_n(u);
                            vertex_t new_v = translation_table.get_n(v);

                            ASSERT(new_u != new_v);

                            flow_network.add(new_u, new_v, 2 * w * distance);
                            continue;
                        }

                        if (is_left_region[v] != is_region_mark) { continue; }
                        if (u < v) { continue; } // no double edges

                        vertex_t new_u = translation_table.get_n(u);
                        vertex_t new_v = translation_table.get_n(v);

                        ASSERT(new_u != new_v);

                        flow_network.add(new_u, new_v, 2 * w * distance);
                    }
                endfor
            }

            // build right region
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

                        flow_network.add(new_u, new_v, 2 * w * distance);
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

        bool cut_is_valid(const graph_t &g,
                          const p_manager_t &p_manager,
                          partition_t left_id,
                          partition_t right_id,
                          std::vector<u8> &is_left) {
            weight_t    left_weight  = p_manager.get_bweight(left_id);
            weight_t    right_weight = p_manager.get_bweight(right_id);
            for (size_t j            = 0; j < left_region_size; ++j) {
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

        bool cut_changes_partition(std::vector<u8> &is_left) {
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

        std::vector<u8> change_boundary(const graph_t &g,
                                        bv_manager_t &bv_manager,
                                        p_manager_t &p_manager,
                                        q_graph_t &q_graph,
                                        std::vector<u8> &is_left,
                                        partition_t left_id,
                                        partition_t right_id) {
            std::vector<u8> changed(is_left.size(), 0);

            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u     = left_region[j];
                vertex_t new_u = translation_table.get_n(u);

                ASSERT(is_left_region[u] == is_region_mark);
                ASSERT(new_u < left_region_size + right_region_size);

                if (is_left[new_u] == 0) {
                    changed[new_u] = 1;
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
                    changed[new_u] = 1;
                    if (bv_manager.is_boundary(u)) {
                        bv_manager.move(g, p_manager, u, right_id, left_id);
                    } else {
                        bv_manager.add_new(g, p_manager, u, left_id);
                    }

                    q_graph.move(g, p_manager, u, right_id, left_id);
                    p_manager.move(u, g.weight(u), right_id, left_id);
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

        void refine_layer(const u64 level,
                                  const u64 max_level,
                                  graph_t& g,
                                  d_oracle_t& d_oracle,
                                  bv_manager_t& bv_manager,
                                  p_manager_t& p_manager,
                                  q_graph_t& q_graph,
                                  size_t layer) override {}
    };
}

#endif //HEIPROMAP_FLOW_BASED_REFINEMENT_H
