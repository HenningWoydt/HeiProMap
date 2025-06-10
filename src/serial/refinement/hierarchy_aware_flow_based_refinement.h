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

#ifndef HIERARCHY_AWARE_FLOW_BASED_REFINEMENT_H
#define HIERARCHY_AWARE_FLOW_BASED_REFINEMENT_H

#include <algorithm>
#include <unordered_set>

#include "../../extern/maxflow-v3.04.src/graph.h"

#include "flow_based_refinement.h"
#include "quotient_graph_refinement.h"
#include "../../commons/indexed_update_heap.h"
#include "../../commons/random_engine.h"
#include "../../commons/statistic_collector.h"
#include "../../commons/utils.h"
#include "../utility/functions.h"
#include "ISerialRefiner.h"
#include "../utility/qap.h"

namespace HeiProMap {
    class HierarchyAwareFlowBasedRefinementConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit HierarchyAwareFlowBasedRefinementConfiguration(const std::string &t_name) : ISerialRefinerConfiguration(t_name) {}

        u64  max_global_iteration       = 1;
        u64  max_local_iteration        = 3;
        f64  alpha                      = 2.0;
        f64  alpha_upper_bound          = 8.0;
        f64  alpha_modifier             = 2.0;
        bool use_closed_vertex_set      = true;
        u64  closed_vertex_sets_repeats = 10;
    };

    class HierarchyAwareFlowBasedRefinement final : public ISerialRefiner {
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

        // array for boundary vertices
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

        AlignedArray<u32> vertex_used;
        u32               vertex_marker = 0;

        AlignedArray<u32> block_used;
        u32               block_marker = 0;

        AlignedArray<Move> moves;
        size_t             moves_size = 0;

        FlowNetwork         flow_network;
        ResidualFlowNetwork residual_flow_network;
        SCCGraph            scc_graph;

        RandomEngine                                         *random_engine    = nullptr;
        const HierarchyAwareFlowBasedRefinementConfiguration *config           = nullptr;
        StatisticCollector                                   *m_stat_collector = nullptr;

    public:
        HierarchyAwareFlowBasedRefinement() = default;

        ~HierarchyAwareFlowBasedRefinement() override = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const f64 t_imbalance,
                        const weight_t t_lmax,
                        const std::vector<partition_t> &t_hierarchy,
                        const std::vector<weight_t> &t_distance,
                        RandomEngine &t_random_engine,
                        const ISerialRefinerConfiguration &i_config,
                        StatisticCollector &t_stat_collect) override {
            m_n         = t_n;
            m_m         = t_m;
            m_k         = t_k;
            m_lmax      = t_lmax;
            m_imbalance = t_imbalance;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;

            random_engine    = &t_random_engine;
            config           = dynamic_cast<const HierarchyAwareFlowBasedRefinementConfiguration *>(&i_config);
            m_stat_collector = &t_stat_collect;

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

            vertex_used.initialize(m_n, 0);
            vertex_marker = 0;

            block_used.initialize(m_k, 0);
            block_marker = 0;

            moves.initialize(m_n);
            moves_size = 0;
        }

        void refine(const u64 level,
                    const u64 max_level,
                    graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph) override {
            for (size_t iteration = 0; iteration < config->max_global_iteration; ++iteration) {
                for (size_t i = 0; i < m_hierarchy.size(); ++i) {
                    refine_layer(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, m_hierarchy.size() - 1 - i);
                }
            }
        }

        void refine_layer(const u64 level,
                          const u64 max_level,
                          graph_t &g,
                          d_oracle_t &d_oracle,
                          bv_manager_t &bv_manager,
                          p_manager_t &p_manager,
                          q_graph_t &q_graph,
                          size_t layer) {
            partition_t n_upper_total_super_blocks = 1;
            partition_t n_total_super_blocks       = m_hierarchy[layer];
            for (size_t i                          = layer + 1; i < m_hierarchy.size(); ++i) { n_upper_total_super_blocks *= m_hierarchy[i]; }
            partition_t ids_per_super_block        = m_k / (n_upper_total_super_blocks * m_hierarchy[layer]);

            for (u64 iteration = 0; iteration < config->max_global_iteration; ++iteration) {
                for (size_t i = 0; i < n_upper_total_super_blocks; ++i) {
                    for (size_t j = 0; j < n_total_super_blocks; ++j) {
                        for (size_t k = j + 1; k < n_total_super_blocks; ++k) {
                            partition_t l_start = i * (n_total_super_blocks * ids_per_super_block) + j * ids_per_super_block;
                            partition_t r_start = i * (n_total_super_blocks * ids_per_super_block) + k * ids_per_super_block;

                            if (ids_per_super_block != 1) {
                                refine_blocks(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, l_start, r_start, ids_per_super_block, m_lmax * ids_per_super_block);
                            }
                        }
                    }
                }
            }
        }

        void refine_blocks(const u64 level,
                           const u64 max_level,
                           graph_t &g,
                           d_oracle_t &d_oracle,
                           bv_manager_t &bv_manager,
                           p_manager_t &p_manager,
                           q_graph_t &q_graph,
                           partition_t l_start,
                           partition_t r_start,
                           partition_t ids_per_super_block,
                           weight_t lmax) {

            f64 alpha             = config->alpha;
            f64 alpha_upper_bound = config->alpha_upper_bound;
            f64 alpha_modifier    = config->alpha_modifier;

            u64 max_local_iteration = config->max_local_iteration;
            u64 iteration           = 0;

            while (iteration < max_local_iteration) {
                iteration += 1;

                // get allowed weight on each side
                weight_t l_lmax = std::ceil((1.0 + (m_imbalance * alpha)) * ((f64) g.weight() / (f64) m_k)) * ids_per_super_block;
                weight_t r_lmax = std::ceil((1.0 + (m_imbalance * alpha)) * ((f64) g.weight() / (f64) m_k)) * ids_per_super_block;

                for (partition_t id = l_start; id < l_start + ids_per_super_block; ++id) { l_lmax -= p_manager.get_bweight(id); }
                for (partition_t id = r_start; id < r_start + ids_per_super_block; ++id) { r_lmax -= p_manager.get_bweight(id); }

                // get boundary vertices
                determine_boundary_vertices(g, bv_manager, p_manager, l_start, r_start, ids_per_super_block);

                // get regions
                determine_regions(g, p_manager, l_start, r_start, ids_per_super_block, l_lmax, r_lmax);

                if (left_region_size + right_region_size == 0) {
                    // if both regions are empty, increase their sizes
                    if (alpha == alpha_upper_bound) { return; }
                    alpha = std::min(alpha_modifier * alpha, alpha_upper_bound);
                    continue;
                }

                // determine penalties for all vertices
                determine_penalties(g, p_manager, d_oracle, l_start, r_start, ids_per_super_block);

                weight_t    l_region_weight = 0;
                weight_t    r_region_weight = 0;
                for (size_t i               = 0; i < left_region_size; ++i) { l_region_weight += g.weight(left_region[i]); }
                for (size_t i               = 0; i < right_region_size; ++i) { r_region_weight += g.weight(right_region[i]); }

                // build a translation table from graph to flow network
                vertex_t    new_u = 0;
                for (size_t i     = 0; i < left_region_size; ++i) { translation_table.add(left_region[i], new_u++); }
                for (size_t i     = 0; i < right_region_size; ++i) { translation_table.add(right_region[i], new_u++); }

                // build flownetwork
                build_flow_network(g, d_oracle, l_start, r_start);

                // solve the flow network
                flow_network.solve();

                std::vector<u8> is_left;
                if (config->use_closed_vertex_set) {
                    // build residual network
                    flow_network.build_residual_network(residual_flow_network);

                    // build scc graph
                    scc_graph.initialize(residual_flow_network, g, translation_table);

                    // reduce the scc graph
                    scc_graph.reduce();

                    // determine best balanced min cut
                    weight_t         l_weight = 0;
                    weight_t         r_weight = 0;
                    for (partition_t id       = l_start; id < l_start + ids_per_super_block; ++id) { l_weight += p_manager.get_bweight(id); }
                    for (partition_t id       = r_start; id < r_start + ids_per_super_block; ++id) { r_weight += p_manager.get_bweight(id); }

                    weight_t left_non_region_weight  = l_weight - l_region_weight;
                    weight_t right_non_region_weight = r_weight - r_region_weight;
                    bool     closure_found           = scc_graph.find_best_closure(left_non_region_weight, right_non_region_weight, lmax, config->closed_vertex_sets_repeats, *random_engine, is_left);

                    if (!closure_found) {
                        if (alpha == 1.0) { return; }
                        alpha = std::max(alpha / alpha_modifier, 1.0);
                        continue;
                    }

                } else {
                    flow_network.get_cut(is_left);

                    if (!cut_is_valid(g, p_manager, l_start, r_start, ids_per_super_block, lmax, is_left)) {
                        if (alpha == 1.0) { return; }
                        alpha = std::max(alpha / alpha_modifier, 1.0);
                        continue;
                    }
                }

                // check if the cut actually changes the partition
                if (get_n_changes(is_left) == 0) {
                    // cut is valid, but does not change anything
                    if (alpha == 1.0) { return; }
                    alpha = std::max(alpha / alpha_modifier, 1.0);
                    continue;
                }

                // cut is valid and changes the partition, increase alpha
                alpha = std::min(alpha * alpha_modifier, alpha_upper_bound);

                change_boundary(g, d_oracle, bv_manager, p_manager, q_graph, is_left, l_start, r_start, ids_per_super_block);

                weight_t         l_weight = 0;
                weight_t         r_weight = 0;
                for (partition_t id       = l_start; id < l_start + ids_per_super_block; ++id) { l_weight += p_manager.get_bweight(id); }
                for (partition_t id       = r_start; id < r_start + ids_per_super_block; ++id) { r_weight += p_manager.get_bweight(id); }

                if(l_weight > lmax){
                    // rebalance left superblock
                    std::cout << "left cant be rebalanced" << std::endl;
                }
                if(r_weight > lmax){
                    // rebalance right superblock
                    std::cout << "right cant be rebalanced" << std::endl;
                }
            }
        }

        void determine_boundary_vertices(const graph_t &g,
                                         bv_manager_t &bv_manager,
                                         const p_manager_t &p_manager,
                                         partition_t l_start,
                                         partition_t r_start,
                                         partition_t ids_per_super_block) {
            left_boundary_size  = 0;
            for (partition_t id = l_start; id < l_start + ids_per_super_block; ++id) {
                forall_bv_id_iu(bv_manager, id, i, u)
                    {
                        forall_guiv(g, u, j, v)
                            {
                                partition_t v_id = p_manager[v];
                                if (r_start <= v_id && v_id < r_start + ids_per_super_block) {
                                    left_boundary[left_boundary_size++] = u;
                                    break;
                                }
                            }
                        endfor
                    }
                endfor
            }

            right_boundary_size = 0;
            for (partition_t id = r_start; id < r_start + ids_per_super_block; ++id) {
                forall_bv_id_iu(bv_manager, id, i, u)
                    {
                        forall_guiv(g, u, j, v)
                            {
                                partition_t v_id = p_manager[v];
                                if (l_start <= v_id && v_id < l_start + ids_per_super_block) {
                                    right_boundary[right_boundary_size++] = u;
                                    break;
                                }
                            }
                        endfor
                    }
                endfor
            }
        }

        void determine_regions(const graph_t &g,
                               const p_manager_t &p_manager,
                               partition_t l_start,
                               partition_t r_start,
                               partition_t ids_per_super_block,
                               weight_t l_lmax,
                               weight_t r_lmax) {
            is_region_mark += 1;
            seen_mark += 2;
            // seen[u] == seen_mark     means u is processed
            // seen[u] == seen_mark - 1 means u is in the queue

            // fill the left boundary in the queue
            size_t queue_idx = 0;
            queue_size = 0;
            for (size_t i = 0; i < left_boundary_size; ++i) {
                vertex_t u = left_boundary[i];
                queue[queue_size++] = u;
                seen[u]             = seen_mark - 1;
            }

            // empty queue in bfs fashion and add vertices in the same superblock
            left_region_size = 0;
            while (queue_idx < queue_size) {
                vertex_t u = queue[queue_idx++];
                if (seen[u] == seen_mark) { continue; }

                if (g.weight(u) <= r_lmax) {
                    r_lmax -= g.weight(u);
                    left_region[left_region_size++] = u;
                    is_left_region[u]               = is_region_mark;
                    forall_guiv(g, u, i, v)
                        {
                            partition_t v_id = p_manager[v];
                            if (!(l_start <= v_id && v_id < l_start + ids_per_super_block)) { continue; }

                            if (seen[v] != seen_mark && seen[v] != seen_mark - 1) {
                                queue[queue_size++] = v;
                                seen[v]             = seen_mark - 1;
                            }
                        }
                    endfor
                }
                seen[u]    = seen_mark;
            }

            // fill the right boundary in the queue
            queue_idx  = 0;
            queue_size = 0;
            for (size_t i = 0; i < right_boundary_size; ++i) {
                vertex_t u = right_boundary[i];
                ASSERT(r_start <= p_manager[u] && p_manager[u] < r_start + ids_per_super_block);
                queue[queue_size++] = u;
                seen[u]             = seen_mark - 1;
            }

            // empty queue in bfs fashion and add vertices in the same superblock
            right_region_size = 0;
            while (queue_idx < queue_size) {
                vertex_t u = queue[queue_idx++];
                if (seen[u] == seen_mark) { continue; }

                if (g.weight(u) <= l_lmax) {
                    l_lmax -= g.weight(u);
                    right_region[right_region_size++] = u;
                    is_right_region[u]                = is_region_mark;
                    forall_guiv(g, u, i, v)
                        {
                            partition_t v_id = p_manager[v];
                            if (!(r_start <= v_id && v_id < r_start + ids_per_super_block)) { continue; }

                            if (seen[v] != seen_mark && seen[v] != seen_mark - 1) {
                                queue[queue_size++] = v;
                                seen[v]             = seen_mark - 1;
                            }
                        }
                    endfor
                }
                seen[u]    = seen_mark;
            }
        }

        void determine_penalties(const graph_t &g,
                                 const p_manager_t &p_manager,
                                 d_oracle_t &d_oracle,
                                 partition_t l_start,
                                 partition_t r_start,
                                 partition_t ids_per_super_block) {
            weight_t         max_l_distance = 0;
            weight_t         max_r_distance = 0;
            for (partition_t id             = 0; id < ids_per_super_block; ++id) {
                max_l_distance = std::max(max_l_distance, d_oracle.get(l_start, l_start + id));
                max_r_distance = std::max(max_r_distance, d_oracle.get(r_start, r_start + id));
            }
            max_l_distance                  = 0;
            max_r_distance                  = 0;

            for (size_t i = 0; i < left_region_size; ++i) {
                vertex_t    u    = left_region[i];
                partition_t u_id = p_manager[u];
                left_penalties[u]  = 0; // penalty for being on the left side
                right_penalties[u] = 0; // penalty for being on the right side

                forall_guivw(g, u, j, v, w)
                    {
                        if (is_left_region[v] == is_region_mark || is_right_region[v] == is_region_mark) { continue; }

                        partition_t v_id             = p_manager[v];
                        bool        in_l_super_block = l_start <= v_id && v_id < l_start + ids_per_super_block;
                        bool        in_r_super_block = r_start <= v_id && v_id < r_start + ids_per_super_block;

                        if (!in_l_super_block && !in_r_super_block) {
                            left_penalties[u] += w * d_oracle.get(l_start, v_id);
                            right_penalties[u] += w * d_oracle.get(r_start, v_id);
                            continue;
                        }

                        left_penalties[u] += w * d_oracle.get(u_id, v_id); // u stays left, add current penalty

                        if (l_start <= v_id && v_id < l_start + ids_per_super_block) {
                            // u moves to the right, but v is in the left
                            right_penalties[u] += w * d_oracle.get(l_start, r_start);
                        } else {
                            // u moves to the right and v is also on the right,
                            right_penalties[u] += w * max_r_distance;
                        }
                    }
                endfor
                left_penalties[u] *= 2;
                right_penalties[u] *= 2;
            }

            for (size_t i = 0; i < right_region_size; ++i) {
                vertex_t    u    = right_region[i];
                partition_t u_id = p_manager[u];
                left_penalties[u]  = 0;
                right_penalties[u] = 0;

                forall_guivw(g, u, j, v, w)
                    {
                        if (is_left_region[v] == is_region_mark || is_right_region[v] == is_region_mark) { continue; }

                        partition_t v_id             = p_manager[v];
                        bool        in_l_super_block = l_start <= v_id && v_id < l_start + ids_per_super_block;
                        bool        in_r_super_block = r_start <= v_id && v_id < r_start + ids_per_super_block;

                        if (!in_l_super_block && !in_r_super_block) {
                            left_penalties[u] += w * d_oracle.get(l_start, v_id);
                            right_penalties[u] += w * d_oracle.get(r_start, v_id);
                            continue;
                        }

                        right_penalties[u] += w * d_oracle.get(u_id, v_id); // u stays right, add current penalty

                        if (r_start <= v_id && v_id < r_start + ids_per_super_block) {
                            // u moves to the left, but v is on the right
                            left_penalties[u] += w * d_oracle.get(r_start, l_start);
                        } else {
                            // u moves to the left and v is on the left
                            left_penalties[u] += w * max_l_distance;
                        }
                    }
                endfor
                left_penalties[u] *= 2;
                right_penalties[u] *= 2;
            }
        }


        void build_flow_network(const graph_t &g,
                                d_oracle_t &d_oracle,
                                partition_t l_start,
                                partition_t r_start) {
            weight_t distance = d_oracle.get(l_start, r_start);

            // build flownetwork
            size_t n = left_region_size + right_region_size;
            flow_network.initialize(n);

            // build the left region
            for (size_t i = 0; i < left_region_size; ++i) {
                vertex_t u = left_region[i];

                forall_guivw(g, u, j, v, w)
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

            // build the right region
            for (size_t i = 0; i < right_region_size; ++i) {
                vertex_t u = right_region[i];

                forall_guivw(g, u, j, v, w)
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


        size_t get_n_changes(std::vector<u8> &is_left) {
            size_t      count = 0;
            for (size_t j     = 0; j < left_region_size; ++j) {
                vertex_t u     = left_region[j];
                vertex_t new_u = translation_table.get_n(u);
                if (is_left[new_u] == 0) {
                    count += 1;
                }
            }
            for (size_t j     = 0; j < right_region_size; ++j) {
                vertex_t u     = right_region[j];
                vertex_t new_u = translation_table.get_n(u);
                if (is_left[new_u] == 1) {
                    count += 1;
                }
            }

            return count;
        }

        bool cut_is_valid(const graph_t &g,
                          const p_manager_t &p_manager,
                          partition_t l_start,
                          partition_t r_start,
                          partition_t ids_per_superblock,
                          weight_t lmax,
                          std::vector<u8> &is_left) {
            weight_t         l_weight = 0;
            weight_t         r_weight = 0;
            for (partition_t id       = 0; id < ids_per_superblock; ++id) {
                l_weight += p_manager.get_bweight(l_start + id);
                r_weight += p_manager.get_bweight(r_start + id);
            }

            for (size_t i = 0; i < left_region_size; ++i) {
                vertex_t u        = left_region[i];
                weight_t u_weight = g.weight(u);
                vertex_t new_u    = translation_table.get_n(u);
                if (is_left[new_u] == 0) {
                    l_weight -= u_weight;
                    r_weight += u_weight;
                }
            }

            for (size_t i = 0; i < right_region_size; ++i) {
                vertex_t u        = right_region[i];
                weight_t u_weight = g.weight(u);
                vertex_t new_u    = translation_table.get_n(u);
                if (is_left[new_u] == 1) {
                    r_weight -= u_weight;
                    l_weight += u_weight;
                }
            }

            return l_weight <= lmax && r_weight <= lmax;
        }

        void change_boundary(graph_t &g,
                             d_oracle_t &d_oracle,
                             bv_manager_t &bv_manager,
                             p_manager_t &p_manager,
                             q_graph_t &q_graph,
                             std::vector<u8> &is_left,
                             partition_t l_start,
                             partition_t r_start,
                             partition_t ids_per_superblock) {
            // move all vertices of the left region
            for (size_t i = 0; i < left_region_size; ++i) {
                vertex_t    u        = left_region[i];
                partition_t u_id     = p_manager[u];
                weight_t    u_weight = g.weight(u);
                vertex_t    new_u    = translation_table.get_n(u);

                ASSERT(l_start <= u_id && u_id < l_start + ids_per_superblock);
                ASSERT(is_left_region[u] == is_region_mark);
                ASSERT(new_u < left_region_size + right_region_size);

                if (is_left[new_u] == 0) {
                    // determine the best right block to move u to, that will not be overloaded
                    partition_t      best_id        = r_start;
                    s64              best_qap_delta = get_u_qap_delta(g, u, u_id, r_start, p_manager, d_oracle);
                    weight_t         best_weight    = p_manager.get_bweight(r_start) + u_weight;
                    for (partition_t id             = r_start; id < r_start + ids_per_superblock; ++id) {
                        // if (p_manager.get_bweight(id) + u_weight > m_lmax) { continue; }
                        s64 qap_delta = get_u_qap_delta(g, u, u_id, id, p_manager, d_oracle);
                        if (qap_delta > best_qap_delta || (qap_delta == best_qap_delta && p_manager.get_bweight(id) + u_weight < best_weight)) {
                            best_id        = id;
                            best_qap_delta = qap_delta;
                            best_weight    = p_manager.get_bweight(id) + u_weight;
                        }
                    }

                    if (bv_manager.is_boundary(u)) {
                        bv_manager.move(g, p_manager, u, u_id, best_id);
                    } else {
                        bv_manager.add_new(g, p_manager, u, best_id);
                    }

                    q_graph.move(g, p_manager, u, u_id, best_id);
                    p_manager.move(u, g.weight(u), u_id, best_id);
                }
            }

            // move vertices of the right region
            for (size_t i = 0; i < right_region_size; ++i) {
                vertex_t    u        = right_region[i];
                partition_t u_id     = p_manager[u];
                weight_t    u_weight = g.weight(u);
                vertex_t    new_u    = translation_table.get_n(u);

                ASSERT(r_start <= u_id && u_id < r_start + ids_per_superblock);
                ASSERT(is_right_region[u] == is_region_mark);
                ASSERT(new_u < left_region_size + right_region_size);

                if (is_left[new_u] == 1) {
                    // determine the best left block to move u to, that will not be overloaded
                    partition_t      best_id        = l_start;
                    s64              best_qap_delta = get_u_qap_delta(g, u, u_id, l_start, p_manager, d_oracle);
                    weight_t         best_weight    = p_manager.get_bweight(l_start) + u_weight;
                    for (partition_t id             = l_start; id < l_start + ids_per_superblock; ++id) {
                        // if (p_manager.get_bweight(id) + u_weight > m_lmax) { continue; }
                        s64 qap_delta = get_u_qap_delta(g, u, u_id, id, p_manager, d_oracle);
                        if (qap_delta > best_qap_delta || (qap_delta == best_qap_delta && p_manager.get_bweight(id) + u_weight < best_weight)) {
                            best_id        = id;
                            best_qap_delta = qap_delta;
                            best_weight    = p_manager.get_bweight(id) + u_weight;
                        }
                    }

                    if (bv_manager.is_boundary(u)) {
                        bv_manager.move(g, p_manager, u, u_id, best_id);
                    } else {
                        bv_manager.add_new(g, p_manager, u, best_id);
                    }

                    q_graph.move(g, p_manager, u, u_id, best_id);
                    p_manager.move(u, g.weight(u), u_id, best_id);
                }
            }
        }

        void rebalance(graph_t &g,
                       d_oracle_t &d_oracle,
                       bv_manager_t &bv_manager,
                       p_manager_t &p_manager,
                       q_graph_t &q_graph,
                       partition_t id_start,
                       partition_t ids_per_superblock) {

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

#endif //HIERARCHY_AWARE_FLOW_BASED_REFINEMENT_H
