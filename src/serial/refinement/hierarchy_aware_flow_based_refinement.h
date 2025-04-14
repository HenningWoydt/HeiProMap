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
#include "../datastructures/functions.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/qap.h"

namespace HeiProMap {
    class HierarchyAwareFlowBasedRefinementConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit HierarchyAwareFlowBasedRefinementConfiguration(const std::string& t_name) : ISerialRefinerConfiguration(t_name) {}
        u64 max_global_iteration       = 1;
        u64 max_local_iteration        = 3;
        f64 alpha                      = 2.0;
        f64 alpha_upper_bound          = 8.0;
        f64 alpha_modifier             = 2.0;
        bool use_closed_vertex_set     = true;
        u64 closed_vertex_sets_repeats = 10;
    };

    class HierarchyAwareFlowBasedRefinement final : public ISerialRefiner {
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

        u32* vertex_used  = nullptr;
        u32 vertex_marker = 0;

        u32* block_used  = nullptr;
        u32 block_marker = 0;

        Move* moves       = nullptr;
        size_t moves_size = 0;

        FlowNetwork flow_network;
        ResidualFlowNetwork residual_flow_network;
        SCCGraph scc_graph;

        RandomEngine* random_engine                                  = nullptr;
        const HierarchyAwareFlowBasedRefinementConfiguration* config = nullptr;
        StatisticCollector* m_stat_collector                         = nullptr;

    public:
        HierarchyAwareFlowBasedRefinement() = default;

        ~HierarchyAwareFlowBasedRefinement() override {
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
            config           = dynamic_cast<const HierarchyAwareFlowBasedRefinementConfiguration*>(&i_config);
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

            vertex_t t_n_64 = round_up_64(t_n);
            vertex_t t_k_64 = round_up_64(t_k);

            vertex_used = (u32*)aligned_alloc(64, t_n_64 * sizeof(u32));
            std::fill_n(vertex_used, t_n_64, 0);
            vertex_marker = 0;

            block_used = (u32*)aligned_alloc(64, t_k_64 * sizeof(u32));
            std::fill_n(block_used, t_k_64, 0);
            block_marker = 0;

            moves      = (Move*)aligned_alloc(64, t_n_64 * sizeof(Move));
            moves_size = 0;
        }

        void refine(const u64 level,
                    const u64 max_level,
                    const graph_t& g,
                    const d_oracle_t& d_oracle,
                    bv_manager_t& bv_manager,
                    p_manager_t& p_manager,
                    q_graph_t& q_graph) override {
            for (size_t iteration = 0; iteration < config->max_global_iteration; ++iteration) {
                for (size_t i = 0; i < m_hierarchy.size() - 1; ++i) {
                    refine_layer(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, m_hierarchy.size() - 1 - i);
                    // rebalance_layer(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, m_hierarchy.size() - 2 - i);
                }
            }
        }

        void refine_layer(const u64 level,
                          const u64 max_level,
                          const graph_t& g,
                          const d_oracle_t& d_oracle,
                          bv_manager_t& bv_manager,
                          p_manager_t& p_manager,
                          q_graph_t& q_graph,
                          size_t layer) {
            partition_t n_total_super_blocks = 1;
            for (size_t i = layer; i < m_hierarchy.size(); ++i) { n_total_super_blocks *= m_hierarchy[i]; }
            partition_t ids_per_super_block = m_k / n_total_super_blocks;

            partition_t n_local_super_blocks = m_hierarchy[layer];
            partition_t n_upper_blocks       = 1;
            for (size_t i = layer + 1; i < m_hierarchy.size(); ++i) { n_upper_blocks *= m_hierarchy[i]; }

            partition_t neighborhood_id_start = 0;
            partition_t neighborhood_id_end   = ids_per_super_block * n_local_super_blocks;

            for (size_t upper_block_id = 0; upper_block_id < n_upper_blocks; ++upper_block_id) {
                refine_neighborhood(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, neighborhood_id_start, neighborhood_id_end, n_local_super_blocks, ids_per_super_block);

                neighborhood_id_start = neighborhood_id_end;
                neighborhood_id_end += ids_per_super_block * n_local_super_blocks;
            }
        }

        void refine_neighborhood(const u64 level,
                                 const u64 max_level,
                                 const graph_t& g,
                                 const d_oracle_t& d_oracle,
                                 bv_manager_t& bv_manager,
                                 p_manager_t& p_manager,
                                 q_graph_t& q_graph,
                                 partition_t neighborhood_id_start,
                                 partition_t neighborhood_id_end,
                                 partition_t n_local_super_blocks,
                                 partition_t ids_per_super_block) {
            f64 alpha             = config->alpha;
            f64 alpha_upper_bound = config->alpha_upper_bound;
            f64 alpha_modifier    = config->alpha_modifier;

            u64 max_local_iteration = config->max_local_iteration;
            u64 iteration           = 0;

            while (iteration < max_local_iteration) {
                ASSERT(max(p_manager.get_bweights()) <= m_lmax);
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

                bool            qap_normal_calculated = false;
                weight_t        qap_normal_change;
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

                        qap_normal_calculated = true;
                        // qap_normal_change     = qap - get_qap(g, p_manager, d_oracle);

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
                    bool     closure_found           = scc_graph.find_best_closure(left_non_region_weight, right_non_region_weight, m_lmax, config->closed_vertex_sets_repeats, *random_engine, is_left);
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

                    if (left_region_size == 5 && right_region_size == 5 && scc_graph.get_n_scc() >=4 && false) {
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
                            vertex_t    u    = left_region[i];
                            weight_t l_penalty = left_penalties[u];
                            weight_t r_penalty = right_penalties[u];
                            std::cout << u << " : " << l_penalty << " -- " << r_penalty << std::endl;
                        }

                        std::cout << "Right-Region Penalties" << std::endl;
                        for (size_t i = 0; i < right_region_size; ++i) {
                            vertex_t    u    = right_region[i];
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
                ASSERT(max(p_manager.get_bweights()) <= m_lmax);
                change_boundary(g, bv_manager, p_manager, q_graph, is_left, left_id, right_id);
                ASSERT(max(p_manager.get_bweights()) <= m_lmax);
                // HEAVYASSERT(assert_correct_boundary(g, p_manager, bv_manager, m_k));

                active_next_round[left_id]  = 1;
                active_next_round[right_id] = 1;
            }
            ASSERT(max(p_manager.get_bweights()) <= m_lmax);
        }

        void rebalance_layer(const u64 level,
                             const u64 max_level,
                             const graph_t& g,
                             const d_oracle_t& d_oracle,
                             bv_manager_t& bv_manager,
                             p_manager_t& p_manager,
                             q_graph_t& q_graph,
                             size_t layer) {
            partition_t n_total_super_blocks = 1;
            for (size_t i = layer; i < m_hierarchy.size(); ++i) { n_total_super_blocks *= m_hierarchy[i]; }
            partition_t ids_per_super_block = m_k / n_total_super_blocks;

            partition_t n_local_super_blocks = m_hierarchy[layer];
            partition_t n_upper_blocks       = 1;
            for (size_t i = layer + 1; i < m_hierarchy.size(); ++i) { n_upper_blocks *= m_hierarchy[i]; }

            partition_t neighborhood_id_start = 0;
            partition_t neighborhood_id_end   = ids_per_super_block * n_local_super_blocks;

            for (size_t upper_block_id = 0; upper_block_id < n_upper_blocks; ++upper_block_id) {
                rebalance_neighborhoods(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, neighborhood_id_start, neighborhood_id_end, n_local_super_blocks, ids_per_super_block);

                neighborhood_id_start = neighborhood_id_end;
                neighborhood_id_end += ids_per_super_block * n_local_super_blocks;
            }
        }

        void rebalance_neighborhoods(const u64 level,
                                     const u64 max_level,
                                     const graph_t& g,
                                     const d_oracle_t& d_oracle,
                                     bv_manager_t& bv_manager,
                                     p_manager_t& p_manager,
                                     q_graph_t& q_graph,
                                     partition_t neighborhood_id_start,
                                     partition_t neighborhood_id_end,
                                     partition_t n_local_super_blocks,
                                     partition_t ids_per_super_block) {
            weight_t blocks_lmax = (weight_t)ids_per_super_block * m_lmax;;
            std::vector<weight_t> blocks_weights(n_local_super_blocks, 0);

            for (partition_t i = 0; i < n_local_super_blocks; ++i) {
                for (partition_t j = 0; j < ids_per_super_block; ++j) {
                    partition_t id = neighborhood_id_start + i * ids_per_super_block + j;
                    blocks_weights[i] += p_manager.get_bweight(id);
                }
            }

            std::vector<weight_t> min_blocks_weights(n_local_super_blocks, 0);
            for (partition_t i = 0; i < n_local_super_blocks; ++i) {
                min_blocks_weights[i] = std::max(blocks_weights[i], blocks_lmax);
            }

            std::vector<size_t> indices(n_local_super_blocks);
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(), [&](size_t i, size_t j) { return blocks_weights[i] > blocks_weights[j]; });

            if (blocks_weights[indices[0]] <= blocks_lmax) { return; }

            while (max(blocks_weights) > blocks_lmax) {
                std::sort(indices.begin(), indices.end(), [&](size_t i, size_t j) { return blocks_weights[i] > blocks_weights[j]; });

                bool move_occurred = false;
                for (size_t i = 0; i < indices.size(); ++i) {
                    partition_t id = indices[i];
                    if (blocks_weights[id] <= blocks_lmax) { continue; }

                    vertex_t vertex_move;
                    partition_t vertex_curr_id;
                    partition_t vertex_new_id;
                    weight_t vertex_weight;
                    s64 vertex_qap_delta = -std::numeric_limits<s64>::max();

                    partition_t start = neighborhood_id_start + id * ids_per_super_block;
                    partition_t end   = neighborhood_id_start + (id + 1) * ids_per_super_block;
                    for (partition_t u_id = start; u_id < end; ++u_id) {
                        forall_bv_id_iu(bv_manager, u_id, k, u)
                        {
                            weight_t u_weight = g.weight(u);

                            block_marker += 1;
                            for (partition_t v_id = neighborhood_id_start; v_id < neighborhood_id_end; ++v_id) {
                                if (u_id == v_id) { continue; }
                                if (block_used[v_id] == block_marker) { continue; }

                                partition_t id1 = (u_id - neighborhood_id_start) / ids_per_super_block;
                                partition_t id2 = (v_id - neighborhood_id_start) / ids_per_super_block;

                                if (blocks_weights[id2] + u_weight >= blocks_weights[id1]) { continue; }

                                s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);

                                if (qap_delta > vertex_qap_delta) {
                                    vertex_qap_delta = qap_delta;
                                    vertex_move      = u;
                                    vertex_curr_id   = u_id;
                                    vertex_new_id    = v_id;
                                    vertex_weight    = u_weight;
                                }

                                block_used[v_id] = block_marker;
                            }
                        }
                        endfor
                    }

                    if (vertex_qap_delta != -std::numeric_limits<s64>::max()) {
                        partition_t id1 = (vertex_curr_id - neighborhood_id_start) / ids_per_super_block;
                        partition_t id2 = (vertex_new_id - neighborhood_id_start) / ids_per_super_block;
                        blocks_weights[id1] -= vertex_weight;
                        blocks_weights[id2] += vertex_weight;

                        min_blocks_weights[id1] = std::min(min_blocks_weights[id1], blocks_weights[id1]);
                        min_blocks_weights[id2] = std::min(min_blocks_weights[id2], blocks_weights[id2]);
                        min_blocks_weights[id1] = std::max(min_blocks_weights[id1], blocks_lmax);
                        min_blocks_weights[id2] = std::max(min_blocks_weights[id2], blocks_lmax);

                        bv_manager.move(g, p_manager, vertex_move, vertex_curr_id, vertex_new_id);
                        q_graph.move(g, p_manager, vertex_move, vertex_curr_id, vertex_new_id);
                        p_manager.move(vertex_move, vertex_weight, vertex_curr_id, vertex_new_id);
                        move_occurred = true;
                        break;
                    }
                }

                if (!move_occurred) {
                    // there is not one available move
                    return;
                }
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
                ASSERT(max(p_manager.get_bweights()) <= m_lmax);
                iteration += 1;

                // get boundary vertices
                determine_boundary_vertices(g, bv_manager, p_manager, left_id, right_id);

                // calc max weight for each bfs
                weight_t lmax             = std::ceil((1.0 + (m_imbalance * alpha)) * ((f64)g.weight() / (f64)m_k));
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
                        weight_t qap            = get_qap(g, p_manager, d_oracle);
                        std::vector<u8> changed = change_boundary(g, bv_manager, p_manager, q_graph, is_left, left_id, right_id);
                        HEAVYASSERT(assert_correct_boundary(g, p_manager, bv_manager, m_k));

                        qap_normal_calculated = true;
                        // qap_normal_change     = qap - get_qap(g, p_manager, d_oracle);

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
                    bool closure_found               = scc_graph.find_best_closure(left_non_region_weight, right_non_region_weight, m_lmax, config->closed_vertex_sets_repeats, *random_engine, is_left);
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
                            vertex_t u       = left_region[i];
                            weight_t u_w     = g.weight(u);
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
                            vertex_t u       = right_region[i];
                            weight_t u_w     = g.weight(u);
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
                ASSERT(max(p_manager.get_bweights()) <= m_lmax);
                change_boundary(g, bv_manager, p_manager, q_graph, is_left, left_id, right_id);
                ASSERT(max(p_manager.get_bweights()) <= m_lmax);
                // HEAVYASSERT(assert_correct_boundary(g, p_manager, bv_manager, m_k));

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
                vertex_t u = left_boundary[i];
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
                vertex_t u = right_boundary[i];
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
                        if (bfs_level[u] < bfs_level[v]) { continue; } // only forward edges allowed

                        weight_t mult = bfs_level[u] == bfs_level[v] ? 1 : 2;

                        vertex_t new_u = translation_table.get_n(u);
                        vertex_t new_v = translation_table.get_n(v);

                        ASSERT(new_u != new_v);

                        flow_network.add(new_u, new_v, mult * w * distance);
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
                        if (bfs_level[u] > bfs_level[v]) { continue; } // only forward edges allowed

                        weight_t mult = bfs_level[u] == bfs_level[v] ? 1 : 2;

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

        std::vector<u8> change_boundary(const graph_t& g,
                                        bv_manager_t& bv_manager,
                                        p_manager_t& p_manager,
                                        q_graph_t& q_graph,
                                        std::vector<u8>& is_left,
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

        void revert_boundary(const graph_t& g,
                             bv_manager_t& bv_manager,
                             p_manager_t& p_manager,
                             q_graph_t& q_graph,
                             std::vector<u8>& changed,
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
