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
        AlignedArray<u8> active_this_round;
        AlignedArray<u8> active_next_round;
        AlignedArray<PairWeight> pairs;
        size_t pairs_size = 0;

        // array for boundary vertices
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

        AlignedArray<u32> vertex_used;
        u32 vertex_marker = 0;

        AlignedArray<u32> block_used;
        u32 block_marker = 0;

        AlignedArray<Move> moves;
        size_t moves_size = 0;

        FlowNetwork flow_network;
        ResidualFlowNetwork residual_flow_network;
        SCCGraph scc_graph;

        RandomEngine* random_engine                                  = nullptr;
        const HierarchyAwareFlowBasedRefinementConfiguration* config = nullptr;
        StatisticCollector* m_stat_collector                         = nullptr;

    public:
        HierarchyAwareFlowBasedRefinement() = default;

        ~HierarchyAwareFlowBasedRefinement() override = default;

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

            // active block scheduling
            active_this_round.initialize(m_k);
            active_next_round.initialize(m_k);
            size_t size = (size_t)m_k * (size_t)m_k;
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
                    graph_t& g,
                    d_oracle_t& d_oracle,
                    bv_manager_t& bv_manager,
                    p_manager_t& p_manager,
                    q_graph_t& q_graph) override {
            for (size_t iteration = 0; iteration < config->max_global_iteration; ++iteration) {
                for (size_t i = m_hierarchy.size() - 1; i < m_hierarchy.size(); ++i) {
                    refine_layer(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, m_hierarchy.size() - 1 - i);
                }
            }
        }

        bool is_connected(const graph_t& g,
                          bv_manager_t& bv_manager,
                          p_manager_t& p_manager,
                          partition_t l_start,
                          partition_t r_start,
                          partition_t ids_per_super_block) {
            for (partition_t id = l_start; id < l_start + ids_per_super_block; ++id) {
                forall_bv_id_iu(bv_manager, id, i, u)
                    {
                        forall_guiv(g, u, j, v)
                            {
                                partition_t v_id = p_manager[v];

                                if (r_start <= v_id && v_id < r_start + ids_per_super_block) {
                                    return true;
                                }
                            }
                        endfor
                    }
                endfor
            }
            return false;
        }

        void refine_layer(const u64 level,
                          const u64 max_level,
                          graph_t& g,
                          d_oracle_t& d_oracle,
                          bv_manager_t& bv_manager,
                          p_manager_t& p_manager,
                          q_graph_t& q_graph,
                          size_t layer) {
            partition_t n_upper_total_super_blocks = 1;
            partition_t n_total_super_blocks       = m_hierarchy[layer];
            for (size_t i = layer + 1; i < m_hierarchy.size(); ++i) { n_upper_total_super_blocks *= m_hierarchy[i]; }
            partition_t ids_per_super_block = m_k / (n_upper_total_super_blocks * m_hierarchy[layer]);

            for (size_t i = 0; i < n_upper_total_super_blocks; ++i) {
                for (size_t j = 0; j < n_total_super_blocks; ++j) {
                    for (size_t k = j + 1; k < n_total_super_blocks; ++k) {
                        partition_t l_start = i * (n_total_super_blocks * ids_per_super_block) + j * ids_per_super_block;
                        partition_t r_start = i * (n_total_super_blocks * ids_per_super_block) + k * ids_per_super_block;

                        refine_blocks(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, l_start, r_start, ids_per_super_block, m_lmax * ids_per_super_block);
                    }
                }
            }
        }

        void refine_blocks(const u64 level,
                           const u64 max_level,
                           graph_t& g,
                           d_oracle_t& d_oracle,
                           bv_manager_t& bv_manager,
                           p_manager_t& p_manager,
                           q_graph_t& q_graph,
                           partition_t l_start,
                           partition_t r_start,
                           partition_t ids_per_super_block,
                           weight_t lmax) {
            s64 qap_before = 0; // get_qap(g, p_manager, d_oracle);

            // get allowed weight on each side
            weight_t l_lmax = lmax;
            weight_t r_lmax = lmax;
            for (partition_t id = l_start; id < l_start + ids_per_super_block; ++id) { l_lmax -= p_manager.get_bweight(id); }
            for (partition_t id = r_start; id < r_start + ids_per_super_block; ++id) { r_lmax -= p_manager.get_bweight(id); }

            // get boundary vertices
            determine_boundary_vertices(g, bv_manager, p_manager, l_start, r_start, ids_per_super_block);

            // get regions
            determine_regions(g, p_manager, l_start, r_start, ids_per_super_block, l_lmax, r_lmax);

            if (left_region_size + right_region_size <= 10) { return; }

            weight_t l_region_weight = 0;
            weight_t r_region_weight = 0;
            for (size_t i = 0; i < left_region_size; ++i) { l_region_weight += g.weight(left_region[i]); }
            for (size_t i = 0; i < right_region_size; ++i) { r_region_weight += g.weight(right_region[i]); }

            // build a translation table from graph to flow network
            vertex_t new_u = 0;
            for (size_t i = 0; i < left_region_size; ++i) { translation_table.add(left_region[i], new_u++); }
            for (size_t i = 0; i < right_region_size; ++i) { translation_table.add(right_region[i], new_u++); }

            // build flownetwork
            build_flow_network(g);

            // solve the flow network
            flow_network.solve();

            std::vector<u8> is_left;
            flow_network.get_cut(is_left);

            vertex_t n_vertices = 0;
            for (partition_t id = 0; id < ids_per_super_block; ++id) {
                n_vertices += p_manager.size(l_start + id) + p_manager.size(r_start + id);
            }
            vertex_t n_vertices_region = left_region_size + right_region_size;
            vertex_t n_changes         = get_n_changes(is_left);

            // std::cout << n_changes << " vertices " << ((f64) n_changes / n_vertices_region)*100.0 <<"% " << ((f64) n_changes / (f64) n_vertices)*100.0 << "% would change partition" << std::endl;

            // check if the cut actually changes the partition
            if (n_changes > 0) {
                // make the changes
                change_boundary(g, d_oracle, bv_manager, p_manager, q_graph, is_left, l_start, r_start, ids_per_super_block);
            }

            s64 qap_after = 0; // get_qap(g, p_manager, d_oracle);

            // std::cout << l_start << " " << r_start << " " << ids_per_super_block << " " << qap_before << " " << qap_after << " " << qap_before - qap_after << std::endl;
        }

        void determine_boundary_vertices(const graph_t& g,
                                         bv_manager_t& bv_manager,
                                         const p_manager_t& p_manager,
                                         partition_t l_start,
                                         partition_t r_start,
                                         partition_t ids_per_super_block) {
            left_boundary_size = 0;
            for (partition_t id = l_start; id < l_start + ids_per_super_block; ++id) {
                forall_bv_id_iu(bv_manager, id, i, u)
                    {
                        forall_guiv(g, u, i, v)
                            {
                                partition_t v_id = p_manager[v];
                                if (r_start <= v_id && v_id < r_start + ids_per_super_block) {
                                    left_boundary[left_boundary_size++] = u;
                                }
                                break;
                            }
                        endfor
                    }
                endfor
            }

            right_boundary_size = 0;
            for (partition_t id = r_start; id < r_start + ids_per_super_block; ++id) {
                forall_bv_id_iu(bv_manager, id, i, u)
                    {
                        forall_guiv(g, u, i, v)
                            {
                                partition_t v_id = p_manager[v];
                                if (l_start <= v_id && v_id < l_start + ids_per_super_block) {
                                    right_boundary[right_boundary_size++] = u;
                                }
                                break;
                            }
                        endfor
                    }
                endfor
            }
        }

        void determine_regions(const graph_t& g,
                               const p_manager_t& p_manager,
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
            queue_size       = 0;
            for (size_t i = 0; i < left_boundary_size; ++i) {
                vertex_t u = left_boundary[i];
                queue[queue_size++] = u;
                seen[u]             = seen_mark - 1;
            }

            // empty queue in bfs fashion and add vertices in the same superblock
            left_region_size = 0;
            while (queue_idx < queue_size) {
                vertex_t u = queue[queue_idx++];

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
                seen[u] = seen_mark;
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
                seen[u] = seen_mark;
            }
        }

        void determine_penalties(const graph_t& g,
                                 const p_manager_t& p_manager,
                                 d_oracle_t& d_oracle,
                                 partition_t l_start,
                                 partition_t r_start,
                                 partition_t ids_per_super_block) {
            for (size_t i = 0; i < left_boundary_size; ++i) {
                vertex_t u         = left_boundary[i];
                partition_t u_id   = p_manager[u];
                left_penalties[u]  = 0;
                right_penalties[u] = 0;

                forall_guivw(g, u, j, v, w)
                    {
                        if (is_left_region[v] == is_region_mark || is_right_region[v] == is_region_mark) { continue; }

                        partition_t v_id = p_manager[v];

                        if (l_start <= v_id && v_id < l_start + ids_per_super_block) {
                            // u is in left superblock and v is in left superblock

                            // if u stays in the left superblock, we pay the current cost
                            left_penalties[u] += w * d_oracle.get(u_id, v_id);

                            // if u goes to the right superblock, we do not know the final distance
                            weight_t max_distance = -1;
                            weight_t min_distance = std::numeric_limits<weight_t>::max();

                            for (partition_t id = r_start; id < r_start + ids_per_super_block; ++id) {
                                weight_t dist = d_oracle.get(id, v_id);
                                max_distance  = std::max(max_distance, dist);
                                min_distance  = std::min(min_distance, dist);
                            }

                            // using max distance is pessimistic and leads to fewer moves
                            // using min distance is optimistic and leads to more movement
                            right_penalties[u] += w * min_distance;
                        } else if (r_start <= v_id && v_id < r_start + ids_per_super_block) {
                            // u is in the left superblock and v is in the right superblock

                            // if u stays in the left superblock, we pay the current cost
                            left_penalties[u] += w * d_oracle.get(u_id, v_id);

                            // if u goes to the right superblock, we do not know the final distance
                            weight_t max_distance = -1;
                            weight_t min_distance = std::numeric_limits<weight_t>::max();

                            for (partition_t id = r_start; id < r_start + ids_per_super_block; ++id) {
                                weight_t dist = d_oracle.get(id, v_id);
                                max_distance  = std::max(max_distance, dist);
                                min_distance  = std::min(min_distance, dist);
                            }

                            // using max distance is pessimistic and leads to fewer moves
                            // using min distance is optimistic and leads to more movement
                            right_penalties[u] += w * min_distance;
                        }
                    }
                endfor
                left_penalties[u] *= 2;
                right_penalties[u] *= 2;
            }

            for (size_t i = 0; i < right_boundary_size; ++i) {
                vertex_t u         = right_boundary[i];
                partition_t u_id   = p_manager[u];
                left_penalties[u]  = 0;
                right_penalties[u] = 0;

                forall_guivw(g, u, j, v, w)
                    {
                        if (is_left_region[v] == is_region_mark || is_right_region[v] == is_region_mark) { continue; }

                        partition_t v_id = p_manager[v];

                        if (l_start <= v_id && v_id < l_start + ids_per_super_block) {
                            // u is in right superblock and v is in right superblock

                            // if u stays in the right superblock, we pay the current cost
                            right_penalties[u] += w * d_oracle.get(u_id, v_id);

                            // if u goes to the left superblock, we do not know the final distance
                            weight_t max_distance = -1;
                            weight_t min_distance = std::numeric_limits<weight_t>::max();

                            for (partition_t id = l_start; id < l_start + ids_per_super_block; ++id) {
                                weight_t dist = d_oracle.get(id, v_id);
                                max_distance  = std::max(max_distance, dist);
                                min_distance  = std::min(min_distance, dist);
                            }

                            // using max distance is pessimistic and leads to fewer moves
                            // using min distance is optimistic and leads to more movement
                            left_penalties[u] += w * min_distance;
                        } else if (r_start <= v_id && v_id < r_start + ids_per_super_block) {
                            // u is in the right superblock and v is in the right superblock

                            // if u stays in the right superblock, we pay the current cost
                            right_penalties[u] += w * d_oracle.get(u_id, v_id);

                            // if u goes to the left superblock, we do not know the final distance
                            weight_t max_distance = -1;
                            weight_t min_distance = std::numeric_limits<weight_t>::max();

                            for (partition_t id = l_start; id < l_start + ids_per_super_block; ++id) {
                                weight_t dist = d_oracle.get(id, v_id);
                                max_distance  = std::max(max_distance, dist);
                                min_distance  = std::min(min_distance, dist);
                            }

                            // using max distance is pessimistic and leads to fewer moves
                            // using min distance is optimistic and leads to more movement
                            left_penalties[u] += w * min_distance;
                        }
                    }
                endfor
                left_penalties[u] *= 2;
                right_penalties[u] *= 2;
            }
        }


        void build_flow_network(const graph_t& g) {
            // build flownetwork
            size_t n = left_region_size + right_region_size;
            flow_network.initialize(n);

            // build the left region
            for (size_t i = 0; i < left_region_size; ++i) {
                vertex_t u = left_region[i];

                forall_guivw(g, u, j, v, w)
                    {
                        if (u < v) { continue; } // no double edges
                        if (is_left_region[v] == is_region_mark || is_right_region[v] == is_region_mark) {
                            vertex_t new_u = translation_table.get_n(u);
                            vertex_t new_v = translation_table.get_n(v);
                            flow_network.add(new_u, new_v, w);
                        }
                    }
                endfor
            }

            // build the right region
            for (size_t i = 0; i < right_region_size; ++i) {
                vertex_t u = right_region[i];

                forall_guivw(g, u, j, v, w)
                    {
                        if (u < v) { continue; } // no double edges
                        if (is_left_region[v] == is_region_mark || is_right_region[v] == is_region_mark) {
                            vertex_t new_u = translation_table.get_n(u);
                            vertex_t new_v = translation_table.get_n(v);
                            flow_network.add(new_u, new_v, w);
                        }
                    }
                endfor
            }
        }


        size_t get_n_changes(std::vector<u8>& is_left) {
            size_t count = 0;
            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u     = left_region[j];
                vertex_t new_u = translation_table.get_n(u);
                if (is_left[new_u] == 0) {
                    count += 1;
                }
            }
            for (size_t j = 0; j < right_region_size; ++j) {
                vertex_t u     = right_region[j];
                vertex_t new_u = translation_table.get_n(u);
                if (is_left[new_u] == 1) {
                    count += 1;
                }
            }

            return count;
        }

        void change_boundary(const graph_t& g,
                             d_oracle_t& d_oracle,
                             bv_manager_t& bv_manager,
                             p_manager_t& p_manager,
                             q_graph_t& q_graph,
                             std::vector<u8>& is_left,
                             partition_t l_start,
                             partition_t r_start,
                             partition_t ids_per_superblock) {
            // move all vertices of the left region
            for (size_t i = 0; i < left_region_size; ++i) {
                vertex_t u        = left_region[i];
                partition_t u_id  = p_manager[u];
                weight_t u_weight = g.weight(u);
                vertex_t new_u    = translation_table.get_n(u);

                ASSERT(l_start <= u_id && u_id < l_start + ids_per_superblock);
                ASSERT(is_left_region[u] == is_region_mark);
                ASSERT(new_u < left_region_size + right_region_size);

                print(is_left);
                if (is_left[new_u] == 0) {
                    // determine the best right block to move u to, that will not be overloaded
                    partition_t best_id  = r_start;
                    s64 best_qap_delta   = get_u_qap_delta(g, u, u_id, r_start, p_manager, d_oracle);
                    weight_t best_weight = p_manager.get_bweight(r_start) + u_weight;
                    for (partition_t id = r_start; id < r_start + ids_per_superblock; ++id) {
                        if (p_manager.get_bweight(id) + u_weight > m_lmax) { continue; }
                        s64 qap_delta = get_u_qap_delta(g, u, u_id, id, p_manager, d_oracle);
                        if (qap_delta > best_qap_delta || (qap_delta == best_qap_delta && p_manager.get_bweight(id) + u_weight < best_weight)) {
                            best_id        = id;
                            best_qap_delta = qap_delta;
                            best_weight    = p_manager.get_bweight(id) + u_weight;
                        }
                    }
                    std::cout << best_qap_delta << std::endl;

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
                vertex_t u        = right_region[i];
                partition_t u_id  = p_manager[u];
                weight_t u_weight = g.weight(u);
                vertex_t new_u    = translation_table.get_n(u);

                ASSERT(r_start <= u_id && u_id < r_start + ids_per_superblock);
                ASSERT(is_right_region[u] == is_region_mark);
                ASSERT(new_u < left_region_size + right_region_size);

                if (is_left[new_u] == 1) {
                    // determine the best left block to move u to, that will not be overloaded
                    partition_t best_id  = l_start;
                    s64 best_qap_delta   = get_u_qap_delta(g, u, u_id, l_start, p_manager, d_oracle);
                    weight_t best_weight = p_manager.get_bweight(l_start) + u_weight;
                    for (partition_t id = l_start; id < l_start + ids_per_superblock; ++id) {
                        if (p_manager.get_bweight(id) + u_weight > m_lmax) { continue; }
                        s64 qap_delta = get_u_qap_delta(g, u, u_id, id, p_manager, d_oracle);
                        if (qap_delta > best_qap_delta || (qap_delta == best_qap_delta && p_manager.get_bweight(id) + u_weight < best_weight)) {
                            best_id        = id;
                            best_qap_delta = qap_delta;
                            best_weight    = p_manager.get_bweight(id) + u_weight;
                        }
                    }

                    // std::cout << best_qap_delta << std::endl;

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
