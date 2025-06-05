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

#ifndef HEIPROMAP_HIERARCHY_AWARE_QUOTIENT_GRAPH_REFINEMENT_H
#define HEIPROMAP_HIERARCHY_AWARE_QUOTIENT_GRAPH_REFINEMENT_H

#include <queue>

#include "ISerialRefiner.h"
#include "../../commons/random_engine.h"
#include "../../commons/statistic_collector.h"
#include "../../commons/utils.h"
#include "../utility/functions.h"
#include "../utility/qap.h"

namespace HeiProMap {
    class HierarchyAwareQuotientGraphRefinementConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit HierarchyAwareQuotientGraphRefinementConfiguration(const std::string& t_name) : ISerialRefinerConfiguration(t_name) {}
    };

    class HierarchyAwareQuotientGraphRefinement final : public ISerialRefiner {
    private:
        vertex_t m_n    = 0;
        vertex_t m_m    = 0;
        partition_t m_k = 0;
        f64 m_imbalance = 0.0;
        weight_t m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;

        AlignedArray<u32> vertex_used;
        u32 vertex_marker = 0;

        AlignedArray<u32> block_used;
        u32 block_marker = 0;

        // priority queues
        std::priority_queue<KWayFMMove> boundary_vertices_sb1;
        std::priority_queue<KWayFMMove> boundary_vertices_sb2;

        // store change
        AlignedArray<Move> moves;
        size_t moves_size      = 0;
        s64 curr_qap_gain      = 0;
        s64 max_qap_gain       = 0;
        s64 curr_edge_cut_gain = 0;
        s64 max_edge_cut_gain  = 0;
        size_t best_idx        = 0;

        RandomEngine* random_engine                                      = nullptr;
        const HierarchyAwareQuotientGraphRefinementConfiguration* config = nullptr;
        StatisticCollector* m_stat_collector                             = nullptr;

    public:
        HierarchyAwareQuotientGraphRefinement() = default;

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
            m_imbalance = t_imbalance;
            m_lmax      = t_lmax;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;

            random_engine    = &t_random_engine;
            config           = dynamic_cast<const HierarchyAwareQuotientGraphRefinementConfiguration*>(&i_config);
            m_stat_collector = &t_stat_collect;

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
            for (size_t rep = 0; rep < 1; ++rep) {
                for (size_t i = 0; i < m_hierarchy.size(); ++i) {
                    refine_layer(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, m_hierarchy.size() - 1 - i);
                }
            }
        }

        void refine_layer(const u64 level,
                          const u64 max_level,
                          const graph_t& g,
                          d_oracle_t& d_oracle,
                          bv_manager_t& bv_manager,
                          p_manager_t& p_manager,
                          q_graph_t& q_graph,
                          size_t layer) {
            for (size_t rep = 0; rep < 3; ++rep) {
                partition_t n_total_super_blocks = 1;
                for (size_t i = layer; i < m_hierarchy.size(); ++i) { n_total_super_blocks *= m_hierarchy[i]; }
                partition_t ids_per_super_block = m_k / n_total_super_blocks;

                for (size_t i = 0; i < n_total_super_blocks; ++i) {
                    for (size_t j = i + 1; j < n_total_super_blocks; ++j) {
                        partition_t n1_start = i * ids_per_super_block;
                        partition_t n2_start = j * ids_per_super_block;

                        if (is_connected(g, bv_manager, p_manager, n1_start, n2_start, ids_per_super_block)) {
                            refine_blocks(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, n1_start, n2_start, ids_per_super_block, m_lmax * ids_per_super_block);

                            weight_t n1_weight = 0;
                            for (partition_t id = n1_start; id < n1_start + ids_per_super_block; ++id) {
                                n1_weight += p_manager.get_bweight(id);
                            }

                            weight_t n2_weight = 0;
                            for (partition_t id = n2_start; id < n2_start + ids_per_super_block; ++id) {
                                n2_weight += p_manager.get_bweight(id);
                            }

                            // std::cout << n1_start << " " << n2_start << " " << ids_per_super_block << " " << n1_weight << " " << n2_weight << " " << m_lmax * ids_per_super_block << std::endl;
                        }
                    }
                }
            }
        }

        bool is_connected(const graph_t& g,
                          bv_manager_t& bv_manager,
                          p_manager_t& p_manager,
                          partition_t n1_start,
                          partition_t n2_start,
                          partition_t ids_per_super_block) {
            for (partition_t id = n1_start; id < n1_start + ids_per_super_block; ++id) {
                forall_bv_id_iu(bv_manager, id, i, u)
                    {
                        forall_guiv(g, u, j, v)
                            {
                                partition_t v_id = p_manager[v];

                                if (n2_start <= v_id && v_id < n2_start + ids_per_super_block) {
                                    return true;
                                }
                            }
                        endfor
                    }
                endfor
            }
            return false;
        }

        bool is_connected(const graph_t& g,
                          p_manager_t& p_manager,
                          vertex_t u,
                          partition_t id_start,
                          partition_t ids_per_super_block) {
            forall_guiv(g, u, j, v)
                {
                    partition_t v_id = p_manager[v];

                    if (id_start <= v_id && v_id < id_start + ids_per_super_block) { return true; }
                }
            endfor

            return false;
        }

        void refine_blocks(const u64 level,
                           const u64 max_level,
                           const graph_t& g,
                           d_oracle_t& d_oracle,
                           bv_manager_t& bv_manager,
                           p_manager_t& p_manager,
                           q_graph_t& q_graph,
                           partition_t n1_start,
                           partition_t n2_start,
                           partition_t ids_per_super_block,
                           weight_t lmax) {
            f64 alpha = 1000.0;
            f64 beta  = std::log(g.get_n());

            /*
            std::cout << "Refining: [";
            for (partition_t id = n1_start; id < n1_start + ids_per_super_block; ++id) {
                std::cout << id << " ";
            }
            std::cout << "] - [";
            for (partition_t id = n2_start; id < n2_start + ids_per_super_block; ++id) {
                std::cout << id << " ";
            }
            std::cout << "] with lmax = " << lmax << std::endl;
            */

            // initialize the heaps
            boundary_vertices_sb1 = std::priority_queue<KWayFMMove>();
            for (partition_t id = n1_start; id < n1_start + ids_per_super_block; ++id) {
                forall_bv_id_iu(bv_manager, id, i, u)
                    {
                        block_marker += 1;
                        forall_guiv(g, u, j, v)
                            {
                                partition_t v_id = p_manager[v];
                                if (block_used[v_id] == block_marker) { continue; }

                                if (n2_start <= v_id && v_id < n2_start + ids_per_super_block) {
                                    s64 delta = get_u_qap_delta(g, u, id, v_id, p_manager, d_oracle);
                                    boundary_vertices_sb1.emplace(u, id, v_id, delta);
                                    block_used[v_id] = block_marker;
                                }
                            }
                        endfor
                    }
                endfor
            }

            boundary_vertices_sb2 = std::priority_queue<KWayFMMove>();
            for (partition_t id = n2_start; id < n2_start + ids_per_super_block; ++id) {
                forall_bv_id_iu(bv_manager, id, i, u)
                    {
                        block_marker += 1;
                        forall_guiv(g, u, j, v)
                            {
                                partition_t v_id = p_manager[v];
                                if (block_used[v_id] == block_marker) { continue; }

                                if (n1_start <= v_id && v_id < n1_start + ids_per_super_block) {
                                    s64 delta = get_u_qap_delta(g, u, id, v_id, p_manager, d_oracle);
                                    boundary_vertices_sb2.emplace(u, id, v_id, delta);
                                    block_used[v_id] = block_marker;
                                }
                            }
                        endfor
                    }
                endfor
            }

            // get the weight of both superblocks
            weight_t weight_sb1 = 0;
            for (partition_t id = n1_start; id < n1_start + ids_per_super_block; ++id) { weight_sb1 += p_manager.get_bweight(id); }
            weight_t weight_sb2 = 0;
            for (partition_t id = n2_start; id < n2_start + ids_per_super_block; ++id) { weight_sb2 += p_manager.get_bweight(id); }

            weight_t max_allowed_weight = std::max(weight_sb1, weight_sb2);
            bool start_was_unbalanced   = weight_sb1 > lmax || weight_sb2 > lmax;

            // start executing moves based on the TopGain method
            vertex_marker += 1;
            moves_size         = 0;
            best_idx           = 0;
            curr_qap_gain      = 0;
            max_qap_gain       = 0;
            size_t max_n_swaps = 100000000;

            f64 steps_since_last_improvement = 0.0;
            f64 qap_gain_mean                = 0.0;
            f64 qap_gain_var                 = 0.0;

            while ((!boundary_vertices_sb1.empty() || !boundary_vertices_sb2.empty()) && moves_size < max_n_swaps) {
                // determine from which block to choose
                bool choose_sb1 = true;
                // 1. if one block is empty, then choose the other one
                if (boundary_vertices_sb1.empty() || boundary_vertices_sb2.empty()) {
                    choose_sb1 = boundary_vertices_sb2.empty();
                } else {
                    // 2. choose the block with greater gain and randomly if even
                    if (boundary_vertices_sb2.top().qap_delta > boundary_vertices_sb1.top().qap_delta) {
                        choose_sb1 = false;
                    } else if (boundary_vertices_sb2.top().qap_delta == boundary_vertices_sb1.top().qap_delta) {
                        choose_sb1 = random_engine->get_f32() < 0.5;
                    }

                    // 3. if one block is overloaded, choose the larger one, if both same sizes, then randomly
                    if (weight_sb1 > lmax && weight_sb1 > weight_sb2) { choose_sb1 = true; }
                    if (weight_sb2 > lmax && weight_sb2 > weight_sb1) { choose_sb1 = false; }
                    if (weight_sb1 > lmax && weight_sb2 > lmax && weight_sb1 == weight_sb2) { choose_sb1 = random_engine->get_f32() < 0.5; }
                }

                // choose the priority queue
                std::priority_queue<KWayFMMove>& boundary_vertices = choose_sb1 ? boundary_vertices_sb1 : boundary_vertices_sb2;

                vertex_t vertex        = boundary_vertices.top().u;
                s64 qap_delta          = boundary_vertices.top().qap_delta;
                weight_t vertex_weight = g.weight(vertex);
                partition_t vertex_id  = boundary_vertices.top().u_id;
                partition_t move_id    = boundary_vertices.top().to_move_id;
                boundary_vertices.pop();

                // if vertex is not in the prev partition, throw away
                if (vertex_used[vertex] == vertex_marker) { continue; }
                if (p_manager[vertex] != vertex_id) { continue; }
                if (choose_sb1 ? !is_connected(g, p_manager, vertex, n2_start, ids_per_super_block) : !is_connected(g, p_manager, vertex, n1_start, ids_per_super_block)) { continue; }
                if (qap_delta != get_u_qap_delta(g, vertex, vertex_id, move_id, p_manager, d_oracle)) { continue; }

                // move the vertex
                moves[moves_size++] = {vertex, vertex_id, move_id};
                curr_qap_gain += qap_delta;
                if (choose_sb1) {
                    weight_sb1 -= vertex_weight;
                    weight_sb2 += vertex_weight;
                } else {
                    weight_sb1 += vertex_weight;
                    weight_sb2 -= vertex_weight;
                }

                bool balanced              = weight_sb1 <= lmax && weight_sb2 <= lmax;
                bool unbalanced_but_better = start_was_unbalanced && weight_sb1 <= max_allowed_weight && weight_sb2 <= max_allowed_weight;
                if (curr_qap_gain >= max_qap_gain && (balanced || unbalanced_but_better)) {
                    best_idx     = moves_size;
                    max_qap_gain = curr_qap_gain;

                    steps_since_last_improvement = 0.0;
                    qap_gain_mean                = 0.0;
                    qap_gain_var                 = 0.0;
                }

                // make move in structures
                p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                vertex_used[vertex] = vertex_marker;

                steps_since_last_improvement += 1.0;
                f64 new_qap_gain_mean = qap_gain_mean + ((f64)qap_delta - qap_gain_mean) / steps_since_last_improvement;
                f64 new_qap_gain_var  = (qap_gain_var + ((f64)qap_delta - qap_gain_mean) * ((f64)qap_delta - new_qap_gain_mean)) / steps_since_last_improvement;

                qap_gain_mean = new_qap_gain_mean;
                qap_gain_var  = new_qap_gain_var;

                if (steps_since_last_improvement > 2.0 && steps_since_last_improvement * qap_gain_mean * qap_gain_mean > alpha * qap_gain_var + beta) { break; }

                // we have to push or update the neighbors that were not moved already
                forall_guiv(g, vertex, i, neighbor)
                    {
                        if (vertex_used[neighbor] == vertex_marker) { continue; }
                        partition_t neighbor_id = p_manager[neighbor];

                        if (choose_sb1) {
                            if (!(n2_start <= neighbor_id && neighbor_id < n2_start + ids_per_super_block)) { continue; }

                            block_marker += 1;
                            forall_guiv(g, neighbor, j, v)
                                {
                                    partition_t v_id = p_manager[v];
                                    if (block_used[v_id] == block_marker) { continue; }

                                    if (n1_start <= v_id && v_id < n1_start + ids_per_super_block) {
                                        s64 delta = get_u_qap_delta(g, neighbor, neighbor_id, v_id, p_manager, d_oracle);
                                        boundary_vertices_sb2.emplace(neighbor, neighbor_id, v_id, delta);
                                        block_used[v_id] = block_marker;
                                    }
                                }
                            endfor
                        } else {
                            if (!(n1_start <= neighbor_id && neighbor_id < n1_start + ids_per_super_block)) { continue; }

                            block_marker += 1;
                            forall_guiv(g, neighbor, j, v)
                                {
                                    partition_t v_id = p_manager[v];
                                    if (block_used[v_id] == block_marker) { continue; }

                                    if (n2_start <= v_id && v_id < n2_start + ids_per_super_block) {
                                        s64 delta = get_u_qap_delta(g, neighbor, neighbor_id, v_id, p_manager, d_oracle);
                                        boundary_vertices_sb1.emplace(neighbor, neighbor_id, v_id, delta);
                                        block_used[v_id] = block_marker;
                                    }
                                }
                            endfor
                        }
                    }
                endfor

                // remove vertex from u if it is not boundary
                while (!boundary_vertices_sb1.empty() && !is_connected(g, p_manager, boundary_vertices_sb1.top().u, n2_start, ids_per_super_block)) { boundary_vertices_sb1.pop(); }

                // remove vertex from v if it is not boundary
                while (!boundary_vertices_sb2.empty() && !is_connected(g, p_manager, boundary_vertices_sb2.top().u, n1_start, ids_per_super_block)) { boundary_vertices_sb2.pop(); }
            }

            // revert all moves in partitioning manager
            for (size_t i = 0; i < moves_size; i++) {
                vertex_t vertex        = moves[moves_size - 1 - i].u;
                weight_t vertex_weight = g.weight(vertex);
                partition_t vertex_id  = moves[moves_size - 1 - i].to_move_id;
                partition_t move_id    = moves[moves_size - 1 - i].u_id;

                ASSERT(p_manager[vertex] == vertex_id);

                p_manager.move(vertex, vertex_weight, vertex_id, move_id);
            }

            // make all moves to best index
            for (size_t i = 0; i < best_idx; ++i) {
                vertex_t vertex        = moves[i].u;
                weight_t vertex_weight = g.weight(vertex);
                partition_t vertex_id  = moves[i].u_id;
                partition_t move_id    = moves[i].to_move_id;

                ASSERT(p_manager[vertex] == vertex_id);

                bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                p_manager.move(vertex, vertex_weight, vertex_id, move_id);
            }
        }

        JSONString get_stats() override { return {}; }
    };
}

#endif //HEIPROMAP_HIERARCHY_AWARE_QUOTIENT_GRAPH_REFINEMENT_H
