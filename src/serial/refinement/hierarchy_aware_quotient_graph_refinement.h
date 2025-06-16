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

        u64 max_iteration = 5;
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

        // priority queues
        IndexedMaxHeap<s64> l_boundary_vertices;
        IndexedMaxHeap<s64> r_boundary_vertices;

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

            // priority queues
            l_boundary_vertices.initialize(m_n);
            r_boundary_vertices.initialize(m_n);

            vertex_used.initialize(m_n, 0);
            vertex_marker = 0;

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
                          graph_t& g,
                          d_oracle_t& d_oracle,
                          bv_manager_t& bv_manager,
                          p_manager_t& p_manager,
                          q_graph_t& q_graph,
                          size_t layer) override {
            partition_t n_upper_total_super_blocks = 1;
            partition_t n_total_super_blocks       = m_hierarchy[layer];
            for (size_t i = layer + 1; i < m_hierarchy.size(); ++i) { n_upper_total_super_blocks *= m_hierarchy[i]; }
            partition_t ids_per_super_block = m_k / (n_upper_total_super_blocks * m_hierarchy[layer]);

            for (u64 iteration = 0; iteration < config->max_iteration; ++iteration) {
                for (size_t i = 0; i < n_upper_total_super_blocks; ++i) {
                    for (size_t j = 0; j < n_total_super_blocks; ++j) {
                        for (size_t k = j + 1; k < n_total_super_blocks; ++k) {
                            partition_t l_start = i * (n_total_super_blocks * ids_per_super_block) + j * ids_per_super_block;
                            partition_t r_start = i * (n_total_super_blocks * ids_per_super_block) + k * ids_per_super_block;

                            if (!is_connected(g, bv_manager, p_manager, l_start, r_start, ids_per_super_block)) { continue; }

                            refine_blocks(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, l_start, r_start, ids_per_super_block, m_lmax * ids_per_super_block);
                        }
                    }
                }
            }
        }

        static bool is_connected(const graph_t& g,
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

                                if (r_start <= v_id && v_id < r_start + ids_per_super_block) { return true; }
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
                           partition_t l_start,
                           partition_t r_start,
                           partition_t ids_per_super_block,
                           weight_t lmax) {
            f64 alpha = 1000.0;
            f64 beta  = std::log(g.get_n());

            // initialize the heaps
            l_boundary_vertices.clear();
            for (partition_t id = l_start; id < l_start + ids_per_super_block; ++id) {
                forall_bv_id_iu(bv_manager, id, i, u)
                    {
                        s64 delta         = 0;
                        bool is_connected = false;
                        forall_guivw(g, u, j, v, w)
                            {
                                partition_t v_id = p_manager[v];
                                if (l_start <= v_id && v_id < l_start + ids_per_super_block) {
                                    delta -= w;
                                }
                                if (r_start <= v_id && v_id < r_start + ids_per_super_block) {
                                    delta += w;
                                    is_connected = true;
                                }
                            }
                        endfor
                        if (is_connected) {
                            l_boundary_vertices.push(u, delta);
                        }
                    }
                endfor
            }

            r_boundary_vertices.clear();
            for (partition_t id = r_start; id < r_start + ids_per_super_block; ++id) {
                forall_bv_id_iu(bv_manager, id, i, u)
                    {
                        s64 delta         = 0;
                        bool is_connected = false;
                        forall_guivw(g, u, j, v, w)
                            {
                                partition_t v_id = p_manager[v];
                                if (l_start <= v_id && v_id < l_start + ids_per_super_block) {
                                    delta += w;
                                    is_connected = true;
                                }

                                if (r_start <= v_id && v_id < r_start + ids_per_super_block) {
                                    delta -= w;
                                }
                            }
                        endfor
                        if (is_connected) {
                            r_boundary_vertices.push(u, delta);
                        }
                    }
                endfor
            }

            // get the weight of both superblocks
            weight_t l_weight = 0;
            weight_t r_weight = 0;
            for (partition_t id = l_start; id < l_start + ids_per_super_block; ++id) { l_weight += p_manager.get_bweight(id); }
            for (partition_t id = r_start; id < r_start + ids_per_super_block; ++id) { r_weight += p_manager.get_bweight(id); }

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

            while ((!l_boundary_vertices.empty() || !r_boundary_vertices.empty()) && moves_size < max_n_swaps) {
                // determine from which block to choose
                bool choose_l = true;
                // 1. if one block is empty, then choose the other one
                if (l_boundary_vertices.empty() || r_boundary_vertices.empty()) {
                    choose_l = r_boundary_vertices.empty();
                } else {
                    // 2. choose the block with greater gain and randomly if even
                    if (r_boundary_vertices.top() > l_boundary_vertices.top()) {
                        choose_l = false;
                    } else if (r_boundary_vertices.top() == l_boundary_vertices.top()) {
                        choose_l = random_engine->get_f32() < 0.5;
                    }

                    // 3. if one block is overloaded, choose the larger one, if both same sizes, then randomly
                    if (l_weight > lmax && l_weight > r_weight) { choose_l = true; }
                    if (r_weight > lmax && r_weight > l_weight) { choose_l = false; }
                    if (l_weight > lmax && r_weight > lmax && l_weight == r_weight) { choose_l = random_engine->get_f32() < 0.5; }
                }

                // choose the priority queue
                IndexedMaxHeap<s64>& boundary_vertices = choose_l ? l_boundary_vertices : r_boundary_vertices;

                vertex_t vertex        = boundary_vertices.top_key();
                s64 qap_delta          = boundary_vertices.top();
                weight_t vertex_weight = g.weight(vertex);
                partition_t vertex_id  = p_manager[vertex];
                partition_t move_id    = choose_l ? r_start : l_start;
                boundary_vertices.pop();

                // greedy way where to put the vertex
                s64 move_qap_delta = get_u_qap_delta(g, vertex, vertex_id, move_id, p_manager, d_oracle);;
                if (choose_l) {
                    for (partition_t id = r_start; id < r_start + ids_per_super_block; ++id) {
                        s64 real_qap_delta = get_u_qap_delta(g, vertex, vertex_id, id, p_manager, d_oracle);
                        if (real_qap_delta > move_qap_delta && p_manager.get_bweight(id) + vertex_weight <= m_lmax) {
                            move_id = id;
                            move_qap_delta = real_qap_delta;
                        }
                    }
                } else {
                    for (partition_t id = l_start; id < l_start + ids_per_super_block; ++id) {
                        s64 real_qap_delta = get_u_qap_delta(g, vertex, vertex_id, id, p_manager, d_oracle);
                        if (real_qap_delta > move_qap_delta && p_manager.get_bweight(id) + vertex_weight <= m_lmax) {
                            move_id = id;
                            move_qap_delta = real_qap_delta;
                        }
                    }
                }

                // move the vertex
                moves[moves_size++] = {vertex, vertex_id, move_id};
                curr_qap_gain += qap_delta;
                if (choose_l) {
                    l_weight -= vertex_weight;
                    r_weight += vertex_weight;
                } else {
                    l_weight += vertex_weight;
                    r_weight -= vertex_weight;
                }

                if (curr_qap_gain >= max_qap_gain && l_weight <= lmax && r_weight <= lmax) {
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

                        bool in_l = l_start <= neighbor_id && neighbor_id < l_start + ids_per_super_block;
                        bool in_r = r_start <= neighbor_id && neighbor_id < r_start + ids_per_super_block;

                        if (!(in_l || in_r)) { continue; }

                        if (in_l) {
                            s64 delta         = 0;
                            bool is_connected = false;
                            forall_guivw(g, neighbor, j, v, w)
                                {
                                    partition_t v_id = p_manager[v];
                                    if (l_start <= v_id && v_id < l_start + ids_per_super_block) {
                                        delta -= w;
                                    }
                                    if (r_start <= v_id && v_id < r_start + ids_per_super_block) {
                                        delta += w;
                                        is_connected = true;
                                    }
                                }
                            endfor
                            if (is_connected) {
                                l_boundary_vertices.push_update(neighbor, delta);
                            }
                        } else {
                            s64 delta         = 0;
                            bool is_connected = false;
                            forall_guivw(g, neighbor, j, v, w)
                                {
                                    partition_t v_id = p_manager[v];
                                    if (l_start <= v_id && v_id < l_start + ids_per_super_block) {
                                        delta += w;
                                        is_connected = true;
                                    }

                                    if (r_start <= v_id && v_id < r_start + ids_per_super_block) {
                                        delta -= w;
                                    }
                                }

                            endfor
                            if (is_connected) {
                                r_boundary_vertices.push_update(neighbor, delta);
                            }
                        }
                    }
                endfor

                // remove vertex from u if it is not boundary
                while (!l_boundary_vertices.empty() && !is_connected(g, p_manager, l_boundary_vertices.top_key(), r_start, ids_per_super_block)) { l_boundary_vertices.pop(); }

                // remove vertex from v if it is not boundary
                while (!r_boundary_vertices.empty() && !is_connected(g, p_manager, r_boundary_vertices.top_key(), l_start, ids_per_super_block)) { r_boundary_vertices.pop(); }
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
