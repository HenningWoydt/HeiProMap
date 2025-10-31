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

#ifndef HEIPROMAP_HIERARCHY_AWARE_MULTI_WAY_FM_REFINEMENT_H
#define HEIPROMAP_HIERARCHY_AWARE_MULTI_WAY_FM_REFINEMENT_H

#include <algorithm>

#include "../utility/utils.h"
#include "../utility/functions.h"
#include "ISerialRefiner.h"
#include "../utility/qap.h"

namespace HeiProMap {
    class HierarchyAwareMultiWayFMConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit HierarchyAwareMultiWayFMConfiguration(const std::string &t_name) : ISerialRefinerConfiguration(t_name) {}

        u64 max_iteration = 1; // how many iterations to run the algorithm at most
        f64 alpha         = 10.0;
        f64 beta          = 1.0;
    };

    inline partition_t get_island_id(const partition_t u_id, const partition_t ids_per_island) {
        return u_id / ids_per_island;
    }

    /**
     * Since the top level of the hierarchy is the most important, try to optimize it the most.
     * Aggregate all partitions of the islands and then try to find moves between the islands instead of individual partitions.
     * If moves between the islands have been found, then try to distribute it onto the individual partitions.
     */
    class HierarchyAwareMultiWayFMRefinement final : public ISerialRefiner {
        vertex_t                 m_n         = 0;
        vertex_t                 m_m         = 0;
        partition_t              m_k         = 0;
        f64                      m_imbalance = 0.0;
        weight_t                 m_lmax      = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t>    m_distance;

        AlignedArray<u32> vertex_used;
        u32               vertex_marker = 0;

        AlignedArray<u32> block_used;
        u32               block_marker = 0;

        // IndexedMaxHeap<KWayFMMove> heap;
        // std::priority_queue<KWayFMMove> heap;

        AlignedArray<Move> moves;
        size_t             moves_size = 0;

        RandomEngine                                *random_engine = nullptr;
        const HierarchyAwareMultiWayFMConfiguration *config        = nullptr;

    public:
        HierarchyAwareMultiWayFMRefinement() = default;

        ~HierarchyAwareMultiWayFMRefinement() override = default;

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

            random_engine = &t_random_engine;
            config        = dynamic_cast<const HierarchyAwareMultiWayFMConfiguration *>(&i_config);

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
            for (size_t iteration = 0; iteration < config->max_iteration; ++iteration) {
                for (size_t i = 0; i < m_hierarchy.size() - 1; ++i) {
                    refine_layer(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, m_hierarchy.size() - 1 - i);
                    // rebalance_layer(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, m_hierarchy.size() - 2 - i);
                }
            }
        }

#define IN_SAME_BLOCK(id1, id2, ids_per_block) (id1 / ids_per_block) == (id2 / ids_per_block)
#define IN_NEIGHBORING_BLOCK(v_id, neighborhood_id_start, neighborhood_id_end) ((neighborhood_id_start <= v_id) && ( v_id < neighborhood_id_end))

        void refine_layer(const u64 level,
                          const u64 max_level,
                          const graph_t &g,
                          d_oracle_t &d_oracle,
                          bv_manager_t &bv_manager,
                          p_manager_t &p_manager,
                          q_graph_t &q_graph,
                          size_t layer) {
            partition_t n_total_super_blocks = 1;
            for (size_t i                    = layer; i < m_hierarchy.size(); ++i) { n_total_super_blocks *= m_hierarchy[i]; }
            partition_t ids_per_super_block  = m_k / n_total_super_blocks;

            partition_t n_local_super_blocks = m_hierarchy[layer];
            partition_t n_upper_blocks       = 1;
            for (size_t i                    = layer + 1; i < m_hierarchy.size(); ++i) { n_upper_blocks *= m_hierarchy[i]; }

            partition_t neighborhood_id_start = 0;
            partition_t neighborhood_id_end   = ids_per_super_block * n_local_super_blocks;

            for (size_t upper_block_id = 0; upper_block_id < n_upper_blocks; ++upper_block_id) {
                refine_neighborhood(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, neighborhood_id_start, neighborhood_id_end, n_local_super_blocks, ids_per_super_block);

                neighborhood_id_start = neighborhood_id_end;
                neighborhood_id_end += ids_per_super_block * n_local_super_blocks;
            }
        }

        void refine_neighborhood([[maybe_unused]] const u64 level,
                                 [[maybe_unused]] const u64 max_level,
                                 const graph_t &g,
                                 d_oracle_t &d_oracle,
                                 bv_manager_t &bv_manager,
                                 p_manager_t &p_manager,
                                 q_graph_t &q_graph,
                                 partition_t neighborhood_id_start,
                                 partition_t neighborhood_id_end,
                                 partition_t n_local_super_blocks,
                                 partition_t ids_per_super_block) {
            f64 alpha = config->alpha;
            f64 beta  = std::log(g.get_n());

            weight_t              blocks_lmax = (weight_t) ids_per_super_block * m_lmax;
            std::vector<weight_t> blocks_weights(n_local_super_blocks, 0);

            for (partition_t i = 0; i < n_local_super_blocks; ++i) {
                for (partition_t j = 0; j < ids_per_super_block; ++j) {
                    partition_t id = neighborhood_id_start + i * ids_per_super_block + j;
                    blocks_weights[i] += p_manager.get_bweight(id);
                }
            }

            std::priority_queue<KWayFMMove> heap;
            // collect all boundary vertices that can be moved between the superblocks
            for (partition_t                u_id = neighborhood_id_start; u_id < neighborhood_id_end; ++u_id) {
                forall_bv_id_iu(bv_manager, u_id, i, u)
                    {
                        weight_t u_weight = g.weight(u);

                        block_marker += 1;
                        forall_guiv(g, u, j, v)
                            {
                                partition_t v_id = p_manager[v];

                                if (IN_SAME_BLOCK(u_id, v_id, ids_per_super_block)) { continue; }
                                if (!IN_NEIGHBORING_BLOCK(v_id, neighborhood_id_start, neighborhood_id_end)) { continue; }
                                if (block_used[v_id] == block_marker) { continue; }
                                if (blocks_weights[(v_id - neighborhood_id_start) / ids_per_super_block] + u_weight > blocks_lmax) { continue; }

                                s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                                heap.emplace(u, u_id, v_id, qap_delta);

                                block_used[v_id] = block_marker;
                            }
                        endfor
                    }
                endfor
            }

            // start moving the vertices on the blocks on the island
            moves_size = 0;
            size_t best_idx      = 0;
            s64    max_qap_gain  = 0;
            s64    curr_qap_gain = 0;

            f64 steps_since_last_improvement = 0.0;
            f64 qap_gain_mean                = 0.0;
            f64 qap_gain_var                 = 0.0;

            vertex_marker += 1;
            while (!heap.empty()) {
                KWayFMMove move = heap.top();
                heap.pop();

                vertex_t    vertex        = move.u;
                partition_t vertex_id     = p_manager[vertex];
                weight_t    vertex_weight = g.weight(vertex);
                partition_t move_id       = move.to_move_id;

                if (vertex_used[vertex] == vertex_marker) { continue; }
                if (blocks_weights[(move_id - neighborhood_id_start) / ids_per_super_block] + vertex_weight > blocks_lmax) { continue; }

                s64 temp_qap_delta = get_u_qap_delta(g, vertex, vertex_id, move_id, p_manager, d_oracle);
                if (temp_qap_delta != move.qap_delta) { continue; }

                moves[moves_size++] = Move(vertex, vertex_id, move_id);
                curr_qap_gain += move.qap_delta;
                if (curr_qap_gain > max_qap_gain) {
                    best_idx     = moves_size;
                    max_qap_gain = curr_qap_gain;

                    steps_since_last_improvement = 0.0;
                    qap_gain_mean                = 0.0;
                    qap_gain_var                 = 0.0;
                }

                // make move in structures
                p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                blocks_weights[(move_id - neighborhood_id_start) / ids_per_super_block] += vertex_weight;
                blocks_weights[(vertex_id - neighborhood_id_start) / ids_per_super_block] -= vertex_weight;
                vertex_used[vertex] = vertex_marker;

                steps_since_last_improvement += 1.0;
                f64 new_qap_gain_mean = qap_gain_mean + ((f64) move.qap_delta - qap_gain_mean) / steps_since_last_improvement;
                f64 new_qap_gain_var  = (qap_gain_var + ((f64) move.qap_delta - qap_gain_mean) * ((f64) move.qap_delta - new_qap_gain_mean)) / steps_since_last_improvement;

                qap_gain_mean = new_qap_gain_mean;
                qap_gain_var  = new_qap_gain_var;

                if (steps_since_last_improvement > 2.0 && steps_since_last_improvement * qap_gain_mean * qap_gain_mean > alpha * qap_gain_var + beta) { break; }

                // we have to push or update the neighbors that were not moved already
                forall_guiv(g, vertex, i, neighbor)
                    {
                        partition_t neighbor_id     = p_manager[neighbor];
                        weight_t    neighbor_weight = g.weight(neighbor);

                        if (!IN_NEIGHBORING_BLOCK(neighbor_id, neighborhood_id_start, neighborhood_id_end)) { continue; }
                        if (vertex_used[neighbor] == vertex_marker) { continue; }
                        if (!is_boundary(g, p_manager, neighbor)) { continue; }

                        block_marker += 1;
                        forall_guiv(g, neighbor, j, v)
                            {
                                partition_t v_id = p_manager[v];

                                if (IN_SAME_BLOCK(neighbor_id, v_id, ids_per_super_block)) { continue; }
                                if (!IN_NEIGHBORING_BLOCK(v_id, neighborhood_id_start, neighborhood_id_end)) { continue; }
                                if (block_used[v_id] == block_marker) { continue; }
                                if (blocks_weights[(v_id - neighborhood_id_start) / ids_per_super_block] + neighbor_weight > blocks_lmax) { continue; }

                                s64 qap_delta = get_u_qap_delta(g, neighbor, neighbor_id, v_id, p_manager, d_oracle);
                                heap.emplace(neighbor, neighbor_id, v_id, qap_delta);
                                block_used[v_id] = block_marker;
                            }
                        endfor
                    }
                endfor
            }

            // revert all moves in partitioning manager
            for (size_t i = 0; i < moves_size; i++) {
                vertex_t    vertex        = moves[moves_size - 1 - i].u;
                weight_t    vertex_weight = g.weight(vertex);
                partition_t vertex_id     = moves[moves_size - 1 - i].to_move_id;
                partition_t move_id       = moves[moves_size - 1 - i].u_id;

                p_manager.move(vertex, vertex_weight, vertex_id, move_id);
            }

            // make all moves to best index
            for (size_t i = 0; i < best_idx; ++i) {
                vertex_t    vertex        = moves[i].u;
                weight_t    vertex_weight = g.weight(vertex);
                partition_t vertex_id     = moves[i].u_id;
                partition_t move_id       = moves[i].to_move_id;

                bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                p_manager.move(vertex, vertex_weight, vertex_id, move_id);
            }
        }

        void rebalance_layer(const u64 level,
                             const u64 max_level,
                             const graph_t &g,
                             d_oracle_t &d_oracle,
                             bv_manager_t &bv_manager,
                             p_manager_t &p_manager,
                             q_graph_t &q_graph,
                             size_t layer) {
            partition_t n_total_super_blocks = 1;
            for (size_t i                    = layer; i < m_hierarchy.size(); ++i) { n_total_super_blocks *= m_hierarchy[i]; }
            partition_t ids_per_super_block  = m_k / n_total_super_blocks;

            partition_t n_local_super_blocks = m_hierarchy[layer];
            partition_t n_upper_blocks       = 1;
            for (size_t i                    = layer + 1; i < m_hierarchy.size(); ++i) { n_upper_blocks *= m_hierarchy[i]; }

            partition_t neighborhood_id_start = 0;
            partition_t neighborhood_id_end   = ids_per_super_block * n_local_super_blocks;

            for (size_t upper_block_id = 0; upper_block_id < n_upper_blocks; ++upper_block_id) {
                rebalance_neighborhoods(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, neighborhood_id_start, neighborhood_id_end, n_local_super_blocks, ids_per_super_block);

                neighborhood_id_start = neighborhood_id_end;
                neighborhood_id_end += ids_per_super_block * n_local_super_blocks;
            }
        }

        void rebalance_neighborhoods([[maybe_unused]] const u64 level,
                                     [[maybe_unused]] const u64 max_level,
                                     const graph_t &g,
                                     d_oracle_t &d_oracle,
                                     bv_manager_t &bv_manager,
                                     p_manager_t &p_manager,
                                     q_graph_t &q_graph,
                                     partition_t neighborhood_id_start,
                                     partition_t neighborhood_id_end,
                                     partition_t n_local_super_blocks,
                                     partition_t ids_per_super_block) {
            weight_t blocks_lmax = (weight_t) ids_per_super_block * m_lmax;;
            std::vector<weight_t> blocks_weights(n_local_super_blocks, 0);

            for (partition_t i = 0; i < n_local_super_blocks; ++i) {
                for (partition_t j = 0; j < ids_per_super_block; ++j) {
                    partition_t id = neighborhood_id_start + i * ids_per_super_block + j;
                    blocks_weights[i] += p_manager.get_bweight(id);
                }
            }

            std::vector<weight_t> min_blocks_weights(n_local_super_blocks, 0);
            for (partition_t      i = 0; i < n_local_super_blocks; ++i) {
                min_blocks_weights[i] = std::max(blocks_weights[i], blocks_lmax);
            }

            std::vector<size_t> indices(n_local_super_blocks);
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(), [&](size_t i, size_t j) { return blocks_weights[i] > blocks_weights[j]; });

            if (blocks_weights[indices[0]] <= blocks_lmax) { return; }

            while (max(blocks_weights) > blocks_lmax) {
                std::sort(indices.begin(), indices.end(), [&](size_t i, size_t j) { return blocks_weights[i] > blocks_weights[j]; });

                bool        move_occurred = false;
                for (size_t i             = 0; i < indices.size(); ++i) {
                    partition_t id = indices[i];
                    if (blocks_weights[id] <= blocks_lmax) { continue; }

                    vertex_t    vertex_move;
                    partition_t vertex_curr_id;
                    partition_t vertex_new_id;
                    weight_t    vertex_weight;
                    s64         vertex_qap_delta = -std::numeric_limits<s64>::max();

                    partition_t      start = neighborhood_id_start + id * ids_per_super_block;
                    partition_t      end   = neighborhood_id_start + (id + 1) * ids_per_super_block;
                    for (partition_t u_id  = start; u_id < end; ++u_id) {
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

        void refine_layer([[maybe_unused]] const u64 level,
                          [[maybe_unused]] const u64 max_level,
                          [[maybe_unused]] graph_t &g,
                          [[maybe_unused]] d_oracle_t &d_oracle,
                          [[maybe_unused]] bv_manager_t &bv_manager,
                          [[maybe_unused]] p_manager_t &p_manager,
                          [[maybe_unused]] q_graph_t &q_graph,
                          [[maybe_unused]] size_t layer) override {}
    };
}

#endif //HEIPROMAP_HIERARCHY_AWARE_MULTI_WAY_FM_REFINEMENT_H
