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

#ifndef HEIPROMAP_GREEDY_KWAY_PARTITIONER_H
#define HEIPROMAP_GREEDY_KWAY_PARTITIONER_H

#include <string>

#include "../serial_definitions_1.h"
#include "../serial_definitions_2.h"
#include "../serial_definitions_3.h"
#include "../../commons/random_engine.h"
#include "../../commons/statistic_collector.h"

namespace HeiProMap {
    class GreedyKWayPartitionerConfiguration {};

    class GreedyKWayPartitioner {
    private:
        AlignedArray<weight_t> weights;
        AlignedArray<weight_t> edge_cut_saved;
        AlignedArray<partition_t> new_partition;

        partition_t m_k = 0;
        std::vector<std::vector<vertex_t>> blocks;

    public:
        void initialize(vertex_t max_n, partition_t k) {
            new_partition.initialize(max_n);
            m_k = k;
        }

        /**
         * Partitions the subgraph of g where the vertices have the corresponding id.
         *
         * @param g The graph to partition.
         * @param p_manager The partition manager.
         * @param id The id of the subgraph to partition.
         * @param k The number of partitions.
         * @param lmax The maximum weight of a partition.
         * @param t_random_engine The random engine.
         * @param i_config The configuration.
         * @param t_stat_collect The statistic collector.
         */
        void partition(const graph_t& g,
                       deep_p_manager_t& p_manager,
                       deep_bv_manager_t& bv_manager,
                       deep_q_graph_t& q_graph,
                       partition_t id,
                       partition_t id_increment,
                       partition_t k,
                       weight_t lmax,
                       s32 hierarchy_level,
                       RandomEngine& t_random_engine,
                       const GreedyKWayPartitionerConfiguration& i_config) {
            weights.initialize(k, 0);
            new_partition.initialize(g.get_n(), k);

            forall_gu(g, u)
                {
                    partition_t u_id = p_manager[u];
                    if (u_id != id) { continue; }

                    edge_cut_saved.initialize(k, 0);
                    forall_guivw(g, u, i, v, w)
                        {
                            partition_t v_id = p_manager[v];
                            if (v_id != id) { continue; } // wrong partition
                            if (new_partition[v] == k) { continue; } // v was not placed yet

                            edge_cut_saved[new_partition[v]] += w;
                        }
                    endfor

                    weight_t u_weight   = g.weight(u);
                    partition_t best_id = 0;
                    for (partition_t new_id = 0; new_id < k; ++new_id) {
                        if (weights[new_id] < weights[best_id]) { best_id = new_id; }
                    }

                    for (partition_t new_id = 0; new_id < k; ++new_id) {
                        if (weights[new_id] + u_weight > lmax) { continue; }

                        if (edge_cut_saved[new_id] > edge_cut_saved[best_id] || (edge_cut_saved[new_id] == edge_cut_saved[best_id] && weights[new_id] < weights[best_id])) {
                            best_id = new_id;
                        }
                    }

                    weights[best_id] += u_weight;
                    new_partition[u] = best_id;
                }
            endfor

            for (vertex_t u = 0; u < g.get_n(); ++u) {
                if (p_manager[u] != id) { continue; }
                if (new_partition[u] == 0) { continue; }
                partition_t move_id = id + id_increment * new_partition[u];

                partition_t u_id  = id;
                weight_t u_weight = g.weight(u);

                bv_manager.move(g, p_manager, u, u_id, move_id);
                q_graph.move(g, p_manager, u, u_id, move_id);
                p_manager.move(u, u_weight, u_id, move_id);
            }

            for (partition_t i = 0; i < k; ++i) {
                partition_t move_id = id + id_increment * i;
                p_manager.set_lmax(move_id, lmax);
                p_manager.set_hierarchy_level(move_id, hierarchy_level - 1);
            }
        }

        void partition_full_balance(const graph_t& g,
                                    deep_p_manager_t& p_manager,
                                    deep_bv_manager_t& bv_manager,
                                    deep_q_graph_t& q_graph,
                                    partition_t id,
                                    partition_t id_increment,
                                    partition_t k,
                                    weight_t lmax,
                                    s32 hierarchy_level,
                                    RandomEngine& t_random_engine,
                                    const GreedyKWayPartitionerConfiguration& i_config) {
            weights.initialize(k, 0);

            for (size_t i = 0; i < blocks[id].size(); ++i) {
                vertex_t u        = blocks[id][i];
                weight_t u_weight = g.weight(u);

                partition_t best_id = 0;
                for (partition_t new_id = 0; new_id < k; ++new_id) {
                    if (weights[new_id] < weights[best_id]) { best_id = new_id; }
                }

                weights[best_id] += u_weight;
                new_partition[u] = best_id;
            }

            for (size_t i = 0; i < blocks[id].size(); ++i) {
                vertex_t u = blocks[id][i];
                if (new_partition[u] == 0) { continue; }
                partition_t move_id = id + id_increment * new_partition[u];

                partition_t u_id  = id;
                weight_t u_weight = g.weight(u);

                bv_manager.move(g, p_manager, u, u_id, move_id);
                q_graph.move(g, p_manager, u, u_id, move_id);
                p_manager.move(u, u_weight, u_id, move_id);
            }

            for (partition_t i = 0; i < k; ++i) {
                partition_t move_id = id + id_increment * i;
                p_manager.set_lmax(move_id, lmax);
                p_manager.set_hierarchy_level(move_id, hierarchy_level - 1);
            }
        }

        void determine_all_blocks(const graph_t& g,
                                  deep_p_manager_t& p_manager) {
            for (partition_t id = 0; id < std::min((partition_t)blocks.size(), m_k); ++id) {
                blocks[id].clear();
            }
            if (blocks.size() < m_k) {
                blocks.resize(m_k);
            }

            forall_gu(g, u) {
                partition_t u_id  = p_manager[u];
                blocks[u_id].push_back(u);
            }
            endfor
        }
    };
}

#endif //HEIPROMAP_GREEDY_KWAY_PARTITIONER_H
