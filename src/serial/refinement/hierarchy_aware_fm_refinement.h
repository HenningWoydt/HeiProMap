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

#ifndef HEIPROMAP_HIERARCHY_AWARE_FM_REFINEMENT_H
#define HEIPROMAP_HIERARCHY_AWARE_FM_REFINEMENT_H

#include <algorithm>
#include <random>

#include "../datastructures/distance_oracle.h"
#include "../datastructures/functions.h"
#include "../interfaces/ISerialBoundaryVertexManager.h"
#include "../interfaces/ISerialQuotientGraph.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/qap.h"
#include "../../commons/utils.h"
#include "k_way_fm_refinement_Faraj20.h"

namespace HeiProMap {
    class HierarchyAwareFMConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit HierarchyAwareFMConfiguration(const std::string &t_name) : ISerialRefinerConfiguration(t_name) {}

        u64 max_iteration = 1; // how many iterations to run the algorithm at most
    };

    inline partition_t get_island_id(partition_t u_id, const std::vector<std::vector<partition_t>> &island_ids) {
        for (partition_t i = 0; i < island_ids.size(); ++i) {
            if (std::find(island_ids[i].begin(), island_ids[i].end(), u_id) != island_ids[i].end()) {
                return i;
            }
        }
        std::cout << "Error u_id " << u_id << " not found!" << std::endl;
        exit(EXIT_FAILURE);
    }

    /**
     * Since the top level of the hierarchy is the most important, try to optimize it the most.
     * Aggregate all partitions of the islands and then try to find moves between the islands instead of individual partitions.
     * If moves between the islands have been found, then try to distribute it onto the individual partitions.
     */
    class HierarchyAwareFMRefinement final : public ISerialRefiner {
        vertex_t                 m_n    = 0;
        vertex_t                 m_m    = 0;
        partition_t              m_k    = 0;
        weight_t                 m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t>    m_distance;

        u32 *vertex_used  = nullptr;
        u32 vertex_marker = 0;

        u32 *block_used  = nullptr;
        u32 block_marker = 0;

        IndexedMaxHeap<KWayFMMove> heap;

        RandomEngine                        *random_engine    = nullptr;
        const HierarchyAwareFMConfiguration *config           = nullptr;
        StatisticCollector                  *m_stat_collector = nullptr;

    public:
        HierarchyAwareFMRefinement() = default;

        ~HierarchyAwareFMRefinement() override {
            free(vertex_used);
            free(block_used);
        }

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
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
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;

            random_engine    = &t_random_engine;
            config           = dynamic_cast<const HierarchyAwareFMConfiguration *>(&i_config);
            m_stat_collector = &t_stat_collect;

            vertex_t t_n_64 = round_up_64(t_n);
            vertex_t t_k_64 = round_up_64(t_k);

            vertex_used = (u32 *) aligned_alloc(64, t_n_64 * sizeof(u32));
            std::fill_n(vertex_used, t_n_64, 0);
            vertex_marker = 0;

            block_used = (u32 *) aligned_alloc(64, t_k_64 * sizeof(u32));
            std::fill_n(block_used, t_k_64, 0);
            block_marker = 0;

            heap.initialize(m_n);
        }

        void refine(const u64 level,
                    const u64 max_level,
                    const graph_t &g,
                    const d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph) override {
            auto sp = std::chrono::high_resolution_clock::now();

            size_t                                n_islands      = m_hierarchy.back();
            size_t                                ids_per_island = m_k / n_islands;
            std::vector<std::vector<partition_t>> island_ids(n_islands, std::vector<partition_t>(ids_per_island));
            std::vector<weight_t>                 islands_weight(n_islands, 0.0);
            weight_t                              island_lmax    = (weight_t) ids_per_island * m_lmax;

            size_t      id = 0;
            for (size_t i  = 0; i < n_islands; i++) {
                for (size_t j = 0; j < ids_per_island; j++) {
                    island_ids[i][j] = id;
                    islands_weight[i] += p_manager.get_bweight(id);
                    id += 1;
                }
            }

            s64 qap_gain = 0;
            for (u64 iteration = 0; iteration < config->max_iteration; ++iteration) {
                // for each island, move vertices to it
                for (partition_t island_id = 0; island_id < n_islands; ++island_id) {
                    auto sp_island = std::chrono::high_resolution_clock::now();
                    heap.clear();

                    // gather boundary vertices of all other islands, and check if moving them to any neighboring block of this island is beneficial
                    vertex_marker += 1;
                    for (auto v_id: island_ids[island_id]) {
                        forall_bv_id_iu(bv_manager, v_id, i, v)
                            {
                                if (vertex_used[v] == vertex_marker) { continue; }
                                forall_guiv(g, v, j, vv)
                                    {
                                        if (vertex_used[vv] == vertex_marker) { continue; }

                                        partition_t vv_id = p_manager[vv];
                                        if (exists(island_ids[island_id], vv_id)) { continue; }

                                        partition_t best_vvv_id;
                                        s64         best_qap_delta = -std::numeric_limits<s64>::max();
                                        block_marker += 1;
                                        forall_guiv(g, vv, k, vvv)
                                            {
                                                partition_t vvv_id = p_manager[vvv];
                                                if (block_used[vvv_id] == block_marker) { continue; }
                                                if (!exists(island_ids[island_id], vvv_id)) { continue; }
                                                s64 qap_delta = get_u_qap_delta(g, vv, vv_id, vvv_id, p_manager, d_oracle);
                                                if (qap_delta > best_qap_delta) {
                                                    best_qap_delta = qap_delta;
                                                    best_vvv_id    = vvv_id;
                                                }
                                                block_used[vvv_id] = block_marker;
                                            }
                                        endfor
                                        if (best_qap_delta != -std::numeric_limits<s64>::max()) {
                                            heap.push_update(vv, {vv, vv_id, best_vvv_id, best_qap_delta});
                                        }
                                        vertex_used[vv] = vertex_marker;
                                    }
                                endfor
                                vertex_used[v] = vertex_marker;
                            }
                        endfor
                    }

                    // start moving the vertices on the blocks on the island
                    std::vector<KWayFMMove> moves;
                    size_t                  best_idx      = 0;
                    s64                     max_qap_gain  = 0;
                    s64                     curr_qap_gain = 0;

                    vertex_marker += 1;
                    while (!heap.empty()) {
                        KWayFMMove move = heap.top();
                        heap.pop();

                        vertex_t    vertex        = move.u;
                        partition_t vertex_id     = p_manager[vertex];
                        weight_t    vertex_weight = g.weight(vertex);
                        partition_t move_id       = move.to_move_id;

                        if (islands_weight[island_id] + vertex_weight > island_lmax) { continue; }

                        moves.push_back(move);
                        curr_qap_gain += move.qap_delta;
                        if (curr_qap_gain > max_qap_gain) {
                            best_idx     = moves.size();
                            max_qap_gain = curr_qap_gain;
                        }

                        // make move in structures
                        p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                        vertex_used[vertex] = vertex_marker;
                        islands_weight[island_id] += vertex_weight;
                        islands_weight[get_island_id(vertex_id, island_ids)] -= vertex_weight;

                        // we have to push or update the neighbors that were not moved already
                        forall_guiv(g, vertex, i, neighbor)
                            {
                                partition_t neighbor_id     = p_manager[neighbor];
                                weight_t    neighbor_weight = g.weight(neighbor);
                                if (vertex_used[neighbor] == vertex_marker) { continue; }
                                if (islands_weight[island_id] + neighbor_weight > island_lmax) { continue; }
                                if (!is_boundary(g, p_manager, neighbor)) { continue; }

                                block_marker += 1;
                                s64         best_qap_delta = -std::numeric_limits<s64>::max();
                                partition_t best_v_id;
                                forall_guiv(g, neighbor, j, v)
                                    {
                                        partition_t v_id     = p_manager[v];
                                        if (block_used[v_id] == block_marker) { continue; }
                                        if (!exists(island_ids[island_id], v_id)) { continue; }

                                        s64 qap_delta = get_u_qap_delta(g, neighbor, neighbor_id, v_id, p_manager, d_oracle);
                                        if (qap_delta > best_qap_delta) {
                                            best_qap_delta = qap_delta;
                                            best_v_id      = v_id;
                                        }
                                        block_used[v_id] = block_marker;
                                    }
                                endfor
                                if (best_qap_delta != -std::numeric_limits<s64>::max()) {
                                    heap.push_update(neighbor, KWayFMMove(neighbor, neighbor_id, best_v_id, best_qap_delta));
                                }
                            }
                        endfor
                    }

                    // revert all moves in partitioning manager
                    for (size_t i = 0; i < moves.size(); i++) {
                        vertex_t    vertex        = moves[moves.size() - 1 - i].u;
                        weight_t    vertex_weight = g.weight(vertex);
                        partition_t vertex_id     = moves[moves.size() - 1 - i].to_move_id;
                        partition_t move_id       = moves[moves.size() - 1 - i].u_id;

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

                    auto ep_island = std::chrono::high_resolution_clock::now();
                    // std::cout << "checked island " << island_id << " in " << get_seconds(sp_island, ep_island) << std::endl;
                    // std::cout << "theoretical initial max qap gain " << initial_max_qap_gain << std::endl;
                    // std::cout << "max qap gain " << max_qap_gain << " at move " << best_idx << " of " << moves.size() << std::endl;
                    // std::cout << "is_overloaded " << p_manager.is_overloaded() << " max weight " << max(p_manager.get_bweights()) << " of " << m_lmax << std::endl;
                    // std::cout << std::endl;

                    qap_gain += max_qap_gain;
                }

                // now we have to rebalance the actual blocks and spend at most qap_gain to do so



            }
            std::cout << "-- Level " << level << " made " << qap_gain*2 << std::endl;
        }

        void refine_exhaustive_3(const u64 level,
                                 const u64 max_level,
                                 const graph_t &g,
                                 const d_oracle_t &d_oracle,
                                 bv_manager_t &bv_manager,
                                 p_manager_t &p_manager,
                                 q_graph_t &q_graph) {
            auto sp = std::chrono::high_resolution_clock::now();

            size_t                                n_islands      = m_hierarchy.back();
            size_t                                ids_per_island = m_k / n_islands;
            std::vector<std::vector<partition_t>> island_ids(n_islands, std::vector<partition_t>(ids_per_island));
            std::vector<weight_t>                 islands_weight(n_islands, 0.0);
            weight_t                              island_lmax    = (weight_t) ids_per_island * m_lmax;

            size_t      id = 0;
            for (size_t i  = 0; i < n_islands; i++) {
                for (size_t j = 0; j < ids_per_island; j++) {
                    island_ids[i][j] = id;
                    islands_weight[i] += p_manager.get_bweight(id);
                    id += 1;
                }
            }

            for (u64 iteration = 0; iteration < config->max_iteration; ++iteration) {
                // for each island, move vertices to it
                for (partition_t island_id = 0; island_id < n_islands; ++island_id) {
                    // gather boundary vertices of all other islands, and check if moving them to any block of this island is beneficial
                    std::vector<vertex_t> curr_boundary;

                    for (auto b_id: island_ids[island_id]) {
                        forall_bv_id_iu(bv_manager, b_id, i, u)
                            {
                                forall_guiv(g, u, j, v)
                                    {
                                        partition_t v_id = p_manager[v];
                                        if (exists(island_ids[island_id], v_id)) { continue; }
                                        curr_boundary.push_back(v);
                                    }
                                endfor
                            }
                        endfor
                    }

                    // make the boundary unique
                    std::sort(curr_boundary.begin(), curr_boundary.end());
                    curr_boundary.erase(std::unique(curr_boundary.begin(), curr_boundary.end()), curr_boundary.end());

                    std::vector<std::vector<partition_t>> move_ids(curr_boundary.size());

                    for (size_t i = 0; i < curr_boundary.size(); ++i) {
                        vertex_t u = curr_boundary[i];
                        forall_guiv(g, u, j, v)
                            {
                                partition_t v_id = p_manager[v];
                                if (exists(island_ids[island_id], v_id)) { move_ids[i].push_back(v_id); }
                            }
                        endfor
                        std::sort(move_ids[i].begin(), move_ids[i].end());
                        move_ids[i].erase(std::unique(move_ids[i].begin(), move_ids[i].end()), move_ids[i].end());
                    }

                    // check the best result if 3 moves at once are possible, all three vertices should be connected
                    u64 count_overload = 0;
                    u64 neg_qap_delta  = 0;
                    u64 pos_qap_delta  = 0;

                    s64 total_qap_delta = 0;
                    vertex_marker += 1;
                    for (size_t i = 0; i < curr_boundary.size(); ++i) {
                        s64         best_qap_delta = -std::numeric_limits<s64>::max();
                        vertex_t    best_v, best_vv, best_vvv;
                        weight_t    best_v_weight, best_vv_weight, best_vvv_weight;
                        partition_t best_v_id, best_vv_id, best_vvv_id;
                        partition_t best_new_v_id, best_new_vv_id, best_new_vvv_id;

                        vertex_t    v        = curr_boundary[i];
                        weight_t    v_weight = g.weight(v);
                        partition_t v_id     = p_manager[v];

                        if (vertex_used[v] == vertex_marker) { continue; }
                        if (!bv_manager.is_boundary(v)) { continue; }
                        if (islands_weight[island_id] + v_weight > island_lmax) {
                            count_overload += 1;
                            continue;
                        }

                        for (size_t j = i + 1; j < curr_boundary.size(); ++j) {
                            vertex_t    vv        = curr_boundary[j];
                            weight_t    vv_weight = g.weight(vv);
                            partition_t vv_id     = p_manager[vv];

                            if (vertex_used[vv] == vertex_marker) { continue; }
                            if (islands_weight[island_id] + v_weight + vv_weight > island_lmax) {
                                count_overload += 1;
                                continue;
                            }

                            for (size_t k = j + 1; k < curr_boundary.size(); ++k) {
                                vertex_t    vvv        = curr_boundary[k];
                                weight_t    vvv_weight = g.weight(vvv);
                                partition_t vvv_id     = p_manager[vvv];

                                if (vertex_used[vvv] == vertex_marker) { continue; }
                                if (islands_weight[island_id] + v_weight + vv_weight + vvv_weight > island_lmax) {
                                    count_overload += 1;
                                    continue;
                                }

                                for (partition_t new_v_id: move_ids[i]) {
                                    for (partition_t new_vv_id: move_ids[j]) {
                                        for (partition_t new_vvv_id: move_ids[k]) {

                                            s64 qap_delta = get_qap_delta(g, v, vv, vvv, v_id, vv_id, vvv_id, new_v_id, new_vv_id, new_vvv_id, p_manager, d_oracle);
                                            if (qap_delta <= 0) { neg_qap_delta += 1; }
                                            if (qap_delta > 0) { pos_qap_delta += 1; }

                                            if (qap_delta > best_qap_delta) {
                                                best_qap_delta = qap_delta;

                                                best_v   = v;
                                                best_vv  = vv;
                                                best_vvv = vvv;

                                                best_v_weight   = v_weight;
                                                best_vv_weight  = vv_weight;
                                                best_vvv_weight = vvv_weight;

                                                best_v_id   = v_id;
                                                best_vv_id  = vv_id;
                                                best_vvv_id = vvv_id;

                                                best_new_v_id   = new_v_id;
                                                best_new_vv_id  = new_vv_id;
                                                best_new_vvv_id = new_vvv_id;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        vertex_used[v] = vertex_marker;

                        if (best_qap_delta > 0) {

                            if (bv_manager.is_boundary(best_v)) {
                                bv_manager.move(g, p_manager, best_v, best_v_id, best_new_v_id);
                            }
                            q_graph.move(g, p_manager, best_v, best_v_id, best_new_v_id);
                            p_manager.move(best_v, best_v_weight, best_v_id, best_new_v_id);

                            if (bv_manager.is_boundary(best_vv)) {
                                bv_manager.move(g, p_manager, best_vv, best_vv_id, best_new_vv_id);
                            }
                            q_graph.move(g, p_manager, best_vv, best_vv_id, best_new_vv_id);
                            p_manager.move(best_vv, best_vv_weight, best_vv_id, best_new_vv_id);

                            if (bv_manager.is_boundary(best_vvv)) {
                                bv_manager.move(g, p_manager, best_vvv, best_vvv_id, best_new_vvv_id);
                            }
                            q_graph.move(g, p_manager, best_vvv, best_vvv_id, best_new_vvv_id);
                            p_manager.move(best_vvv, best_vvv_weight, best_vvv_id, best_new_vvv_id);

                            vertex_used[best_v]   = vertex_marker;
                            vertex_used[best_vv]  = vertex_marker;
                            vertex_used[best_vvv] = vertex_marker;

                            islands_weight[island_id] += best_v_weight + best_vv_weight + best_vvv_weight;
                            islands_weight[get_island_id(best_v_id, island_ids)] -= best_v_weight;
                            islands_weight[get_island_id(best_vv_id, island_ids)] -= best_vv_weight;
                            islands_weight[get_island_id(best_vvv_id, island_ids)] -= best_vvv_weight;

                            total_qap_delta += best_qap_delta;
                        }
                    }

                    std::cout << "total qap delta = " << total_qap_delta << " overload = " << count_overload << " neg_qap " << neg_qap_delta << " pos_qap " << pos_qap_delta << std::endl;
                }
            }
            std::cout << std::endl;
        }

        JSONString get_stats()
        override { return {}; };

    };
}

#endif //HEIPROMAP_HIERARCHY_AWARE_FM_REFINEMENT_H
