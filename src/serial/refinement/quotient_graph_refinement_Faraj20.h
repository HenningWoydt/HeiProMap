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

#ifndef HEIPROMAP_QUOTIENT_GRAPH_REFINEMENT_FARAJ20_H
#define HEIPROMAP_QUOTIENT_GRAPH_REFINEMENT_FARAJ20_H

#include <algorithm>
#include <random>

#include "../datastructures/distance_oracle.h"
#include "../datastructures/functions.h"
#include "../datastructures/indexed_max_heap.h"
#include "../interfaces/ISerialQuotientGraph.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/assert_state.h"
#include "../utility/qap.h"
#include "../utility/utils.h"

namespace HeiProMap {
    struct QuotientGraphRefinementFaraj20Configuration {
        u64 max_iteration         = 1; // how many iterations to run the algorithm at most
        u64 max_moves_without_max = 5000; // stops if after 5000 moves no new maximum has been discovered
    };

    struct Pair {
        partition_t id1;
        partition_t id2;
    };

    /**
     * Executes quotient graph refinement as described in
     * > Marcelo Fonseca Faraj, Alexander van der Grinten, Henning Meyerhenke, Jesper Larsson Träff, and Christian Schulz.
     * > High-quality Hierarchical Process Mapping.
     * > In 18th International Symposium on Experimental Algorithms, SEA 2020, June 16-18, 2020, Catania, Italy, volume 160 of LIPIcs, pages 4:1–4:15.
     * > Schloss Dagstuhl - Leibniz-Zentrum für Informatik, 2020.
     *
     * This includes the TopGain and active scheduling scheme from
     * > Sanders, P., Schulz, C. (2011).
     * > Engineering Multilevel Graph Partitioning Algorithms.
     * > In: Demetrescu, C., Halldórsson, M.M. (eds) Algorithms – ESA 2011. ESA 2011. Lecture Notes in Computer Science, vol 6942. Springer, Berlin, Heidelberg
     *
     */
    class QuotientGraphRefinementFaraj20 final : public ISerialRefiner {
    private:
        vertex_t m_n    = 0;
        vertex_t m_m    = 0;
        partition_t m_k = 0;
        weight_t m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;
        u64 m_seed = 0;

        // priority queues
        IndexedMaxHeap<s64> boundary_vertices_u;
        IndexedMaxHeap<s64> boundary_vertices_v;

        // store change
        vertex_t* moves   = nullptr;
        size_t moves_size = 0;
        s64 curr_qap_gain = 0;
        s64 max_qap_gain  = 0;
        size_t best_idx   = 0;

        // active block scheduling
        u8* active_this_round = nullptr;
        u8* active_next_round = nullptr;
        Pair* pairs           = nullptr;
        size_t pairs_size     = 0;

        // store which vertices have been moved
        u32* vertex_used = nullptr;
        u32 vertex_mark  = 0;

        std::mt19937 gen;
        std::uniform_real_distribution<float> dis;

    public:
        QuotientGraphRefinementFaraj20() = default;

        ~QuotientGraphRefinementFaraj20() override {
            free(vertex_used);
            free(moves);
            free(active_this_round);
            free(active_next_round);
            free(pairs);
        }

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const weight_t t_lmax,
                        const std::vector<partition_t>& t_hierarchy,
                        const std::vector<weight_t>& t_distance,
                        const u64 t_seed) override {
            vertex_t t_n_64        = round_up_64(t_n);
            partition_t t_k_64     = round_up_64(t_k);
            partition_t t_k_t_k_64 = round_up_64(t_k * t_k);

            m_n         = t_n;
            m_m         = t_m;
            m_k         = t_k;
            m_lmax      = t_lmax;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;
            m_seed      = t_seed;

            // priority queues
            boundary_vertices_u.initialize(m_n);
            boundary_vertices_v.initialize(m_n);

            vertex_mark = 0;
            vertex_used = (u32*)aligned_alloc(64, t_n_64 * sizeof(u32));
            std::fill_n(vertex_used, t_n_64, vertex_mark);

            moves      = (vertex_t*)aligned_alloc(64, t_n_64 * sizeof(vertex_t));
            moves_size = 0;

            // active block scheduling
            active_this_round = (u8*)aligned_alloc(64, t_k_64 * sizeof(u8));
            active_next_round = (u8*)aligned_alloc(64, t_k_64 * sizeof(u8));
            pairs             = (Pair*)aligned_alloc(64, t_k_t_k_64 * sizeof(Pair));
            pairs_size        = 0;

            gen.seed(m_seed);
            dis = std::uniform_real_distribution<float>(0.0f, 1.0f);
        }

        template <typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle, typename TSerialQuotientGraph>
        void refine(QuotientGraphRefinementFaraj20Configuration& config,
                    [[maybe_unused]] TSerialGraph& g,
                    [[maybe_unused]] TSerialActiveVertexManager& av_manager,
                    [[maybe_unused]] TSerialBoundaryVertexManager& bv_manager,
                    [[maybe_unused]] TSerialPartitionManager& p_manager,
                    [[maybe_unused]] TSerialDistanceOracle& d_oracle,
                    [[maybe_unused]] TSerialQuotientGraph& q_graph) {
            static_assert(std::is_base_of_v<ISerialGraph, TSerialGraph>, "TSerialGraph must inherit from ISerialGraph");
            static_assert(std::is_base_of_v<ISerialActiveVertexManager, TSerialActiveVertexManager>, "TSerialActiveVertexManager must inherit from ISerialActiveVertexManager");
            static_assert(std::is_base_of_v<ISerialBoundaryVertexManager, TSerialBoundaryVertexManager>, "TSerialBoundaryVertexManager must inherit from ISerialBoundaryVertexManager");
            static_assert(std::is_base_of_v<ISerialPartitionManager, TSerialPartitionManager>, "TSerialPartitionManager must inherit from ISerialPartitionManager");
            static_assert(std::is_base_of_v<ISerialDistanceOracle, TSerialDistanceOracle>, "TSerialDistanceOracle must inherit from ISerialDistanceOracle");
            static_assert(std::is_base_of_v<ISerialQuotientGraph, TSerialQuotientGraph>, "TSerialQuotientGraph must inherit from ISerialQuotientGraph");

            std::fill_n(active_this_round, m_k, 1);
            std::fill_n(active_next_round, m_k, 0);

            bool one_pair_active = true;
            for (u64 iteration = 0; iteration < config.max_iteration && one_pair_active; ++iteration) {
                one_pair_active = false;

                // determine all pairs in the quotient graph
                pairs_size = 0;
                for (partition_t u_id = 0; u_id < m_k; ++u_id) {
                    for (partition_t v_id = u_id + 1; v_id < m_k; ++v_id) {
                        if (q_graph.has_edge(u_id, v_id) && (active_this_round[u_id] || active_this_round[v_id])) {
                            pairs[pairs_size++] = {u_id, v_id};
                        }
                    }
                }

                // for each pair do fm refinement
                for (size_t j = 0; j < pairs_size; ++j) {
                    auto [u_id, v_id] = pairs[j];

                    // add all boundary vertices with gain
                    boundary_vertices_u.clear();
                    for (vertex_t u : bv_manager[u_id]) {
                        for (const auto [v, w] : g[u]) {
                            if (p_manager[v] == v_id) {
                                // u is connected to block v_id
                                s64 qap_delta_u = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                                boundary_vertices_u.push(u, qap_delta_u);
                                break;
                            }
                        }
                    }

                    boundary_vertices_v.clear();
                    for (vertex_t v : bv_manager[v_id]) {
                        for (const auto [u, w] : g[v]) {
                            if (p_manager[u] == u_id) {
                                // v is connected to block u_id
                                s64 qap_delta = get_u_qap_delta(g, v, v_id, u_id, p_manager, d_oracle);
                                boundary_vertices_v.push(v, qap_delta);
                                break;
                            }
                        }
                    }


                    // start executing moves based on the TopGain method
                    vertex_mark += 1;
                    moves_size                   = 0;
                    best_idx                     = 0;
                    curr_qap_gain                = 0;
                    max_qap_gain                 = 0;
                    u32 moves_since_last_maximum = 0;
                    while (!boundary_vertices_u.empty() || !boundary_vertices_v.empty()) {

                        // remove vertex from u if it is not boundary
                        if (!boundary_vertices_u.empty() && !is_boundary(g, p_manager, boundary_vertices_u.top_key())) {
                            boundary_vertices_u.pop();
                            continue;
                        }

                        // remove vertex from v if it is not boundary
                        if (!boundary_vertices_v.empty() && !is_boundary(g, p_manager, boundary_vertices_v.top_key())) {
                            boundary_vertices_v.pop();
                            continue;
                        }


                        // determine from which block to choose
                        bool choose_u = true;
                        // 1. if one block is empty, then choose the other one
                        if (boundary_vertices_u.empty() || boundary_vertices_v.empty()) {
                            if (boundary_vertices_u.empty()) {
                                choose_u = false;
                            }
                        } else {
                            // 2. choose the block with greater gain and randomly if even
                            if (boundary_vertices_v.top() > boundary_vertices_u.top()) {
                                choose_u = false;
                            } else if (boundary_vertices_v.top() == boundary_vertices_u.top()) {
                                choose_u = dis(gen) < 0.5;
                            }

                            // 3. if one block is overloaded, choose the larger one, if both same sizes, then randomly
                            weight_t u_id_weight = p_manager.get_bweight(u_id);
                            weight_t v_id_weight = p_manager.get_bweight(v_id);

                            if (u_id_weight > m_lmax && u_id_weight > v_id_weight) { choose_u = true; }
                            if (v_id_weight > m_lmax && v_id_weight > u_id_weight) { choose_u = false; }
                            if (u_id_weight > m_lmax && v_id_weight > m_lmax && u_id_weight == v_id_weight) { choose_u = dis(gen) < 0.5; }
                        }

                        // choose the priority queue
                        IndexedMaxHeap<s64>& boundary_vertices = choose_u ? boundary_vertices_u : boundary_vertices_v;
                        vertex_t vertex                        = boundary_vertices.top_key();
                        partition_t vertex_id                  = p_manager[vertex];
                        partition_t move_id                    = choose_u ? v_id : u_id;
                        weight_t vertex_weight                 = g.get_weight(vertex);
                        s64 qap_delta                          = boundary_vertices.top();

                        boundary_vertices.pop();

                        // move the vertex
                        moves[moves_size++] = vertex;
                        curr_qap_gain += qap_delta;
                        moves_since_last_maximum += 1;
                        if (curr_qap_gain > max_qap_gain) {
                            best_idx                 = moves_size;
                            max_qap_gain             = curr_qap_gain;
                            moves_since_last_maximum = 0;
                        }
                        vertex_used[vertex] = vertex_mark;

                        // make move in structures
                        p_manager.move(vertex, vertex_weight, vertex_id, move_id);

                        // break search if too many moves without improvement
                        if (moves_since_last_maximum > config.max_moves_without_max) { break; }

                        // we have to push or update the neighbors that were not moved already
                        for (const auto [neighbor, w] : g[vertex]) {
                            partition_t neighbor_id = p_manager[neighbor];
                            if (vertex_used[neighbor] == vertex_mark || !(neighbor_id == u_id || neighbor_id == v_id) || !is_boundary(g, p_manager, vertex)) {
                                continue;
                            }

                            partition_t new_id = neighbor_id == vertex_id ? move_id : vertex_id;

                            s64 new_qap_delta = get_u_qap_delta(g, neighbor, neighbor_id, new_id, p_manager, d_oracle);

                            if (neighbor_id == u_id) {
                                boundary_vertices_u.push_update(neighbor, new_qap_delta);
                            } else {
                                boundary_vertices_v.push_update(neighbor, new_qap_delta);
                            }
                        }
                    }

                    // revert all moves in partitioning manager
                    for (size_t i = 0; i < moves_size; i++) {
                        vertex_t vertex        = moves[moves_size - 1 - i];
                        weight_t vertex_weight = g.get_weight(vertex);
                        partition_t vertex_id  = p_manager[vertex];
                        partition_t move_id    = u_id == vertex_id ? v_id : u_id;

                        p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                    }

                    // make all moves to best index
                    for (size_t i = 0; i < best_idx; ++i) {
                        vertex_t vertex        = moves[i];
                        weight_t vertex_weight = g.get_weight(vertex);
                        partition_t vertex_id  = p_manager[vertex];
                        partition_t move_id    = u_id == vertex_id ? v_id : u_id;

                        bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                        q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                        p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                    }

                    if (max_qap_gain > 0) {
                        active_next_round[u_id] = 1;
                        active_next_round[v_id] = 1;
                        one_pair_active         = true;
                    }
                }
                std::swap(active_this_round, active_next_round);
                std::fill_n(active_next_round, m_k, 0);
            }
        }
    };
}

#endif //HEIPROMAP_QUOTIENT_GRAPH_REFINEMENT_FARAJ20_H
