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

#ifndef HEIPROMAP_NEGATIVE_CYCLE_DETECTION_H
#define HEIPROMAP_NEGATIVE_CYCLE_DETECTION_H

#include <cmath>
#include <limits>
#include <vector>

#include <omp.h>

#include "../definitions.h"
#include "../datastructures/block_conn.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/distance_oracle.h"
#include "../datastructures/partition_manager.h"
#include "../datastructures/quotient_graph.h"
#include "../utility/aligned_array.h"
#include "../utility/profiler.h"
#include "../utility/qap.h"
#include "../utility/random_engine.h"

namespace HeiProMap {
    class NegativeCycleConfiguration final {
    public:
        NegativeCycleConfiguration() = default;

        explicit NegativeCycleConfiguration(const std::string &t_name) {
            name = t_name;
        }

        std::string name = "Negative Cycle Refinement";
        bool enabled = false;
        u64 random_tries = 10;
        u64 max_iterations = 10;
        u64 max_path_length = 8;
        u64 threshold = 10;
    };

    class NegativeCycleRefinement final {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
        u64 m_threads = 1;

        struct Move {
            vertex_t u;
            partition_t source_id;
            partition_t target_id;
            weight_t qap_delta;
        };

        std::vector<int> m_invalid_at_depth;
        std::vector<u64> m_invalid_by_move_id;
        std::vector<u64> m_active_move_id_at_depth;
        u64 m_next_move_id = 1;

        std::vector<bool> m_block_in_use;

        std::vector<std::vector<int> > m_all_to_all_dists;
        bool m_distances_invalid = true;

        std::vector<weight_t> m_top_qap_prefix_sums;

        struct MoveGraph {
            std::vector<std::vector<Move> > adj_list;
        };

        MoveGraph move_graph;

        NegativeCycleConfiguration config;
        RandomEngine random_engine;

    public:
        NegativeCycleRefinement() = default;

        ~NegativeCycleRefinement() = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const u64 t_threads,
                        const u64 seed,
                        const NegativeCycleConfiguration &i_config) {
            m_n = t_n;
            m_m = t_m;
            m_k = t_k;
            m_threads = t_threads;

            config = i_config;
            random_engine = RandomEngine(seed);

            m_invalid_at_depth.assign(m_n, -1);
            m_invalid_by_move_id.assign(m_n, 0);
            m_active_move_id_at_depth.assign(config.max_path_length + 1, 0);
            m_next_move_id = 1;

            m_block_in_use.assign(m_k, false);
            m_all_to_all_dists.assign(m_k, std::vector<int>(m_k, -1));
            m_distances_invalid = true;
            m_top_qap_prefix_sums.assign(config.max_path_length + 1, 0);
            move_graph.adj_list.resize(m_k);
        }

        void refine(graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph,
                    block_conn_t &block_conn,
                    const AlignedArray<weight_t> &lmax_constraints,
                    bool uniform_v_weights,
                    bool uniform_e_weights) {
            if (uniform_v_weights && uniform_e_weights) refine_impl<true, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
            else if (uniform_v_weights) refine_impl<true, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
            else if (uniform_e_weights) refine_impl<false, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
            else refine_impl<false, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn, lmax_constraints);
        }

    private:
        void compute_all_to_all_distances(q_graph_t &q_graph) {
            if (!m_distances_invalid) return;

            for (partition_t start_id = 0; start_id < m_k; ++start_id) {
                auto &dists = m_all_to_all_dists[start_id];
                std::fill(dists.begin(), dists.end(), -1);
                std::vector<partition_t> queue;
                queue.reserve(m_k);

                dists[start_id] = 0;
                queue.push_back(start_id);

                size_t head = 0;
                while (head < queue.size()) {
                    partition_t u = queue[head++];
                    int d = dists[u];

                    q_graph.for_each_neighbor(u, [&](const partition_t v, const weight_t) {
                        if (dists[v] == -1) {
                            dists[v] = d + 1;
                            queue.push_back(v);
                        }
                    });
                }
            }
            m_distances_invalid = false;
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        bool find_cycle(partition_t curr_id,
                        partition_t start_id,
                        std::vector<Move> &path,
                        weight_t &total_qap_delta,
                        weight_t min_qap_delta,
                        weight_t max_qap_delta,
                        graph_t &g,
                        p_manager_t &p_manager,
                        std::vector<weight_t> &curr_bweights,
                        const AlignedArray<weight_t> &lmax_constraints) {
            if (path.size() >= config.max_path_length) return false;

            // Pruning: if we cannot achieve a positive total_qap_delta even with the best remaining moves
            weight_t remaining_steps = config.max_path_length - path.size();
            if (total_qap_delta + m_top_qap_prefix_sums[remaining_steps] <= 0) return false;

            auto &moves = move_graph.adj_list[curr_id];
            
            m_active_move_id_at_depth[path.size()] = m_next_move_id++;

            for (size_t i = 0; i < std::min<size_t>(moves.size(), config.threshold); ++i) {
                const auto &m = moves[i];

                int d = m_invalid_at_depth[m.u];
                if (d != -1 && d < (int)path.size() && m_invalid_by_move_id[m.u] == m_active_move_id_at_depth[d]) continue;
                if (m_block_in_use[m.target_id]) continue;

                // Reachability Pruning: if target_id is start_id, it closes the cycle.
                // If not, we check if we can reach start_id from target_id within remaining steps.
                if (m.target_id != start_id) {
                    if (m_all_to_all_dists[m.target_id][start_id] == -1 || m_all_to_all_dists[m.target_id][start_id] > (int) remaining_steps - 1) {
                        continue;
                    }
                }

                weight_t u_w = t_uniform_v_weights ? 1 : g.v_weights[m.u];

                // Check balance: curr_bweight + entering_w - u_w <= lmax
                weight_t entering_w = 0;
                if (!path.empty()) {
                    entering_w = t_uniform_v_weights ? 1 : g.v_weights[path.back().u];
                }
                if (curr_bweights[curr_id] + entering_w - u_w > lmax_constraints[curr_id]) continue;

                if (m.target_id == start_id) {
                    // Check balance for the start block: start_bweight + u_w - first_w <= lmax
                    weight_t first_w = t_uniform_v_weights ? 1 : g.v_weights[path.empty() ? m.u : path[0].u];
                    if (curr_bweights[start_id] + u_w - first_w <= lmax_constraints[start_id]) {
                        if (total_qap_delta + m.qap_delta > 0) {
                            path.push_back(m);
                            total_qap_delta += m.qap_delta;

                            // Mark the final move as invalid so it can be used for filtering
                            m_invalid_at_depth[m.u] = path.size() - 1;
                            m_invalid_by_move_id[m.u] = m_active_move_id_at_depth[path.size() - 1];
                            for (u64 j = g.neighborhoods[m.u]; j < g.neighborhoods[m.u + 1]; ++j) {
                                vertex_t v = g.edges_v[j];
                                m_invalid_at_depth[v] = path.size() - 1;
                                m_invalid_by_move_id[v] = m_active_move_id_at_depth[path.size() - 1];
                            }

                            return true;
                        }
                    }
                }

                // Forward
                path.push_back(m);
                total_qap_delta += m.qap_delta;
                m_block_in_use[m.target_id] = true;
                
                m_invalid_at_depth[m.u] = path.size() - 1;
                m_invalid_by_move_id[m.u] = m_active_move_id_at_depth[path.size() - 1];
                for (u64 j = g.neighborhoods[m.u]; j < g.neighborhoods[m.u + 1]; ++j) {
                    vertex_t v = g.edges_v[j];
                    m_invalid_at_depth[v] = path.size() - 1;
                    m_invalid_by_move_id[v] = m_active_move_id_at_depth[path.size() - 1];
                }

                if (find_cycle<t_uniform_v_weights, t_uniform_e_weights>(m.target_id, start_id, path, total_qap_delta, min_qap_delta, max_qap_delta, g, p_manager, curr_bweights, lmax_constraints)) {
                    return true;
                }

                // Backtrack: NO DECREMENT
                m_block_in_use[m.target_id] = false;
                total_qap_delta -= m.qap_delta;
                path.pop_back();
            }

            return false;
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void refine_impl(graph_t &g,
                         d_oracle_t &d_oracle,
                         bv_manager_t &bv_manager,
                         p_manager_t &p_manager,
                         q_graph_t &q_graph,
                         block_conn_t &block_conn,
                         const AlignedArray<weight_t> &lmax_constraints) {
            for (size_t global_i = 0; global_i < config.max_iterations; global_i++) {
                HEIPROMAP_PROFILE_SCOPE("refinement", "NegativeCycleRefinement", "reset_graph");
                for (partition_t i = 0; i < m_k; ++i) move_graph.adj_list[i].clear();

                HEIPROMAP_PROFILE_SCOPE("refinement", "NegativeCycleRefinement", "find_all_moves");
                // Find all potential moves
                // For efficiency, only consider boundary vertices
                weight_t min_qap_delta = std::numeric_limits<weight_t>::max();
                weight_t max_qap_delta = -std::numeric_limits<weight_t>::max();
                std::vector<weight_t> all_qaps;
                for (partition_t b_id = 0; b_id < m_k; ++b_id) {
                    for (size_t i = 0; i < bv_manager.size(b_id); ++i) {
                        vertex_t u = bv_manager.get(b_id, i);
                        partition_t u_id = p_manager[u];

                        for (u64 j = block_conn.start(u); j < block_conn.end(u); ++j) {
                            partition_t target_id = block_conn.get_id(j);
                            if (target_id == u_id) continue;

                            weight_t qap_delta = get_u_qap_delta_t<t_uniform_e_weights>(g, u, u_id, target_id, p_manager, d_oracle, block_conn);

                            if (qap_delta > -100) {
                            move_graph.adj_list[u_id].push_back({u, u_id, target_id, qap_delta});
                            min_qap_delta = std::min(min_qap_delta, qap_delta);
                            max_qap_delta = std::max(max_qap_delta, qap_delta);
                            all_qaps.push_back(qap_delta);
                            }
                        }
                    }
                }

                HEIPROMAP_PROFILE_SCOPE("refinement", "NegativeCycleRefinement", "sort_top_k");
                // Precompute top-K prefix sums for tighter pruning
                std::sort(all_qaps.begin(), all_qaps.end(), std::greater<weight_t>());
                m_top_qap_prefix_sums[0] = 0;
                for (size_t i = 1; i <= config.max_path_length; ++i) {
                    if (i <= all_qaps.size()) {
                        m_top_qap_prefix_sums[i] = m_top_qap_prefix_sums[i - 1] + all_qaps[i - 1];
                    } else {
                        m_top_qap_prefix_sums[i] = m_top_qap_prefix_sums[i - 1];
                    }
                }

                HEIPROMAP_PROFILE_SCOPE("refinement", "NegativeCycleRefinement", "sort_graph");
                for (partition_t i = 0; i < m_k; ++i) {
                    std::sort(move_graph.adj_list[i].begin(), move_graph.adj_list[i].end(), [](const Move &a, const Move &b) {
                        return a.qap_delta > b.qap_delta;
                    });
                }

                std::vector<weight_t> curr_bweights(p_manager.k);
                for (partition_t i = 0; i < m_k; ++i) curr_bweights[i] = p_manager.get_bweight(i);

                std::vector<Move> path;
                u64 random_tries = config.random_tries;
                std::fill(m_block_in_use.begin(), m_block_in_use.end(), false);
                std::fill(m_active_move_id_at_depth.begin(), m_active_move_id_at_depth.end(), 0);

                for (u64 try_i = 0; try_i < random_tries; ++try_i) {
                    partition_t start_id = random_engine.get_u64() % m_k;
                    if (move_graph.adj_list[start_id].empty()) continue;

                    path.clear();
                    HEIPROMAP_PROFILE_SCOPE("refinement", "NegativeCycleRefinement", "bfs");
                    compute_all_to_all_distances(q_graph);

                    HEIPROMAP_PROFILE_SCOPE("refinement", "NegativeCycleRefinement", "find_cycle");
                    weight_t total_qap_delta = 0;
                    if (find_cycle<t_uniform_v_weights, t_uniform_e_weights>(start_id, start_id, path, total_qap_delta, min_qap_delta, max_qap_delta, g, p_manager, curr_bweights, lmax_constraints)) {
                        m_distances_invalid = true;
                        // apply all moves
                        // std::cout << "path_size: " << path.size() << " qap: " << total_qap_delta << std::endl;
                        HEIPROMAP_PROFILE_SCOPE("refinement", "NegativeCycleRefinement", "apply_path");
                        for (const auto &m: path) {
                            partition_t u_id = p_manager[m.u];
                            weight_t u_weight = t_uniform_v_weights ? 1 : g.v_weights[m.u];
                            bv_manager.move(g, p_manager, m.u, u_id, m.target_id);
                            q_graph.move(g, p_manager, m.u, u_id, m.target_id);
                            block_conn.move(g, m.u, u_id, m.target_id);
                            p_manager.move_serial(m.u, u_weight, u_id, m.target_id);
                            curr_bweights[u_id] -= u_weight;
                            curr_bweights[m.target_id] += u_weight;
                        }

                        // Remove all moves that conflict with the found cycle
                        for (partition_t i = 0; i < m_k; ++i) {
                            auto &adj = move_graph.adj_list[i];
                            adj.erase(std::remove_if(adj.begin(), adj.end(), [&](const Move &m) {
                                int d = m_invalid_at_depth[m.u];
                                return d != -1 && d < (int)path.size() && m_invalid_by_move_id[m.u] == m_active_move_id_at_depth[d];
                            }), adj.end());
                        }

                        // State Cleanup: Reset invalid counts and block in use for the next search
                        std::fill(m_block_in_use.begin(), m_block_in_use.end(), false);
                        std::fill(m_active_move_id_at_depth.begin(), m_active_move_id_at_depth.end(), 0);
                    }
                }
            }
        }
    };
}

#endif //HEIPROMAP_NEGATIVE_CYCLE_DETECTION_H
