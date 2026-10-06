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

#ifndef HEIPROMAP_SIMPLE_QUOTIENT_GRAPH_REFINEMENT_H
#define HEIPROMAP_SIMPLE_QUOTIENT_GRAPH_REFINEMENT_H

#include <cmath>
#include <limits>
#include <vector>

#include <omp.h>

#include "../definitions.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/csr_graph.h"
#include "../distance_oracles/matrix_distance_oracle.h"
#include "../datastructures/partition_manager.h"
#include "../datastructures/distance_1_matching.h"
#include "../datastructures/active_block_scheduling.h"
#include "../datastructures/block_adjacency.h"
#include "../utility/aligned_array.h"
#include "../utility/profiler.h"
#include "../utility/random_engine.h"
#include "../utility/indexed_max_heap.h"
#include "simple_label_propagation_refinement.h"

namespace HeiProMap {
    class SimpleQuotientGraphRefinementConfiguration {
    public:
        explicit SimpleQuotientGraphRefinementConfiguration(std::string t_name) : name(std::move(t_name)) {}

        std::string name;
        bool enabled = false;
        u64 max_iteration = 1;
        u64 min_n_steps = 3;
        f64 alpha = 5.0;
        bool use_preemptive_exit = true;
        bool use_edge_cut = true;
    };

    template<bool LARGE_K>
    class SimpleQuotientGraphRefinement {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
        u64 m_threads = 1;

        std::vector<RandomEngine> rnd_engines;
        SimpleQuotientGraphRefinementConfiguration config = SimpleQuotientGraphRefinementConfiguration("simple_qg");

        ActiveBlockScheduling active_block_scheduling;
        Distance1Matching<LARGE_K> d1_matcher;
        BlockAdjacency block_adjacency;

        AlignedArray<u32> vertex_used;
        u32 global_vertex_mark = 0;

        // Dense ID mapping: shared across threads (matching ensures disjoint vertex sets)
        AlignedArray<vertex_t> m_dense_id;
        AlignedArray<u32> m_dense_epoch;
        u32 m_dense_global_epoch = 0;

        std::vector<AlignedArray<vertex_t>> thread_moves;
        std::vector<IndexedMaxHeap<weight_t>> thread_heap_u;
        std::vector<IndexedMaxHeap<weight_t>> thread_heap_v;
        std::vector<AlignedArray<std::pair<size_t, weight_t>>> thread_seed_u;
        std::vector<AlignedArray<std::pair<size_t, weight_t>>> thread_seed_v;
        std::vector<AlignedArray<vertex_t>> thread_dense_to_vertex;

        template<bool t_uniform_e_weights, typename DistanceOracleT>
        static inline __attribute__((always_inline)) weight_t compute_qap_delta(const graph_t &g,
                                                                                const vertex_t u,
                                                                                const partition_t old_id,
                                                                                const partition_t new_id,
                                                                                const p_manager_t &p_manager,
                                                                                DistanceOracleT &d_oracle) {
            weight_t qap_delta = 0;
            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t v = g.edges_v[i];
                const partition_t v_id = p_manager[v];
                if constexpr (t_uniform_e_weights) {
                    qap_delta += d_oracle.get(v_id, old_id) - d_oracle.get(v_id, new_id);
                } else {
                    const weight_t w = g.edges_w[i];
                    qap_delta += (d_oracle.get(v_id, old_id) - d_oracle.get(v_id, new_id)) * w;
                }
            }
            return qap_delta;
        }

        template<bool t_uniform_e_weights>
        static inline __attribute__((always_inline)) weight_t compute_edge_cut_delta(const graph_t &g,
                                                                                     const vertex_t u,
                                                                                     const partition_t old_id,
                                                                                     const partition_t new_id,
                                                                                     const p_manager_t &p_manager) {
            weight_t delta = 0;
            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t v = g.edges_v[i];
                const partition_t v_id = p_manager[v];
                if constexpr (t_uniform_e_weights) {
                    delta -= (v_id != new_id);
                    delta += (v_id != old_id);
                } else {
                    const weight_t w = g.edges_w[i];
                    delta -= (v_id != new_id) * w;
                    delta += (v_id != old_id) * w;
                }
            }
            return delta;
        }

        template<bool t_uniform_e_weights, bool t_use_edge_cut, typename DistanceOracleT>
        static inline weight_t compute_gain_and_connected(const graph_t &g,
                                                          const vertex_t u,
                                                          const partition_t old_id,
                                                          const partition_t new_id,
                                                          bool &is_connected,
                                                          const p_manager_t &p_manager,
                                                          DistanceOracleT &d_oracle) {
            is_connected = false;
            weight_t delta = 0;
            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                const vertex_t v = g.edges_v[i];
                const partition_t v_id = p_manager[v];
                is_connected |= (v_id == new_id);
                if constexpr (t_use_edge_cut) {
                    if constexpr (t_uniform_e_weights) {
                        delta -= (v_id != new_id);
                        delta += (v_id != old_id);
                    } else {
                        const weight_t w = g.edges_w[i];
                        delta -= (v_id != new_id) * w;
                        delta += (v_id != old_id) * w;
                    }
                } else {
                    if constexpr (t_uniform_e_weights) {
                        delta += d_oracle.get(v_id, old_id) - d_oracle.get(v_id, new_id);
                    } else {
                        const weight_t w = g.edges_w[i];
                        delta += (d_oracle.get(v_id, old_id) - d_oracle.get(v_id, new_id)) * w;
                    }
                }
            }
            return delta;
        }

        static inline bool is_connected_to(const graph_t &g,
                                            const p_manager_t &p_manager,
                                            const vertex_t u,
                                            const partition_t id) {
            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                if (p_manager[g.edges_v[i]] == id) return true;
            }
            return false;
        }

    public:
        SimpleQuotientGraphRefinement() = default;
        ~SimpleQuotientGraphRefinement() = default;

        size_t heap_bytes() const {
            size_t bytes = active_block_scheduling.heap_bytes() + d1_matcher.heap_bytes()
                         + block_adjacency.heap_bytes() + vertex_used.heap_bytes()
                         + m_dense_id.heap_bytes() + m_dense_epoch.heap_bytes();
            for (size_t i = 0; i < thread_moves.size(); ++i) {
                bytes += thread_moves[i].heap_bytes();
                bytes += thread_dense_to_vertex[i].heap_bytes();
            }
            for (size_t i = 0; i < thread_heap_u.size(); ++i) {
                bytes += thread_heap_u[i].heap_bytes() + thread_heap_v[i].heap_bytes();
                bytes += thread_seed_u[i].heap_bytes();
                bytes += thread_seed_v[i].heap_bytes();
            }
            return bytes;
        }

        void set_pool(MemoryPool *pool) {
            active_block_scheduling.set_pool(pool);
            d1_matcher.set_pool(pool);
            block_adjacency.set_pool(pool);
            vertex_used.set_pool(pool);
            m_dense_id.set_pool(pool);
            m_dense_epoch.set_pool(pool);
        }

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const u64 t_threads,
                        const u64 t_seed,
                        const SimpleQuotientGraphRefinementConfiguration &t_config,
                        MemoryPool *pool = nullptr) {
            set_pool(pool);
            HEIPROMAP_PROFILE_SCOPE("misc", "SimpleQGRefinement", "initialize");

            m_n = t_n;
            m_m = t_m;
            m_k = t_k;
            m_threads = t_threads;
            config = t_config;

            rnd_engines.resize(m_threads);
            for (u64 t = 0; t < m_threads; ++t) {
                rnd_engines[t] = RandomEngine(t_seed + t);
            }

            active_block_scheduling.initialize(m_k, pool);
            d1_matcher.initialize(m_k, pool);
            block_adjacency.initialize(m_k, m_threads, pool);

            global_vertex_mark = 0;
            vertex_used.initialize(m_n, 0);

            m_dense_id.initialize(m_n);
            m_dense_epoch.initialize(m_n, 0);
            m_dense_global_epoch = 0;

            thread_moves.resize(m_threads);
            thread_heap_u.resize(m_threads);
            thread_heap_v.resize(m_threads);
            for (u64 t = 0; t < m_threads; ++t) {
                thread_moves[t].set_pool(pool);
                thread_heap_u[t].set_pool(pool);
                thread_heap_v[t].set_pool(pool);
            }
            thread_seed_u.resize(m_threads);
            thread_seed_v.resize(m_threads);
            thread_dense_to_vertex.resize(m_threads);
            for (u64 t = 0; t < m_threads; ++t) {
                thread_seed_u[t].set_pool(pool);
                thread_seed_v[t].set_pool(pool);
                thread_dense_to_vertex[t].set_pool(pool);
            }
        }

        template<typename DistanceOracleT>
        void refine(graph_t &g,
                    DistanceOracleT &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager) {
            if (g.uniform_v_weights && g.uniform_e_weights) {
                refine_impl<true, true>(g, d_oracle, bv_manager, p_manager);
            } else if (g.uniform_v_weights) {
                refine_impl<true, false>(g, d_oracle, bv_manager, p_manager);
            } else if (g.uniform_e_weights) {
                refine_impl<false, true>(g, d_oracle, bv_manager, p_manager);
            } else {
                refine_impl<false, false>(g, d_oracle, bv_manager, p_manager);
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights, typename DistanceOracleT>
        void refine_impl(graph_t &g,
                         DistanceOracleT &d_oracle,
                         bv_manager_t &bv_manager,
                         p_manager_t &p_manager) {
            active_block_scheduling.reset(m_k);

            std::vector<std::vector<std::pair<partition_t, partition_t>>> static_matchings;

            for (u64 iteration = 0; iteration < config.max_iteration; ++iteration) {
                HEIPROMAP_PROFILE_SCOPE("refinement", "SimpleQGRefinement", "recompute_adjacency");
                block_adjacency.compute(g, p_manager, bv_manager, m_threads);

                d1_matcher.reset_used_edges();
                bool positive_move_occurred = false;

                {
                    HEIPROMAP_PROFILE_SCOPE("refinement", "SimpleQGRefinement", "matching");
                    d1_matcher.compute_static_matchings(block_adjacency, p_manager, active_block_scheduling.active_this_round, static_matchings, 512 * 512, 8);
                }

                if (static_matchings.empty()) break;

                for (size_t round = 0; round < static_matchings.size(); ++round) {
                    const auto &matching = static_matchings[round];
                    if (matching.empty()) continue;

                    u32 base_mark = global_vertex_mark + 1;
                    global_vertex_mark += static_cast<u32>(matching.size());
                    m_dense_global_epoch++;
                    if (m_dense_global_epoch == 0) {
                        m_dense_epoch.initialize(m_n, 0);
                        m_dense_global_epoch = 1;
                    }

                    HEIPROMAP_PROFILE_SCOPE("refinement", "SimpleQGRefinement", "fm_pass");

                    #pragma omp parallel for num_threads(m_threads) schedule(dynamic)
                    for (size_t j = 0; j < matching.size(); ++j) {
                        const partition_t u_id = matching[j].first;
                        const partition_t v_id = matching[j].second;
                        const u64 tid = omp_get_thread_num();
                        const u32 mark = base_mark + static_cast<u32>(j);

                        const bool last_level_pair = config.use_edge_cut && d_oracle.last_level_pair(u_id, v_id);

                        bool improved;
                        if (last_level_pair) {
                            improved = refine_pair<t_uniform_v_weights, t_uniform_e_weights, true>(
                                g, d_oracle, bv_manager, p_manager,
                                u_id, v_id, mark,
                                thread_moves[tid], thread_heap_u[tid], thread_heap_v[tid],
                                thread_seed_u[tid], thread_seed_v[tid],
                                thread_dense_to_vertex[tid], rnd_engines[tid]);
                        } else {
                            improved = refine_pair<t_uniform_v_weights, t_uniform_e_weights, false>(
                                g, d_oracle, bv_manager, p_manager,
                                u_id, v_id, mark,
                                thread_moves[tid], thread_heap_u[tid], thread_heap_v[tid],
                                thread_seed_u[tid], thread_seed_v[tid],
                                thread_dense_to_vertex[tid], rnd_engines[tid]);
                        }

                        if (improved) {
                            active_block_scheduling.activate(u_id, v_id);
                            #pragma omp atomic write
                            positive_move_occurred = true;
                        }
                    }
                }

                if (!positive_move_occurred) break;

                active_block_scheduling.next_round();
            }
        }

    private:
        vertex_t get_or_assign_dense(const vertex_t v, AlignedArray<vertex_t> &dense_to_vertex,
                                     vertex_t &dense_counter) {
            if (m_dense_epoch[v] == m_dense_global_epoch) {
                return m_dense_id[v];
            }
            vertex_t did = dense_counter++;
            m_dense_id[v] = did;
            m_dense_epoch[v] = m_dense_global_epoch;
            dense_to_vertex[did] = v;
            return did;
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights, bool t_use_edge_cut, typename DistanceOracleT>
        bool refine_pair(graph_t &g,
                         DistanceOracleT &d_oracle,
                         bv_manager_t &bv_manager,
                         p_manager_t &p_manager,
                         const partition_t u_id,
                         const partition_t v_id,
                         const u32 mark,
                         AlignedArray<vertex_t> &moves,
                         IndexedMaxHeap<weight_t> &heap_u,
                         IndexedMaxHeap<weight_t> &heap_v,
                         AlignedArray<std::pair<size_t, weight_t>> &seed_u,
                         AlignedArray<std::pair<size_t, weight_t>> &seed_v,
                         AlignedArray<vertex_t> &dense_to_vertex,
                         RandomEngine &rng) {
            const f64 alpha = config.alpha * (f64) d_oracle.get(u_id, v_id);
            const f64 beta = std::log(g.n) * (f64) d_oracle.get(u_id, v_id);

            const size_t estimated_capacity = p_manager.n_vertices[u_id] + p_manager.n_vertices[v_id];
            heap_u.ensure_capacity(estimated_capacity);
            heap_v.ensure_capacity(estimated_capacity);

            vertex_t dense_counter = 0;
            if (dense_to_vertex.size() < estimated_capacity) {
                dense_to_vertex.initialize(estimated_capacity);
            }
            if (seed_u.size() < estimated_capacity) {
                seed_u.initialize(estimated_capacity);
            }
            if (seed_v.size() < estimated_capacity) {
                seed_v.initialize(estimated_capacity);
            }

            // Seed heaps: batch insert with O(n) heapify
            vertex_t seed_u_size = 0;
            vertex_t seed_v_size = 0;
            for (const vertex_t u : bv_manager.boundary(u_id)) {
                bool connected;
                weight_t delta = compute_gain_and_connected<t_uniform_e_weights, t_use_edge_cut>(
                    g, u, u_id, v_id, connected, p_manager, d_oracle);
                if (connected) {
                    vertex_t did = get_or_assign_dense(u, dense_to_vertex, dense_counter);
                    seed_u[seed_u_size++] = {did, delta};
                }
            }

            for (const vertex_t v : bv_manager.boundary(v_id)) {
                bool connected;
                weight_t delta = compute_gain_and_connected<t_uniform_e_weights, t_use_edge_cut>(
                    g, v, v_id, u_id, connected, p_manager, d_oracle);
                if (connected) {
                    vertex_t did = get_or_assign_dense(v, dense_to_vertex, dense_counter);
                    seed_v[seed_v_size++] = {did, delta};
                }
            }

            heap_u.push_many_heapify(seed_u.get_ptr(), seed_u_size);
            heap_v.push_many_heapify(seed_v.get_ptr(), seed_v_size);
            const u64 n_init_moves = seed_u_size + seed_v_size;

            if (n_init_moves == 0) return false;

            weight_t curr_qap_gain = 0;
            weight_t max_qap_gain = 0;
            size_t best_idx = 0;

            u64 steps_since_last_improvement = 0;
            f64 qap_gain_mean = 0.0;
            f64 qap_gain_var = 0.0;

            const weight_t u_id_weight_initial = p_manager.get_bweight(u_id);
            const weight_t v_id_weight_initial = p_manager.get_bweight(v_id);
            bool best_is_balanced = (u_id_weight_initial <= p_manager.lmax[u_id] && v_id_weight_initial <= p_manager.lmax[v_id]);

            vertex_t moves_size = 0;
            if (moves.size() < estimated_capacity) {
                moves.initialize(estimated_capacity);
            }

            // FM exploration with lazy connectivity check
            while ((!heap_u.empty() || !heap_v.empty()) && moves_size < n_init_moves) {
                while (!heap_u.empty() && !is_connected_to(g, p_manager, dense_to_vertex[heap_u.top_key()], v_id)) { heap_u.pop(); }
                while (!heap_v.empty() && !is_connected_to(g, p_manager, dense_to_vertex[heap_v.top_key()], u_id)) { heap_v.pop(); }
                if (heap_u.empty() && heap_v.empty()) break;

                bool choose_u = true;

                if (heap_u.empty() || heap_v.empty()) {
                    choose_u = heap_v.empty();
                } else {
                    if (heap_v.top() > heap_u.top()) {
                        choose_u = false;
                    } else if (heap_v.top() == heap_u.top()) {
                        choose_u = rng.get_f32() < 0.5f;
                    }

                    const weight_t u_w = p_manager.get_bweight(u_id);
                    const weight_t v_w = p_manager.get_bweight(v_id);
                    if (u_w > p_manager.lmax[u_id] && u_w > v_w) { choose_u = true; }
                    if (v_w > p_manager.lmax[v_id] && v_w > u_w) { choose_u = false; }
                    if (u_w > p_manager.lmax[u_id] && v_w > p_manager.lmax[v_id] && u_w == v_w) {
                        choose_u = rng.get_f32() < 0.5f;
                    }
                }

                IndexedMaxHeap<weight_t> &heap = choose_u ? heap_u : heap_v;

                const vertex_t vertex = dense_to_vertex[heap.top_key()];
                const weight_t qap_delta = heap.top();
                const weight_t vertex_weight = t_uniform_v_weights ? 1 : g.v_weights[vertex];
                const partition_t vertex_id = choose_u ? u_id : v_id;
                const partition_t move_id = choose_u ? v_id : u_id;
                heap.pop();

                moves[moves_size++] = vertex;
                curr_qap_gain += qap_delta;

                const weight_t current_move_w = p_manager.get_bweight(move_id) + vertex_weight;
                const weight_t current_vertex_w = p_manager.get_bweight(vertex_id) - vertex_weight;
                const bool current_is_balanced = (current_move_w <= p_manager.lmax[move_id] && current_vertex_w <= p_manager.lmax[vertex_id]);

                bool update_best = false;
                if (current_is_balanced) {
                    if (!best_is_balanced || curr_qap_gain >= max_qap_gain) {
                        update_best = true;
                    }
                } else if (!best_is_balanced) {
                    const weight_t move_id_weight_initial = move_id == u_id ? u_id_weight_initial : v_id_weight_initial;
                    if (curr_qap_gain >= max_qap_gain && current_move_w <= std::max(move_id_weight_initial, p_manager.lmax[move_id])) {
                        update_best = true;
                    }
                }

                if (update_best) {
                    best_idx = moves_size;
                    max_qap_gain = curr_qap_gain;
                    best_is_balanced = current_is_balanced;
                    steps_since_last_improvement = 0;
                    qap_gain_mean = 0.0;
                    qap_gain_var = 0.0;
                }

                p_manager.move_serial(vertex, vertex_weight, vertex_id, move_id);
                vertex_used[vertex] = mark;

                steps_since_last_improvement++;
                const f64 new_mean = qap_gain_mean + ((f64) qap_delta - qap_gain_mean) / (f64) steps_since_last_improvement;
                const f64 new_var = (qap_gain_var + ((f64) qap_delta - qap_gain_mean) * ((f64) qap_delta - new_mean)) / (f64) steps_since_last_improvement;
                qap_gain_mean = new_mean;
                qap_gain_var = new_var;

                if (config.use_preemptive_exit) {
                    if (steps_since_last_improvement > config.min_n_steps &&
                        (f64) steps_since_last_improvement * qap_gain_mean * qap_gain_mean > alpha * qap_gain_var + beta) {
                        break;
                    }
                }

                // Update neighbor gains
                for (size_t i = g.neighborhoods[vertex]; i < g.neighborhoods[vertex + 1]; ++i) {
                    const vertex_t neighbor = g.edges_v[i];
                    if (vertex_used[neighbor] == mark) continue;

                    const partition_t neighbor_id = p_manager[neighbor];
                    if (neighbor_id != u_id && neighbor_id != v_id) continue;

                    const partition_t new_id = neighbor_id == vertex_id ? move_id : vertex_id;

                    bool is_connected;
                    weight_t new_delta = compute_gain_and_connected<t_uniform_e_weights, t_use_edge_cut>(
                        g, neighbor, neighbor_id, new_id, is_connected, p_manager, d_oracle);

                    if (!is_connected) continue;

                    if (dense_counter >= dense_to_vertex.size()) {
                        AlignedArray<vertex_t> new_dtv(dense_to_vertex.get_pool());
                        new_dtv.initialize(dense_counter * 2);
                        std::memcpy(new_dtv.get_ptr(), dense_to_vertex.get_ptr(), dense_counter * sizeof(vertex_t));
                        dense_to_vertex = std::move(new_dtv);
                    }
                    vertex_t ndid = get_or_assign_dense(neighbor, dense_to_vertex, dense_counter);
                    if (ndid >= heap_u.capacity()) {
                        heap_u.grow(dense_counter);
                        heap_v.grow(dense_counter);
                    }

                    if (neighbor_id == u_id) {
                        heap_u.push_update(ndid, new_delta);
                    } else {
                        heap_v.push_update(ndid, new_delta);
                    }
                }
            }

            // Revert all moves
            for (size_t i = 0; i < moves_size; i++) {
                const vertex_t vertex = moves[moves_size - 1 - i];
                const weight_t vertex_weight = t_uniform_v_weights ? 1 : g.v_weights[vertex];
                const partition_t vertex_id = p_manager[vertex];
                const partition_t move_id = u_id == vertex_id ? v_id : u_id;
                vertex_used[vertex] = global_vertex_mark - 1;
                p_manager.move_serial(vertex, vertex_weight, vertex_id, move_id);
            }

            // Apply moves up to best_idx
            for (size_t i = 0; i < best_idx; ++i) {
                const vertex_t vertex = moves[i];
                const weight_t vertex_weight = t_uniform_v_weights ? 1 : g.v_weights[vertex];
                const partition_t vertex_id = p_manager[vertex];
                const partition_t move_id = u_id == vertex_id ? v_id : u_id;
                vertex_used[vertex] = mark;

                bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                p_manager.move_serial(vertex, vertex_weight, vertex_id, move_id);
            }

            return max_qap_gain > 0;
        }
    };
}

#endif //HEIPROMAP_SIMPLE_QUOTIENT_GRAPH_REFINEMENT_H
