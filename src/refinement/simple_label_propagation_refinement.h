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

#ifndef HEIPROMAP_SIMPLE_LABEL_PROPAGATION_REFINEMENT_H
#define HEIPROMAP_SIMPLE_LABEL_PROPAGATION_REFINEMENT_H

#include <cmath>
#include <limits>
#include <vector>

#include <omp.h>

#include "../definitions.h"
#include "../datastructures/block_conn.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/csr_graph.h"
#include "../distance_oracles/distance_oracle.h"
#include "../datastructures/partition_manager.h"
#include "../datastructures/distance_1_matching.h"
#include "../datastructures/active_block_scheduling.h"
#include "../utility/aligned_array.h"
#include "../utility/profiler.h"
#include "../utility/random_engine.h"

namespace HeiProMap {
    class SimpleLabelPropagationConfiguration {
    public:
        explicit SimpleLabelPropagationConfiguration(std::string t_name) : name(std::move(t_name)) {
        }

        std::string name;
        bool enabled = false;
        u64 max_iteration = 25;

        bool force_parallel_alg = false;
        bool use_edge_cut = true;
        bool use_active_scheduling = true;
        bool use_static_matchings = true;

        u64 max_matching_rounds = 512*512;
        u64 min_matching_threshold = 8;

        bool enable_q_graph = true;
        bool enable_block_conn = true;
    };

    class BlockAdjacency {
    public:
        struct Entry {
            partition_t target = 0;
            mutable u32 used_epoch = 0;
        };

    private:
        partition_t m_k = 0;
        std::vector<std::vector<Entry>> m_adj;
        std::vector<AlignedArray<u32>> m_seen;
        std::vector<u32> m_epoch;

    public:
        void initialize(const partition_t k) {
            m_k = k;
            m_adj.clear();
            m_adj.resize(m_k);
        }

        void initialize(const partition_t k, const u64 num_threads) {
            m_k = k;
            m_adj.clear();
            m_adj.resize(m_k);
            m_seen.resize(num_threads);
            m_epoch.assign(num_threads, 0);
            for (u64 t = 0; t < num_threads; ++t) {
                m_seen[t].initialize(m_k, 0);
            }
        }

        template<typename GraphT, typename PartitionManagerT, typename BoundaryVertexManagerT>
        void compute(const GraphT &g,
                     const PartitionManagerT &p_manager,
                     const BoundaryVertexManagerT &bv_manager,
                     const u64 num_threads) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "BlockAdjacency", "compute");

            if (m_seen.size() < num_threads) {
                m_seen.resize(num_threads);
                m_epoch.resize(num_threads, 0);
                for (u64 t = 0; t < num_threads; ++t) {
                    if (m_seen[t].size() < m_k) m_seen[t].initialize(m_k, 0);
                }
            }

            #pragma omp parallel num_threads(num_threads)
            {
                const u64 tid = omp_get_thread_num();
                auto &seen = m_seen[tid];
                u32 &epoch = m_epoch[tid];

                #pragma omp for schedule(static)
                for (partition_t id = 0; id < m_k; ++id) {
                    auto &adj = m_adj[id];
                    if (bv_manager.size(id) == 0) {
                        adj.clear();
                        continue;
                    }

                    epoch++;
                    if (epoch == 0) { seen.fill(0); epoch = 1; }

                    adj.clear();
                    for (const vertex_t u : bv_manager.boundary(id)) {
                        for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                            const partition_t v_id = p_manager[g.edges_v[i]];
                            if (id != v_id && seen[v_id] != epoch) {
                                seen[v_id] = epoch;
                                adj.push_back(Entry{v_id});
                            }
                        }
                    }
                }
            }
        }

        const std::vector<Entry> &neighbors(const partition_t x) const { return m_adj[x]; }
        size_t degree(const partition_t x) const { return m_adj[x].size(); }

        template<typename F>
        void for_each_neighbor(const partition_t x, F &&f) const {
            for (const auto &e : m_adj[x]) {
                f(e.target, weight_t(1));
            }
        }

        size_t edge_index(const partition_t u_id, const partition_t v_id) const {
            const partition_t min_part = std::min(u_id, v_id);
            const partition_t max_part = std::max(u_id, v_id);
            return static_cast<size_t>(min_part) * static_cast<size_t>(m_k) + static_cast<size_t>(max_part);
        }
    };

    template<bool LARGE_K>
    class SimpleLabelPropagationRefinement {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
        u64 m_threads = 1;

        std::vector<RandomEngine> rnd_engines;
        SimpleLabelPropagationConfiguration config = SimpleLabelPropagationConfiguration("simple_label_propagation");

        ActiveBlockScheduling active_block_scheduling;
        Distance1Matching<LARGE_K> d1_matcher;
        BlockAdjacency block_adjacency;

        static inline __attribute__((always_inline)) void apply_move(graph_t &g,
                                                                     bv_manager_t &bv_manager,
                                                                     p_manager_t &p_manager,
                                                                     const vertex_t u,
                                                                     const weight_t u_weight,
                                                                     const partition_t from_id,
                                                                     const partition_t to_id) {
            bv_manager.move(g, p_manager, u, from_id, to_id);
            p_manager.move_serial(u, u_weight, from_id, to_id);
        }

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

    public:
        SimpleLabelPropagationRefinement() = default;

        ~SimpleLabelPropagationRefinement() = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const u64 t_threads,
                        const u64 t_seed,
                        const SimpleLabelPropagationConfiguration &t_config) {
            m_n = t_n;
            m_m = t_m;
            m_k = t_k;
            m_threads = t_threads;

            config = t_config;

            rnd_engines.resize(m_threads);
            for (u64 t = 0; t < m_threads; ++t) {
                rnd_engines[t] = RandomEngine(t_seed + t);
            }

            active_block_scheduling.initialize(m_k);
            d1_matcher.initialize(m_k);
            block_adjacency.initialize(m_k, m_threads);
        }

        template<typename DistanceOracleT, typename QGraphT>
        void refine(graph_t &g,
                    DistanceOracleT &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    QGraphT &q_graph,
                    block_conn_t &block_conn) {
            if (g.uniform_v_weights && g.uniform_e_weights) {
                refine_impl<true, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn);
            } else if (g.uniform_v_weights) {
                refine_impl<true, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn);
            } else if (g.uniform_e_weights) {
                refine_impl<false, true>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn);
            } else {
                refine_impl<false, false>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn);
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights, typename DistanceOracleT, typename QGraphT>
        void refine_impl(graph_t &g,
                         DistanceOracleT &d_oracle,
                         bv_manager_t &bv_manager,
                         p_manager_t &p_manager,
                         QGraphT &q_graph,
                         block_conn_t &block_conn) {
            if (config.force_parallel_alg || m_threads > 1) {
                refine_impl_parallel<t_uniform_v_weights, t_uniform_e_weights>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn);
            } else {
                refine_impl_serial<t_uniform_v_weights, t_uniform_e_weights>(g, d_oracle, bv_manager, p_manager, q_graph, block_conn);
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights, typename DistanceOracleT, typename QGraphT>
        void refine_impl_parallel(graph_t &g,
                                  DistanceOracleT &d_oracle,
                                  bv_manager_t &bv_manager,
                                  p_manager_t &p_manager,
                                  QGraphT &q_graph,
                                  block_conn_t &block_conn) {
            active_block_scheduling.reset(m_k);

            std::vector<std::vector<std::pair<partition_t, partition_t> > > static_matchings;
            std::vector<std::pair<partition_t, partition_t> > matching;

            std::vector<std::vector<vertex_t>> thread_local_boundary(m_threads);

            for (u64 iteration = 0; iteration < config.max_iteration; ++iteration) {
                HEIPROMAP_PROFILE_SCOPE("refinement", "SimpleLPRefinement", "recompute_adjacency");
                block_adjacency.compute(g, p_manager, bv_manager, m_threads);

                d1_matcher.reset_used_edges();
                bool positive_move_occurred = false;

                if (config.use_static_matchings) {
                    HEIPROMAP_PROFILE_SCOPE("refinement", "SimpleLPRefinement", "matching");
                    d1_matcher.compute_static_matchings(block_adjacency, p_manager, active_block_scheduling.active_this_round, static_matchings, config.max_matching_rounds, config.min_matching_threshold);

                    if (static_matchings.empty()) {
                        break;
                    }

                    for (size_t i = 0; i < static_matchings.size(); ++i) {
                        matching.assign(static_matchings[i].begin(), static_matchings[i].end());

                        if (matching.empty()) {
                            continue;
                        }

                        HEIPROMAP_PROFILE_SCOPE("refinement", "SimpleLPRefinement", "process");

                        #pragma omp parallel for num_threads(m_threads) schedule(dynamic)
                        for (size_t j = 0; j < matching.size(); ++j) {
                            const partition_t A = matching[j].first;
                            const partition_t B = matching[j].second;

                            const u64 tid = omp_get_thread_num();
                            RandomEngine &rng = rnd_engines[tid];

                            const size_t size_A = bv_manager.size(A);
                            const size_t size_B = bv_manager.size(B);
                            auto &local_boundary = thread_local_boundary[tid];
                            local_boundary.clear();
                            local_boundary.insert(local_boundary.end(), bv_manager.boundary(A).begin(), bv_manager.boundary(A).end());
                            local_boundary.insert(local_boundary.end(), bv_manager.boundary(B).begin(), bv_manager.boundary(B).end());

                            const bool last_level_pair = config.use_edge_cut && d_oracle.last_level_pair(A, B);
                            bool moved_in_pair = false;

                            for (const vertex_t u: local_boundary) {
                                const partition_t u_id = p_manager[u];
                                if (u_id != A && u_id != B) { continue; }

                                const partition_t target_id = (u_id == A) ? B : A;
                                const weight_t u_weight = t_uniform_v_weights ? 1 : g.v_weights[u];

                                if (p_manager.get_bweight(target_id) + u_weight > p_manager.lmax[target_id]) { continue; }

                                weight_t qap_delta;
                                if (last_level_pair) {
                                    qap_delta = compute_edge_cut_delta<t_uniform_e_weights>(g, u, u_id, target_id, p_manager);
                                } else {
                                    qap_delta = compute_qap_delta<t_uniform_e_weights>(g, u, u_id, target_id, p_manager, d_oracle);
                                }

                                if (qap_delta > 0 || (qap_delta == 0 && rng.get_f32() < 0.5f)) {
                                    apply_move(g, bv_manager, p_manager, u, u_weight, u_id, target_id);
                                    if (qap_delta > 0) {
                                        moved_in_pair = true;
                                    }
                                }
                            }

                            if (moved_in_pair) {
                                if (config.use_active_scheduling) {
                                    active_block_scheduling.activate(A, B);
                                }
                                #pragma omp atomic write
                                positive_move_occurred = true;
                            }
                        }
                    }
                } else {
                    std::vector<std::pair<partition_t, partition_t> > dyn_matching;

                    HEIPROMAP_PROFILE_SCOPE("refinement", "SimpleLPRefinement", "matching");
                    d1_matcher.find_matching(q_graph, p_manager, active_block_scheduling.active_this_round, dyn_matching);

                    if (dyn_matching.empty()) {
                        break;
                    }

                    bool found_matching = true;
                    while (found_matching) {
                        HEIPROMAP_PROFILE_SCOPE("refinement", "SimpleLPRefinement", "process");

                        #pragma omp parallel for num_threads(m_threads) schedule(dynamic)
                        for (size_t i = 0; i < dyn_matching.size(); ++i) {
                            const partition_t A = dyn_matching[i].first;
                            const partition_t B = dyn_matching[i].second;

                            const u64 tid = omp_get_thread_num();
                            RandomEngine &rng = rnd_engines[tid];

                            const size_t size_A = bv_manager.size(A);
                            const size_t size_B = bv_manager.size(B);
                            auto &local_boundary = thread_local_boundary[tid];
                            local_boundary.clear();
                            local_boundary.insert(local_boundary.end(), bv_manager.boundary(A).begin(), bv_manager.boundary(A).end());
                            local_boundary.insert(local_boundary.end(), bv_manager.boundary(B).begin(), bv_manager.boundary(B).end());

                            const bool last_level_pair = config.use_edge_cut && d_oracle.last_level_pair(A, B);
                            bool moved_in_pair = false;

                            for (const vertex_t u: local_boundary) {
                                const partition_t u_id = p_manager[u];
                                if (u_id != A && u_id != B) { continue; }

                                const partition_t target_id = (u_id == A) ? B : A;
                                const weight_t u_weight = t_uniform_v_weights ? 1 : g.v_weights[u];

                                if (p_manager.get_bweight(target_id) + u_weight > p_manager.lmax[target_id]) { continue; }

                                weight_t qap_delta;
                                if (last_level_pair) {
                                    qap_delta = compute_edge_cut_delta<t_uniform_e_weights>(g, u, u_id, target_id, p_manager);
                                } else {
                                    qap_delta = compute_qap_delta<t_uniform_e_weights>(g, u, u_id, target_id, p_manager, d_oracle);
                                }

                                if (qap_delta > 0 || (qap_delta == 0 && rng.get_f32() < 0.5f)) {
                                    apply_move(g, bv_manager, p_manager, u, u_weight, u_id, target_id);
                                    if (qap_delta > 0) {
                                        moved_in_pair = true;
                                    }
                                }
                            }

                            if (moved_in_pair) {
                                if (config.use_active_scheduling) {
                                    active_block_scheduling.activate(A, B);
                                }
                                #pragma omp atomic write
                                positive_move_occurred = true;
                            }
                        }

                        if constexpr (LARGE_K) {
                            if (dyn_matching.size() < config.min_matching_threshold) {
                                break;
                            }
                        }

                        HEIPROMAP_PROFILE_SCOPE("refinement", "SimpleLPRefinement", "matching");
                        found_matching = d1_matcher.find_matching(block_adjacency, p_manager, active_block_scheduling.active_this_round, dyn_matching);
                    }
                }

                if (!positive_move_occurred) {
                    break;
                }

                if (config.use_active_scheduling) {
                    active_block_scheduling.next_round();
                }
            }

            HEIPROMAP_PROFILE_SCOPE("refinement", "SimpleLPRefinement", "recompute_final");
            if (config.enable_q_graph) {
                q_graph.compute_from_scratch(g, p_manager, bv_manager, m_threads);
            }
            if (config.enable_block_conn) {
                block_conn.compute_from_scratch(g, p_manager, m_threads);
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights, typename DistanceOracleT, typename QGraphT>
        void refine_impl_serial(graph_t &g,
                                DistanceOracleT &d_oracle,
                                bv_manager_t &bv_manager,
                                p_manager_t &p_manager,
                                QGraphT &q_graph,
                                block_conn_t &block_conn) {
            RandomEngine &random_engine = rnd_engines[0];

            bool positive_move_occurred = true;
            for (u64 iteration = 0; iteration < config.max_iteration && positive_move_occurred; ++iteration) {
                positive_move_occurred = false;

                HEIPROMAP_PROFILE_SCOPE("refinement", "SimpleLPRefinement", "process_vertices");
                for (vertex_t u = 0; u < g.n; ++u) {
                    if (!bv_manager.is_boundary(u)) { continue; }

                    const weight_t u_weight = t_uniform_v_weights ? 1 : g.v_weights[u];
                    const partition_t u_id = p_manager[u];

                    partition_t best_id = NO_ID;
                    weight_t best_qap_delta = -std::numeric_limits<weight_t>::max();
                    f32 counter = 0.0f;

                    // Find neighboring blocks on the fly
                    partition_t neighbor_blocks[256];
                    size_t n_neighbor_blocks = 0;

                    for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                        const partition_t v_id = p_manager[g.edges_v[i]];
                        if (v_id == u_id) { continue; }

                        bool found = false;
                        for (size_t j = 0; j < n_neighbor_blocks; ++j) {
                            if (neighbor_blocks[j] == v_id) {
                                found = true;
                                break;
                            }
                        }
                        if (!found && n_neighbor_blocks < 256) {
                            neighbor_blocks[n_neighbor_blocks++] = v_id;
                        }
                    }

                    for (size_t bi = 0; bi < n_neighbor_blocks; ++bi) {
                        const partition_t id = neighbor_blocks[bi];
                        const weight_t v_id_weight = p_manager.get_bweight(id);

                        if (v_id_weight + u_weight > p_manager.lmax[id]) { continue; }

                        const weight_t qap_delta = compute_qap_delta<t_uniform_e_weights>(g, u, u_id, id, p_manager, d_oracle);
                        if (qap_delta > best_qap_delta) {
                            best_id = id;
                            best_qap_delta = qap_delta;
                            counter = 1.0f;
                        } else if (qap_delta == best_qap_delta) {
                            counter += 1.0f;
                            if (random_engine.get_f32() < 1.0f / counter) {
                                best_id = id;
                            }
                        }
                    }

                    if (best_qap_delta > 0 || (best_qap_delta == 0 && random_engine.get_f32() < 0.5f)) {
                        apply_move(g, bv_manager, p_manager, u, u_weight, u_id, best_id);
                        positive_move_occurred |= (best_qap_delta > 0);
                    }
                }
            }

            HEIPROMAP_PROFILE_SCOPE("refinement", "SimpleLPRefinement", "recompute_final");
            if (config.enable_q_graph) {
                q_graph.compute_from_scratch(g, p_manager, bv_manager, m_threads);
            }
            if (config.enable_block_conn) {
                block_conn.compute_from_scratch(g, p_manager, m_threads);
            }
        }
    };
}

#endif //HEIPROMAP_SIMPLE_LABEL_PROPAGATION_REFINEMENT_H
