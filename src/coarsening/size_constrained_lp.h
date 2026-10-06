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

#ifndef HEIPROMAP_SIZE_CONSTRAINED_LP_H
#define HEIPROMAP_SIZE_CONSTRAINED_LP_H

#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

#include <omp.h>

#include "../definitions.h"
#include "../datastructures/csr_graph.h"
#include "../distance_oracles/matrix_distance_oracle.h"
#include "../datastructures/partition_manager.h"
#include "../utility/aligned_array.h"
#include "../utility/mapping.h"
#include "../utility/profiler.h"
#include "../utility/random_engine.h"
#include "../utility/small_map.h"
#include "../utility/utils.h"
#include "edge_rating.h"

namespace HeiProMap {

    // =========================================================================
    // Configuration
    // =========================================================================

    class SizeConstrainedLPConfiguration {
    public:
        u64 max_rounds = 5;
        f64 min_threshold = 0.05;
        f64 f = 64;
        EdgeRatingFunction rating_function = EdgeRatingFunction::WEIGHT;
        bool use_degree_ordering = true;
        bool force_parallel_alg = false;
    };

    // =========================================================================
    // Size Constrained Label Propagation
    // =========================================================================

    class SizeConstrainedLP {
        // graph dimensions
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;

        // vertex ordering
        AlignedArray<vertex_t> flat_vertices;
        AlignedArray<vertex_t> bucket_sizes;
        AlignedArray<vertex_t> bucket_offsets;
        std::vector<std::vector<vertex_t>> thread_bucket_counts;
        std::vector<std::vector<vertex_t>> thread_bucket_offsets_local;

        // cluster state
        AlignedArray<weight_t> cluster_weights;
        AlignedArray<vertex_t> cluster_count;

        // active vertex tracking
        AlignedArray<u8> active;
        AlignedArray<u8> active_next;

        // remapping
        AlignedArray<vertex_t> remap;
        AlignedArray<vertex_t> singletons;

        // config
        SizeConstrainedLPConfiguration config;
        RandomEngine random_engine;

    public:
        size_t heap_bytes() const {
            return flat_vertices.heap_bytes() + bucket_sizes.heap_bytes()
                 + bucket_offsets.heap_bytes() + cluster_weights.heap_bytes()
                 + cluster_count.heap_bytes() + active.heap_bytes()
                 + active_next.heap_bytes() + remap.heap_bytes()
                 + singletons.heap_bytes();
        }

        void release_memory() {
            flat_vertices.free_memory();
            bucket_sizes.free_memory();
            bucket_offsets.free_memory();
            cluster_weights.free_memory();
            cluster_count.free_memory();
            active.free_memory();
            active_next.free_memory();
            remap.free_memory();
            singletons.free_memory();
        }

        void set_pool(MemoryPool *pool) {
            flat_vertices.set_pool(pool);
            bucket_sizes.set_pool(pool);
            bucket_offsets.set_pool(pool);
            cluster_weights.set_pool(pool);
            cluster_count.set_pool(pool);
            active.set_pool(pool);
            active_next.set_pool(pool);
            remap.set_pool(pool);
            singletons.set_pool(pool);
        }

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const u64 t_seed,
                        const SizeConstrainedLPConfiguration &t_config,
                        MemoryPool *pool = nullptr) {
            set_pool(pool);
            HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "initialize");

            m_n = t_n;
            m_m = t_m;
            m_k = t_k;

            config = t_config;
            random_engine = RandomEngine(t_seed);
        }

        // =====================================================================
        // Public interface
        // =====================================================================

        void cluster(const size_t level,
                     const graph_t &g,
                     const p_manager_t &p_manager,
                     Mapping &mapping,
                     const f64 imbalance,
                     const u64 threads,
                     const weight_t lmax) {
            dispatch_with_templates(level, g, p_manager, mapping, imbalance, threads, lmax);
            release_memory();
        }

        void cluster(const size_t level,
                     const graph_t &g,
                     const p_manager_t &p_manager,
                     Mapping &mapping,
                     const f64 imbalance,
                     const u64 threads) {
            const weight_t lmax = std::ceil((1.0 + imbalance) * (static_cast<f64>(g.g_weight) / static_cast<f64>(p_manager.k)));
            cluster(level, g, p_manager, mapping, imbalance, threads, lmax);
        }

    private:
        // =====================================================================
        // Template dispatch
        // =====================================================================

        void dispatch_with_templates(const size_t level,
                                     const graph_t &g,
                                     const p_manager_t &p_manager,
                                     Mapping &mapping,
                                     const f64 imbalance,
                                     const u64 threads,
                                     const weight_t lmax) {
            auto dispatch = [&]<EdgeRatingFunction R>() {
                if (g.uniform_v_weights && g.uniform_e_weights) {
                    cluster_templated<true, true, R>(level, g, p_manager, mapping, imbalance, threads, lmax);
                } else if (g.uniform_v_weights) {
                    cluster_templated<true, false, R>(level, g, p_manager, mapping, imbalance, threads, lmax);
                } else if (g.uniform_e_weights) {
                    cluster_templated<false, true, R>(level, g, p_manager, mapping, imbalance, threads, lmax);
                } else {
                    cluster_templated<false, false, R>(level, g, p_manager, mapping, imbalance, threads, lmax);
                }
            };

            switch (config.rating_function) {
                case EdgeRatingFunction::WEIGHT:          dispatch.template operator()<EdgeRatingFunction::WEIGHT>(); break;
                case EdgeRatingFunction::EXPANSION:       dispatch.template operator()<EdgeRatingFunction::EXPANSION>(); break;
                case EdgeRatingFunction::EXPANSIONSTAR:   dispatch.template operator()<EdgeRatingFunction::EXPANSIONSTAR>(); break;
                case EdgeRatingFunction::EXPANSIONSTARSTAR: dispatch.template operator()<EdgeRatingFunction::EXPANSIONSTARSTAR>(); break;
                case EdgeRatingFunction::INNEROUTER:      dispatch.template operator()<EdgeRatingFunction::INNEROUTER>(); break;
            }
        }

        // =====================================================================
        // Setup helpers
        // =====================================================================

        void init_vertex_order(const graph_t &g, const size_t B, const u64 threads) {
            HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "flat_vertices");
            flat_vertices.initialize(g.n);

            if (!config.use_degree_ordering) {
                #pragma omp parallel for num_threads(threads)
                for (vertex_t u = 0; u < g.n; ++u) {
                    flat_vertices[u] = u;
                }
                return;
            }

            bucket_sizes.initialize(B, 0);
            bucket_offsets.initialize(B);

            if (thread_bucket_counts.size() < threads) {
                thread_bucket_counts.resize(threads);
                thread_bucket_offsets_local.resize(threads);
            }

            // parallel bucket count
            #pragma omp parallel num_threads(threads)
            {
                const u64 tid = omp_get_thread_num();
                thread_bucket_counts[tid].assign(B, 0);
                auto &local_counts = thread_bucket_counts[tid];
                #pragma omp for schedule(static)
                for (vertex_t u = 0; u < g.n; ++u) {
                    const size_t d = g.deg(u);
                    const size_t b = (d == 0) ? 0 : floor_log2(d);
                    local_counts[b]++;
                }
            }
            for (u64 t = 0; t < threads; ++t) {
                for (size_t b = 0; b < B; ++b) {
                    bucket_sizes[b] += thread_bucket_counts[t][b];
                }
            }

            bucket_offsets[0] = 0;
            for (size_t i = 1; i < B; ++i) {
                bucket_offsets[i] = bucket_offsets[i - 1] + bucket_sizes[i - 1];
            }

            // per-thread scatter offsets
            for (size_t b = 0; b < B; ++b) {
                vertex_t off = bucket_offsets[b];
                for (u64 t = 0; t < threads; ++t) {
                    thread_bucket_offsets_local[t].resize(B);
                    thread_bucket_offsets_local[t][b] = off;
                    off += thread_bucket_counts[t][b];
                }
            }

            #pragma omp parallel num_threads(threads)
            {
                const u64 tid = omp_get_thread_num();
                auto &local_offsets = thread_bucket_offsets_local[tid];
                #pragma omp for schedule(static)
                for (vertex_t u = 0; u < g.n; ++u) {
                    const size_t d = g.deg(u);
                    const size_t b = (d == 0) ? 0 : floor_log2(d);
                    flat_vertices[local_offsets[b]++] = u;
                }
            }

            // recompute bucket_offsets for shuffle phase
            bucket_offsets[0] = 0;
            for (size_t i = 1; i < B; ++i) {
                bucket_offsets[i] = bucket_offsets[i - 1] + bucket_sizes[i - 1];
            }
        }

        void init_cluster_state(const graph_t &g, Mapping &mapping, const u64 threads) {
            HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "cluster_weights");
            cluster_weights.initialize(g.n);
            cluster_count.initialize(g.n);
            #pragma omp parallel for num_threads(threads)
            for (vertex_t u = 0; u < g.n; ++u) {
                mapping.set(u, u);
                cluster_weights[u] = g.v_weights[u];
                cluster_count[u] = 1;
            }
        }

        void init_active(const graph_t &g, const u64 threads) {
            HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "active");
            active.initialize(g.n);
            active_next.initialize(g.n);
            #pragma omp parallel for num_threads(threads)
            for (vertex_t u = 0; u < g.n; ++u) {
                active[u] = 1;
                active_next[u] = 1;
            }
        }

        // =====================================================================
        // Shuffle
        // =====================================================================

        void shuffle_buckets(const graph_t &g, const size_t B, const u64 threads) {
            HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "shuffle_buckets");
            u64 shuffle_seed = random_engine.get_u64();

            constexpr size_t MIN_CHUNK = 1024;
            std::vector<std::pair<size_t, size_t>> chunks;
            for (size_t i = 0; i < B - 1; ++i) {
                const size_t beg = bucket_offsets[i];
                const size_t end = bucket_offsets[i + 1];
                const size_t bucket_size = end - beg;
                if (bucket_size <= MIN_CHUNK) {
                    chunks.emplace_back(beg, end);
                } else {
                    const size_t n_chunks = std::min(threads, (bucket_size + MIN_CHUNK - 1) / MIN_CHUNK);
                    const size_t chunk_size = bucket_size / n_chunks;
                    for (size_t c = 0; c < n_chunks; ++c) {
                        const size_t c_beg = beg + c * chunk_size;
                        const size_t c_end = (c + 1 == n_chunks) ? end : c_beg + chunk_size;
                        chunks.emplace_back(c_beg, c_end);
                    }
                }
            }

            #pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
            for (size_t i = 0; i < chunks.size(); ++i) {
                std::mt19937 local_gen(shuffle_seed + i);
                fast_shuffle_unchecked(flat_vertices.get_ptr() + chunks[i].first,
                                      flat_vertices.get_ptr() + chunks[i].second, local_gen);
            }
        }

        // =====================================================================
        // Matching phases
        // =====================================================================

        template<bool t_uniform_v_weights, bool t_uniform_e_weights, EdgeRatingFunction t_rating_function>
        void merge_when_identity([[maybe_unused]] const size_t level,
                                 const graph_t &g,
                                 const p_manager_t &p_manager,
                                 Mapping &mapping,
                                 const weight_t max_w,
                                 const u64 threads) {
            HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "merge_when_identity");

            std::vector<u8> used(g.n, 0);

            #pragma omp parallel for num_threads(threads) schedule(static)
            for (vertex_t a = 0; a < g.n; ++a) {
                if (used[a] == 1) { continue; }
                const partition_t a_id = p_manager[a];
                const weight_t a_w = t_uniform_v_weights ? 1 : g.v_weights[a];

                vertex_t best_b = static_cast<vertex_t>(-1);
                weight_t best_ew = 0;

                for (size_t i = g.neighborhoods[a]; i < g.neighborhoods[a + 1]; ++i) {
                    const vertex_t b = g.edges_v[i];
                    if (used[b] == 1 || p_manager[b] != a_id) { continue; }

                    const weight_t b_w = t_uniform_v_weights ? 1 : g.v_weights[b];
                    if (a_w + b_w <= max_w) {
                        const weight_t ew = t_uniform_e_weights ? 1 : g.edges_w[i];
                        if (ew > best_ew || best_b == static_cast<vertex_t>(-1)) {
                            best_ew = ew;
                            best_b = b;
                        }
                    }
                }

                if (best_b != static_cast<vertex_t>(-1)) {
                    u8 expected_a = 0;
                    if (__atomic_compare_exchange_n(&used[a], &expected_a, 1, false, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE)) {
                        u8 expected_b = 0;
                        if (__atomic_compare_exchange_n(&used[best_b], &expected_b, 1, false, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE)) {
                            mapping.set(a, best_b);
                        } else {
                            __atomic_store_n(&used[a], 0, __ATOMIC_RELEASE);
                        }
                    }
                }
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights, EdgeRatingFunction t_rating_function>
        void merge_singletons([[maybe_unused]] const size_t level,
                              const graph_t &g,
                              const p_manager_t &p_manager,
                              Mapping &mapping,
                              const weight_t max_w,
                              const u64 threads) {
            HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "merge_singletons");

            // parallel singleton collection
            std::vector<std::vector<vertex_t>> thread_singletons(threads);
            #pragma omp parallel num_threads(threads)
            {
                const u64 tid = omp_get_thread_num();
                auto &local = thread_singletons[tid];
                local.clear();
                #pragma omp for schedule(static)
                for (vertex_t u = 0; u < g.n; ++u) {
                    if (cluster_count[mapping.get(u)] == 1) {
                        local.push_back(u);
                    }
                }
            }
            vertex_t singletons_size = 0;
            for (u64 t = 0; t < threads; ++t) singletons_size += thread_singletons[t].size();

            singletons.initialize(singletons_size);
            vertex_t offset = 0;
            for (u64 t = 0; t < threads; ++t) {
                std::memcpy(singletons.get_ptr() + offset, thread_singletons[t].data(), thread_singletons[t].size() * sizeof(vertex_t));
                offset += thread_singletons[t].size();
            }

            if (singletons_size == 0) { return; }

            // LP over singletons
            #pragma omp parallel num_threads(threads)
            {
                FlatMap<vertex_t, f32> flat_map;
                flat_map.reserve(128);

                #pragma omp for schedule(static)
                for (size_t i = 0; i < singletons_size; ++i) {
                    const vertex_t u = singletons[i];
                    const vertex_t cur_id = mapping.get(u);

                    if (cluster_count[cur_id] != 1) { continue; }

                    const partition_t u_id = p_manager[u];
                    const weight_t u_w = t_uniform_v_weights ? 1 : g.v_weights[u];
                    f32 current_id_w = 0;

                    flat_map.clear();
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const vertex_t v = g.edges_v[j];
                        if (u_id != p_manager[v]) { continue; }

                        const weight_t v_w = t_uniform_v_weights ? 1 : g.v_weights[v];
                        const weight_t ew = t_uniform_e_weights ? 1 : g.edges_w[j];
                        const f32 edge_rating = compute_edge_rating<t_uniform_v_weights, t_uniform_e_weights, t_rating_function>(g, u, v, u_w, v_w, ew);

                        const vertex_t id = mapping.get(v);
                        if (id == cur_id) {
                            current_id_w += edge_rating;
                        } else if (u_w + cluster_weights[id] <= max_w) {
                            flat_map.add(id, edge_rating);
                        }
                    }

                    vertex_t best_id = cur_id;
                    f32 best_weight = current_id_w;
                    for (const auto &[id, w]: flat_map) {
                        if (w > best_weight && u_w + cluster_weights[id] <= max_w) {
                            best_weight = w;
                            best_id = id;
                        }
                    }

                    if (best_id == cur_id) { continue; }

                    mapping.set(u, best_id);
                    #pragma omp atomic
                    cluster_weights[best_id] += u_w;
                    #pragma omp atomic
                    cluster_weights[cur_id] -= u_w;
                    #pragma omp atomic
                    cluster_count[best_id] += 1;
                    #pragma omp atomic
                    cluster_count[cur_id] -= 1;
                }
            }
        }

        // =====================================================================
        // Core LP iterations
        // =====================================================================

        template<bool t_uniform_v_weights, bool t_uniform_e_weights, EdgeRatingFunction t_rating_function>
        u64 cluster_threaded(const graph_t &g,
                             const p_manager_t &p_manager,
                             Mapping &mapping,
                             const weight_t max_w,
                             const u64 round,
                             const u64 threads) {
            HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "cluster_threaded");
            u64 n_moved = 0;

            #pragma omp parallel num_threads(threads)
            {
                FlatMap<vertex_t, f32> flat_map;
                flat_map.reserve(128);

                #pragma omp for schedule(static) reduction(+:n_moved)
                for (size_t i = 0; i < g.n; ++i) {
                    const vertex_t u = flat_vertices[i];
                    if (active[u] == 0) { continue; }

                    const weight_t u_w = t_uniform_v_weights ? 1 : g.v_weights[u];
                    const partition_t u_id = p_manager[u];
                    const vertex_t current_id = mapping.get(u);
                    f32 current_id_w = 0;

                    vertex_t best_id = current_id;
                    f32 best_weight = 0;

                    flat_map.clear();
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const vertex_t v = g.edges_v[j];
                        if (u_id != p_manager[v]) { continue; }

                        const weight_t v_w = t_uniform_v_weights ? 1 : g.v_weights[v];
                        const weight_t ew = t_uniform_e_weights ? 1 : g.edges_w[j];
                        const f32 edge_rating = compute_edge_rating<t_uniform_v_weights, t_uniform_e_weights, t_rating_function>(g, u, v, u_w, v_w, ew);

                        const vertex_t id = mapping.get(v);
                        if (id == current_id) {
                            current_id_w += edge_rating;
                        } else if (u_w + cluster_weights[id] <= max_w) {
                            const f32 new_w = flat_map.add_and_ret(id, edge_rating);
                            if (new_w > best_weight) {
                                best_weight = new_w;
                                best_id = id;
                            }
                        }
                    }

                    if (current_id_w >= best_weight) {
                        best_id = current_id;
                    }

                    if (best_id != current_id) {
                        mapping.set(u, best_id);
                        #pragma omp atomic
                        cluster_weights[best_id] += u_w;
                        #pragma omp atomic
                        cluster_weights[current_id] -= u_w;

                        n_moved += 1;
                        if (round > 0) {
                            for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                                const vertex_t v = g.edges_v[j];
                                #pragma omp atomic write
                                active_next[v] = 1;
                            }
                        }
                    }
                }
            }
            return n_moved;
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights, EdgeRatingFunction t_rating_function>
        u64 cluster_serial(const graph_t &g,
                           const p_manager_t &p_manager,
                           Mapping &mapping,
                           const weight_t max_w,
                           const u64 round) {
            HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "cluster_serial");
            u64 n_moved = 0;
            FlatMap<vertex_t, f32> flat_map;
            flat_map.reserve(128);

            for (size_t i = 0; i < g.n; ++i) {
                const vertex_t u = flat_vertices[i];
                if (active[u] == 0) { continue; }

                const weight_t u_w = t_uniform_v_weights ? 1 : g.v_weights[u];
                const partition_t u_id = p_manager[u];
                const vertex_t current_id = mapping.get(u);
                f32 current_id_w = 0;

                vertex_t best_id = current_id;
                f32 best_weight = 0;

                flat_map.clear();
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const vertex_t v = g.edges_v[j];
                    if (u_id != p_manager[v]) { continue; }

                    const weight_t v_w = t_uniform_v_weights ? 1 : g.v_weights[v];
                    const weight_t ew = t_uniform_e_weights ? 1 : g.edges_w[j];
                    const f32 edge_rating = compute_edge_rating<t_uniform_v_weights, t_uniform_e_weights, t_rating_function>(g, u, v, u_w, v_w, ew);

                    const vertex_t id = mapping.get(v);
                    if (id == current_id) {
                        current_id_w += edge_rating;
                    } else if (u_w + cluster_weights[id] <= max_w) {
                        const f32 new_w = flat_map.add_and_ret(id, edge_rating);
                        if (new_w > best_weight) {
                            best_weight = new_w;
                            best_id = id;
                        }
                    }
                }

                if (current_id_w >= best_weight) {
                    best_weight = current_id_w;
                    best_id = current_id;
                }

                if (best_id != current_id) {
                    mapping.set(u, best_id);
                    cluster_weights[best_id] += u_w;
                    cluster_weights[current_id] -= u_w;

                    n_moved += 1;
                    if (round > 0) {
                        for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                            const vertex_t v = g.edges_v[j];
                            active_next[v] = 1;
                        }
                    }
                }
            }
            return n_moved;
        }

        // =====================================================================
        // Remapping
        // =====================================================================

        void remap_clusters(const graph_t &g, Mapping &mapping, const u64 threads) {
            HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "calc_map");
            remap.initialize(g.n, 0);

            #pragma omp parallel for num_threads(threads)
            for (vertex_t u = 0; u < g.n; ++u) {
                remap[mapping.get(u)] = 1;
            }

            // parallel prefix sum
            std::vector<vertex_t> thread_counts(threads, 0);
            #pragma omp parallel num_threads(threads)
            {
                const u64 tid = omp_get_thread_num();
                const vertex_t chunk = (g.n + threads - 1) / threads;
                const vertex_t start = std::min(g.n, static_cast<vertex_t>(tid * chunk));
                const vertex_t end = std::min(g.n, static_cast<vertex_t>(start + chunk));
                vertex_t local_count = 0;
                for (vertex_t u = start; u < end; ++u) {
                    if (remap[u] == 1) local_count++;
                }
                thread_counts[tid] = local_count;

                #pragma omp barrier

                #pragma omp single
                {
                    vertex_t sum = 0;
                    for (u64 i = 0; i < threads; ++i) {
                        vertex_t c = thread_counts[i];
                        thread_counts[i] = sum;
                        sum += c;
                    }
                    mapping.set_coarse_n(sum);
                }

                vertex_t local_id = thread_counts[tid];
                for (vertex_t u = start; u < end; ++u) {
                    if (remap[u] == 1) {
                        remap[u] = local_id++;
                    } else {
                        remap[u] = m_n;
                    }
                }
            }

            #pragma omp parallel for num_threads(threads)
            for (vertex_t u = 0; u < g.n; ++u) {
                mapping.set(u, remap[mapping.get(u)]);
            }
        }

        // =====================================================================
        // Main algorithm
        // =====================================================================

        template<bool t_uniform_v_weights, bool t_uniform_e_weights, EdgeRatingFunction t_rating_function>
        void cluster_templated([[maybe_unused]] const size_t level,
                               const graph_t &g,
                               const p_manager_t &p_manager,
                               Mapping &mapping,
                               [[maybe_unused]] const f64 imbalance,
                               const u64 threads,
                               const weight_t lmax) {
            mapping.initialize(g.n);

            // compute max vertex weight and degree
            weight_t max_v_w = 0;
            vertex_t max_deg = 0;
            HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "max");

            #pragma omp parallel for num_threads(threads) reduction(max:max_v_w,max_deg)
            for (vertex_t u = 0; u < g.n; ++u) {
                max_v_w = std::max(max_v_w, g.v_weights[u]);
                max_deg = std::max(max_deg, g.deg(u));
            }

            const weight_t W = std::ceil(static_cast<f64>(lmax) / config.f);
            const weight_t max_w = std::max(max_v_w, W);
            const size_t B = (max_deg == 0) ? 1 : (floor_log2(max_deg) + 1);

            // setup
            init_vertex_order(g, B, threads);
            init_cluster_state(g, mapping, threads);
            init_active(g, threads);

            // LP rounds
            for (u64 round = 0; round < config.max_rounds; ++round) {
                if (config.use_degree_ordering) {
                    shuffle_buckets(g, B, threads);
                }

                u64 n_moved;
                if (config.force_parallel_alg || threads > 1) {
                    n_moved = cluster_threaded<t_uniform_v_weights, t_uniform_e_weights, t_rating_function>(g, p_manager, mapping, max_w, round, threads);
                } else {
                    n_moved = cluster_serial<t_uniform_v_weights, t_uniform_e_weights, t_rating_function>(g, p_manager, mapping, max_w, round);
                }

                if (static_cast<f64>(n_moved) < static_cast<f64>(g.n) * config.min_threshold) {
                    break;
                }

                HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "swap_active");
                std::swap(active, active_next);
                #pragma omp parallel for num_threads(threads)
                for (vertex_t u = 0; u < g.n; ++u) {
                    active_next[u] = 0;
                }
            }

            // recompute cluster_count
            #pragma omp parallel for num_threads(threads)
            for (vertex_t u = 0; u < g.n; ++u) {
                cluster_count[u] = 0;
            }
            #pragma omp parallel for num_threads(threads)
            for (vertex_t u = 0; u < g.n; ++u) {
                #pragma omp atomic
                cluster_count[mapping.get(u)] += 1;
            }

            merge_singletons<t_uniform_v_weights, t_uniform_e_weights, t_rating_function>(level, g, p_manager, mapping, max_w, threads);

            // check for identity mapping
            bool ident_mapping = true;
            HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "is_identity_mapping");
            for (vertex_t u = 0; u < g.n; ++u) {
                if (u != mapping.get(u)) {
                    ident_mapping = false;
                    break;
                }
            }

            if (ident_mapping) {
                merge_when_identity<t_uniform_v_weights, t_uniform_e_weights, t_rating_function>(level, g, p_manager, mapping, max_w, threads);
            }

            remap_clusters(g, mapping, threads);
        }
    };
}

#endif //HEIPROMAP_SIZE_CONSTRAINED_LP_H
