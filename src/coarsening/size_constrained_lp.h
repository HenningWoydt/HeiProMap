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
#include "../distance_oracles/distance_oracle.h"
#include "../datastructures/partition_manager.h"
#include "../utility/aligned_array.h"
#include "../utility/mapping.h"
#include "../utility/profiler.h"
#include "../utility/random_engine.h"
#include "../utility/small_map.h"
#include "../utility/utils.h"

namespace HeiProMap {
    class SizeConstrainedLPConfiguration {
    public:
        u64 max_rounds = 5;
        f64 min_threshold = 0.05;
        f64 f = 64;
        EdgeRatingFunction rating_function = EdgeRatingFunction::WEIGHT;
        bool use_degree_ordering = true;
        bool force_parallel_alg = false;
    };

    class SizeConstrainedLP {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;

        AlignedArray<vertex_t> flat_vertices;
        AlignedArray<vertex_t> bucket_sizes;
        AlignedArray<vertex_t> bucket_offsets;
        AlignedArray<weight_t> cluster_weights;
        AlignedArray<vertex_t> cluster_count;
        AlignedArray<u8> active;
        AlignedArray<u8> active_next;
        AlignedArray<vertex_t> remap;
        AlignedArray<vertex_t> singletons;

        SizeConstrainedLPConfiguration config;
        RandomEngine random_engine;

    public:
        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const u64 t_seed,
                        const SizeConstrainedLPConfiguration &t_config) {
            HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "initialize");

            m_n = t_n;
            m_m = t_m;
            m_k = t_k;

            config = t_config;
            random_engine = RandomEngine(t_seed);
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights, EdgeRatingFunction t_rating_function>
        [[nodiscard]] static inline __attribute__((always_inline)) f32 compute_edge_rating(const graph_t &g,
                                                                                           const vertex_t u,
                                                                                           const vertex_t v,
                                                                                           const weight_t u_w,
                                                                                           const weight_t v_w,
                                                                                           const weight_t ew) noexcept {
            if constexpr (t_rating_function == EdgeRatingFunction::WEIGHT) {
                return static_cast<f32>(ew);
            } else if constexpr (t_rating_function == EdgeRatingFunction::EXPANSION) {
                return static_cast<f32>(ew) / static_cast<f32>(u_w + v_w);
            } else if constexpr (t_rating_function == EdgeRatingFunction::EXPANSIONSTAR) {
                return static_cast<f32>(ew) / static_cast<f32>(u_w * v_w);
            } else if constexpr (t_rating_function == EdgeRatingFunction::EXPANSIONSTARSTAR) {
                return static_cast<f32>(ew * ew) / static_cast<f32>(u_w * v_w);
            } else if constexpr (t_rating_function == EdgeRatingFunction::INNEROUTER) {
                weight_t out_v = 0;
                for (u64 i_edge = g.neighborhoods[v]; i_edge < g.neighborhoods[v + 1]; ++i_edge) {
                    out_v += g.edges_w[i_edge];
                }
                weight_t out_u = 0;
                for (u64 i_edge = g.neighborhoods[u]; i_edge < g.neighborhoods[u + 1]; ++i_edge) {
                    out_u += g.edges_w[i_edge];
                }
                return static_cast<f32>(ew) / static_cast<f32>(out_v + out_u - (2 * ew));
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights, EdgeRatingFunction t_rating_function>
        void merge_when_identity([[maybe_unused]] const size_t level,
                                 const graph_t &g,
                                 const p_manager_t &p_manager,
                                 Mapping &mapping,
                                 const weight_t max_w) {
            HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "merge_when_identity");

            std::vector<u8> used(g.n, 0);
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
                    used[a] = 1;
                    used[best_b] = 1;
                    mapping.set(a, best_b);
                }
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights, EdgeRatingFunction t_rating_function>
        void merge_singletons([[maybe_unused]] const size_t level,
                              const graph_t &g,
                              const p_manager_t &p_manager,
                              Mapping &mapping,
                              const weight_t max_w) {
            HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "merge_singletons");

            vertex_t singletons_size = 0;
            singletons.initialize(g.n);

            for (vertex_t u = 0; u < g.n; ++u) {
                const vertex_t id = mapping.get(u);
                if (cluster_count[id] == 1) {
                    singletons[singletons_size++] = u;
                }
            }

            if (singletons_size == 0) { return; }

            FlatMap<vertex_t, f32> flat_map;
            flat_map.reserve(128);

            for (size_t i = 0; i < singletons_size; ++i) {
                const vertex_t u = singletons[i];
                const vertex_t cur_id = mapping.get(u);

                // Might not be a singleton anymore if merged into earlier in this loop
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

                // Apply merge: u moves from cur_id -> best_id
                mapping.set(u, best_id);

                cluster_weights[best_id] += u_w;
                cluster_count[best_id] += 1;

                cluster_weights[cur_id] -= u_w;
                cluster_count[cur_id] -= 1; // becomes 0
            }
        }

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

                #pragma omp for schedule(dynamic) reduction(+:n_moved)
                for (size_t i = 0; i < g.n; ++i) {
                    const vertex_t u = flat_vertices[i];
                    if (active[u] == 0) { continue; }

                    const weight_t u_w = t_uniform_v_weights ? 1 : g.v_weights[u];
                    const partition_t u_id = p_manager[u];
                    const vertex_t current_id = mapping.get(u);
                    f32 current_id_w = 0;

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
                            flat_map.add(id, edge_rating);
                        }
                    }

                    vertex_t best_id = current_id;
                    f32 best_weight = current_id_w;
                    for (const auto &[id, w]: flat_map) {
                        if (w > best_weight) {
                            best_weight = w;
                            best_id = id;
                        }
                    }

                    if (best_id != current_id) {
                        mapping.set(u, best_id);
                        #pragma omp atomic
                        cluster_weights[best_id] += u_w;
                        #pragma omp atomic
                        cluster_weights[current_id] -= u_w;
                        #pragma omp atomic
                        cluster_count[best_id] += 1;
                        #pragma omp atomic
                        cluster_count[current_id] -= 1;

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
                    cluster_count[best_id] += 1;
                    cluster_count[current_id] -= 1;

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

        void cluster(const size_t level,
                     const graph_t &g,
                     const p_manager_t &p_manager,
                     Mapping &mapping,
                     const f64 imbalance,
                     const u64 threads,
                     const weight_t lmax) {
            auto dispatch_with_rating = [&](auto rating_func_const) {
                constexpr EdgeRatingFunction rating_func = rating_func_const;
                if (g.uniform_v_weights && g.uniform_e_weights) {
                    cluster_templated<true, true, rating_func>(level, g, p_manager, mapping, imbalance, threads, lmax);
                } else if (g.uniform_v_weights) {
                    cluster_templated<true, false, rating_func>(level, g, p_manager, mapping, imbalance, threads, lmax);
                } else if (g.uniform_e_weights) {
                    cluster_templated<false, true, rating_func>(level, g, p_manager, mapping, imbalance, threads, lmax);
                } else {
                    cluster_templated<false, false, rating_func>(level, g, p_manager, mapping, imbalance, threads, lmax);
                }
            };

            switch (config.rating_function) {
                case EdgeRatingFunction::WEIGHT:
                    dispatch_with_rating(std::integral_constant<EdgeRatingFunction, EdgeRatingFunction::WEIGHT>{});
                    break;
                case EdgeRatingFunction::EXPANSION:
                    dispatch_with_rating(std::integral_constant<EdgeRatingFunction, EdgeRatingFunction::EXPANSION>{});
                    break;
                case EdgeRatingFunction::EXPANSIONSTAR:
                    dispatch_with_rating(std::integral_constant<EdgeRatingFunction, EdgeRatingFunction::EXPANSIONSTAR>{});
                    break;
                case EdgeRatingFunction::EXPANSIONSTARSTAR:
                    dispatch_with_rating(std::integral_constant<EdgeRatingFunction, EdgeRatingFunction::EXPANSIONSTARSTAR>{});
                    break;
                case EdgeRatingFunction::INNEROUTER:
                    dispatch_with_rating(std::integral_constant<EdgeRatingFunction, EdgeRatingFunction::INNEROUTER>{});
                    break;
            }
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

        template<bool t_uniform_v_weights, bool t_uniform_e_weights, EdgeRatingFunction t_rating_function>
        void cluster_templated([[maybe_unused]] const size_t level,
                               const graph_t &g,
                               const p_manager_t &p_manager,
                               Mapping &mapping,
                               [[maybe_unused]] const f64 imbalance,
                               const u64 threads,
                               const weight_t lmax) {
            mapping.initialize(g.n);

            // Determine maximum vertex weight and degree
            weight_t max_v_w = 0;
            vertex_t max_deg = 0;
            HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "max");
            for (vertex_t u = 0; u < g.n; ++u) {
                max_v_w = std::max(max_v_w, g.v_weights[u]);
                max_deg = std::max(max_deg, g.deg(u));
            }

            const weight_t W = std::ceil(static_cast<f64>(lmax) / config.f);
            const weight_t max_w = std::max(max_v_w, W);
            const size_t B = (max_deg == 0) ? 1 : (floor_log2(max_deg) + 1);

            // Setup vertex visit order
            HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "flat_vertices");
            flat_vertices.initialize(g.n);
            if (config.use_degree_ordering) {
                bucket_sizes.initialize(B, 0);
                bucket_offsets.initialize(B);

                for (vertex_t u = 0; u < g.n; ++u) {
                    const size_t d = g.deg(u);
                    const size_t b = (d == 0) ? 0 : floor_log2(d);
                    bucket_sizes[b]++;
                }

                bucket_offsets[0] = 0;
                for (size_t i = 1; i < B; ++i) {
                    bucket_offsets[i] = bucket_offsets[i - 1] + bucket_sizes[i - 1];
                }

                for (vertex_t u = 0; u < g.n; ++u) {
                    const size_t d = g.deg(u);
                    const size_t b = (d == 0) ? 0 : floor_log2(d);
                    flat_vertices[bucket_offsets[b]++] = u;
                }
            } else {
                for (vertex_t u = 0; u < g.n; ++u) {
                    flat_vertices[u] = u;
                }
            }

            // Setup cluster weights and counts (each vertex starts in its own cluster)
            HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "cluster_weights");
            cluster_weights.initialize(g.n);
            cluster_count.initialize(g.n);
            for (vertex_t u = 0; u < g.n; ++u) {
                mapping.set(u, u);
                cluster_weights[u] = g.v_weights[u];
                cluster_count[u] = 1;
            }

            // Setup active tracking
            HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "active");
            active.initialize(g.n, 1);
            active_next.initialize(g.n, 1);

            for (u64 round = 0; round < config.max_rounds; ++round) {
                u64 n_moved = 0;

                if (config.use_degree_ordering) {
                    HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "shuffle_buckets");
                    for (size_t i = 0; i < B - 1; ++i) {
                        const size_t beg = bucket_offsets[i];
                        const size_t end = bucket_offsets[i + 1];
                        fast_shuffle_unchecked(flat_vertices.get_ptr() + beg, flat_vertices.get_ptr() + end, random_engine.generator);
                    }
                }

                if (config.force_parallel_alg || threads > 1) {
                    n_moved = cluster_threaded<t_uniform_v_weights, t_uniform_e_weights, t_rating_function>(g, p_manager, mapping, max_w, round, threads);
                } else {
                    n_moved = cluster_serial<t_uniform_v_weights, t_uniform_e_weights, t_rating_function>(g, p_manager, mapping, max_w, round);
                }

                if (static_cast<f64>(n_moved) < static_cast<f64>(g.n) * config.min_threshold) {
                    break;
                }

                // Swap active buffers for next round
                HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "swap_active");
                std::swap(active, active_next);
                active_next.initialize(g.n, 0);
            }

            merge_singletons<t_uniform_v_weights, t_uniform_e_weights, t_rating_function>(level, g, p_manager, mapping, max_w);

            // Check if coarsening resulted in an identity mapping
            bool ident_mapping = true;

            HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "is_identity_mapping");
            for (vertex_t u = 0; u < g.n; ++u) {
                if (u != mapping.get(u)) {
                    ident_mapping = false;
                    break;
                }
            }

            if (ident_mapping) {
                merge_when_identity<t_uniform_v_weights, t_uniform_e_weights, t_rating_function>(level, g, p_manager, mapping, max_w);
            }

            // Map clusters to a continuous range [0, coarse_n)
            HEIPROMAP_PROFILE_SCOPE("coarsening", "SizeConstrainedLP", "calc_map");
            remap.initialize(g.n, m_n);
            vertex_t new_id = 0;
            for (vertex_t u = 0; u < g.n; ++u) {
                const vertex_t id = mapping.get(u);
                if (remap[id] == m_n) {
                    remap[id] = new_id++;
                }
            }
            mapping.set_coarse_n(new_id);

            for (vertex_t u = 0; u < g.n; ++u) {
                mapping.set(u, remap[mapping.get(u)]);
            }
        }
    };
}

#endif //HEIPROMAP_SIZE_CONSTRAINED_LP_H
