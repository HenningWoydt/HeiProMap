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

#include <map>

#include "../definitions.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/distance_oracle.h"
#include "../datastructures/partition_manager.h"
#include "../utility/random_engine.h"
#include "../utility/small_map.h"
#include "../utility/mapping.h"

namespace HeiProMap {
    class SizeConstrainedLPConfiguration {
    public:
        u64 max_rounds = 5;
        f64 min_threshold = 0.05;
        f64 f = 8;
        weight_t multiplier = 2;
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

        const SizeConstrainedLPConfiguration *config = nullptr;
        RandomEngine random_engine = RandomEngine(0);

    public:
        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const u64 seed,
                        const SizeConstrainedLPConfiguration &i_config) {
            ScopedTimer _t("coarsening", "SizeConstrainedLP", "initialize");

            m_n = t_n;
            m_m = t_m;
            m_k = t_k;

            config = &i_config;
            random_engine = RandomEngine(seed);
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void merge_when_identity([[maybe_unused]] const size_t level,
                                 const graph_t &g,
                                 [[maybe_unused]] const p_manager_t &p_manager,
                                 Mapping &mapping,
                                 weight_t max_w) {
            ScopedTimer _t("coarsening", "SizeConstrainedLP", "merge_when_identity");

            // 2) collect active cluster ids
            std::vector<vertex_t> ids;
            std::vector<u32> used(g.n, 0);
            ids.reserve(g.n);
            for (vertex_t id = 0; id < g.n; ++id) {
                ids.push_back(id);
            }

            // 3) sort by cluster weight ascending (merge small into something feasible)
            std::sort(ids.begin(), ids.end(), [&](vertex_t a, vertex_t b) { return g.v_weights[a] < g.v_weights[b]; });

            // 4) greedy merging:
            // for each small cluster a, find the smallest cluster b (b != a) s.t. cw[a]+cw[b] <= max_w
            // and move all members of a into b.
            for (size_t ia = 0; ia < ids.size(); ++ia) {
                vertex_t a = ids[ia];
                if (used[a] == 1) { continue; }

                // find best target b
                vertex_t best_b = (vertex_t) -1;

                for (size_t ib = 0; ib < ids.size(); ++ib) {
                    vertex_t b = ids[ids.size() - ib - 1];
                    if (b == a) continue;
                    if (used[b] == 1) { continue; }

                    if (p_manager[a] != p_manager[b]) { continue; }

                    weight_t sum = g.v_weights[a] + g.v_weights[b];
                    best_b = b;
                    if (sum <= max_w) {
                        best_b = b;
                        break;
                    }
                }

                if (best_b == (vertex_t) -1) {
                    // No feasible merge target for this cluster under max_w.
                    // In identity-mapping case this happens only if max_w < 2*min_vertex_weight etc.
                    continue;
                }

                used[a] = 1;
                used[best_b] = 1;

                mapping.set(a, best_b);
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void merge_singletons([[maybe_unused]] const size_t level,
                              const graph_t &g,
                              [[maybe_unused]] const p_manager_t &p_manager,
                              Mapping &mapping,
                              weight_t max_w) {
            ScopedTimer _t("coarsening", "SizeConstrainedLP", "merge_singletons");

            // 2) collect singleton vertices
            vertex_t singletons_size = 0;
            singletons.initialize(g.n);

            for (vertex_t u = 0; u < g.n; ++u) {
                vertex_t id = mapping.get(u);
                if (cluster_count[id] == 1) {
                    singletons[singletons_size] = u;
                    singletons_size += 1;
                }
            }

            if (singletons_size == 0) { return; }

            // 3) For each singleton u: choose best neighbor cluster by summed edge weight
            FlatMap<vertex_t, weight_t> flat_map;
            flat_map.reserve(128);

            for (size_t i = 0; i < singletons_size; ++i) {
                vertex_t u = singletons[i];
                vertex_t cur_id = mapping.get(u);
                partition_t u_id = p_manager[u];

                // might not be a singleton anymore if we already merged it earlier
                if (cluster_count[cur_id] != 1) { continue; }

                weight_t u_w = g.v_weights[u];
                vertex_t current_id = mapping.get(u);
                weight_t current_id_w = 0;

                flat_map.clear();
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const vertex_t v = g.edges_v[j];
                    const weight_t w = g.edges_w[j];
                    partition_t v_id = p_manager[v];
                    if (u_id != v_id) { continue; }

                    vertex_t id = mapping.get(v);
                    if (id == current_id) {
                        current_id_w += w;
                    } else {
                        if (u_w + cluster_weights[id] > max_w) { continue; }
                        flat_map.add(id, w);
                    }
                }

                vertex_t best_id = current_id;
                weight_t best_weight = current_id_w;
                for (auto [id, w]: flat_map) {
                    // check constraint only for candidates != current_id (same logic as your loop)
                    if (w > best_weight && u_w + cluster_weights[id] <= max_w) {
                        best_weight = w;
                        best_id = id;
                    }
                }

                if (best_id == cur_id) { continue; }

                // 4) apply merge: u moves from cur_id -> best_id
                mapping.set(u, best_id);

                cluster_weights[best_id] += u_w;
                cluster_count[best_id] += 1;

                cluster_weights[cur_id] -= u_w;
                cluster_count[cur_id] -= 1; // becomes 0
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void cluster([[maybe_unused]] const size_t level,
                     const graph_t &g,
                     [[maybe_unused]] const p_manager_t &p_manager,
                     Mapping &mapping,
                     f64 imbalance,
                     u64 threads) {
            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            weight_t max_w = 0; // (weight_t) ((f64) m_l_max / config->f);
            vertex_t max_deg = 0;
            // get max w and max deg
            {
                ScopedTimer _t("coarsening", "SizeConstrainedLP", "max");

                // determine the maximum allowed cluster weight
                for (vertex_t u = 0; u < g.n; ++u) {
                    max_w = std::max(max_w, g.v_weights[u]);
                    max_deg = std::max(max_deg, g.deg(u));
                }
                max_w *= config->multiplier;
                max_w = std::min(max_w, lmax);
            }
            const size_t B = (max_deg == 0) ? 1 : (floor_log2(max_deg) + 1);
            // get all vertices
            {
                ScopedTimer _t("coarsening", "SizeConstrainedLP", "flat_vertices");

                flat_vertices.initialize(g.n);
                bucket_sizes.initialize(B, 0);
                bucket_offsets.initialize(B);

                for (vertex_t u = 0; u < g.n; ++u) {
                    size_t d = g.deg(u);
                    size_t b = (d == 0) ? 0 : floor_log2(d);
                    bucket_sizes[b]++;
                }

                bucket_offsets[0] = 0;
                for (size_t i = 1; i < B; ++i) { bucket_offsets[i] = bucket_offsets[i - 1] + bucket_sizes[i - 1]; }

                for (vertex_t u = 0; u < g.n; ++u) {
                    size_t d = g.deg(u);
                    size_t b = (d == 0) ? 0 : floor_log2(d);
                    flat_vertices[bucket_offsets[b]] = u;
                    bucket_offsets[b] += 1;
                }
            }
            // setup cluster weights
            {
                ScopedTimer _t("coarsening", "SizeConstrainedLP", "cluster_weights");

                // set each vertex to its own id
                cluster_weights.initialize(g.n);
                cluster_count.initialize(g.n);
                for (vertex_t u = 0; u < g.n; ++u) {
                    mapping.set(u, u);
                    cluster_weights[u] = g.v_weights[u];
                    cluster_count[u] = 1;
                }
            }
            // setup active
            {
                ScopedTimer _t_active("coarsening", "SizeConstrainedLP", "active");

                active.initialize(g.n, 1);
                active_next.initialize(g.n, 1);
            }

            u64 n_moved = 0;

            for (u64 round = 0; round < config->max_rounds; ++round) {
                n_moved = 0;
                //
                {
                    ScopedTimer _t("coarsening", "SizeConstrainedLP", "shuffle_buckets");

                    for (size_t i = 0; i < B - 1; ++i) {
                        [[maybe_unused]] size_t beg = bucket_offsets[i];
                        [[maybe_unused]] size_t end = bucket_offsets[i + 1];

                        fast_shuffle_unchecked(flat_vertices.get_ptr() + beg, flat_vertices.get_ptr() + end, random_engine.generator);
                    }
                }
                // run clustering
                if (threads > 1) {
                    ScopedTimer _t("coarsening", "SizeConstrainedLP", "cluster_threaded");
                    #pragma omp parallel num_threads(threads)
                    {
                        FlatMap<vertex_t, weight_t> flat_map;
                        flat_map.reserve(128);

                        #pragma omp for schedule(dynamic) reduction(+:n_moved)
                        for (size_t i = 0; i < g.n; ++i) {
                            vertex_t u = flat_vertices[i];
                            if (active[u] == 0) { continue; }

                            weight_t u_w = t_uniform_v_weights ? 1 : g.v_weights[u];
                            partition_t u_id = p_manager[u];
                            vertex_t current_id = mapping.get(u);
                            weight_t current_id_w = 0;

                            flat_map.clear();
                            for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                                const vertex_t v = g.edges_v[j];
                                const weight_t w = g.edges_w[j];
                                partition_t v_id = p_manager[v];
                                if (u_id != v_id) { continue; }

                                vertex_t id = mapping.get(v);
                                if (id == current_id) {
                                    current_id_w += w;
                                } else {
                                    if (u_w + cluster_weights[id] > max_w) { continue; }
                                    flat_map.add(id, w);
                                }
                            }

                            vertex_t best_id = current_id;
                            weight_t best_weight = current_id_w;
                            for (auto [id, w]: flat_map) {
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
                                        const vertex_t v = g.edges_v[j]; {
                                            #pragma omp atomic write
                                            active_next[v] = 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                } else {
                    ScopedTimer _t("coarsening", "SizeConstrainedLP", "cluster_serial");
                    FlatMap<vertex_t, weight_t> flat_map;
                    flat_map.reserve(128);

                    for (size_t i = 0; i < g.n; ++i) {
                        vertex_t u = flat_vertices[i];
                        if (active[u] == 0) { continue; }

                        weight_t u_w = t_uniform_v_weights ? 1 : g.v_weights[u];
                        partition_t u_id = p_manager[u];

                        vertex_t current_id = mapping.get(u);
                        weight_t current_id_w = 0;

                        vertex_t best_id = current_id;
                        weight_t best_weight = 0;

                        flat_map.clear();
                        if constexpr (t_uniform_e_weights) {
                            for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                                const vertex_t v = g.edges_v[j];
                                partition_t v_id = p_manager[v];
                                if (u_id != v_id) { continue; }

                                vertex_t id = mapping.get(v);
                                if (id == current_id) {
                                    current_id_w += 1;
                                } else {
                                    if (u_w + cluster_weights[id] > max_w) { continue; }

                                    weight_t new_w = flat_map.add_and_ret(id, 1);
                                    if (new_w > best_weight) {
                                        best_weight = new_w;
                                        best_id = id;
                                    }
                                }
                            }
                        } else {
                            for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                                const vertex_t v = g.edges_v[j];
                                const weight_t w = g.edges_w[j];
                                partition_t v_id = p_manager[v];
                                if (u_id != v_id) { continue; }

                                vertex_t id = mapping.get(v);
                                if (id == current_id) {
                                    current_id_w += w;
                                } else {
                                    if (u_w + cluster_weights[id] > max_w) { continue; }

                                    weight_t new_w = flat_map.add_and_ret(id, w);
                                    if (new_w > best_weight) {
                                        best_weight = new_w;
                                        best_id = id;
                                    }
                                }
                            }
                        }

                        if (current_id_w > best_weight) {
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
                }
                // swap active
                {
                    ScopedTimer _t("coarsening", "SizeConstrainedLP", "swap_active");

                    std::swap(active, active_next);
                    active_next.initialize(g.n, 0);
                }

                if ((f64) n_moved < (f64) g.n * config->min_threshold) {
                    break;
                }
            }

            merge_singletons<t_uniform_v_weights, t_uniform_e_weights>(level, g, p_manager, mapping, max_w);

            bool ident_mapping = true;
            // is identity mapping
            {
                ScopedTimer _t("coarsening", "SizeConstrainedLP", "is_identity_mapping");

                for (vertex_t u = 0; u < g.n; ++u) {
                    ident_mapping &= u == mapping.get(u);
                }
            }
            if (ident_mapping) {
                merge_when_identity<t_uniform_v_weights, t_uniform_e_weights>(level, g, p_manager, mapping, max_w);
            }

            // map to a continuous range
            {
                ScopedTimer _t("coarsening", "SizeConstrainedLP", "calc_map");

                // mapping starts at 0 and increments
                remap.initialize(g.n, m_n);
                vertex_t new_id = 0;
                for (vertex_t u = 0; u < g.n; ++u) {
                    const vertex_t id = mapping.get(u);
                    if (remap[id] == m_n) {
                        remap[id] = new_id;
                        new_id += 1;
                    }
                }
                mapping.set_coarse_n(new_id);

                // remap
                for (vertex_t u = 0; u < g.n; ++u) {
                    const vertex_t id = mapping.get(u);
                    const vertex_t map_id = remap[id];
                    mapping.set(u, map_id);
                }
            }
        }
    };
}

#endif //HEIPROMAP_SIZE_CONSTRAINED_LP_H
