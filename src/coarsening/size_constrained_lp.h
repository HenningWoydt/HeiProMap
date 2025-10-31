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
#include "../definitions_1.h"
#include "../definitions_2.h"
#include "../utility/random_engine.h"
#include "../utility/small_map.h"
#include "../utility/mapping.h"

namespace HeiProMap {

    class SizeConstrainedLPConfiguration {
    public:
        u64 max_rounds = 5;
        f64 min_threshold = 0.05;
        f64 f = 64;
    };

    class SizeConstrainedLP {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
        weight_t m_l_max = 0;

        const SizeConstrainedLPConfiguration *config = nullptr;
        RandomEngine *random_engine = nullptr;

    public:
        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const weight_t t_l_max,
                        RandomEngine &t_random_engine,
                        const SizeConstrainedLPConfiguration &i_config) {
            ScopedTimer _t("io", "SizeConstrainedLP", "initialize");

            m_n = t_n;
            m_m = t_m;
            m_k = t_k;
            m_l_max = t_l_max;

            config = dynamic_cast<const SizeConstrainedLPConfiguration *>(&i_config);
            random_engine = &t_random_engine;
        }

        void cluster([[maybe_unused]] const size_t level,
                     const graph_t &g,
                     [[maybe_unused]] const p_manager_t &p_manager,
                     Mapping &mapping) {
            ScopedTimer _t_max("coarsening", "SizeConstrainedLP", "max");

            // determine the maximum allowed cluster weight
            weight_t max_w = (weight_t) ((f64) m_l_max / config->f);
            vertex_t max_deg = 0;
            forall_gu(g, u)
                {
                    max_w = std::max(max_w, g.weight(u));
                    max_deg = std::max(max_deg, g.size(u));
                }
            endfor

            _t_max.stop();
            ScopedTimer _t_flat_vertices("coarsening", "SizeConstrainedLP", "flat_vertices");

            // get list of all vertices
            AlignedArray<vertex_t> flat_vertices;
            flat_vertices.initialize(g.get_n());

            _t_flat_vertices.stop();
            ScopedTimer _t_bucket_sizes("coarsening", "SizeConstrainedLP", "bucket_sizes");

            const size_t B = (max_deg == 0) ? 1 : (floor_log2(max_deg) + 1);
            AlignedArray<vertex_t> bucket_sizes;
            bucket_sizes.initialize(B, 0);
            forall_gu(g, u)
                {
                    size_t d = g.size(u);
                    size_t b = (d == 0) ? 0 : floor_log2(d);
                    bucket_sizes[b]++;
                }
            endfor

            _t_bucket_sizes.stop();
            ScopedTimer _t_bucket_offsets("coarsening", "SizeConstrainedLP", "bucket_offsets");

            AlignedArray<vertex_t> bucket_offsets;
            bucket_offsets.initialize(B);
            bucket_offsets[0] = 0;
            for (size_t i = 1; i < B; ++i) { bucket_offsets[i] = bucket_offsets[i - 1] + bucket_sizes[i - 1]; }

            _t_bucket_offsets.stop();
            ScopedTimer _t_fill_flat_vertices("coarsening", "SizeConstrainedLP", "fill_flat_vertices");

            forall_gu(g, u)
                {
                    size_t d = g.size(u);
                    size_t b = (d == 0) ? 0 : floor_log2(d);
                    flat_vertices[bucket_offsets[b]] = u;
                    bucket_offsets[b] += 1;
                }
            endfor

            _t_fill_flat_vertices.stop();
            ScopedTimer _t_cluster_weights("coarsening", "SizeConstrainedLP", "cluster_weights");

            // set each vertex to its own id
            AlignedArray<weight_t> cluster_weights;
            cluster_weights.initialize(g.get_n());
            forall_gu(g, u)
                {
                    mapping.set_u(u, u);
                    cluster_weights[u] = g.weight(u);
                }
            endfor

            _t_cluster_weights.stop();
            ScopedTimer _t_active("coarsening", "SizeConstrainedLP", "active");

            AlignedArray<u8> active;
            active.initialize(g.get_n(), 1);
            AlignedArray<u8> active_next;
            active_next.initialize(g.get_n(), 0);

            _t_active.stop();
            ScopedTimer _t_flat_maps("coarsening", "SizeConstrainedLP", "flat_maps");


            FlatMap<vertex_t, weight_t> flat_map(128);
            u64 n_moved = 0;

            _t_flat_maps.stop();

            for (u64 round = 0; round < config->max_rounds; ++round) {
                ScopedTimer _t_reset_n_moved("coarsening", "SizeConstrainedLP", "reset_n_moved");

                n_moved = 0;

                _t_reset_n_moved.stop();
                ScopedTimer _t_cluster("coarsening", "SizeConstrainedLP", "cluster");

                for (size_t i = 0; i < g.get_n(); ++i) {
                    vertex_t u = flat_vertices[i];
                    if (active[u] == 0) { continue; }

                    weight_t u_w = g.weight(u);
                    vertex_t current_id = mapping.get_map_u(u);

                    flat_map.clear();
                    forall_guivw(g, u, j, v, w) {
                            vertex_t id = mapping.get_map_u(v);
                            flat_map[id] += w;
                        }
                    endfor

                    vertex_t best_id = current_id;
                    weight_t best_weight = 0;
                    for (auto [id, w_sum]: flat_map) {
                        weight_t cluster_weight = cluster_weights[id];
                        if (id != current_id && u_w + cluster_weight > max_w) { continue; }
                        if (w_sum > best_weight) {
                            best_weight = w_sum;
                            best_id = id;
                        }
                    }

                    // If you intend to move u, update mapping and cluster weights here.
                    if (best_id != current_id) {
                        // remove from current, add to new (keep these atomic or protected if parallel)
                        mapping.set_u(u, best_id);
                        cluster_weights[best_id] += u_w;
                        cluster_weights[current_id] -= u_w;

                        n_moved += 1;
                        if (round > 0) {
                            forall_guiv(g, u, j, v) {
                                    active_next[v] = 1;
                                }
                            endfor
                        }
                    }
                }

                _t_cluster.stop();
                ScopedTimer _t_swap_active("coarsening", "SizeConstrainedLP", "swap_active");

                if (round > 0) {
                    std::swap(active, active_next);
                    active_next.initialize(g.get_n(), 0);
                }

                _t_swap_active.stop();

                if ((f64) n_moved < (f64) g.get_n() * config->min_threshold) {
                    break;
                }
            }

            ScopedTimer _t_calc_map("coarsening", "SizeConstrainedLP", "calc_map");

            // mapping starts at 0 and increments
            AlignedArray<vertex_t> remap;
            remap.initialize(g.get_n(), m_n);
            vertex_t new_id = 0;
            forall_gu(g, u)
                {
                    const vertex_t id = mapping.get_map_u(u);
                    if (remap[id] == m_n) {
                        remap[id] = new_id;
                        new_id += 1;
                    }
                }
            endfor
            mapping.set_coarse_n(new_id);

            _t_calc_map.stop();
            ScopedTimer _t_remap("coarsening", "SizeConstrainedLP", "remap");

            // remap
            forall_gu(g, u)
                {
                    const vertex_t id = mapping.get_map_u(u);
                    const vertex_t map_id = remap[id];
                    mapping.set_u(u, map_id);
                }
            endfor
        }
    };
}

#endif //HEIPROMAP_SIZE_CONSTRAINED_LP_H
