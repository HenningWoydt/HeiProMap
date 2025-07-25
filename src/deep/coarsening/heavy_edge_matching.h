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

#ifndef HEIPROMAP_HEAVY_EDGE_MATCHING_H
#define HEIPROMAP_HEAVY_EDGE_MATCHING_H

#include <memory>
#include <vector>

#include "../../commons/definitions.h"
#include "../../commons/random_engine.h"
#include "../../serial/serial_definitions_1.h"
#include "../deep_definitions_1.h"
#include "../deep_definitions_2.h"
#include "../deep_definitions_3.h"

namespace HeiProMap {
    inline f32 edge_noise(vertex_t u, vertex_t v) {
        // Make (u, v) order-independent
        vertex_t a = u < v ? u : v;
        vertex_t b = u < v ? v : u;

        // Combine the two into one value
        u64 key = ((u64) a << 32) | (u64) b;

        // Simple 64-bit hash (xorshift-based)
        key ^= (key >> 33);
        key *= 0xff51afd7ed558ccdULL;
        key ^= (key >> 33);
        key *= 0xc4ceb9fe1a85ec53ULL;
        key ^= (key >> 33);

        // Map to [0, 1)
        f32 noise = (key & 0xFFFFFF) / (f32) 0x1000000; // 24-bit mantissa
        return noise * 0.01f; // small noise factor, e.g. 1%
    }

    inline u64 hash_neighborhood(const deep_graph_t &g,
                                 const vertex_t u) {
        u64 hash = 0;
        forall_guiv(g, u, i, v) {
                hash ^= (v * 11400714819323198485ull) ^ (v >> 32);
            }
        endfor
        return hash;
    }

    class ParallelHeavyEdgeMatchingConfiguration {
    public:
        u64 max_hem_iterations = 5;
        u64 max_twin_iterations = 5;
        u64 max_relative_iterations = 5;
    };

    class ParallelHeavyEdgeMatching {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
        weight_t m_l_max = 0;
        u64 m_threads = 1;

        const ParallelHeavyEdgeMatchingConfiguration *config = nullptr;
        RandomEngine *random_engine = nullptr;

    public:
        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const weight_t t_l_max,
                        const u64 t_threads,
                        RandomEngine &t_random_engine,
                        const ParallelHeavyEdgeMatchingConfiguration &i_config) {
            m_n = t_n;
            m_m = t_m;
            m_k = t_k;
            m_l_max = t_l_max;
            m_threads = t_threads;

            config = dynamic_cast<const ParallelHeavyEdgeMatchingConfiguration *>(&i_config);
            random_engine = &t_random_engine;
        }

        void match(const size_t level,
                   const deep_graph_t &g,
                   const deep_p_manager_t &p_manager,
                   Matching &matching) {
            std::vector<vertex_t> preferred(g.get_n());
            std::vector<f32> rating(g.get_n());

            heavy_edge_matching(level, g, p_manager, matching, preferred, rating);

            if ((f64) matching.size() * 2 / (f64) g.get_n() < 0.75) {
                twin_matching(level, g, p_manager, matching, preferred, rating);
            }

            if ((f64) matching.size() * 2 / (f64) g.get_n() < 0.75) {
                relative_matching(level, g, p_manager, matching, preferred, rating);
            }
        }

        void heavy_edge_matching(const size_t level,
                                 const deep_graph_t &g,
                                 const deep_p_manager_t &p_manager,
                                 Matching &matching,
                                 std::vector<vertex_t> &preferred,
                                 std::vector<f32> &rating) {
            for (size_t iteration = 0; iteration < config->max_hem_iterations; ++iteration) {
                std::iota(preferred.begin(), preferred.end(), 0);
                std::fill(rating.begin(), rating.end(), -1.0f);

#pragma omp parallel for num_threads(m_threads) default(none) shared(g, matching, preferred, rating) schedule(static)
                for (vertex_t u = 0; u < g.get_n(); ++u) {
                    if (matching.is_matched(u)) { continue; }
                    forall_guivw(g, u, i, v, w) {
                            if (matching.is_matched(v)) { continue; }
                            if (g.weight(u) + g.weight(v) > m_l_max) { continue; }

                            f32 r = (f32) (w * w) / (f32) (g.weight(u) * g.weight(v));
                            r += edge_noise(u, v);

                            if (r > rating[u]) {
                                rating[u] = r;
                                preferred[u] = v;
                            }
                        }
                    endfor
                }

#pragma omp parallel for num_threads(m_threads) default(none) shared(g, matching, preferred, rating) schedule(static)
                for (vertex_t u = 0; u < g.get_n(); ++u) {
                    vertex_t v = preferred[u];
                    vertex_t pref_v = preferred[v];

                    if (u == pref_v && u < v) {
                        matching.add(u, v);
                    }
                }
            }
        }

        void twin_matching(const size_t level,
                           const deep_graph_t &g,
                           const deep_p_manager_t &p_manager,
                           Matching &matching,
                           std::vector<vertex_t> &preferred,
                           std::vector<f32> &rating) {
            std::vector<std::pair<vertex_t, u64> > hashes(g.get_n());

#pragma omp parallel for num_threads(m_threads) default(none) shared(g, hashes) schedule(static)
            for (vertex_t u = 0; u < g.get_n(); ++u) {
                hashes[u] = {u, hash_neighborhood(g, u)};
            }
            std::sort(hashes.begin(), hashes.end());

            for (size_t iteration = 0; iteration < config->max_twin_iterations; ++iteration) {
                std::iota(preferred.begin(), preferred.end(), 0);
                std::fill(rating.begin(), rating.end(), -1.0f);

#pragma omp parallel for num_threads(m_threads) default(none) shared(g, hashes, matching, preferred, rating) schedule(static)
                for (u64 i = 0; i < hashes.size(); ++i) {
                    vertex_t u = hashes[i].first;
                    u64 hash = hashes[i].second;
                    if (matching.is_matched(u)) { continue; }

                    // search left
                    int idx = (int) i - 1;
                    while (idx >= 0 && hashes[idx].second == hash) {
                        vertex_t v = hashes[idx].first;
                        if (matching.is_matched(v)) {
                            idx -= 1;
                            continue;
                        }
                        if (g.weight(u) + g.weight(v) > m_l_max) {
                            idx -= 1;
                            continue;
                        }

                        f32 r = edge_noise(u, v);
                        if (r > rating[u]) {
                            rating[u] = r;
                            preferred[u] = v;
                        }

                        idx -= 1;
                    }

                    // search right
                    idx = (int) i + 1;
                    while (idx < (int) g.get_n() && hashes[idx].second == hash) {
                        vertex_t v = hashes[idx].first;
                        if (matching.is_matched(v)) {
                            idx += 1;
                            continue;
                        }
                        if (g.weight(u) + g.weight(v) > m_l_max) {
                            idx += 1;
                            continue;
                        }

                        f32 r = edge_noise(u, v);
                        if (r > rating[u]) {
                            rating[u] = r;
                            preferred[u] = v;
                        }

                        idx += 1;
                    }
                }

#pragma omp parallel for num_threads(m_threads) default(none) shared(g, hashes, matching, preferred, rating) schedule(static)
                for (u64 i = 0; i < hashes.size(); ++i) {
                    vertex_t u = hashes[i].first;

                    vertex_t v = preferred[u];
                    vertex_t pref_v = preferred[v];

                    if (u == pref_v && u < v) {
                        matching.add(u, v);
                    }
                }
            }
        }

        void relative_matching(const size_t level,
                               const deep_graph_t &g,
                               const deep_p_manager_t &p_manager,
                               Matching &matching,
                               std::vector<vertex_t> &preferred,
                               std::vector<f32> &rating) {
            for (size_t iteration = 0; iteration < config->max_relative_iterations; ++iteration) {
                std::iota(preferred.begin(), preferred.end(), 0);
                std::fill(rating.begin(), rating.end(), -1.0f);

#pragma omp parallel for num_threads(m_threads) default(none) shared(g, matching, preferred, rating) schedule(static)
                for (vertex_t u = 0; u < g.get_n(); ++u) {
                    if (matching.is_matched(u)) { continue; }

                    forall_guivw(g, u, i, mid_v, mid_w) {
                            vertex_t mid_deg = g.size(mid_v);
                            if (mid_deg > 10) { continue; } // avoid matchmaker with high degree

                            forall_guivw(g, mid_v, j, v, w) {
                                    if (u == v) { continue; }
                                    if (matching.is_matched(v)) { continue; }
                                    if (g.weight(u) + g.weight(v) > m_l_max) { continue; }

                                    f32 r = (f32) ((mid_w + w) * (mid_w + w)) / (f32) (g.weight(u) * g.weight(v));
                                    r += edge_noise(u, v);

                                    if (r > rating[u]) {
                                        rating[u] = r;
                                        preferred[u] = v;
                                    }
                                }
                            endfor
                        }
                    endfor
                }

#pragma omp parallel for num_threads(m_threads) default(none) shared(g, matching, preferred, rating) schedule(static)
                for (vertex_t u = 0; u < g.get_n(); ++u) {
                    vertex_t v = preferred[u];
                    vertex_t pref_v = preferred[v];

                    if (u == pref_v && u < v) {
                        matching.add(u, v);
                    }
                }
            }
        }
    };
}

#endif //HEIPROMAP_HEAVY_EDGE_MATCHING_H
