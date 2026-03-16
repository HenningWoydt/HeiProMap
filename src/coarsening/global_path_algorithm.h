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

#ifndef HEIPROMAP_GLOBAL_PATH_ALGORITHM_H
#define HEIPROMAP_GLOBAL_PATH_ALGORITHM_H

#include <algorithm>
#include <chrono>
#include <queue>
#include <vector>

#include "../definitions.h"
#include "../utility/JSON_utils.h"
#include "../utility/random_engine.h"

namespace HeiProMap {
    struct Neighbors {
        vertex_t n1;
        vertex_t n2;
        f32 w1;
        f32 w2;
    };

#define IS_ENDPOINT(n, u) n.n2 == u
#define IS_NOT_ENDPOINT(n, u) n.n2 != u
#define IS_ONE_ENDPOINT(n, u) n.n1 != u && n.n2 == u
#define IS_UNMATCHED(n, u) n.n1 == u

    class GlobalPathAlgorithmConfiguration {
    public:
        size_t random_level = 4;
    };

    /**
     * Computes a matching based on the Global Path Algorithm from
     * > Jens Maue and Peter Sanders.
     * > Engineering Algorithms for Approximate Weighted Matching.
     * > Experimental Algorithms, 6th International Workshop, WEA 2007, Rome, Italy, June 6-8, 2007, Proceedings.
     */
    class GlobalPathAlgorithmMatcher {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;

        const GlobalPathAlgorithmConfiguration *config = nullptr;
        RandomEngine *random_engine = nullptr;

        AlignedArray<Neighbors> m_neighbors;
        AlignedArray<u32> path_id;
        AlignedArray<u32> path_length;

        AlignedArray<EdgeUVW> edges;
        size_t edges_size = 0;

        // for DP
        AlignedArray<f32> dp_w;
        AlignedArray<s64> dp_m;
        AlignedArray<u8> dp_take;
        AlignedArray<vertex_t> dp_edges;

        Matching dp_cycle_matches1;
        Matching dp_cycle_matches2;

    public:
        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        RandomEngine &t_random_engine,
                        const GlobalPathAlgorithmConfiguration &i_config) {
            ScopedTimer _t("io", "GlobalPathAlgorithmMatcher", "initialize");

            m_n = t_n;
            m_m = t_m;
            m_k = t_k;

            config = dynamic_cast<const GlobalPathAlgorithmConfiguration *>(&i_config);
            random_engine = &t_random_engine;

            m_neighbors.initialize(m_n);
            path_id.initialize(m_n);
            path_length.initialize(m_n);

            edges.initialize(m_m);

            dp_w.initialize(m_n);
            dp_m.initialize(m_n);
            dp_take.initialize(m_n);
            dp_edges.initialize(m_n);

            dp_cycle_matches1.initialize(m_n);
            dp_cycle_matches2.initialize(m_n);
        }

        template<typename PartitionManagerT>
        void match(const size_t level,
                   const graph_t &g,
                   const PartitionManagerT &p_manager,
                   Mapping &mapping,
                   f64 imbalance) {
            ScopedTimer _t_match("coarsening", "GlobalPathAlgorithmMatcher", "match");

            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            Matching matching;
            matching.initialize(g.n);

            if (level < config->random_level) {
                // use a random matching
                random_matching(level, g, matching, lmax);

                matching.set_translation();
                mapping.set_coarse_n(matching.get_n_coarse_nodes());
                for (vertex_t u = 0; u < matching.get_n(); ++u) {
                    mapping.set(u, matching.get_n(u));
                }

                return;
            }

            compute_ratings(g, p_manager, lmax);

            std::sort(edges.get_ptr(), edges.get_ptr() + edges_size, std::greater<>());

            for (vertex_t u = 0; u < g.n; ++u) {
                m_neighbors[u].n1 = u;
                m_neighbors[u].n2 = u;
            }

            for (size_t i = 0; i < edges_size; ++i) {
                auto [u, v, w] = edges[i];

                if (IS_NOT_ENDPOINT(m_neighbors[u], u) || IS_NOT_ENDPOINT(m_neighbors[v], v)) {
                    // u or v is not an endpoint
                    continue;
                }

                bool u_unmatched = IS_UNMATCHED(m_neighbors[u], u);
                bool v_unmatched = IS_UNMATCHED(m_neighbors[v], v);

                if (u_unmatched && v_unmatched) {
                    // both are unmatched, only one new path of length 1
                    m_neighbors[u].n1 = v;
                    m_neighbors[u].w1 = w;
                    m_neighbors[v].n1 = u;
                    m_neighbors[v].w1 = w;
                    path_id[u] = u;
                    path_id[v] = u;
                    path_length[u] = 1;
                    continue;
                }

                u32 u_id = path_id[u];
                u32 v_id = path_id[v];

                if (u_unmatched || v_unmatched) {
                    if (v_unmatched) {
                        std::swap(u, v);
                        std::swap(u_id, v_id);
                    }
                    // only one unmatched, enlarge path
                    m_neighbors[u].n1 = v;
                    m_neighbors[u].w1 = w;
                    m_neighbors[v].n2 = u;
                    m_neighbors[v].w2 = w;
                    path_id[u] = v_id;
                    path_length[v_id] += 1;
                    continue;
                }

                // cycle
                if (u_id == v_id) {
                    if (path_length[u_id] & 1) {
                        // same path and odd length size, close the cycle
                        path_length[u_id] += 1; // increase path length

                        // for u set v as a neighbor
                        m_neighbors[u].n2 = v;
                        m_neighbors[u].w2 = w;

                        // for v set u as a neighbor
                        m_neighbors[v].n2 = u;
                        m_neighbors[v].w2 = w;

                        // solve the cycle
                        solve_cycle(g, u, path_length[u_id], matching);
                        path_length[u_id] = 0;
                    }
                    continue;
                }

                // two paths, both u and v connect larger paths

                // for u set v as a neighbor
                m_neighbors[u].n2 = v;
                m_neighbors[u].w2 = w;

                // for v set u as a neighbor
                m_neighbors[v].n2 = u;
                m_neighbors[v].w2 = w;

                vertex_t v1 = v;
                vertex_t v2 = u;
                u32 id1 = v_id;
                u32 id2 = u_id;
                if (path_length[u_id] > path_length[v_id]) {
                    std::swap(v1, v2);
                    std::swap(id1, id2);
                }

                path_length[id1] += 1 + path_length[id2];
                while (m_neighbors[v2].n2 != v2) {
                    path_id[v2] = id1;
                    vertex_t temp_last_vertex = v1;
                    v1 = v2;
                    v2 = m_neighbors[v2].n1 == temp_last_vertex ? m_neighbors[v2].n2 : m_neighbors[v2].n1;
                }
                path_id[v2] = id1;
            }

            // process all paths
            forall_gu(g, u)
                {
                    if (IS_ONE_ENDPOINT(m_neighbors[u], u) && path_length[path_id[u]] > 0) {
                        solve_path(g, u, path_length[path_id[u]], matching);
                        path_length[path_id[u]] = 0;
                    }
                }
            endfor

            _t_match.stop();

            if ((f64) matching.size() * 2 < 0.75 * (f64) g.n) {
                two_hop_degree_one(level, g, p_manager, matching, imbalance);
            }
            if ((f64) matching.size() * 2 < 0.75 * (f64) g.n) {
                two_hop_twins(level, g, p_manager, matching, imbalance);
            }
            if ((f64) matching.size() * 2 < 0.75 * (f64) g.n) {
                two_hop_matchmaker(level, g, p_manager, matching, imbalance);
            }

            /*
#if ASSERT_ENABLED
            for (size_t i = 0; i < matching.size(); ++i) {
                const auto& [u, v] = matching[i];
                ASSERT(u != v);
            }
#endif

#if ASSERT_ENABLED
            std::vector<u8> hit(g.get_n(), 0);
            for (size_t i = 0; i < matching.size(); ++i) {
                const auto& [u, v] = matching[i];
                hit[u] += 1;
                hit[v] += 1;

                ASSERT(hit[u] == 1);
                ASSERT(hit[v] == 1);
            }
#endif
             */

            ScopedTimer _t("coarsening", "GlobalPathAlgorithmMatcher", "get_mapping");
            matching.set_translation();
            mapping.set_coarse_n(matching.get_n_coarse_nodes());
            for (vertex_t u = 0; u < matching.get_n(); ++u) {
                mapping.set(u, matching.get_n(u));
            }
        }

        template<typename PartitionManagerT>
        void compute_ratings(const graph_t &g, const PartitionManagerT &p_manager, weight_t lmax) {
            edges_size = 0;
            std::vector<EdgeUVW> local_edges;

            forall_gu(g, u)
                {
                    weight_t u_w = g.v_weights[u];

                    forall_guivw(g, u, j, v, w)
                        {
                            if (u > v) { continue; }
                            if (p_manager[u] != p_manager[v]) { continue; }
                            weight_t v_w = g.v_weights[v];

                            // if (u_w > 1.5*av_manager.get_n_active() / 20.0 * m_k) { continue; }
                            // if (v_w > 1.5*av_manager.get_n_active() / 20.0 * m_k) { continue; }

                            if (u_w + v_w > lmax) { continue; }

                            f32 edge_rating;

                            // edge_rating = ((f32) w) / (g.size(u) * g.size(v));
                            // edge_rating = (f32) w / (f32) (u_w * v_w);
                            edge_rating = ((f32) (w * w)) / ((f32) (u_w * v_w));
                            // edge_rating = ((f32) (w * w)) / ((f32) (u_w + v_w));
                            // edge_rating = ((f32) (w * w * w)) / ((f32) (u_w * v_w));
                            // edge_rating = (f32) w / (f32) (u_w * v_w * u_w * v_w);
                            // edge_rating = (f32)w;

                            // edges[edges_size++] = {u, v, edge_rating};
                            local_edges.push_back({u, v, edge_rating});
                        }
                    endfor
                }
            endfor

            // Critical section or lock-free append
            for (const auto &e: local_edges) {
                edges[edges_size++] = e;
            }
        }

        f32 solve_path_length_1([[maybe_unused]] const graph_t &g,
                                const vertex_t u,
                                Matching &matching) {
            vertex_t uu = u;
            vertex_t vv = m_neighbors[u].n1;
            f32 w = m_neighbors[u].w1;

            matching.add(uu, vv);

            return w;
        }

        f32 solve_path_length_2([[maybe_unused]] const graph_t &g,
                                const vertex_t u,
                                Matching &matching) {
            vertex_t v1 = u;
            vertex_t v2 = m_neighbors[u].n1;
            vertex_t v3;
            f32 w1 = m_neighbors[u].w1;
            f32 w2;

            if (m_neighbors[v2].n1 == v1) {
                v3 = m_neighbors[v2].n2;
                w2 = m_neighbors[v2].w2;
            } else {
                v3 = m_neighbors[v2].n1;
                w2 = m_neighbors[v2].w1;
            }

            vertex_t uu, vv;
            f32 w;
            if (w1 > w2) {
                uu = v1;
                vv = v2;
                w = w1;
            } else {
                uu = v2;
                vv = v3;
                w = w2;
            }

            matching.add(uu, vv);
            return w;
        }

        f32 solve_path(const graph_t &g,
                       const vertex_t u,
                       const u32 length,
                       Matching &matching) {
            // special case of length 1
            if (length == 1) {
                return solve_path_length_1(g, u, matching);
            }

            // special case of length 2
            if (length == 2) {
                return solve_path_length_2(g, u, matching);
            }

            vertex_t v1, v2, v3;
            f32 w;

            s64 i = 0;

            // first edge of the path
            v1 = u;
            dp_edges[i] = v1;
            if (m_neighbors[v1].n1 == v1) {
                v2 = m_neighbors[u].n2;
                w = m_neighbors[u].w2;
            } else {
                v2 = m_neighbors[u].n1;
                w = m_neighbors[u].w1;
            }
            dp_edges[i + 1] = v2; // save edge

            // init dp
            dp_w[i] = w;
            dp_m[i] = -1;
            dp_take[i] = 1;
            i += 1;

            // second edge of the path
            if (m_neighbors[v2].n1 == v1) {
                v3 = m_neighbors[v2].n2;
                w = m_neighbors[v2].w2;
            } else {
                v3 = m_neighbors[v2].n1;
                w = m_neighbors[v2].w1;
            }
            dp_edges[i + 1] = v3; // save edge

            // init dp
            if (w > dp_w[i - 1]) {
                dp_w[i] = w;
                dp_m[i] = -1;
                dp_take[i] = 1;
            } else {
                dp_w[i] = dp_w[i - 1];
                dp_m[i] = 0;
                dp_take[i] = 0;
            }
            i += 1;

            // all other edges of the path
            v1 = v2;
            v2 = v3;
            while (m_neighbors[v2].n1 != v2 && m_neighbors[v2].n2 != v2) {
                if (m_neighbors[v2].n1 == v1) {
                    v3 = m_neighbors[v2].n2;
                    w = m_neighbors[v2].w2;
                } else {
                    v3 = m_neighbors[v2].n1;
                    w = m_neighbors[v2].w1;
                }
                dp_edges[i + 1] = v3; // save edge

                // dp
                if (w + dp_w[i - 2] > dp_w[i - 1]) {
                    dp_w[i] = w + dp_w[i - 2];
                    dp_m[i] = i - 2;
                    dp_take[i] = 1;
                } else {
                    dp_w[i] = dp_w[i - 1];
                    dp_m[i] = i - 1;
                    dp_take[i] = 0;
                }

                v1 = v2;
                v2 = v3;
                i += 1;
            }

            s64 idx = i - 1;
            while (idx != -1) {
                if (dp_take[idx]) {
                    vertex_t uu = dp_edges[idx];
                    vertex_t vv = dp_edges[idx + 1];

                    matching.add(uu, vv);
                }
                idx = dp_m[idx];
            }
            return dp_w[i - 1];
        }

        void solve_cycle(const graph_t &g,
                         const vertex_t u,
                         const u32 length,
                         [[maybe_unused]] Matching &matching) {
            vertex_t n1 = m_neighbors[u].n1;
            vertex_t n2 = m_neighbors[u].n2;
            Neighbors original_u = m_neighbors[u];
            Neighbors original_n1 = m_neighbors[original_u.n1];
            Neighbors original_n2 = m_neighbors[original_u.n2];

            f32 matching_weight1 = 0.0;
            f32 matching_weight2 = 0.0;
            dp_cycle_matches1.clear();
            dp_cycle_matches2.clear();

            // cut connection between u and n1, n2 should point to self,
            m_neighbors[u].n1 = original_u.n2;
            m_neighbors[u].w1 = original_u.w2;
            m_neighbors[u].n2 = u;

            if (m_neighbors[n1].n1 == u) {
                m_neighbors[n1].n1 = original_n1.n2;
                m_neighbors[n1].w1 = original_n1.w2;
                m_neighbors[n1].n2 = n1;
            } else {
                m_neighbors[n1].n2 = n1;
            }
            matching_weight1 = solve_path(g, u, length - 1, dp_cycle_matches1);
            m_neighbors[u] = original_u;
            m_neighbors[n1] = original_n1;

            // cut connection between u and n2
            m_neighbors[u].n2 = u;

            if (m_neighbors[n2].n1 == u) {
                m_neighbors[n2].n1 = original_n2.n2;
                m_neighbors[n2].w1 = original_n2.w2;
                m_neighbors[n2].n2 = n2;
            } else {
                m_neighbors[n2].n2 = n2;
            }
            matching_weight2 = solve_path(g, u, length - 1, dp_cycle_matches2);
            m_neighbors[u] = original_u;
            m_neighbors[n2] = original_n2;

            Matching *dp_cycle_matches = &dp_cycle_matches1;
            if (matching_weight2 > matching_weight1) {
                dp_cycle_matches = &dp_cycle_matches2;
            }

            for (size_t i = 0; i < dp_cycle_matches->size(); ++i) {
                // matching.add((*dp_cycle_matches)[i].u, (*dp_cycle_matches)[i].v);
            }
        }

        void random_matching([[maybe_unused]] const size_t level,
                             const graph_t &g,
                             Matching &matching,
                             weight_t lmax) {
            std::vector<u8> is_matched(g.n, 0);

            forall_gu(g, u)
                {
                    if (is_matched[u]) { continue; }
                    weight_t u_w = g.v_weights[u];
                    forall_guiv(g, u, j, v)
                        {
                            if (is_matched[v]) { continue; }
                            weight_t v_w = g.v_weights[v];

                            if (u_w + v_w > lmax) { continue; }

                            is_matched[u] = 1;
                            is_matched[v] = 1;

                            matching.add(u, v);
                            break;
                        }
                    endfor
                }
            endfor

            /*
            #if ASSERT_ENABLED
                        for (size_t i = 0; i < matching.size(); ++i) {
                            const auto& [u, v] = matching[i];
                            ASSERT(u != v);
                        }
            #endif

            #if ASSERT_ENABLED
                        std::vector<u8> hit(g.get_n(), 0);
                        for (size_t i = 0; i < matching.size(); ++i) {
                            const auto& [u, v] = matching[i];
                            hit[u] += 1;
                            hit[v] += 1;

                            ASSERT(hit[u] == 1);
                            ASSERT(hit[v] == 1);
                        }
            #endif
             */
        }

        static inline u64 splitmix64(u64 x) {
            x += 0x9e3779b97f4a7c15ULL;
            x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
            x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
            return x ^ (x >> 31);
        }

        static f64 small_noise(vertex_t u, vertex_t v) {
            // order-independent: (u,v) and (v,u) give the same result
            u64 a = static_cast<u64>(std::min(u, v));
            u64 b = static_cast<u64>(std::max(u, v));

            // combine the pair into one 64-bit value, then hash it
            u64 key = a;
            key = key * 0x9e3779b97f4a7c15ULL + b;
            u64 h = splitmix64(key);

            // map to [0, 1)
            constexpr f64 inv = 1.0 / static_cast<f64>(std::numeric_limits<u64>::max());
            f64 x = static_cast<f64>(h) * inv;

            // tiny positive perturbation
            constexpr f64 eps = 1e-12;
            return eps * x;
        }

        template<typename PartitionManagerT>
        void two_hop_degree_one(const size_t level,
                                const graph_t &g,
                                const PartitionManagerT &p_manager,
                                Matching &matching,
                                f64 imbalance) {
            ScopedTimer _t("coarsening", "GlobalPathAlgorithmMatcher", "two_hop_degree_one");

            std::vector<vertex_t> preferred(g.n);
            std::iota(preferred.begin(), preferred.end(), 0);

            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            forall_gu(g, u)
                {
                    if (g.deg(u) != 1) { continue; }
                    if (matching.is_matched(u)) { continue; }

                    vertex_t preferred_vertex = u;
                    f64 best_rating = -std::numeric_limits<f64>::max();

                    vertex_t middle_vertex = g.edges_v[g.neighborhoods[u]];
                    weight_t middle_w = g.edges_w[g.neighborhoods[u]];

                    forall_guivw(g, middle_vertex, i, v, w)
                        {
                            if (u == v) { continue; }
                            if (g.deg(v) != 1) { continue; }
                            if (matching.is_matched(v)) { continue; }
                            if (g.v_weights[u] + g.v_weights[v] > lmax) { continue; }

                            f64 rating = (f64) (w + middle_w) + small_noise(u, v);

                            if (rating > best_rating) {
                                best_rating = rating;
                                preferred_vertex = v;
                            }
                        }
                    endfor
                    preferred[u] = preferred_vertex;
                }
            endfor

            forall_gu(g, u)
                {
                    if (g.deg(u) != 1) { continue; }
                    if (matching.is_matched(u)) { continue; }

                    vertex_t v = preferred[u];
                    if (u == v) { continue; }

                    if (preferred[v] == u && u < v) {
                        matching.add(u, v);
                    }
                }
            endfor
        }

        static inline uint64_t hash_combine_u64(uint64_t a, uint64_t b) {
            return splitmix64(a ^ (b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2)));
        }

        static inline uint64_t hash_edge(vertex_t v, weight_t w) {
            uint64_t hv = splitmix64(static_cast<uint64_t>(v));
            uint64_t hw = splitmix64(static_cast<uint64_t>(w));
            return hash_combine_u64(hv, hw);
        }

        static inline uint64_t neighborhood_hash(const graph_t &g, vertex_t u) {
            uint64_t x = splitmix64(static_cast<uint64_t>(g.deg(u)));
            uint64_t s1 = 0;
            uint64_t s2 = 0;

            forall_guivw(g, u, i, v, w)
                {
                    uint64_t he = hash_edge(v, w);
                    x ^= he;
                    s1 += he;
                    s2 += splitmix64(he);
                }
            endfor

            return hash_combine_u64(x, hash_combine_u64(s1, s2));
        }

        static inline bool same_neighborhood(const graph_t &g, vertex_t u, vertex_t v) {
            if (g.deg(u) != g.deg(v)) {
                return false;
            }

            std::vector<std::pair<vertex_t, weight_t> > nu;
            std::vector<std::pair<vertex_t, weight_t> > nv;
            nu.reserve(g.deg(u));
            nv.reserve(g.deg(v));

            forall_guivw(g, u, i, x, w)
                {
                    nu.emplace_back(x, w);
                }
            endfor

            forall_guivw(g, v, i, x, w)
                {
                    nv.emplace_back(x, w);
                }
            endfor

            std::sort(nu.begin(), nu.end());
            std::sort(nv.begin(), nv.end());

            return nu == nv;
        }

        template<typename PartitionManagerT>
        void two_hop_twins(const size_t level,
                           const graph_t &g,
                           const PartitionManagerT &p_manager,
                           Matching &matching,
                           f64 imbalance) {
            ScopedTimer _t("coarsening", "GlobalPathAlgorithmMatcher", "two_hop_twins");

            struct Candidate {
                uint64_t hash;
                vertex_t u;
            };

            std::vector<Candidate> candidates;
            candidates.reserve(g.n);

            weight_t lmax = std::ceil((1.0 + imbalance) * (static_cast<f64>(g.g_weight) / static_cast<f64>(p_manager.k)));

            forall_gu(g, u)
                {
                    if (matching.is_matched(u)) { continue; }
                    candidates.push_back({neighborhood_hash(g, u), u});
                }
            endfor

            std::sort(candidates.begin(), candidates.end(), [](const Candidate &a, const Candidate &b) {
                if (a.hash != b.hash) return a.hash < b.hash;
                return a.u < b.u;
            });

            size_t begin = 0;
            while (begin < candidates.size()) {
                size_t end = begin + 1;
                while (end < candidates.size() && candidates[end].hash == candidates[begin].hash) {
                    ++end;
                }

                // Split hash bucket into exact-equality subgroups
                std::vector<std::vector<vertex_t> > exact_groups;

                for (size_t i = begin; i < end; ++i) {
                    vertex_t u = candidates[i].u;
                    if (matching.is_matched(u)) { continue; }

                    bool placed = false;
                    for (auto &group: exact_groups) {
                        if (same_neighborhood(g, u, group.front())) {
                            group.push_back(u);
                            placed = true;
                            break;
                        }
                    }

                    if (!placed) {
                        exact_groups.push_back({u});
                    }
                }

                // Match greedily inside each exact group
                for (auto &group: exact_groups) {
                    std::sort(group.begin(), group.end(), [&](vertex_t a, vertex_t b) {
                        if (g.v_weights[a] != g.v_weights[b]) {
                            return g.v_weights[a] < g.v_weights[b];
                        }
                        return a < b;
                    });

                    for (size_t i = 0; i + 1 < group.size();) {
                        vertex_t u = group[i];
                        if (matching.is_matched(u)) {
                            ++i;
                            continue;
                        }

                        size_t j = i + 1;
                        while (j < group.size()) {
                            vertex_t v = group[j];
                            if (!matching.is_matched(v) &&
                                g.v_weights[u] + g.v_weights[v] <= lmax) {
                                matching.add(u, v);
                                break;
                            }
                            ++j;
                        }

                        ++i;
                    }
                }

                begin = end;
            }
        }

        template<typename PartitionManagerT>
        void two_hop_matchmaker(const size_t level,
                                const graph_t &g,
                                const PartitionManagerT &p_manager,
                                Matching &matching,
                                f64 imbalance) {
            ScopedTimer _t("coarsening", "GlobalPathAlgorithmMatcher", "two_hop_matchmaker");

            std::vector<vertex_t> preferred(g.n);
            std::iota(preferred.begin(), preferred.end(), 0);

            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));

            forall_gu(g, u)
                {
                    if (matching.is_matched(u)) { continue; }

                    vertex_t preferred_vertex = u;
                    f64 best_rating = -std::numeric_limits<f64>::max();

                    forall_guivw(g, u, j, middle_vertex, middle_w)
                        {
                            forall_guivw(g, middle_vertex, i, v, w)
                                {
                                    if (u == v) { continue; }
                                    if (matching.is_matched(v)) { continue; }
                                    if (g.v_weights[u] + g.v_weights[v] > lmax) { continue; }

                                    f64 rating = (f64) (w + middle_w) + small_noise(u, v);

                                    if (rating > best_rating) {
                                        best_rating = rating;
                                        preferred_vertex = v;
                                    }
                                }
                            endfor
                            preferred[u] = preferred_vertex;
                        }
                    endfor
                }
            endfor

            forall_gu(g, u)
                {
                    if (matching.is_matched(u)) { continue; }

                    vertex_t v = preferred[u];
                    if (u == v) { continue; }

                    if (preferred[v] == u && u < v) {
                        matching.add(u, v);
                    }
                }
            endfor
        }
    };
}

#endif //HEIPROMAP_GLOBAL_PATH_ALGORITHM_H
