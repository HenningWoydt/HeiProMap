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
#include <numeric>
#include <set>
#include <vector>

#include "../../definitions.h"
#include "../../macros.h"
#include "../interfaces/ISerialMatcher.h"

namespace HeiProMap {
    class Path {
    public:
        std::vector<vertex_t> vertices;
        std::vector<f64> weights;

        explicit Path(vertex_t v) { vertices.push_back(v); }

        size_t size() const { return weights.size(); }

        void add(vertex_t v, f64 w) {
            vertices.push_back(v);
            weights.push_back(w);
        }

        void add(Path& p, f64 weight) {
            ASSERT(&p != this);
            for (vertex_t v : p.vertices) {
                vertices.push_back(v);
            }
            weights.push_back(weight);
            for (f64 w : p.weights) {
                weights.push_back(w);
            }
        }

        void reverse() {
            std::reverse(vertices.begin(), vertices.end());
            std::reverse(weights.begin(), weights.end());
        }

        vertex_t front() {
            return vertices[0];
        }

        vertex_t back() {
            return vertices.back();
        }

        void clear() {
            vertices.clear();
            weights.clear();
        }

        bool empty() const {
            return weights.empty();
        }

        bool is_cycle() {
            return vertices.size() > 1 && vertices.front() == vertices.back();
        }

        Path copy() const {
            Path p(0);
            p.vertices = vertices;
            p.weights  = weights;
            return p;
        }

        void pop_front() {
            vertices.erase(vertices.begin());
            weights.erase(weights.begin());
        }

        void pop_back() {
            vertices.pop_back();
            weights.pop_back();
        }
    };

    class DPSolver {
    public:
        void solve(Path& p, std::vector<EdgeUV>& matches) {
            f64 matching_weight = 0.0;
            if (p.is_cycle()) {
                solve_cycle(p, matches);
            } else {
                solve_line(p, matches, matching_weight);
            }
        }

    private:
        std::vector<f64> w;
        std::vector<u32> m;
        std::vector<bool> take;

        void solve_cycle(Path& p, std::vector<EdgeUV>& matches) {
            // remove e1
            f64 matching_weight_1 = 0.0;
            std::vector<EdgeUV> matches_1;
            Path p1 = p.copy();
            p1.pop_front();
            solve_line(p1, matches_1, matching_weight_1);

            // remove e2
            f64 matching_weight_2 = 0.0;
            std::vector<EdgeUV> matches_2;
            Path p2 = p.copy();
            p2.pop_front();
            solve_line(p2, matches_2, matching_weight_2);

            if (matching_weight_1 > matching_weight_2) {
                for (const auto& e : matches_1) {
                    matches.emplace_back(e);
                }
            } else {
                for (const auto& e : matches_2) {
                    matches.emplace_back(e);
                }
            }
        }

        void solve_line(Path& p, std::vector<EdgeUV>& matches, f64& matching_weight) {
            if (p.size() == 1) {
                matches.emplace_back(p.vertices[0], p.vertices[1]);
                return;
            }

            w.resize(p.size());
            m.resize(p.size());
            take.resize(p.size());
            w[0]    = 0;
            m[0]    = 0;
            take[0] = false;

            w[1]    = p.weights[0];
            m[1]    = 1;
            take[1] = true;

            for (size_t i = 2; i < p.size(); ++i) {
                if (p.weights[i] + w[i - 2] > w[i - 1]) {
                    w[i]    = p.weights[i] + w[i - 2];
                    m[i]    = i - 2;
                    take[i] = true;
                } else {
                    w[i]    = w[i - 1];
                    m[i]    = i - 1;
                    take[i] = false;
                }
            }

            u32 idx         = m[p.size() - 1];
            matching_weight = w[p.size() - 1];
            while (idx > 1) {
                if (take[idx]) {
                    matches.emplace_back(p.vertices[idx - 1], p.vertices[idx]);
                }
                idx = m[idx];
            }
            if (idx == 1) {
                matches.emplace_back(p.vertices[idx - 1], p.vertices[idx]);
            }
        }
    };

    class GlobalPathAlgorithmMatcherNew final : public ISerialMatcher {
        vertex_t m_n = 0;
        vertex_t m_m = 0;

        std::vector<vertex_t> neighbor1;
        std::vector<vertex_t> neighbor2;
        std::vector<f64> weight1;
        std::vector<f64> weight2;
        std::vector<u32> path_id;
        std::vector<u32> path_length;

        std::vector<vertex_t> cycle_vertices;

        weight_t l_max = 0;

        std::vector<EdgeUVW> edges;

        // for DP
        std::vector<f64> w;
        std::vector<u32> m;
        std::vector<bool> take;
        std::vector<EdgeUV> dp_edges;
        std::vector<f64> dp_weights;

        std::vector<EdgeUV> dp_cycle_matches1;
        std::vector<EdgeUV> dp_cycle_matches2;

    public:
        GlobalPathAlgorithmMatcherNew() = default;

        void initialize(const vertex_t n, const vertex_t m, const weight_t t_l_max) override {
            m_n = n;
            m_m = m;

            neighbor1.resize(n);
            neighbor2.resize(n);
            weight1.resize(n);
            weight2.resize(n);
            path_id.resize(n);
            path_length.resize(n);

            edges.reserve(m);
            l_max = t_l_max;
        }

        template <typename TSerialGraph, typename TSerialActiveVertexManager>
        void match(TSerialGraph& g,
                   TSerialActiveVertexManager& av_manager,
                   std::vector<EdgeUV>& matches) {
            edges.clear();
            for (vertex_t u : av_manager) {
                for (const auto& [v, w] : g[u]) {
                    // f64 edge_rating = (((f64) w) * ((f64) w)) / (g.get_weight(u) * g.get_weight(v));
                    // f64 edge_rating = ((f64) w) / (g.size(u) * g.size(v));
                    f64 edge_rating = (f64)w / (g.get_weight(u) * g.get_weight(v));
                    edges.emplace_back(u, v, edge_rating);
                }
            }
            std::sort(edges.begin(), edges.end(), std::greater<EdgeUVW>());

            std::iota(neighbor1.begin(), neighbor1.end(), 0);
            std::iota(neighbor2.begin(), neighbor2.end(), 0);
            std::iota(path_id.begin(), path_id.end(), 0);
            std::fill(path_length.begin(), path_length.end(), 0);

            cycle_vertices.clear();
            for (const auto& [u, v, w] : edges) {
                bool u_is_endpoint = neighbor1[u] == u || neighbor2[u] == u;
                bool v_is_endpoint = neighbor1[v] == v || neighbor2[v] == v;

                if (!u_is_endpoint || !v_is_endpoint) {
                    // one of the vertices is not endpoint
                    continue;
                }

                // cycle
                if (path_id[u] == path_id[v] && path_length[path_id[u]] & 1) {
                    // same path and odd length size, close the cycle
                    path_length[path_id[u]] += 1; // increase path length

                    // for u set v as a neighbor
                    if (neighbor1[u] == u) {
                        neighbor1[u] = v;
                        weight1[u]   = w;
                    } else {
                        neighbor2[u] = v;
                        weight2[u]   = w;
                    }

                    // for v set u as a neighbor
                    if (neighbor1[v] == v) {
                        neighbor1[v] = u;
                        weight1[v]   = w;
                    } else {
                        neighbor2[v] = u;
                        weight2[v]   = w;
                    }
                    cycle_vertices.push_back(u);
                    continue;
                }
                if (path_id[u] == path_id[v]) { continue; } // not applicable

                ASSERT(path_id[u] != path_id[v]);

                // two paths

                // for u set v as a neighbor
                if (neighbor1[u] == u) {
                    neighbor1[u] = v;
                    weight1[u]   = w;
                } else {
                    neighbor2[u] = v;
                    weight2[u]   = w;
                }

                // for v set u as a neighbor
                if (neighbor1[v] == v) {
                    neighbor1[v] = u;
                    weight1[v]   = w;
                } else {
                    neighbor2[v] = u;
                    weight2[v]   = w;
                }

                if (path_length[path_id[u]] <= path_length[path_id[v]]) {
                    // pull the path of u into the path of v
                    path_length[path_id[v]] += 1 + path_length[path_id[u]];

                    // set all path_ids of vertices on the u-path to path_id[v]
                    vertex_t last_vertex = v;
                    vertex_t curr_vertex = u;
                    path_id[curr_vertex] = path_id[v];
                    while (neighbor1[curr_vertex] != curr_vertex && neighbor2[curr_vertex] != curr_vertex) {
                        path_id[curr_vertex]      = path_id[v];
                        vertex_t temp_last_vertex = last_vertex;
                        last_vertex               = curr_vertex;
                        curr_vertex               = neighbor1[curr_vertex] == temp_last_vertex ? neighbor2[curr_vertex] : neighbor1[curr_vertex];
                    }
                    path_id[curr_vertex] = path_id[v];
                } else {
                    // pull the path of v into the path of u
                    path_length[path_id[u]] += 1 + path_length[path_id[v]];

                    // set all path_ids of vertices on the v-path to path_id[u]
                    vertex_t last_vertex = u;
                    vertex_t curr_vertex = v;
                    path_id[curr_vertex] = path_id[u];
                    while (neighbor1[curr_vertex] != curr_vertex && neighbor2[curr_vertex] != curr_vertex) {
                        path_id[curr_vertex]      = path_id[u];
                        vertex_t temp_last_vertex = last_vertex;
                        last_vertex               = curr_vertex;
                        curr_vertex               = neighbor1[curr_vertex] == temp_last_vertex ? neighbor2[curr_vertex] : neighbor1[curr_vertex];
                    }
                    path_id[curr_vertex] = path_id[u];
                }
            }

            // process all paths
            size_t n_paths       = 0;
            size_t paths_lengths = 0;
            for (vertex_t u : av_manager) {
                bool is_endpoint = (neighbor1[u] == u && neighbor2[u] != u) || (neighbor1[u] != u && neighbor2[u] == u);

                if (is_endpoint) {
                    n_paths += 1;
                    paths_lengths += path_length[path_id[u]];
                    solve_path(u, matches);
                }
            }

            // process all cycles
            size_t n_cycles        = 0;
            size_t cycless_lengths = 0;
            for (vertex_t u : cycle_vertices) {
                n_cycles += 1;
                cycless_lengths += path_length[path_id[u]];
                solve_cycle(u, matches);
            }


#if ASSERT_ENABLED
            for (const EdgeUV& e : matches) {
                ASSERT(e.u != e.v);
                ASSERT(av_manager.is_active(e.u));
                ASSERT(av_manager.is_active(e.v));
            }
#endif

#if ASSERT_ENABLED
            std::vector<u8> hit(g.get_n(), 0);
            for (const auto& e : matches) {
                hit[e.u] += 1;
                hit[e.v] += 1;

                ASSERT(hit[e.u] != 2);
                ASSERT(hit[e.v] != 2);
            }
#endif
        }

        bool assert_no_duplicate_vertex(const std::vector<EdgeUV>& dp_edges) {
            if (dp_edges.empty()) { return true; }
            std::vector<vertex_t> vertex_set;
            vertex_set.push_back(dp_edges[0].u);
            vertex_set.push_back(dp_edges[0].v);
            for (size_t i = 1; i < dp_edges.size(); i++) {
                vertex_set.push_back(dp_edges[i].v);
            }
            std::sort(vertex_set.begin(), vertex_set.end());
            return no_duplicates_sorted(vertex_set);
        }

        f64 solve_path(vertex_t u, std::vector<EdgeUV>& matches) {
            // save the path into better format
            dp_edges.clear();
            dp_weights.clear();

            ASSERT(!(neighbor1[u] == u && neighbor2[u] == u));

            if (neighbor1[u] == u) {
                dp_edges.emplace_back(u, neighbor2[u]);
                dp_weights.emplace_back(weight2[u]);
            } else {
                dp_edges.emplace_back(u, neighbor1[u]);
                dp_weights.emplace_back(weight1[u]);
            }

            vertex_t last_vertex = u;
            vertex_t curr_vertex = neighbor1[u] == u ? neighbor2[u] : neighbor1[u];
            while (neighbor1[curr_vertex] != curr_vertex && neighbor2[curr_vertex] != curr_vertex) {
                if (neighbor1[curr_vertex] == last_vertex) {
                    dp_edges.emplace_back(curr_vertex, neighbor2[curr_vertex]);
                    dp_weights.emplace_back(weight2[curr_vertex]);

                    last_vertex = curr_vertex;
                    curr_vertex = neighbor2[curr_vertex];
                } else {
                    dp_edges.emplace_back(curr_vertex, neighbor1[curr_vertex]);
                    dp_weights.emplace_back(weight1[curr_vertex]);

                    last_vertex = curr_vertex;
                    curr_vertex = neighbor1[curr_vertex];
                }
            }
            ASSERT(assert_no_duplicate_vertex(dp_edges));

            if (dp_edges.size() == 1) {
                matches.emplace_back(dp_edges[0].u, dp_edges[0].v);
                // destroy the endpoints
                vertex_t v1   = dp_edges[0].u;
                vertex_t v2   = dp_edges.back().v;
                neighbor1[v1] = -1;
                neighbor2[v1] = -1;
                neighbor1[v2] = -1;
                neighbor2[v2] = -1;
                return dp_weights[0];
            }

            w.resize(dp_edges.size());
            m.resize(dp_edges.size());
            take.resize(dp_edges.size());

            w[0]    = 0;
            m[0]    = 0;
            take[0] = false;

            w[1]    = dp_weights[0];
            m[1]    = 1;
            take[1] = true;

            for (size_t i = 2; i < dp_edges.size(); ++i) {
                if (dp_weights[i] + w[i - 2] > w[i - 1]) {
                    w[i]    = dp_weights[i] + w[i - 2];
                    m[i]    = i - 2;
                    take[i] = true;
                } else {
                    w[i]    = w[i - 1];
                    m[i]    = i - 1;
                    take[i] = false;
                }
            }

            u32 idx             = m[dp_edges.size() - 1];
            f64 matching_weight = 0.0;
            while (idx > 1) {
                if (take[idx]) {
                    matches.push_back(dp_edges[idx]);
                    matching_weight += dp_weights[idx];
                }
                idx = m[idx];
            }
            if (idx == 1) {
                matches.push_back(dp_edges[idx]);
                matching_weight += dp_weights[idx];
            }

            // destroy the endpoints
            vertex_t v1   = dp_edges[0].u;
            vertex_t v2   = dp_edges.back().v;
            neighbor1[v1] = -1;
            neighbor2[v1] = -1;
            neighbor1[v2] = -1;
            neighbor2[v2] = -1;
            return matching_weight;
        }

        void solve_cycle(vertex_t u, std::vector<EdgeUV>& matches) {
            vertex_t left  = neighbor1[u];
            vertex_t right = neighbor2[u];

            vertex_t left_left   = neighbor1[left];
            vertex_t left_right  = neighbor2[left];
            vertex_t right_left  = neighbor1[right];
            vertex_t right_right = neighbor2[right];

            f64 matching_weight1 = 0.0;
            f64 matching_weight2 = 0.0;
            dp_cycle_matches1.clear();
            dp_cycle_matches2.clear();

            // remove left edge and solve
            neighbor1[u] = u;
            if (neighbor1[left] == u) { neighbor1[left] = left; } else { neighbor2[left] = left; }
            matching_weight1 = solve_path(u, dp_cycle_matches1);
            neighbor1[u]     = left;
            neighbor2[u]     = right;
            neighbor1[left]  = left_left;
            neighbor2[left]  = left_right;

            // remove right edge and solve
            neighbor2[u] = u;
            if (neighbor1[right] == u) { neighbor1[right] = right; } else { neighbor2[right] = right; }
            matching_weight2 = solve_path(u, dp_cycle_matches2);
            neighbor1[u]     = left;
            neighbor2[u]     = right;
            neighbor1[right] = right_left;
            neighbor2[right] = right_right;

            if (matching_weight1 >= matching_weight2) {
                for (auto& e : dp_cycle_matches1) {
                    matches.emplace_back(e);
                }
            } else {
                for (auto& e : dp_cycle_matches2) {
                    matches.emplace_back(e);
                }
            }
        }
    };

    /**
     * Computes a matching based on the Global Path Algorithm from
     * > Jens Maue and Peter Sanders.
     * > Engineering Algorithms for Approximate Weighted Matching.
     * > Experimental Algorithms, 6th International Workshop, WEA 2007, Rome, Italy, June 6-8, 2007, Proceedings.
     */
    class GlobalPathAlgorithmMatcher final : public ISerialMatcher {
        vertex_t m_n = 0;
        vertex_t m_m = 0;

        std::vector<u8> is_endpoint;
        std::vector<u32> path_id;

        weight_t l_max = 0;

        std::vector<EdgeUVW> edges;

        std::vector<Path> paths;

        GlobalPathAlgorithmMatcherNew new_matcher;

    public:
        GlobalPathAlgorithmMatcher() = default;

        void initialize(const vertex_t n, const vertex_t m, const weight_t t_l_max) override {
            m_n = n;
            m_m = m;

            edges.reserve(m);
            is_endpoint.resize(n);
            path_id.resize(n);
            l_max = t_l_max;

            new_matcher.initialize(n, m, t_l_max);
        }

        template <typename TSerialGraph, typename TSerialActiveVertexManager>
        void match(TSerialGraph& g,
                   TSerialActiveVertexManager& av_manager,
                   std::vector<EdgeUV>& matches) {
            // new_matcher.match(g, av_manager, matches);
            // return;

            edges.clear();
            for (vertex_t u : av_manager) {
                for (const auto& [v, w] : g[u]) {
                    // f64 edge_rating = (((f64) w) * ((f64) w)) / (g.get_weight(u) * g.get_weight(v));
                    // f64 edge_rating = ((f64) w) / (g.size(u) * g.size(v));
                    f64 edge_rating = (f64)w / (g.get_weight(u) * g.get_weight(v));
                    edges.emplace_back(u, v, edge_rating);
                }
            }
            std::sort(edges.begin(), edges.end(), std::greater<EdgeUVW>());

            std::cout << "no_duplicates " << no_duplicates(edges) << std::endl;

            std::fill(is_endpoint.begin(), is_endpoint.end(), 1);
            std::iota(path_id.begin(), path_id.end(), 0);
            paths.clear();
            for (vertex_t u = 0; u < is_endpoint.size(); ++u) {
                paths.emplace_back(u);
            }

            for (const auto& [u, v, w] : edges) {
                if (!is_endpoint[u] || !is_endpoint[v]) {
                    continue; // not applicable
                }

                // cycle
                if (path_id[u] == path_id[v] && (paths[path_id[u]].size() & 1)) {
                    // same path and odd length size, close the cycle
                    is_endpoint[u] = false;
                    is_endpoint[v] = false;
                    if (paths[path_id[u]].back() == v) {
                        paths[path_id[u]].reverse();
                    }
                    paths[path_id[u]].add(v, w);
                    continue;
                }
                if (path_id[u] == path_id[v]) { continue; } // not applicable

                ASSERT(path_id[u] != path_id[v]);

                // two paths
                vertex_t first_vertex, second_vertex;
                if (paths[path_id[u]].front() == u || paths[path_id[u]].back() == u) {
                    first_vertex  = u;
                    second_vertex = v;
                } else {
                    first_vertex  = v;
                    second_vertex = u;
                }

                // correct order
                if (paths[path_id[first_vertex]].front() == first_vertex) { paths[path_id[first_vertex]].reverse(); }
                if (paths[path_id[second_vertex]].back() == second_vertex) { paths[path_id[second_vertex]].reverse(); }

                // add the paths together
                paths[path_id[first_vertex]].add(paths[path_id[second_vertex]], w);
                paths[path_id[second_vertex]].clear();

                vertex_t p_id        = path_id[first_vertex];
                vertex_t last_vertex = paths[p_id].back();

                path_id[last_vertex] = path_id[first_vertex];

                is_endpoint[first_vertex]                         = false;
                is_endpoint[second_vertex]                        = false;
                is_endpoint[paths[path_id[first_vertex]].front()] = true;
                is_endpoint[paths[path_id[first_vertex]].back()]  = true;
            }

            // use dynamic programming to compute maximum weight matchings
            DPSolver dp_solver;
            for (Path& p : paths) {
                if (!p.empty()) {
                    dp_solver.solve(p, matches);
                }
            }


#if ASSERT_ENABLED
            for (const EdgeUV& e : matches) {
                ASSERT(e.u != e.v);
                ASSERT(av_manager.is_active(e.u));
                ASSERT(av_manager.is_active(e.v));
            }
#endif

#if ASSERT_ENABLED
            std::vector<u8> hit(g.get_n(), 0);
            for (const auto& e : matches) {
                hit[e.u] += 1;
                hit[e.v] += 1;

                if (hit[e.u] == 2) {
                    ASSERT(false);
                }
                if (hit[e.v] == 2) {
                    ASSERT(false);
                }
            }
#endif
        }
    };
}

#endif //HEIPROMAP_GLOBAL_PATH_ALGORITHM_H
