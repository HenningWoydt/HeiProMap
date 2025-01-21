#ifndef HEIPROMAP_GLOBAL_PATH_ALGORITHM_H
#define HEIPROMAP_GLOBAL_PATH_ALGORITHM_H

#include <vector>
#include <algorithm>
#include <numeric>
#include <deque>

#include "../../definitions.h"
#include "../../macros.h"
#include "../interfaces/ISerialMatcher.h"

namespace HeiProMap {
    class Path {
    public:
        std::vector<vertex_t> vertices;
        std::vector<f64>      weights;

        explicit Path(vertex_t v) { vertices.push_back(v); }

        size_t size() const { return weights.size(); }

        void add(vertex_t v, f64 w) {
            vertices.push_back(v);
            weights.push_back(w);
        }

        void add(Path &p, f64 weight) {
            ASSERT(&p != this);
            for (vertex_t v: p.vertices) {
                vertices.push_back(v);
            }
            weights.push_back(weight);
            for (f64 w: p.weights) {
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
        void solve(Path &p, std::vector<EdgeUV> &matches) {
            f64 matching_weight = 0.0;
            if (p.is_cycle()) {
                solve_cycle(p, matches);
            } else {
                solve_line(p, matches, matching_weight);
            }
        }

    private:
        std::vector<f64>  w;
        std::vector<u32>  m;
        std::vector<bool> take;

        void solve_cycle(Path &p, std::vector<EdgeUV> &matches) {
            // remove e1
            f64                 matching_weight_1 = 0.0;
            std::vector<EdgeUV> matches_1;
            Path                p1                = p.copy();
            p1.pop_front();
            solve_line(p1, matches_1, matching_weight_1);

            // remove e2
            f64                 matching_weight_2 = 0.0;
            std::vector<EdgeUV> matches_2;
            Path                p2                = p.copy();
            p2.pop_front();
            solve_line(p2, matches_2, matching_weight_2);

            if (matching_weight_1 > matching_weight_2) {
                for (const auto &e: matches_1) {
                    matches.emplace_back(e);
                }
            } else {
                for (const auto &e: matches_2) {
                    matches.emplace_back(e);
                }
            }
        }

        void solve_line(Path &p, std::vector<EdgeUV> &matches, f64 &matching_weight) {
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
                if (p.weights[1] + w[i - 2] > w[i - 1]) {
                    w[i]    = p.weights[1] + w[i - 2];
                    m[i]    = i - 2;
                    take[i] = true;
                } else {
                    w[i]    = w[i - 2];
                    m[i]    = i - 1;
                    take[i] = false;
                }
            }

            u32 idx = m[p.size() - 1];
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

    /**
     * Computes a matching based on the Global Path Algorithm from
     * > Jens Maue and Peter Sanders.
     * > Engineering Algorithms for Approximate Weighted Matching.
     * > Experimental Algorithms, 6th International Workshop, WEA 2007, Rome, Italy, June 6-8, 2007, Proceedings.
     */
    class GlobalPathAlgorithmMatcher final : public ISerialMatcher {
        std::vector<bool> is_endpoint;
        std::vector<u32>  path_id;

        weight_t l_max = 0;

        std::vector<EdgeUVW> edges;

        std::vector<Path> paths;

    public:
        GlobalPathAlgorithmMatcher() = default;

        void initialize(const size_t n, const weight_t t_l_max) override {
            is_endpoint.resize(n);
            path_id.resize(n);
            l_max = t_l_max;
        }

        template<typename TSerialGraph, typename TSerialActiveVertexManager>
        void match(TSerialGraph &g,
                   TSerialActiveVertexManager &av_manager,
                   std::vector<EdgeUV> &matches) {
            edges.clear();
            for (vertex_t u: av_manager) {
                for (const auto &[v, w]: g[u]) {
                    f64 edge_rating = (((f64) w) * ((f64) w)) / (g.get_weight(u) * g.get_weight(v));
                    edges.push_back({u, v, edge_rating});
                }
            }
            std::sort(edges.begin(), edges.end());

            std::fill(is_endpoint.begin(), is_endpoint.end(), true);
            std::iota(path_id.begin(), path_id.end(), 0);
            paths.clear();
            for (vertex_t u = 0; u < is_endpoint.size(); ++u) {
                paths.emplace_back(u);
            }

            for (const auto &[u, v, w]: edges) {
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
                } else if (path_id[u] == path_id[v]) {
                    continue; // not applicable
                }

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
            DPSolver  dp_solver;
            for (Path &p: paths) {
                if (!p.empty()) {
                    dp_solver.solve(p, matches);
                }
            }


#if ASSERT_ENABLED
            for (const EdgeUV &e: matches) {
                ASSERT(e.u != e.v);
                ASSERT(av_manager.is_active(e.u));
                ASSERT(av_manager.is_active(e.v));
            }
#endif

#if ASSERT_ENABLED
            std::vector<u8> hit(g.get_n(), 0);
            for (const auto &e: matches) {
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
