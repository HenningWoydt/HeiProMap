#ifndef SERIALPROCESSMAPPING_GREEDY_EDGE_MATCHER_H
#define SERIALPROCESSMAPPING_GREEDY_EDGE_MATCHER_H

/*

#include "../datastructures/graph.h"
#include "../datastructures/iterators/active_vertex_iterator.h"

namespace HeiProMap {

    class GreedyEdgeMatcher {
    private:
        Graph *p_g = nullptr;
        std::vector<s32> marker;

        std::vector<EdgeUVW> edges;
        std::vector<EdgeUVW> edges_help;
        std::vector<weight_t> help;

    public:
        GreedyEdgeMatcher() = default;

        void initialize(Graph *t_g){
            p_g = t_g;
            marker.resize((*p_g).get_n(), -1);
        }

        void match(std::vector<Edge> &matches, s32 level) {
            Graph &g = *p_g;
            ASSERT(marker.size() == g.get_n());

            edges.clear();
            edges.reserve(g.get_m());

            for (ActiveVertexIterator avi(g); avi.not_end(); avi.next()) {
                vertex_t u = avi.get();
                ASSERT(g.get_vertex_state(u) == 1);

                for (EdgeW &e: g[u]) {
                    ASSERT(u != e.v);
                    if (u < e.v) {
                        f64 rating = (f64) e.w / (f64) (g[u].size() * g[e.v].size());

                        edges.emplace_back(u, e.v, rating);
                    }
                }
            }
            std::sort(edges.begin(), edges.end(), std::greater<>());

            // use counting sort
            // counting_sort(edges, edges_help, help, min_w, max_w);

            matches.clear();
            for (EdgeUVW &e: edges) {
                if (marker[e.u] != level && marker[e.v] != level) {
                    marker[e.u] = level;
                    marker[e.v] = level;
                    if (g[e.u].size() > g[e.v].size()) {
                        matches.emplace_back(e.u, e.v);
                    } else {
                        matches.emplace_back(e.v, e.u);
                    }
                }
            }

            for (const Edge &e: matches) {
                ASSERT(e.u != e.v);
                ASSERT(g.get_vertex_state(e.u) == 1);
                ASSERT(g.get_vertex_state(e.v) == 1);
            }

            for (u64 i = 0; i < matches.size(); ++i) {
                for (u64 j = i + 1; j < matches.size(); ++j) {
                    bool not_used_twice = true;
                    not_used_twice &= matches[i].u != matches[j].u;
                    not_used_twice &= matches[i].u != matches[j].v;
                    not_used_twice &= matches[i].v != matches[j].u;
                    not_used_twice &= matches[i].v != matches[j].v;
                    ASSERT(not_used_twice);
                }
            }

        }
    };

}

 */

#endif //SERIALPROCESSMAPPING_GREEDY_EDGE_MATCHER_H
