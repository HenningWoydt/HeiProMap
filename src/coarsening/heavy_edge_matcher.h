#ifndef SERIALPROCESSMAPPING_HEAVY_EDGE_MATCHER_H
#define SERIALPROCESSMAPPING_HEAVY_EDGE_MATCHER_H

#include "../datastructures/graph.h"
#include "../datastructures/iterators/active_vertex_iterator.h"

namespace SPM {

    class HeavyEdgeMatcher {
    private:
        Graph *p_g = nullptr;
        std::vector<s32> marker;

    public:
        HeavyEdgeMatcher() = default;

        void initialize(Graph *t_g){
            p_g = t_g;
            marker.resize((*p_g).get_n(), -1);
        }

        void match(std::vector<Edge> &matches, s32 level) {
            ASSERT(p_g != nullptr);
            Graph &g = *p_g;
            ASSERT(marker.size() == g.get_n());

            matches.clear();

            // first check degree 1 vertices
            for (ActiveVertexIterator avi(g); avi.not_end(); avi.next()) {
                vertex_t u = avi.get();
                ASSERT(g.get_vertex_state(u) == 1);
                if (marker[u] != level && g[u].size() == 1) {
                    if (marker[g[u][0].v] != level) {
                        marker[u] = level;
                        marker[g[u][0].v] = level;

                        matches.emplace_back(g[u][0].v, u); // pull u into v
                    }
                }
            }

            // check all other vertices
            for (ActiveVertexIterator avi(g); avi.not_end(); avi.next()) {
                vertex_t u = avi.get();
                ASSERT(g.get_vertex_state(u) == 1);

                if (marker[u] != level && g[u].size() > 1) {

                    bool found = false;
                    size_t best_idx;
                    weight_t max_weight;

                    for (size_t i = 0; i < g[u].size(); ++i) {
                        const EdgeW &e = g[u][i];
                        ASSERT(u != e.v);
                        ASSERT(g.get_vertex_state(e.v) == 1);
                        if (marker[e.v] != level) {
                            if (!found || e.w > max_weight) {
                                best_idx = i;
                                max_weight = e.w;
                                found = true;
                            }
                        }
                    }

                    if (found) {
                        const EdgeW &e = g[u][best_idx];
                        marker[u] = level;
                        marker[e.v] = level;

                        if (g[u].size() > g[e.v].size()) {
                            matches.emplace_back(u, e.v);
                        } else {
                            matches.emplace_back(e.v, u);
                        }
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
                    // ASSERT(not_used_twice);
                }
            }

        }
    };
}

#endif //SERIALPROCESSMAPPING_HEAVY_EDGE_MATCHER_H
