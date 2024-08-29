#ifndef SERIALPROCESSMAPPING_SIMPLE_EDGE_MATCHER_H
#define SERIALPROCESSMAPPING_SIMPLE_EDGE_MATCHER_H

#include "../datastructures/graph.h"

/*

namespace HeiProMap {
    class SimpleEdgeMatcher {
    public:
        void match(Graph &g, std::vector<Edge> &matches, std::vector<s32> &marker, s32 level) {
            ASSERT(marker.size() == g.get_n());

            matches.clear();

            for (ActiveVertexIterator avi(g); avi.not_end(); avi.next()) {
                vertex_t u = avi.get();
                ASSERT(g.get_vertex_state(u) == 1);

                if (marker[u] != level) {
                    for (EdgeW &e: g[u]) {
                        ASSERT(u != e.v);
                        ASSERT(g.get_vertex_state(e.v) == 1);
                        if (marker[e.v] != level) {
                            marker[u] = level;
                            marker[e.v] = level;

                            if (g[u].size() > g[e.v].size()) {
                                matches.emplace_back(u, e.v);
                            } else {
                                matches.emplace_back(e.v, u);
                            }

                            break;
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

 */

#endif //SERIALPROCESSMAPPING_SIMPLE_EDGE_MATCHER_H
