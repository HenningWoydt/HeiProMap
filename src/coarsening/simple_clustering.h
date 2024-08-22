#ifndef SERIALPROCESSMAPPING_SIMPLE_CLUSTERING_H
#define SERIALPROCESSMAPPING_SIMPLE_CLUSTERING_H

/*

#include "../datastructures/graph.h"
#include "../datastructures/iterators/active_vertex_iterator.h"

namespace HeiProMap {
    class SimpleClustering {
    private:
        Graph *p_g = nullptr;
        std::vector<s32> marker;

    public:
        SimpleClustering() = default;

        void initialize(Graph *t_g){
            p_g = t_g;
            marker.resize((*p_g).get_n(), -1);
        }

        void match(std::vector<Edge> &matches, s32 level) {
            Graph &g = (*p_g);
            ASSERT(marker.size() == g.get_n());

            matches.clear();

            for (ActiveVertexIterator avi(g); avi.not_end(); avi.next()) {
                vertex_t u = avi.get();
                ASSERT(g.get_vertex_state(u) == 1);

                bool can_cluster = false;
                bool deg_1_neighbour = false;
                if (marker[u] != level) {
                    can_cluster = true;
                    deg_1_neighbour |= (g[u].size() == 1);
                    for (EdgeW &e: g[u]) {
                        ASSERT(u != e.v);
                        ASSERT(g.get_vertex_state(e.v) == 1);
                        deg_1_neighbour |= (g[u].size() == 1);
                        if (marker[e.v] == level) {
                            can_cluster = false;
                            break;
                        }
                    }
                }

                if(can_cluster && deg_1_neighbour){
                    marker[u] = level;
                    for (EdgeW &e: g[u]) {
                        marker[e.v] = level;
                        matches.emplace_back(u, e.v);
                    }
                }
            }

            for (ActiveVertexIterator avi(g); avi.not_end(); avi.next()) {
                vertex_t u = avi.get();
                ASSERT(g.get_vertex_state(u) == 1);

                bool can_cluster = false;
                if (marker[u] != level) {
                    can_cluster = true;
                    for (EdgeW &e: g[u]) {
                        ASSERT(u != e.v);
                        ASSERT(g.get_vertex_state(e.v) == 1);
                        if (marker[e.v] == level) {
                            can_cluster = false;
                            break;
                        }
                    }
                }

                if(can_cluster){
                    marker[u] = level;
                    for (EdgeW &e: g[u]) {
                        marker[e.v] = level;
                        matches.emplace_back(u, e.v);
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

#endif //SERIALPROCESSMAPPING_SIMPLE_CLUSTERING_H
