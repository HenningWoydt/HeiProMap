#ifndef HEIDELBERGPROCESSMAPPING_GREEDY_EDGE_MATCHER_H
#define HEIDELBERGPROCESSMAPPING_GREEDY_EDGE_MATCHER_H

#include <algorithm>
#include <vector>

#include "../../definitions.h"
#include "../../macros.h"
#include "../interfaces/ISerialMatcher.h"

namespace HeiProMap {
    template <typename TSerialGraph, typename TSerialActiveVertexManager>
    class GreedyEdgeMatcher final : public ISerialMatcher<TSerialGraph, TSerialActiveVertexManager> {
        u32 mark = 0;
        std::vector<u32> used;

    public:
        GreedyEdgeMatcher() = default;

        void initialize(const size_t n) override {
            mark = 0;
            used.resize(n, 0);
        }

        void match(const TSerialGraph& g,
                   TSerialActiveVertexManager& av_manager,
                   std::vector<EdgeUV>& matches) override {
            std::vector<EdgeUVW> edges;
            edges.reserve(g.get_m());
            for (vertex_t u : av_manager) {
                ASSERT(av_manager.is_active(u));

                for (size_t i = 0; i < g.size(u); ++i) {
                    const vertex_t v  = g.neighbor(u, i);
                    const weight_t ew = g.get_weight(u, i);
                    const f64 rating  = (f64)ew / (f64)(g.size(u) * g.size(v));
                    edges.emplace_back(u, v, rating);
                }
            }
            std::sort(edges.begin(), edges.end(), std::greater<>());

            mark += 1;
            matches.clear();
            for (const auto& e : edges) {
                if (used[e.u] != mark && used[e.v] != mark) {
                    used[e.u] = mark;
                    used[e.v] = mark;
                    if (g.size(e.u) > g.size(e.v)) {
                        matches.emplace_back(e.u, e.v);
                    } else {
                        matches.emplace_back(e.v, e.u);
                    }
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
            for (auto& e : matches) {
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

#endif //HEIDELBERGPROCESSMAPPING_GREEDY_EDGE_MATCHER_H
