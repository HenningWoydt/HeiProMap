#ifndef HEIDELBERGPROCESSMAPPING_GREEDY_EDGE_MATCHER_H
#define HEIDELBERGPROCESSMAPPING_GREEDY_EDGE_MATCHER_H

#include "../interfaces/ISerialMatcher.h"

namespace HeiProMap {

    template<typename TSerialGraph, typename TSerialActiveVertexManager>
    class GreedyEdgeMatcher final : public ISerialMatcher<TSerialGraph, TSerialActiveVertexManager> {
        TSerialGraph *p_g = nullptr;
        TSerialActiveVertexManager *p_av_manager = nullptr;

        u32 mark = 0;
        std::vector<u32> used;

    public:
        GreedyEdgeMatcher() = default;

        void initialize(TSerialGraph *t_p_g,
                        TSerialActiveVertexManager *t_p_av_manager) override {
            ASSERT(t_p_g != nullptr);
            ASSERT(t_p_av_manager != nullptr);

            p_g = t_p_g;
            p_av_manager = t_p_av_manager;

            mark = 0;
            used.resize(p_g->get_n(), 0);
        }

        void match(std::vector<EdgeUV> &matches) override {
            ASSERT(p_g != nullptr);
            ASSERT(p_av_manager != nullptr);
            ASSERT(used.size() == p_g->get_n());

            std::vector<EdgeUVW> edges;
            edges.reserve(p_g->get_m());
            for (p_av_manager->reset_iterator(); p_av_manager->available(); p_av_manager->next()) {
                const vertex_t u = p_av_manager->get();
                ASSERT(p_av_manager->is_active(u));

                for (size_t i = 0; i < p_g->size(u); ++i) {
                    const vertex_t v = p_g->neighbor(u, i);
                    const weight_t ew = p_g->get_weight(u, i);
                    const f64 rating = (f64) ew / (f64) (p_g->size(u) * p_g->size(v));
                    edges.emplace_back(u, v, rating);
                }
            }
            std::sort(edges.begin(), edges.end(), std::greater<>());

            mark += 1;
            matches.clear();
            for (const auto &e: edges) {
                if (used[e.u] != mark && used[e.v] != mark) {
                    used[e.u] = mark;
                    used[e.v] = mark;
                    if (p_g->size(e.u) > p_g->size(e.v)) {
                        matches.emplace_back(e.u, e.v);
                    } else {
                        matches.emplace_back(e.v, e.u);
                    }
                }
            }

#if ASSERT_ENABLED
            for (const EdgeUV &e: matches) {
                ASSERT(e.u != e.v);
                ASSERT(p_av_manager->is_active(e.u));
                ASSERT(p_av_manager->is_active(e.v));
            }
#endif

#if ASSERT_ENABLED
            std::vector<u8> hit(p_g->get_n(), 0);
            for (auto & e : matches) {
                hit[e.u] += 1;
                hit[e.v] += 1;

                if(hit[e.u] == 2){
                    ASSERT(false);
                }
                if(hit[e.v] == 2){
                    ASSERT(false);
                }
            }
#endif

        }
    };
}

#endif //HEIDELBERGPROCESSMAPPING_GREEDY_EDGE_MATCHER_H
