#ifndef HEIDELBERGPROCESSMAPPING_HEAVY_EDGE_MATCHER_H
#define HEIDELBERGPROCESSMAPPING_HEAVY_EDGE_MATCHER_H

#include "../interfaces/ISerialMatcher.h"

namespace HeiProMap {

    template<typename TSerialGraph, typename TSerialActiveVertexManager>
    class HeavyEdgeMatcher final : public ISerialMatcher<TSerialGraph, TSerialActiveVertexManager> {
        TSerialGraph *p_g = nullptr;
        TSerialActiveVertexManager *p_av_manager = nullptr;

        u32 mark = 0;
        std::vector<u32> used;

    public:
        HeavyEdgeMatcher() = default;

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

            mark += 1;
            matches.clear();

            // first check vertices with degree 1
            for (p_av_manager->reset_iterator(); p_av_manager->available(); p_av_manager->next()) {
                vertex_t u = p_av_manager->get();
                ASSERT(p_av_manager->is_active(u));

                if (used[u] != mark && p_g->size(u) == 1) {
                    vertex_t v = p_g->neighbor(u, 0);
                    if (used[v] != mark) {
                        used[u] = mark;
                        used[v] = mark;

                        matches.emplace_back(v, u); // pull u into v
                    }
                }
            }

            // check all other vertices
            for (p_av_manager->reset_iterator(); p_av_manager->available(); p_av_manager->next()) {
                vertex_t u = p_av_manager->get();
                ASSERT(p_av_manager->is_active(u));

                if (used[u] != mark) {
                    size_t best_idx;
                    weight_t max_weight = 0;

                    for (size_t i = 0; i < p_g->size(u); ++i) {
                        vertex_t v = p_g->neighbor(u, i);
                        weight_t ew = p_g->get_weight(u, i);
                        ASSERT(u != v);
                        ASSERT(p_av_manager->is_active(v));
                        if (used[v] != mark) {
                            if (ew > max_weight) {
                                best_idx = i;
                                max_weight = ew;
                            }
                        }
                    }

                    if (max_weight != 0) {
                        vertex_t v = p_g->neighbor(u, best_idx);
                        used[u] = mark;
                        used[v] = mark;

                        if (p_g->size(u) > p_g->size(v)) {
                            matches.emplace_back(u, v);
                        } else {
                            matches.emplace_back(v, u);
                        }
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

#endif //HEIDELBERGPROCESSMAPPING_HEAVY_EDGE_MATCHER_H
