#ifndef SERIALPROCESSMAPPING_HEAVY_EDGE_MATCHER_H
#define SERIALPROCESSMAPPING_HEAVY_EDGE_MATCHER_H

#include "../interfaces/IGraph.h"
#include "../interfaces/IMatcher.h"

namespace HeiProMap {

    class HeavyEdgeMatcher : public IMatcher {
    private:
        IGraph *m_p_g = nullptr;
        IActiveVertexManager *m_p_av_manager = nullptr;

        u32 m_mark = 0;
        std::vector<u32> m_used;

    public:
        HeavyEdgeMatcher() = default;

        void initialize(IGraph *t_p_g,
                        IActiveVertexManager *t_p_av_manager) final {
            ASSERT(t_p_g != nullptr);
            ASSERT(t_p_av_manager != nullptr);

            m_p_g = t_p_g;
            m_p_av_manager = t_p_av_manager;

            m_mark = 0;
            m_used.resize(m_p_g->get_n(), 0);
        }

        void match(std::vector<Edge> &matches) final {
            ASSERT(m_p_g != nullptr);
            ASSERT(m_p_av_manager != nullptr);
            ASSERT(m_used.size() == m_p_g->get_n());

            m_mark += 1;
            matches.clear();

            // first check vertices with degree 1
            for (m_p_av_manager->reset_iterator(); m_p_av_manager->available(); m_p_av_manager->next()) {
                vertex_t u = m_p_av_manager->get();
                ASSERT(m_p_av_manager->is_active(u));

                if (m_used[u] != m_mark && m_p_g->n_neighbors(u) == 1) {
                    const EdgeW &e = m_p_g->neighbor(u, 0);
                    if (m_used[e.v] != m_mark) {
                        m_used[u] = m_mark;
                        m_used[e.v] = m_mark;

                        matches.emplace_back(e.v, u); // pull u into v
                    }
                }
            }

            // check all other vertices
            for (m_p_av_manager->reset_iterator(); m_p_av_manager->available(); m_p_av_manager->next()) {
                vertex_t u = m_p_av_manager->get();
                ASSERT(m_p_av_manager->is_active(u));

                if (m_used[u] != m_mark) {
                    size_t best_idx;
                    weight_t max_weight = 0;

                    for (size_t i = 0; i < m_p_g->n_neighbors(u); ++i) {
                        const EdgeW &e = m_p_g->neighbor(u, i);
                        ASSERT(u != e.v);
                        ASSERT(m_p_av_manager->is_active(e.v));
                        if (m_used[e.v] != m_mark) {
                            if (e.w > max_weight) {
                                best_idx = i;
                                max_weight = e.w;
                            }
                        }
                    }

                    if (max_weight != 0) {
                        const EdgeW &e = m_p_g->neighbor(u, best_idx);
                        m_used[u] = m_mark;
                        m_used[e.v] = m_mark;

                        if (m_p_g->n_neighbors(u) > m_p_g->n_neighbors(e.v)) {
                            matches.emplace_back(u, e.v);
                        } else {
                            matches.emplace_back(e.v, u);
                        }
                    }
                }
            }

#if ASSERT_ENABLED
            for (const Edge &e: matches) {
                ASSERT(e.u != e.v);
                ASSERT(m_p_av_manager->is_active(e.u));
                ASSERT(m_p_av_manager->is_active(e.v));
            }
#endif

#if ASSERT_ENABLED
            std::vector<u8> hit(m_p_g->get_n(), 0);
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

#endif //SERIALPROCESSMAPPING_HEAVY_EDGE_MATCHER_H
