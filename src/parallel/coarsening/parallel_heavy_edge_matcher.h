#ifndef HEIDELBERGPROCESSMAPPING_PARALLEL_HEAVY_EDGE_MATCHER_H
#define HEIDELBERGPROCESSMAPPING_PARALLEL_HEAVY_EDGE_MATCHER_H

#include "../interfaces/IParallelMatcher.h"

namespace HeiProMap {

    class ParallelHeavyEdgeMatcher : public IParallelMatcher {
    private:
        IParallelGraph *m_p_g = nullptr;
        IParallelActiveVertexManager *m_p_av_manager = nullptr;

        u64 m_n_threads = 1;

        u32 m_mark = 0;
        std::vector<u32> m_used;

    public:
        ParallelHeavyEdgeMatcher() = default;

        void initialize(IParallelGraph *t_p_g,
                        IParallelActiveVertexManager *t_p_av_manager,
                        u64 n_threads) final {
            ASSERT(t_p_g != nullptr);
            ASSERT(t_p_av_manager != nullptr);

            m_p_g = t_p_g;
            m_p_av_manager = t_p_av_manager;

            m_n_threads = n_threads;

            m_mark = 0;
            m_used.resize(m_p_g->get_n(), 0);
        }

        void match(std::vector<EdgeUV> &matches) final {
            ASSERT(m_p_g != nullptr);
            ASSERT(m_p_av_manager != nullptr);
            ASSERT(m_used.size() == m_p_g->get_n());

            m_mark += 1;
            matches.clear();

            // first check vertices with degree 1
            for (m_p_av_manager->reset_iterator(); m_p_av_manager->available(); m_p_av_manager->next()) {
                vertex_t u = m_p_av_manager->get();
                ASSERT(m_p_av_manager->is_active(u));

                if (m_used[u] != m_mark && m_p_g->size(u) == 1) {
                    vertex_t v = m_p_g->neighbor(u, 0);
                    if (m_used[v] != m_mark) {
                        m_used[u] = m_mark;
                        m_used[v] = m_mark;

                        matches.emplace_back(v, u); // pull u into v
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

                    for (size_t i = 0; i < m_p_g->size(u); ++i) {
                        vertex_t v = m_p_g->neighbor(u, i);
                        weight_t ew = m_p_g->get_weight(u, i);
                        ASSERT(u != v);
                        ASSERT(m_p_av_manager->is_active(v));
                        if (m_used[v] != m_mark) {
                            if (ew > max_weight) {
                                best_idx = i;
                                max_weight = ew;
                            }
                        }
                    }

                    if (max_weight != 0) {
                        vertex_t v = m_p_g->neighbor(u, best_idx);
                        m_used[u] = m_mark;
                        m_used[v] = m_mark;

                        if (m_p_g->size(u) > m_p_g->size(v)) {
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

#endif //HEIDELBERGPROCESSMAPPING_PARALLEL_HEAVY_EDGE_MATCHER_H
