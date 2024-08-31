#ifndef HEIDELBERGPROCESSMAPPING_PARALLEL_QUOTIENT_GRAPH_H
#define HEIDELBERGPROCESSMAPPING_PARALLEL_QUOTIENT_GRAPH_H


#include "../interfaces/IParallelQuotientGraph.h"

namespace HeiProMap {

    template<typename TParallelGraph,
             typename TParallelActiveVertexManager,
             typename TParallelPartitionManager,
             typename TParallelDistanceOracle>
    class ParallelQuotientGraph : public IParallelQuotientGraph<TParallelGraph, TParallelActiveVertexManager, TParallelPartitionManager, TParallelDistanceOracle> {
    private:
        TParallelGraph *m_p_g = nullptr;
        TParallelPartitionManager *m_p_p_manager = nullptr;
        TParallelDistanceOracle *m_p_d_oracle = nullptr;

        partition_t m_k = 0;
        u64 m_n_threads = 1;

        std::vector<weight_t> m_adj_mtx;

        std::vector<s32> m_used;
        std::vector<s32> m_used_edge;
        s32 m_mark = -1;

    public:
        void initialize(TParallelGraph *t_p_g,
                        TParallelPartitionManager *t_p_p_manager,
                        TParallelDistanceOracle *t_p_d_oracle,
                        partition_t t_k,
                        u64 t_n_threads) final {
            m_p_g = t_p_g;
            m_p_p_manager = t_p_p_manager;
            m_p_d_oracle = t_p_d_oracle;

            m_k = t_k;
            m_n_threads = t_n_threads;

            m_adj_mtx.resize(m_k * m_k, 0);

            m_used.resize(m_k, -1);
            m_used_edge.resize(m_k * m_k, -1);
        }

        void add_edge(partition_t u, partition_t v, weight_t w) final {
            partition_t min = std::min(u, v);
            partition_t max = std::max(u, v);
            m_adj_mtx[min * m_k + max] += w;
        }

        void remove_edge(partition_t u, partition_t v, weight_t w) final {
            partition_t min = std::min(u, v);
            partition_t max = std::max(u, v);
            m_adj_mtx[min * m_k + max] -= w;
        }

        bool has_edge(partition_t u, partition_t v) final {
            partition_t min = std::min(u, v);
            partition_t max = std::max(u, v);
            return m_adj_mtx[min * m_k + max] > 0;
        }

        void move(vertex_t u, partition_t old_id, partition_t new_id) final {
            for(size_t i = 0; i < m_p_g->size(u); ++i){
                vertex_t v = m_p_g->neighbor(u, i);
                partition_t v_id = (*m_p_p_manager)[v];
                weight_t w = m_p_g->get_weight(u, i);

                m_adj_mtx[std::min(old_id, v_id) * m_k + std::max(old_id, v_id)] -= w;
                m_adj_mtx[std::min(new_id, v_id) * m_k + std::max(new_id, v_id)] += w;
            }
        }

        void get_matching_hierarchy(std::vector<std::vector<QGraphUV>> &matching_hierarchy) final {
            matching_hierarchy.clear();

            std::fill(m_used.begin(), m_used.end(), 0);
            std::fill(m_used_edge.begin(), m_used_edge.end(), 0);
            m_mark = 1;

            bool found_match = true;
            while(found_match){
                m_mark += 1;
                matching_hierarchy.emplace_back();
                found_match = false;

                for(partition_t id1 = 0; id1 < m_k; ++id1){
                    for(partition_t id2 = id1 + 1; id2 < m_k; ++id2){
                        if(m_used[id2] == m_mark || m_used[id1] == m_mark || m_used_edge[id1 * m_k + id2] != 0) { continue; }

                        if(has_edge(id1, id2)){
                            matching_hierarchy.back().emplace_back(id1, id2);

                            // mark the edge as used
                            m_used_edge[id1 * m_k + id2] = m_mark;

                            // mark the neighborhood of id1 and id2 as used
                            for(partition_t id_temp = 0; id_temp < m_k; ++id_temp){
                                if(has_edge(id1, id_temp)){ m_used[id_temp] = m_mark; }
                            }
                            for(partition_t id_temp = 0; id_temp < m_k; ++id_temp){
                                if(has_edge(id2, id_temp)){ m_used[id_temp] = m_mark; }
                            }
                            found_match = true;
                            break;
                        }
                    }
                }
            }

            matching_hierarchy.pop_back();
        }
    };

}


#endif //HEIDELBERGPROCESSMAPPING_PARALLEL_QUOTIENT_GRAPH_H
