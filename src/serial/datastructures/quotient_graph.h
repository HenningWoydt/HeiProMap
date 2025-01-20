#ifndef HEIDELBERGPROCESSMAPPING_QUOTIENT_GRAPH_H
#define HEIDELBERGPROCESSMAPPING_QUOTIENT_GRAPH_H

#include "distance_oracle.h"
#include "../../definitions.h"
#include "../interfaces/ISerialQuotientGraph.h"
#include "../utility/utils.h"

namespace HeiProMap {
    class QuotientGraph final : public ISerialQuotientGraph {
    private:
        partition_t k = 0;

        std::vector<weight_t> m_adj_mtx;

    public:
        QuotientGraph() = default;

        void initialize(partition_t t_k) override {
            k = t_k;

            m_adj_mtx.resize(k * k, 0);
        }

        void add_edge(partition_t u, partition_t v, weight_t w) override {
            partition_t min = std::min(u, v);
            partition_t max = std::max(u, v);
            m_adj_mtx[min * k + max] += w;
        }

        void remove_edge(partition_t u, partition_t v, weight_t w) override {
            partition_t min = std::min(u, v);
            partition_t max = std::max(u, v);
            m_adj_mtx[min * k + max] -= w;
        }

        bool has_edge(partition_t u, partition_t v) override {
            partition_t min = std::min(u, v);
            partition_t max = std::max(u, v);
            return m_adj_mtx[min * k + max] > 0;
        }

        template<typename TSerialGraph, typename TSerialPartitionManager>
        void move(TSerialGraph &g, TSerialPartitionManager &p_manager, vertex_t u, partition_t old_id, partition_t new_id) {
            ASSERT(new_id < k);
            ASSERT(new_id != old_id);

            // first remove all edges from old_id
            for (size_t i = 0; i < g.size(u); ++i) {
                vertex_t    v    = g.neighbor(u, i);
                weight_t    w    = g.get_weight(u, i);
                partition_t v_id = p_manager[v];

                // remove old edge, if existed
                if(old_id != v_id) {
                    remove_edge(old_id, v_id, w);
                }
                // add new edge, if has to exist
                if(new_id != v_id){
                    add_edge(new_id, v_id, w);
                }
            }
        }
    };
}

#endif //HEIDELBERGPROCESSMAPPING_QUOTIENT_GRAPH_H
