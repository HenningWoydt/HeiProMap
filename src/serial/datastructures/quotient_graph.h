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

        void add_edge(partition_t u, partition_t v, weight_t w) {
            partition_t min = std::min(u, v);
            partition_t max = std::max(u, v);
            m_adj_mtx[min * k + max] += w;
        }

        void remove_edge(partition_t u, partition_t v, weight_t w) {
            partition_t min = std::min(u, v);
            partition_t max = std::max(u, v);
            m_adj_mtx[min * k + max] -= w;
        }

        bool has_edge(partition_t u, partition_t v) {
            partition_t min = std::min(u, v);
            partition_t max = std::max(u, v);
            return m_adj_mtx[min * k + max] > 0;
        }
    };
}

#endif //HEIDELBERGPROCESSMAPPING_QUOTIENT_GRAPH_H
