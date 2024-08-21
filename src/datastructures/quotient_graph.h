#ifndef SERIALPROCESSMAPPING_QUOTIENT_GRAPH_H
#define SERIALPROCESSMAPPING_QUOTIENT_GRAPH_H

#include "../utility/definitions.h"
#include "../utility/macros.h"
#include "../utility/utils.h"
#include "graph.h"
#include "../utility/qap.h"
#include "iterators/active_vertex_iterator.h"
#include "distance_oracle.h"
#include "boundary_vertex_manager.h"

namespace SPM {

    class QuotientGraph {
    private:
        partition_t k = 0;

        std::vector<weight_t> m_adj_mtx;

    public:
        QuotientGraph() = default;

        void initialize(partition_t t_k) {
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

#endif //SERIALPROCESSMAPPING_QUOTIENT_GRAPH_H
