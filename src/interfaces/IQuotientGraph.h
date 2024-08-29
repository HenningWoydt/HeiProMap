#ifndef HEIDELBERGPROCESSMAPPING_IQUOTIENTGRAPH_H
#define HEIDELBERGPROCESSMAPPING_IQUOTIENTGRAPH_H

#include "../definitions.h"
#include "../macros.h"

namespace HeiProMap {

    class IQuotientGraph {
    public:
        virtual void add_edge(partition_t u, partition_t v, weight_t w) = 0;
        virtual void remove_edge(partition_t u, partition_t v, weight_t w) = 0;
        virtual bool has_edge(partition_t u, partition_t v) = 0;

        virtual void get_matching_hierarchy(std::vector<std::vector<QGraphUV>> &matching_hierarchy) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IQUOTIENTGRAPH_H
