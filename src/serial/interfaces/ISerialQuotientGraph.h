#ifndef HEIDELBERGPROCESSMAPPING_ISERIALQUOTIENTGRAPH_H
#define HEIDELBERGPROCESSMAPPING_ISERIALQUOTIENTGRAPH_H

#include "../../definitions.h"

namespace HeiProMap {
    class ISerialQuotientGraph {
    public:
        virtual ~ISerialQuotientGraph() = default;
        virtual void initialize(partition_t k) = 0;
        virtual void add_edge(partition_t u, partition_t v, weight_t w) = 0;
        virtual void remove_edge(partition_t u, partition_t v, weight_t w)  = 0;
        virtual bool has_edge(partition_t u, partition_t v) = 0;

        template<typename TSerialGraph, typename TSerialPartitionManager>
        void move(TSerialGraph &g, TSerialPartitionManager &p_manager, vertex_t u, partition_t old_id, partition_t new_id) {}
    };
}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALQUOTIENTGRAPH_H
