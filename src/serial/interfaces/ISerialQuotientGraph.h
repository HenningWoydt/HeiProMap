#ifndef HEIDELBERGPROCESSMAPPING_ISERIALQUOTIENTGRAPH_H
#define HEIDELBERGPROCESSMAPPING_ISERIALQUOTIENTGRAPH_H

#include "../../definitions.h"

namespace HeiProMap {
    class ISerialQuotientGraph {
    public:
        virtual ~ISerialQuotientGraph() = default;
        virtual void initialize(partition_t k) = 0;
    };
}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALQUOTIENTGRAPH_H
