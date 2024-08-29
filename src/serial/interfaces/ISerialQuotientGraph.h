#ifndef HEIDELBERGPROCESSMAPPING_ISERIALQUOTIENTGRAPH_H
#define HEIDELBERGPROCESSMAPPING_ISERIALQUOTIENTGRAPH_H

#include "../../interfaces/IQuotientGraph.h"

namespace HeiProMap {

    class ISerialQuotientGraph : public IQuotientGraph {
    public:
        virtual void initialize(partition_t k) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALQUOTIENTGRAPH_H
