#ifndef HEIDELBERGPROCESSMAPPING_IPARALLELQUOTIENTGRAPH_H
#define HEIDELBERGPROCESSMAPPING_IPARALLELQUOTIENTGRAPH_H

#include "../../interfaces/IQuotientGraph.h"

namespace HeiProMap {

    class IParallelQuotientGraph : public IQuotientGraph {
    public:
        virtual void initialize(partition_t k,
                                u64 n_threads) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IPARALLELQUOTIENTGRAPH_H
