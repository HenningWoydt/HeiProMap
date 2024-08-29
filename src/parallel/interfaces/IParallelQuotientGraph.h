#ifndef HEIDELBERGPROCESSMAPPING_IPARALLELQUOTIENTGRAPH_H
#define HEIDELBERGPROCESSMAPPING_IPARALLELQUOTIENTGRAPH_H

#include "../../interfaces/IQuotientGraph.h"

namespace HeiProMap {

    class IParallelQuotientGraph : public IQuotientGraph {
    public:
        virtual void initialize(IParallelGraph *t_p_g,
                                IParallelPartitionManager *t_p_p_manager,
                                IParallelDistanceOracle *t_p_d_oracle,
                                partition_t k,
                                u64 n_threads) = 0;

        virtual void move(vertex_t u, partition_t old_id, partition_t new_id) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IPARALLELQUOTIENTGRAPH_H
