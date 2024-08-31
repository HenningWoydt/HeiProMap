#ifndef HEIDELBERGPROCESSMAPPING_IPARALLELQUOTIENTGRAPH_H
#define HEIDELBERGPROCESSMAPPING_IPARALLELQUOTIENTGRAPH_H

#include "../../interfaces/IQuotientGraph.h"

namespace HeiProMap {

    template<typename TParallelGraph,
             typename TParallelActiveVertexManager,
             typename TParallelPartitionManager,
             typename TParallelDistanceOracle>
    class IParallelQuotientGraph : public IQuotientGraph {
        static_assert(std::is_base_of<IParallelGraph, TParallelGraph>::value, "TParallelGraph must inherit from IParallelGraph");
        static_assert(std::is_base_of<IParallelPartitionManager<TParallelGraph, TParallelActiveVertexManager>, TParallelPartitionManager>::value, "TParallelPartitionManager must inherit from IParallelPartitionManager");
        static_assert(std::is_base_of<IParallelDistanceOracle, TParallelDistanceOracle>::value, "TParallelDistanceOracle must inherit from IParallelDistanceOracle");
    public:
        virtual void initialize(TParallelGraph *t_p_g,
                                TParallelPartitionManager *t_p_p_manager,
                                TParallelDistanceOracle *t_p_d_oracle,
                                partition_t k,
                                u64 n_threads) = 0;

        virtual void move(vertex_t u, partition_t old_id, partition_t new_id) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IPARALLELQUOTIENTGRAPH_H
