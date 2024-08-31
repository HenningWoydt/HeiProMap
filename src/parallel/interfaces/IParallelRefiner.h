#ifndef HEIDELBERGPROCESSMAPPING_IPARALLELREFINER_H
#define HEIDELBERGPROCESSMAPPING_IPARALLELREFINER_H

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../../interfaces/IRefiner.h"
#include "IParallelActiveVertexManager.h"
#include "IParallelBoundaryVertexManager.h"
#include "IParallelDistanceOracle.h"
#include "IParallelQuotientGraph.h"
#include "../datastructures/parallel_static_csr_graph.h"
#include "../datastructures/parallel_distance_oracle.h"

namespace HeiProMap {

    template<typename TParallelGraph,
             typename TParallelActiveVertexManager,
             typename TParallelBoundaryVertexManager,
             typename TParallelPartitionManager,
             typename TParallelDistanceOracle,
             typename TParallelQuotientGraph>
    class IParallelRefiner : public IRefiner {
        static_assert(std::is_base_of<IParallelGraph, TParallelGraph>::value, "TParallelGraph must inherit from IParallelGraph");
        static_assert(std::is_base_of<IParallelActiveVertexManager<TParallelGraph>, TParallelActiveVertexManager>::value, "TParallelActiveVertexManager must inherit from IParallelActiveVertexManager");
        static_assert(std::is_base_of<IParallelBoundaryVertexManager<TParallelGraph, TParallelActiveVertexManager, TParallelPartitionManager>, TParallelBoundaryVertexManager>::value, "TParallelBoundaryVertexManager must inherit from IParallelBoundaryVertexManager");
        static_assert(std::is_base_of<IParallelPartitionManager<TParallelGraph, TParallelActiveVertexManager>, TParallelPartitionManager>::value, "TParallelPartitionManager must inherit from IParallelPartitionManager");
        static_assert(std::is_base_of<IParallelDistanceOracle, TParallelDistanceOracle>::value, "TParallelDistanceOracle must inherit from IParallelDistanceOracle");
        static_assert(std::is_base_of<IParallelQuotientGraph<TParallelGraph, TParallelActiveVertexManager, TParallelPartitionManager, TParallelDistanceOracle>, TParallelQuotientGraph>::value, "TParallelQuotientGraph must inherit from IParallelQuotientGraph");
    public:
        // initialization
        virtual void initialize(TParallelGraph *t_p_g,
                                TParallelActiveVertexManager *t_p_av_manager,
                                TParallelBoundaryVertexManager *t_p_bv_manager,
                                TParallelPartitionManager *t_p_p_manger,
                                TParallelDistanceOracle *t_p_d_oracle,
                                TParallelQuotientGraph *t_p_qgraph,
                                std::vector<partition_t> &t_hierarchy,
                                std::vector<weight_t> &t_distance,
                                weight_t t_lmax,
                                u64 n_threads) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IPARALLELREFINER_H
