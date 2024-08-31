#ifndef HEIDELBERGPROCESSMAPPING_IPARALLELPARTITIONER_H
#define HEIDELBERGPROCESSMAPPING_IPARALLELPARTITIONER_H

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../../interfaces/IPartitioner.h"
#include "IParallelGraph.h"
#include "IParallelActiveVertexManager.h"
#include "IParallelBoundaryVertexManager.h"

namespace HeiProMap {

    template<typename TParallelGraph, typename TParallelActiveVertexManager, typename TParallelBoundaryVertexManager, typename TParallelPartitionManager>
    class IParallelPartitioner : public IPartitioner {
        static_assert(std::is_base_of<IParallelGraph, TParallelGraph>::value, "TParallelGraph must inherit from IParallelGraph");
        static_assert(std::is_base_of<IParallelActiveVertexManager<TParallelGraph>, TParallelActiveVertexManager>::value, "TParallelActiveVertexManager must inherit from IParallelActiveVertexManager");
        static_assert(std::is_base_of<IParallelBoundaryVertexManager<TParallelGraph, TParallelActiveVertexManager, TParallelPartitionManager>, TParallelBoundaryVertexManager>::value, "TParallelBoundaryVertexManager must inherit from IParallelBoundaryVertexManager");
        static_assert(std::is_base_of<IParallelPartitionManager<TParallelGraph, TParallelActiveVertexManager>, TParallelPartitionManager>::value, "TParallelActiveVertexManager must inherit from IParallelActiveVertexManager");
    public:
        // initialization
        virtual void initialize(TParallelGraph *t_p_g,
                                TParallelActiveVertexManager *t_p_av_manager,
                                TParallelBoundaryVertexManager *t_p_bv_manager,
                                TParallelPartitionManager *t_p_p_manager,
                                std::vector<partition_t> &t_hierarchy,
                                std::vector<weight_t> &t_distance,
                                f64 t_imbalance,
                                u64 n_threads) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IPARALLELPARTITIONER_H
