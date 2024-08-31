#ifndef HEIDELBERGPROCESSMAPPING_ISERIALPARTITIONER_H
#define HEIDELBERGPROCESSMAPPING_ISERIALPARTITIONER_H

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../../interfaces/IPartitioner.h"
#include "ISerialGraph.h"
#include "ISerialActiveVertexManager.h"
#include "ISerialBoundaryVertexManager.h"

namespace HeiProMap {

    template<typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager>
    class ISerialPartitioner : public IPartitioner {
        static_assert(std::is_base_of<ISerialGraph, TSerialGraph>::value, "TSerialGraph must inherit from ISerialGraph");
        static_assert(std::is_base_of<ISerialActiveVertexManager<TSerialGraph>, TSerialActiveVertexManager>::value, "TSerialActiveVertexManager must inherit from ISerialActiveVertexManager");
        static_assert(std::is_base_of<ISerialBoundaryVertexManager<TSerialGraph, TSerialActiveVertexManager, TSerialPartitionManager>, TSerialBoundaryVertexManager>::value, "TSerialBoundaryVertexManager must inherit from ISerialBoundaryVertexManager");
        static_assert(std::is_base_of<ISerialPartitionManager<TSerialGraph, TSerialActiveVertexManager>, TSerialPartitionManager>::value, "TSerialPartitionManager must inherit from ISerialPartitionManager");
    public:
        // initialization
        virtual void initialize(TSerialGraph *t_p_g,
                                TSerialActiveVertexManager *t_p_av_manager,
                                TSerialBoundaryVertexManager *t_p_bv_manager,
                                TSerialPartitionManager *t_p_p_manager,
                                std::vector<partition_t> &t_hierarchy,
                                std::vector<weight_t> &t_distance,
                                f64 t_imbalance) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALPARTITIONER_H
