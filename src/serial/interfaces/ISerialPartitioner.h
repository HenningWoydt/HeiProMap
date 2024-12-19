#ifndef HEIDELBERGPROCESSMAPPING_ISERIALPARTITIONER_H
#define HEIDELBERGPROCESSMAPPING_ISERIALPARTITIONER_H

#include <vector>

#include "ISerialActiveVertexManager.h"
#include "ISerialBoundaryVertexManager.h"
#include "ISerialGraph.h"

namespace HeiProMap {
    template <typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialPartitionManager>
    class ISerialPartitioner {
        static_assert(std::is_base_of_v<ISerialGraph, TSerialGraph>, "TSerialGraph must inherit from ISerialGraph");
        static_assert(std::is_base_of_v<ISerialActiveVertexManager, TSerialActiveVertexManager>, "TSerialActiveVertexManager must inherit from ISerialActiveVertexManager");
        static_assert(std::is_base_of_v<ISerialPartitionManager, TSerialPartitionManager>, "TSerialPartitionManager must inherit from ISerialPartitionManager");

    public:
        virtual ~ISerialPartitioner() = default;
        virtual void partition(TSerialGraph& g,
                               TSerialActiveVertexManager& av_manager,
                               TSerialPartitionManager& p_manager,
                               const std::vector<partition_t>& hierarchy,
                               const std::vector<weight_t>& distance,
                               f64 t_imbalance) = 0;
    };
}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALPARTITIONER_H
