#ifndef HEIDELBERGPROCESSMAPPING_ISERIALPARTITIONER_H
#define HEIDELBERGPROCESSMAPPING_ISERIALPARTITIONER_H

#include <vector>

#include "ISerialActiveVertexManager.h"
#include "ISerialBoundaryVertexManager.h"
#include "ISerialGraph.h"

namespace HeiProMap {
    class ISerialPartitioner {
    public:
        virtual ~ISerialPartitioner() = default;

        template<typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialPartitionManager>
        void partition(TSerialGraph &g,
                       TSerialActiveVertexManager &av_manager,
                       TSerialPartitionManager &p_manager,
                       const std::vector<partition_t> &hierarchy,
                       const std::vector<weight_t> &distance,
                       f64 t_imbalance) {}
    };
}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALPARTITIONER_H
