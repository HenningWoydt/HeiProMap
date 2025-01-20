#ifndef HEIDELBERGPROCESSMAPPING_ISERIALREFINER_H
#define HEIDELBERGPROCESSMAPPING_ISERIALREFINER_H

#include <vector>

#include "ISerialActiveVertexManager.h"
#include "ISerialBoundaryVertexManager.h"
#include "ISerialDistanceOracle.h"

namespace HeiProMap {
    class ISerialRefiner {
    public:
        virtual ~ISerialRefiner() = default;
        // initialization
        virtual void initialize(vertex_t n,
                                std::vector<partition_t>& t_hierarchy,
                                std::vector<weight_t>& t_distance,
                                weight_t t_lmax) = 0;

        template <typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle>
        void refine(TSerialGraph& g,
                    TSerialActiveVertexManager& av_manager,
                    TSerialBoundaryVertexManager& bv_manager,
                    TSerialPartitionManager& p_manager,
                    TSerialDistanceOracle& d_oracle) {}
    };
}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALREFINER_H
