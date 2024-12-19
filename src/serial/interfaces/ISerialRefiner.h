#ifndef HEIDELBERGPROCESSMAPPING_ISERIALREFINER_H
#define HEIDELBERGPROCESSMAPPING_ISERIALREFINER_H

#include <vector>

#include "ISerialActiveVertexManager.h"
#include "ISerialBoundaryVertexManager.h"
#include "ISerialDistanceOracle.h"

namespace HeiProMap {
    template <typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle>
    class ISerialRefiner {
        static_assert(std::is_base_of_v<ISerialGraph, TSerialGraph>, "TSerialGraph must inherit from ISerialGraph");
        static_assert(std::is_base_of_v<ISerialActiveVertexManager, TSerialActiveVertexManager>, "TSerialActiveVertexManager must inherit from ISerialActiveVertexManager");
        static_assert(std::is_base_of_v<ISerialBoundaryVertexManager<TSerialGraph, TSerialActiveVertexManager, TSerialPartitionManager>, TSerialBoundaryVertexManager>, "TSerialBoundaryVertexManager must inherit from ISerialBoundaryVertexManager");
        static_assert(std::is_base_of_v<ISerialPartitionManager, TSerialPartitionManager>, "TSerialPartitionManager must inherit from ISerialPartitionManager");
        static_assert(std::is_base_of_v<ISerialDistanceOracle, TSerialDistanceOracle>, "TSerialDistanceOracle must inherit from ISerialDistanceOracle");

    public:
        virtual ~ISerialRefiner() = default;
        // initialization
        virtual void initialize(vertex_t n,
                                std::vector<partition_t>& t_hierarchy,
                                std::vector<weight_t>& t_distance,
                                weight_t t_lmax) = 0;
    };
}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALREFINER_H
