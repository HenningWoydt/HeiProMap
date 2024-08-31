#ifndef HEIDELBERGPROCESSMAPPING_ISERIALREFINER_H
#define HEIDELBERGPROCESSMAPPING_ISERIALREFINER_H

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../../interfaces/IRefiner.h"
#include "ISerialActiveVertexManager.h"
#include "ISerialBoundaryVertexManager.h"
#include "ISerialDistanceOracle.h"

namespace HeiProMap {

    template<typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle>
    class ISerialRefiner : public IRefiner {
        static_assert(std::is_base_of<ISerialGraph, TSerialGraph>::value, "TSerialGraph must inherit from ISerialGraph");
        static_assert(std::is_base_of<ISerialActiveVertexManager<TSerialGraph>, TSerialActiveVertexManager>::value, "TSerialActiveVertexManager must inherit from ISerialActiveVertexManager");
        static_assert(std::is_base_of<ISerialBoundaryVertexManager<TSerialGraph, TSerialActiveVertexManager, TSerialPartitionManager>, TSerialBoundaryVertexManager>::value, "TSerialBoundaryVertexManager must inherit from ISerialBoundaryVertexManager");
        static_assert(std::is_base_of<ISerialPartitionManager<TSerialGraph, TSerialActiveVertexManager>, TSerialPartitionManager>::value, "TSerialPartitionManager must inherit from ISerialPartitionManager");
        static_assert(std::is_base_of<ISerialDistanceOracle, TSerialDistanceOracle>::value, "TSerialDistanceOracle must inherit from ISerialDistanceOracle");
    public:
        // initialization
        virtual void initialize(TSerialGraph *t_p_g,
                                TSerialActiveVertexManager *t_p_av_manager,
                                TSerialBoundaryVertexManager *t_p_bv_manager,
                                TSerialPartitionManager *t_p_p_manger,
                                TSerialDistanceOracle *t_p_d_oracle,
                                std::vector<partition_t> &t_hierarchy,
                                std::vector<weight_t> &t_distance,
                                weight_t t_lmax) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALREFINER_H
