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

    class ISerialRefiner : public IRefiner {
    public:
        // initialization
        virtual void initialize(ISerialGraph *t_p_g,
                                ISerialActiveVertexManager *t_p_av_manager,
                                ISerialBoundaryVertexManager *t_p_bv_manager,
                                ISerialPartitionManager *t_p_p_manger,
                                ISerialDistanceOracle *t_p_d_oracle,
                                std::vector<partition_t> &t_hierarchy,
                                std::vector<weight_t> &t_distance,
                                weight_t t_lmax) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALREFINER_H
