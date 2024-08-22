#ifndef HEIDELBERGPROCESSMAPPING_IREFINER_H
#define HEIDELBERGPROCESSMAPPING_IREFINER_H

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../utility/definitions.h"
#include "../utility/utils.h"
#include "../utility/macros.h"

#include "IGraph.h"
#include "IActiveVertexManager.h"
#include "IPartitionManager.h"
#include "IDistanceOracle.h"
#include "IBoundaryVertexManager.h"

namespace HeiProMap {

    class IRefiner {
    public:
        // initialization
        virtual void initialize(IGraph *t_p_g,
                                IActiveVertexManager *t_p_av_manager,
                                IBoundaryVertexManager *t_p_bv_manager,
                                IPartitionManager *t_p_p_manger,
                                IDistanceOracle *t_p_d_oracle,
                                std::vector<partition_t> &t_hierarchy,
                                std::vector<weight_t> &t_distance,
                                weight_t t_lmax) = 0;

        // refine
        virtual void refine() = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IREFINER_H
