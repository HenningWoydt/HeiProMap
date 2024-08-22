#ifndef HEIDELBERGPROCESSMAPPING_IPARTITIONER_H
#define HEIDELBERGPROCESSMAPPING_IPARTITIONER_H

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
#include "IBoundaryVertexManager.h"

namespace HeiProMap {

    class IPartitioner {
    public:
        // initialization
        virtual void initialize(IGraph *t_p_g,
                                IActiveVertexManager *t_p_av_manager,
                                IBoundaryVertexManager *t_p_bv_manager,
                                IPartitionManager *t_p_p_manager,
                                std::vector<partition_t> &t_hierarchy,
                                std::vector<weight_t> &t_distance,
                                f64 t_imbalance) = 0;

        // partition
        virtual void partition() = 0;


    };

}

#endif //HEIDELBERGPROCESSMAPPING_IPARTITIONER_H
