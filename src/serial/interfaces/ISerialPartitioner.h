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

    class ISerialPartitioner : public IPartitioner {
    public:
        // initialization
        virtual void initialize(ISerialGraph *t_p_g,
                                ISerialActiveVertexManager *t_p_av_manager,
                                ISerialBoundaryVertexManager *t_p_bv_manager,
                                ISerialPartitionManager *t_p_p_manager,
                                std::vector<partition_t> &t_hierarchy,
                                std::vector<weight_t> &t_distance,
                                f64 t_imbalance) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALPARTITIONER_H
