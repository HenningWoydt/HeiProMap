#ifndef HEIDELBERGPROCESSMAPPING_ISERIALPARTITIONMANAGER_H
#define HEIDELBERGPROCESSMAPPING_ISERIALPARTITIONMANAGER_H

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../../definitions.h"
#include "../../macros.h"
#include "../../serial/utility/utils.h"
#include "ISerialGraph.h"
#include "ISerialActiveVertexManager.h"

namespace HeiProMap {

    class ISerialPartitionManager : public IPartitionManager{
    public:
        // initialization
        virtual void initialize(ISerialGraph *t_p_g,
                                ISerialActiveVertexManager *t_p_av_manager,
                                partition_t k) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALPARTITIONMANAGER_H
