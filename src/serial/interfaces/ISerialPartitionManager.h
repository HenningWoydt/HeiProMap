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

    template<typename TSerialGraph, typename TSerialActiveVertexManager>
    class ISerialPartitionManager : public IPartitionManager{
        static_assert(std::is_base_of<ISerialGraph, TSerialGraph>::value, "TSerialGraph must inherit from ISerialGraph");
        static_assert(std::is_base_of<ISerialActiveVertexManager<TSerialGraph>, TSerialActiveVertexManager>::value, "TSerialActiveVertexManager must inherit from ISerialActiveVertexManager");
    public:
        // initialization
        virtual void initialize(TSerialGraph *t_p_g,
                                TSerialActiveVertexManager *t_p_av_manager,
                                partition_t k) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALPARTITIONMANAGER_H
