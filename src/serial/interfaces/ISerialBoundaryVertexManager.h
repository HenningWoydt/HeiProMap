#ifndef HEIDELBERGPROCESSMAPPING_ISERIALBOUNDARYVERTEXMANAGER_H
#define HEIDELBERGPROCESSMAPPING_ISERIALBOUNDARYVERTEXMANAGER_H

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../../interfaces/IBoundaryVertexManager.h"
#include "ISerialGraph.h"
#include "ISerialActiveVertexManager.h"
#include "ISerialPartitionManager.h"

namespace HeiProMap {

    template<typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialPartitionManager>
    class ISerialBoundaryVertexManager : public IBoundaryVertexManager {
        static_assert(std::is_base_of<ISerialGraph, TSerialGraph>::value, "TSerialGraph must inherit from ISerialGraph");
        static_assert(std::is_base_of<ISerialActiveVertexManager<TSerialGraph>, TSerialActiveVertexManager>::value, "TSerialActiveVertexManager must inherit from ISerialActiveVertexManager");
        static_assert(std::is_base_of<ISerialPartitionManager<TSerialGraph, TSerialActiveVertexManager>, TSerialPartitionManager>::value, "TSerialPartitionManager must inherit from ISerialPartitionManager");
    public:
        // initialization
        virtual void initialize(TSerialGraph *t_p_g,
                                TSerialActiveVertexManager *t_p_av_manager,
                                TSerialPartitionManager *t_p_p_manger,
                                partition_t t_k) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALBOUNDARYVERTEXMANAGER_H
