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

    class ISerialBoundaryVertexManager : public IBoundaryVertexManager {
    public:
        // initialization
        virtual void initialize(ISerialGraph *t_p_g,
                                ISerialActiveVertexManager *t_p_av_manager,
                                ISerialPartitionManager *t_p_p_manger,
                                partition_t t_k) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALBOUNDARYVERTEXMANAGER_H
