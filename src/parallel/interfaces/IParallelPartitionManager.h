#ifndef HEIDELBERGPROCESSMAPPING_IPARALLELPARTITIONMANAGER_H
#define HEIDELBERGPROCESSMAPPING_IPARALLELPARTITIONMANAGER_H

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../../definitions.h"
#include "../../macros.h"
#include "../../serial/utility/utils.h"

#include "../../interfaces/IGraph.h"
#include "../../interfaces/IActiveVertexManager.h"
#include "IParallelGraph.h"
#include "../../interfaces/IPartitionManager.h"
#include "IParallelActiveVertexManager.h"

namespace HeiProMap {

    class IParallelPartitionManager : public IPartitionManager {
    public:
        // initialization
        virtual void initialize(IParallelGraph *t_p_g,
                                IParallelActiveVertexManager *t_p_av_manager,
                                partition_t k,
                                u64 n_threads) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IPARALLELPARTITIONMANAGER_H
