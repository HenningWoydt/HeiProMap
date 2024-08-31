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

    template<typename TParallelGraph, typename TParallelActiveVertexManager>
    class IParallelPartitionManager : public IPartitionManager {
        static_assert(std::is_base_of<IParallelGraph, TParallelGraph>::value, "TParallelGraph must inherit from IParallelGraph");
        static_assert(std::is_base_of<IParallelActiveVertexManager<TParallelGraph>, TParallelActiveVertexManager>::value, "TParallelActiveVertexManager must inherit from IParallelActiveVertexManager");
    public:
        // initialization
        virtual void initialize(TParallelGraph *t_p_g,
                                TParallelActiveVertexManager *t_p_av_manager,
                                partition_t k,
                                u64 n_threads) = 0;

        virtual bool is_boundary(vertex_t u) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IPARALLELPARTITIONMANAGER_H
