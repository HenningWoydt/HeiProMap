#ifndef HEIDELBERGPROCESSMAPPING_IPARALLELBOUNDARYVERTEXMANAGER_H
#define HEIDELBERGPROCESSMAPPING_IPARALLELBOUNDARYVERTEXMANAGER_H

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../../interfaces/IBoundaryVertexManager.h"
#include "IParallelGraph.h"
#include "IParallelActiveVertexManager.h"
#include "IParallelPartitionManager.h"

namespace HeiProMap {

    template<typename TParallelGraph, typename TParallelActiveVertexManager, typename TParallelPartitionManager>
    class IParallelBoundaryVertexManager : public IBoundaryVertexManager {
        static_assert(std::is_base_of<IParallelGraph, TParallelGraph>::value, "TParallelGraph must inherit from IParallelGraph");
        static_assert(std::is_base_of<IParallelActiveVertexManager<TParallelGraph>, TParallelActiveVertexManager>::value, "TParallelActiveVertexManager must inherit from IParallelActiveVertexManager");
        static_assert(std::is_base_of<IParallelPartitionManager<TParallelGraph, TParallelActiveVertexManager>, TParallelPartitionManager>::value, "TParallelActiveVertexManager must inherit from IParallelActiveVertexManager");
    public:
        // initialization
        virtual void initialize(TParallelGraph *t_p_g,
                                TParallelActiveVertexManager *t_p_av_manager,
                                TParallelPartitionManager *t_p_p_manger,
                                partition_t t_k,
                                u64 n_threads) = 0;

        virtual vertex_t get_vertex(size_t idx) const = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IPARALLELBOUNDARYVERTEXMANAGER_H
