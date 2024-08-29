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

    class IParallelBoundaryVertexManager : public IBoundaryVertexManager {
    public:
        // initialization
        virtual void initialize(IParallelGraph *t_p_g,
                                IParallelActiveVertexManager *t_p_av_manager,
                                IParallelPartitionManager *t_p_p_manger,
                                partition_t t_k,
                                u64 n_threads) = 0;

        virtual vertex_t get_vertex(size_t idx) const = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IPARALLELBOUNDARYVERTEXMANAGER_H
