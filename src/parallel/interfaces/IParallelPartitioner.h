#ifndef HEIDELBERGPROCESSMAPPING_IPARALLELPARTITIONER_H
#define HEIDELBERGPROCESSMAPPING_IPARALLELPARTITIONER_H

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../../interfaces/IPartitioner.h"
#include "IParallelGraph.h"
#include "IParallelActiveVertexManager.h"
#include "IParallelBoundaryVertexManager.h"

namespace HeiProMap {

    class IParallelPartitioner : public IPartitioner {
    public:
        // initialization
        virtual void initialize(IParallelGraph *t_p_g,
                                IParallelActiveVertexManager *t_p_av_manager,
                                IParallelBoundaryVertexManager *t_p_bv_manager,
                                IParallelPartitionManager *t_p_p_manager,
                                std::vector<partition_t> &t_hierarchy,
                                std::vector<weight_t> &t_distance,
                                f64 t_imbalance,
                                u64 n_threads) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IPARALLELPARTITIONER_H
