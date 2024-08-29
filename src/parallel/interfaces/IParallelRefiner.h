#ifndef HEIDELBERGPROCESSMAPPING_IPARALLELREFINER_H
#define HEIDELBERGPROCESSMAPPING_IPARALLELREFINER_H

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../../interfaces/IRefiner.h"
#include "IParallelActiveVertexManager.h"
#include "IParallelBoundaryVertexManager.h"
#include "IParallelDistanceOracle.h"
#include "IParallelQuotientGraph.h"

namespace HeiProMap {

    class IParallelRefiner : public IRefiner {
    public:
        // initialization
        virtual void initialize(IParallelGraph *t_p_g,
                                IParallelActiveVertexManager *t_p_av_manager,
                                IParallelBoundaryVertexManager *t_p_bv_manager,
                                IParallelPartitionManager *t_p_p_manger,
                                IParallelDistanceOracle *t_p_d_oracle,
                                IParallelQuotientGraph *t_p_qgraph,
                                std::vector<partition_t> &t_hierarchy,
                                std::vector<weight_t> &t_distance,
                                weight_t t_lmax,
                                u64 n_threads) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IPARALLELREFINER_H
