#ifndef HEIDELBERGPROCESSMAPPING_IPARALLELDISTANCEORACLE_H
#define HEIDELBERGPROCESSMAPPING_IPARALLELDISTANCEORACLE_H

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../../interfaces/IDistanceOracle.h"

namespace HeiProMap {

    class IParallelDistanceOracle : public IDistanceOracle {
    public:
        // initialization
        virtual void initialize(std::vector<partition_t> &t_hierarchy,
                                std::vector<weight_t> &t_distance,
                                u64 n_threads) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IPARALLELDISTANCEORACLE_H
