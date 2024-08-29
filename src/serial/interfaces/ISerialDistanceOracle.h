#ifndef HEIDELBERGPROCESSMAPPING_ISERIALDISTANCEORACLE_H
#define HEIDELBERGPROCESSMAPPING_ISERIALDISTANCEORACLE_H

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../../interfaces/IDistanceOracle.h"

namespace HeiProMap {

    class ISerialDistanceOracle : public IDistanceOracle {
    public:
        // initialization
        virtual void initialize(std::vector<partition_t> &t_hierarchy,
                                std::vector<weight_t> &t_distance) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALDISTANCEORACLE_H
