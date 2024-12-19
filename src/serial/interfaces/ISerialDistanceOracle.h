#ifndef HEIDELBERGPROCESSMAPPING_ISERIALDISTANCEORACLE_H
#define HEIDELBERGPROCESSMAPPING_ISERIALDISTANCEORACLE_H

#include <regex>
#include <vector>

namespace HeiProMap {
    class ISerialDistanceOracle {
    public:
        virtual ~ISerialDistanceOracle() = default;
        virtual void initialize(std::vector<partition_t>& t_hierarchy,
                                std::vector<weight_t>& t_distance) = 0;
    };
}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALDISTANCEORACLE_H
