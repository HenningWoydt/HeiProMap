#ifndef HEIDELBERGPROCESSMAPPING_IREFINER_H
#define HEIDELBERGPROCESSMAPPING_IREFINER_H

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

namespace HeiProMap {

    class IRefiner {
    public:
        // refine
        virtual void refine() = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IREFINER_H
