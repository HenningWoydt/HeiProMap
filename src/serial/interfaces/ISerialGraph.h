#ifndef HEIDELBERGPROCESSMAPPING_ISERIALGRAPH_H
#define HEIDELBERGPROCESSMAPPING_ISERIALGRAPH_H

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../../definitions.h"
#include "../../macros.h"
#include "../../interfaces/IGraph.h"
#include "../utility/utils.h"

namespace HeiProMap {

    class ISerialGraph : public IGraph {
    public:
        // initialization
        virtual void initialize(const std::string &graph_in) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALGRAPH_H
