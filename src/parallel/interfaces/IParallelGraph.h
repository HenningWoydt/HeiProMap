#ifndef HEIDELBERGPROCESSMAPPING_IPARALLELGRAPH_H
#define HEIDELBERGPROCESSMAPPING_IPARALLELGRAPH_H

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../../interfaces/IGraph.h"

namespace HeiProMap {

    class IParallelGraph : public IGraph {
    public:
        // initialization
        virtual void initialize(const std::string &graph_in, u64 n_threads) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IPARALLELGRAPH_H
