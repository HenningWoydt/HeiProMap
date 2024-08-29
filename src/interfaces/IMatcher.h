#ifndef HEIDELBERGPROCESSMAPPING_IMATCHER_H
#define HEIDELBERGPROCESSMAPPING_IMATCHER_H

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../definitions.h"
#include "../macros.h"
#include "../serial/utility/utils.h"
#include "IGraph.h"
#include "IActiveVertexManager.h"

namespace HeiProMap {

    class IMatcher {
    public:
        // matching
        virtual void match(std::vector<EdgeUV> &t_matches) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IMATCHER_H
