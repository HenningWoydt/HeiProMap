#ifndef HEIDELBERGPROCESSMAPPING_IMATCHER_H
#define HEIDELBERGPROCESSMAPPING_IMATCHER_H

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../utility/definitions.h"
#include "../utility/utils.h"
#include "../utility/macros.h"
#include "IGraph.h"
#include "IActiveVertexManager.h"

namespace HeiProMap {

    class IMatcher {
    public:
        // initialize
        virtual void initialize(IGraph *t_p_g,
                                IActiveVertexManager *t_p_av_manager) = 0;

        // matching
        virtual void match(std::vector<Edge> &t_matches) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IMATCHER_H
