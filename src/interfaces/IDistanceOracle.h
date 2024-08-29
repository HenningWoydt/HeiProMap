#ifndef HEIDELBERGPROCESSMAPPING_IDISTANCEORACLE_H
#define HEIDELBERGPROCESSMAPPING_IDISTANCEORACLE_H

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
#include "IPartitionManager.h"

namespace HeiProMap {

    class IDistanceOracle {
    public:
        // distances
        virtual weight_t get(vertex_t u_id, vertex_t v_id) const = 0;

        virtual partition_t get_h(vertex_t u_id, vertex_t v_id) const = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IDISTANCEORACLE_H
