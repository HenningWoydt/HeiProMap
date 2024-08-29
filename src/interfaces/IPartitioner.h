#ifndef HEIDELBERGPROCESSMAPPING_IPARTITIONER_H
#define HEIDELBERGPROCESSMAPPING_IPARTITIONER_H

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
#include "IBoundaryVertexManager.h"

namespace HeiProMap {

    class IPartitioner {
    public:
        // partition
        virtual void partition() = 0;


    };

}

#endif //HEIDELBERGPROCESSMAPPING_IPARTITIONER_H
