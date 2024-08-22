#ifndef HEIDELBERGPROCESSMAPPING_IPARTITIONMANAGER_H
#define HEIDELBERGPROCESSMAPPING_IPARTITIONMANAGER_H

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

    class IPartitionManager {
    public:
        // initialization
        virtual void initialize(IGraph *t_p_g,
                                IActiveVertexManager *t_p_av_manager,
                                partition_t k) = 0;

        // read
        virtual const partition_t &operator[](vertex_t u) const = 0;

        // write
        virtual void set(vertex_t u, partition_t id) = 0;

        virtual void move(vertex_t u, partition_t old_id, partition_t new_id) = 0;

        // weights
        virtual weight_t get_bweight(partition_t id) const = 0;

        virtual std::vector<weight_t> get_bweights() const = 0;

        // uncoarsing
        virtual void uncontract(vertex_t u, vertex_t v) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IPARTITIONMANAGER_H
