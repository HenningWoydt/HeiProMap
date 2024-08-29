#ifndef HEIDELBERGPROCESSMAPPING_IBOUNDARYVERTEXMANAGER_H
#define HEIDELBERGPROCESSMAPPING_IBOUNDARYVERTEXMANAGER_H

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

    class IBoundaryVertexManager {
    public:
        // add
        virtual vertex_t get_n_boundary() const = 0;

        virtual void insert(vertex_t u, partition_t id) = 0;

        virtual void move(vertex_t u, partition_t old_id, partition_t id) = 0;

        // check
        virtual bool is_boundary(vertex_t u) const = 0;

        // uncontract
        virtual void uncontract(vertex_t u, vertex_t v) = 0;

        // iteration
        virtual void reset_iterator() = 0;

        virtual vertex_t get() = 0;

        virtual void next() = 0;

        virtual bool available() = 0;

        // block iteration
        virtual void reset_iterator(partition_t id) = 0;

        virtual vertex_t get(partition_t id) = 0;

        virtual void next(partition_t id) = 0;

        virtual bool available(partition_t id) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IBOUNDARYVERTEXMANAGER_H
