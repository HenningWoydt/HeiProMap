#ifndef HEIDELBERGPROCESSMAPPING_IGRAPH_H
#define HEIDELBERGPROCESSMAPPING_IGRAPH_H

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../definitions.h"
#include "../macros.h"
#include "../serial/utility/utils.h"

namespace HeiProMap {

    class IGraph {
    public:
        virtual ~IGraph() = default;

        // graph properties
        virtual vertex_t get_n() const = 0;
        virtual vertex_t get_m() const = 0;
        virtual weight_t get_weight() const = 0;

        // vertex properties
        virtual weight_t get_weight(vertex_t u) const = 0;
        virtual size_t size(vertex_t u) const = 0;

        // vertex neighbor
        virtual vertex_t neighbor(vertex_t u, size_t idx) const = 0;
        virtual weight_t get_weight(vertex_t u, size_t idx) const = 0;

        // edge manipulation
        virtual bool edge_exists(vertex_t u, vertex_t v) const = 0;

        // coarsing and uncoarsing
        virtual void contract(vertex_t u, vertex_t v) = 0;
        virtual void uncontract(vertex_t u, vertex_t v) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IGRAPH_H
