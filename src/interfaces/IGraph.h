#ifndef HEIDELBERGPROCESSMAPPING_IGRAPH_H
#define HEIDELBERGPROCESSMAPPING_IGRAPH_H

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../utility/definitions.h"
#include "../utility/utils.h"
#include "../utility/macros.h"

namespace HeiProMap {

    class IGraph {
    public:
        // initialization
        virtual void initialize(const std::string &graph_in, u64 n_threads) = 0;

        // graph properties
        virtual vertex_t get_n() const = 0;

        virtual vertex_t get_m() const = 0;

        virtual weight_t get_weight() const = 0;

        // vertex properties
        virtual void set_weight(vertex_t u, weight_t w) = 0;

        virtual weight_t get_weight(vertex_t u) const = 0;

        virtual size_t n_neighbors(vertex_t u) const = 0;

        virtual EdgeW &neighbor(vertex_t u, size_t idx) = 0;

        virtual const EdgeW &neighbor(vertex_t u, size_t idx) const = 0;

        // edge manipulation
        virtual void add_edge(vertex_t u, vertex_t v, weight_t w) = 0;

        virtual bool edge_exists(vertex_t u, vertex_t v) const = 0;

        // coarsing and uncoarsing
        virtual void contract(vertex_t u, vertex_t v) = 0;

        virtual void uncontract(vertex_t u, vertex_t v) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IGRAPH_H
