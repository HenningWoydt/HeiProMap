#ifndef HEIDELBERGPROCESSMAPPING_IACTIVEVERTEXMANAGER_H
#define HEIDELBERGPROCESSMAPPING_IACTIVEVERTEXMANAGER_H

#include "IGraph.h"

namespace HeiProMap {

    class IActiveVertexManager {
    public:
        // initialize
        virtual void initialize(IGraph *t_p_g) = 0;

        // active vertex manipulation
        virtual vertex_t get_n_active() const = 0;

        virtual void activate_vertex(vertex_t u) = 0;

        virtual void disable_vertex(vertex_t u) = 0;

        virtual bool is_active(vertex_t u) const = 0;

        virtual bool is_disabled(vertex_t u) const = 0;

        virtual bool get_state(vertex_t u) const = 0;

        // coarsing and uncoarsing
        virtual void contract(vertex_t u, vertex_t v) = 0;

        virtual void uncontract(vertex_t u, vertex_t v) = 0;

        // iteration
        virtual void reset_iterator() = 0;

        virtual vertex_t get() = 0;

        virtual void next() = 0;

        virtual bool available() = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IACTIVEVERTEXMANAGER_H
