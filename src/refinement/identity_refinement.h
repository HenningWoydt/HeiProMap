#ifndef SERIALPROCESSMAPPING_IDENTITY_REFINEMENT_H
#define SERIALPROCESSMAPPING_IDENTITY_REFINEMENT_H

#include "../utility/definitions.h"
#include "../utility/macros.h"
#include "../utility/utils.h"
#include "../datastructures/graph.h"

namespace SPM {

    class IdentityRefinement {
    private:

    public:
        IdentityRefinement() = default;

        void refine(std::vector<vertex_t> &partition,
                    std::vector<u64> &pweights) {
            // do nothing
        }
    };
}

#endif //SERIALPROCESSMAPPING_IDENTITY_REFINEMENT_H
