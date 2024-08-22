#ifndef SERIALPROCESSMAPPING_IDENTITY_REFINEMENT_H
#define SERIALPROCESSMAPPING_IDENTITY_REFINEMENT_H

#include "../utility/definitions.h"
#include "../utility/macros.h"
#include "../utility/utils.h"
#include "../datastructures/graph.h"
#include "../interfaces/IRefiner.h"

namespace HeiProMap {

    class IdentityRefinement : public IRefiner {
    public:
        // initialization
        void initialize([[maybe_unused]] IGraph *t_p_g,
                        [[maybe_unused]] IActiveVertexManager *t_p_av_manager,
                        [[maybe_unused]] IBoundaryVertexManager *t_p_bv_manager,
                        [[maybe_unused]] IPartitionManager *t_p_p_manger,
                        [[maybe_unused]] IDistanceOracle *t_p_d_oracle,
                        [[maybe_unused]] std::vector<partition_t> &t_hierarchy,
                        [[maybe_unused]] std::vector<weight_t> &t_distance,
                        [[maybe_unused]] weight_t t_lmax) final {}

        // refine
        void refine() final {}
    };
}

#endif //SERIALPROCESSMAPPING_IDENTITY_REFINEMENT_H
