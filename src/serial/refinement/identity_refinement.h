#ifndef HEIDELBERGPROCESSMAPPING_IDENTITY_REFINEMENT_H
#define HEIDELBERGPROCESSMAPPING_IDENTITY_REFINEMENT_H

#include "../../definitions.h"
#include "../../interfaces/IRefiner.h"
#include "../interfaces/ISerialRefiner.h"

namespace HeiProMap {

    template<typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle>
    class IdentityRefinement final : public ISerialRefiner<TSerialGraph, TSerialActiveVertexManager, TSerialBoundaryVertexManager, TSerialPartitionManager, TSerialDistanceOracle> {
    public:
        // initialization
        void initialize([[maybe_unused]] TSerialGraph *t_p_g,
                        [[maybe_unused]] TSerialActiveVertexManager *t_p_av_manager,
                        [[maybe_unused]] TSerialBoundaryVertexManager *t_p_bv_manager,
                        [[maybe_unused]] TSerialPartitionManager *t_p_p_manger,
                        [[maybe_unused]] TSerialDistanceOracle *t_p_d_oracle,
                        [[maybe_unused]] std::vector<partition_t> &t_hierarchy,
                        [[maybe_unused]] std::vector<weight_t> &t_distance,
                        [[maybe_unused]] weight_t t_lmax) override {}

        // refine
        void refine() override {}
    };
}

#endif //HEIDELBERGPROCESSMAPPING_IDENTITY_REFINEMENT_H
