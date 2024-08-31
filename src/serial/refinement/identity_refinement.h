#ifndef SERIALPROCESSMAPPING_IDENTITY_REFINEMENT_H
#define SERIALPROCESSMAPPING_IDENTITY_REFINEMENT_H

#include "../../definitions.h"
#include "../../macros.h"
#include "../utility/utils.h"
#include "../datastructures/graph.h"
#include "../../interfaces/IRefiner.h"
#include "../interfaces/ISerialRefiner.h"
#include "../interfaces/ISerialActiveVertexManager.h"
#include "../interfaces/ISerialBoundaryVertexManager.h"
#include "../interfaces/ISerialDistanceOracle.h"

namespace HeiProMap {

    template<typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle>
    class IdentityRefinement : public ISerialRefiner<TSerialGraph, TSerialActiveVertexManager, TSerialBoundaryVertexManager, TSerialPartitionManager, TSerialDistanceOracle> {
    public:
        // initialization
        void initialize([[maybe_unused]] TSerialGraph *t_p_g,
                        [[maybe_unused]] TSerialActiveVertexManager *t_p_av_manager,
                        [[maybe_unused]] TSerialBoundaryVertexManager *t_p_bv_manager,
                        [[maybe_unused]] TSerialPartitionManager *t_p_p_manger,
                        [[maybe_unused]] TSerialDistanceOracle *t_p_d_oracle,
                        [[maybe_unused]] std::vector<partition_t> &t_hierarchy,
                        [[maybe_unused]] std::vector<weight_t> &t_distance,
                        [[maybe_unused]] weight_t t_lmax) final {}

        // refine
        void refine() final {}
    };
}

#endif //SERIALPROCESSMAPPING_IDENTITY_REFINEMENT_H
