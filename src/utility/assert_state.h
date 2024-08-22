#ifndef HEIDELBERGPROCESSMAPPING_ASSERT_STATE_H
#define HEIDELBERGPROCESSMAPPING_ASSERT_STATE_H

#include "../interfaces/IGraph.h"
#include "../interfaces/IActiveVertexManager.h"
#include "../interfaces/IPartitionManager.h"
#include "../interfaces/IBoundaryVertexManager.h"

namespace HeiProMap {

    bool assert_state_pre_partitioning(IGraph &g,
                                       IActiveVertexManager &av_manager);

    bool assert_state_after_partitioning(IGraph &g,
                                         IActiveVertexManager &av_manager,
                                         IPartitionManager &p_manager,
                                         IBoundaryVertexManager &bv_manager,
                                         partition_t k);

};


#endif //HEIDELBERGPROCESSMAPPING_ASSERT_STATE_H
