#ifndef SERIALPROCESSMAPPING_QAP_H
#define SERIALPROCESSMAPPING_QAP_H

#include "definitions.h"
#include "macros.h"
#include "utils.h"
#include "../datastructures/graph.h"
#include "../datastructures/translation_table.h"
#include "../datastructures/distance_oracle.h"
#include "../datastructures/partition_manager.h"

namespace HeiProMap {
    weight_t get_qap(IGraph &g,
                     IActiveVertexManager &av_manager,
                     IPartitionManager &p_manager,
                     IDistanceOracle &d_oracle);

    weight_t get_u_qap(IGraph &g,
                       vertex_t u,
                       IPartitionManager &p_manager,
                       IDistanceOracle &d_oracle);

    weight_t get_u_qap(IGraph &g,
                       vertex_t u,
                       partition_t id,
                       IPartitionManager &p_manager,
                       IDistanceOracle &d_oracle);
}

#endif //SERIALPROCESSMAPPING_QAP_H
