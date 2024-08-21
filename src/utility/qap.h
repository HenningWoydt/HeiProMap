#ifndef SERIALPROCESSMAPPING_QAP_H
#define SERIALPROCESSMAPPING_QAP_H

#include "definitions.h"
#include "macros.h"
#include "utils.h"
#include "../datastructures/graph.h"
#include "../datastructures/iterators/active_vertex_iterator.h"
#include "../datastructures/translation_table.h"
#include "../datastructures/distance_oracle.h"
#include "../datastructures/partition_manager.h"

namespace SPM {
    void determine_loc(u64 u_id,
                       u64 k,
                       const std::vector<u64> &hierarchy,
                       std::vector<u64> &u_loc);

    u64 determine_distance(u64 u_id,
                           u64 v_id,
                           u64 k,
                           const std::vector<u64> &hierarchy,
                           const std::vector<u64> &distance,
                           std::vector<u64> &u_loc,
                           std::vector<u64> &v_loc);

    u64 get_qap(Graph &g,
                PartitionManager &pm,
                std::vector<u64> &hierarchy,
                std::vector<u64> &distance);

    u64 get_qap(Graph &g,
                vertex_t u,
                std::vector<vertex_t> &partition,
                std::vector<u64> &hierarchy,
                u64 k,
                std::vector<u64> &distance,
                std::vector<u64> &u_loc,
                std::vector<u64> &v_loc);

    s64 get_u_qap(Graph &g,
                  vertex_t u,
                  PartitionManager &pm,
                  DistanceOracle &dist_o);

    s64 get_u_qap(const Graph &g,
                  vertex_t u,
                  partition_t new_u_id,
                  const PartitionManager &pm,
                  const DistanceOracle &dist_o);

    void get_u_qap(Graph &g,
                   vertex_t u,
                   const std::vector<partition_t> &moves_to_check,
                   std::vector<u64> &qap_deltas,
                   PartitionManager &pm,
                   DistanceOracle &dist_o);

    void get_pweights(Graph &g,
                      PartitionManager &pm,
                                  std::vector<u64> &pweights);

    u64 get_lmax(f64 imbalance, u64 k, u64 g_weight);
}

#endif //SERIALPROCESSMAPPING_QAP_H
