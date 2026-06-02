#include <gtest/gtest.h>
#include <algorithm>
#include <numeric>
#include "src/partitioning/global_multisection.h"
#include "src/datastructures/partition_manager.h"

namespace HeiProMap {

// Dummy definition to satisfy the linker for this test
void heipa_multisection_partition_wrapper(graph_t &, partition_t, f64, u64, AlignedArray<partition_t> &, GlobalMultisectionMode, u64) {
    abort(); // Should not be called in this test
}

TEST(GlobalMultisectionRefinementTest, RefinementWorks) {
    // Create a simple line graph: 0-1-2-3
    vertex_t n = 4;
    vertex_t m = 6;
    weight_t v_weights[] = {1, 1, 1, 1};
    size_t neighborhoods[] = {0, 1, 3, 5, 6};
    vertex_t edges_v[] = {1, 0, 2, 1, 3, 2};
    weight_t edges_w[] = {1, 1, 1, 1, 1, 1};

    CSRGraph g;
    g.initialize(n, m, v_weights, neighborhoods, edges_w, edges_v);

    std::vector<partition_t> hierarchy = {4}; // 1 level, split into 4
    std::vector<weight_t> distance = {0}; 
    f64 imbalance = 0.03;
    u64 seed = 42;

    GlobalMultisectionConfiguration config;
    config.mode = GLOBAL_MULTISECTION_KAFFPA_FAST;
    config.kappa = 1;
    config.refine = true; // Overall enable option
    config.label_propagation_config.enabled = true;
    config.label_propagation_config.max_iteration = 2;
    config.quotient_graph_refinement_config.enabled = true;
    config.quotient_graph_refinement_config.max_iteration = 1;
    config.flow_based_refinement_config.enabled = true;
    config.flow_based_refinement_config.max_global_iteration = 1;

    PartitionManager pm;
    pm.initialize(n, 4, 4);

    // This should work and use the refinement
    GlobalMultisectionPartitioner::partition(g, pm, hierarchy, distance, imbalance, config, seed);

    for (vertex_t u = 0; u < n; ++u) {
        EXPECT_LT(pm[u], 4);
    }
}

} // namespace HeiProMap
