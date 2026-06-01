#include <gtest/gtest.h>
#include <vector>
#include "src/partitioning/recursive_bisection.h"
#include "src/datastructures/csr_graph.h"
#include "src/datastructures/partition_manager.h"

namespace HeiProMap {

TEST(RecursiveBisectionTest, CycleGraphBisection) {
    // Create a cycle graph: 0-1-2-3-0
    vertex_t n = 4;
    vertex_t m = 8;
    weight_t v_weights[] = {1, 1, 1, 1};
    size_t neighborhoods[] = {0, 2, 4, 6, 8};
    vertex_t edges_v[] = {1, 3, 0, 2, 1, 3, 0, 2};
    weight_t edges_w[] = {1, 1, 1, 1, 1, 1, 1, 1};

    CSRGraph g;
    g.initialize(n, m, v_weights, neighborhoods, edges_w, edges_v);

    PartitionManager pm;
    pm.initialize(n, 2, 4);

    RecursiveBisectionPartitioner partitioner;
    RecursiveBisectionConfiguration config;
    config.use_full_refine = true;
    config.kappa = 1;
    config.lp_config.enabled = true;

    // Run partition into 2 blocks with 0% imbalance to force 2-2 split
    partitioner.partition(g, pm, 2, 42, 1, 0.0, BisectionMethod::BFS, config);

    EXPECT_EQ(pm.get_bweight(0), 2);
    EXPECT_EQ(pm.get_bweight(1), 2);
}

TEST(RecursiveBisectionTest, PartitionK4) {
    // Create a path graph: 0-1-2-3-4-5-6-7
    vertex_t n = 8;
    vertex_t m = 14;
    std::vector<weight_t> v_weights(n, 1);
    std::vector<size_t> neighborhoods = {0, 1, 3, 5, 7, 9, 11, 13, 14};
    std::vector<vertex_t> edges_v = {1, 0, 2, 1, 3, 2, 4, 3, 5, 4, 6, 5, 7, 6};
    std::vector<weight_t> edges_w(m, 1);

    CSRGraph g;
    g.initialize(n, m, v_weights.data(), neighborhoods.data(), edges_w.data(), edges_v.data());

    PartitionManager pm;
    pm.initialize(n, 4, 8);

    RecursiveBisectionPartitioner partitioner;
    RecursiveBisectionConfiguration config;
    config.use_full_refine = true;
    
    // Use 0% imbalance to force 2-2-2-2 split
    partitioner.partition(g, pm, 4, 42, 1, 0.0, BisectionMethod::HYBRID, config);

    EXPECT_EQ(pm.get_bweight(0), 2);
    EXPECT_EQ(pm.get_bweight(1), 2);
    EXPECT_EQ(pm.get_bweight(2), 2);
    EXPECT_EQ(pm.get_bweight(3), 2);
}

} // namespace HeiProMap
