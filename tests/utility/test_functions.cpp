#include <gtest/gtest.h>
#include "src/utility/functions.h"
#include "src/datastructures/csr_graph.h"
#include "src/datastructures/partition_manager.h"

namespace HeiProMap {

CSRGraph create_functions_test_graph() {
    CSRGraph g;
    g.n = 3;
    g.m = 6; // 3 undirected edges
    g.neighborhoods.initialize(4);
    g.neighborhoods[0] = 0; g.neighborhoods[1] = 2; g.neighborhoods[2] = 4; g.neighborhoods[3] = 6;
    g.edges_v.initialize(6);
    g.edges_v[0] = 1; g.edges_v[1] = 2; // 0 -> {1, 2}
    g.edges_v[2] = 0; g.edges_v[3] = 2; // 1 -> {0, 2}
    g.edges_v[4] = 0; g.edges_v[5] = 1; // 2 -> {0, 1}
    return g;
}

TEST(FunctionsTest, IsBoundary) {
    auto g = create_functions_test_graph();
    PartitionManager pm;
    pm.initialize(g.n, 2, 3);
    pm.partition[0] = 0;
    pm.partition[1] = 1; // Neighbor of 2 is in a different partition
    pm.partition[2] = 0;

    EXPECT_TRUE(is_boundary(g, pm, 0));
    EXPECT_TRUE(is_boundary(g, pm, 1));
    EXPECT_TRUE(is_boundary(g, pm, 2));

    pm.partition[1] = 0;
    EXPECT_FALSE(is_boundary(g, pm, 1));
}

TEST(FunctionsTest, IsConnectedTo) {
    auto g = create_functions_test_graph();
    PartitionManager pm;
    pm.initialize(g.n, 2, 3);
    pm.partition[0] = 0;
    pm.partition[1] = 1;
    pm.partition[2] = 0;

    EXPECT_TRUE(is_connected_to(g, pm, 0, 1));
    EXPECT_TRUE(is_connected_to(g, pm, 0, 0));
    EXPECT_TRUE(is_connected_to(g, pm, 2, 1));
}

} // namespace HeiProMap
