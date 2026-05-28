#include <gtest/gtest.h>
#include "src/utility/assert_state.h"
#include "src/datastructures/csr_graph.h"
#include "src/datastructures/partition_manager.h"

namespace HeiProMap {

// Basic graph for assertion tests
CSRGraph create_assert_test_graph() {
    CSRGraph g;
    g.n = 3;
    g.m = 4;
    g.v_weights.initialize(3, 1);
    g.neighborhoods.initialize(4);
    g.neighborhoods[0] = 0; g.neighborhoods[1] = 2; g.neighborhoods[2] = 4; g.neighborhoods[3] = 4;
    g.edges_v.initialize(4);
    g.edges_v[0] = 1; g.edges_v[1] = 2;
    g.edges_v[2] = 0; g.edges_v[3] = 1;
    return g;
}

#if HEAVYASSERT_ENABLED || ASSERT_ENABLED
TEST(AssertStateTest, CSRStructure) {
    auto g = create_assert_test_graph();
    EXPECT_TRUE(assert_csr_structure(g));
    
    // Introduce an error
    g.neighborhoods[1] = 3;
    EXPECT_DEATH(assert_csr_structure(g), "");
}

TEST(AssertStateTest, NoSelfLoops) {
    auto g = create_assert_test_graph();
    EXPECT_TRUE(assert_no_self_loops(g));
    
    // Introduce a self-loop
    g.edges_v[0] = 0;
    EXPECT_DEATH(assert_no_self_loops(g), "");
}

TEST(AssertStateTest, NoDoubleEdges) {
    auto g = create_assert_test_graph();
    EXPECT_TRUE(assert_no_double_edges(g));
    
    // Introduce a double edge
    g.edges_v[1] = 1;
    EXPECT_DEATH(assert_no_double_edges(g), "");
}

TEST(AssertStateTest, BWeights) {
    auto g = create_assert_test_graph();
    g.g_weight = 3;
    PartitionManager p_manager;
    p_manager.initialize(g.n, 2, g.g_weight);
    p_manager.partition[0] = 0;
    p_manager.partition[1] = 0;
    p_manager.partition[2] = 1;
    p_manager.recalculate_weights(g);

    EXPECT_TRUE(assert_bweights(g, p_manager, 2));

    // Introduce error
    p_manager.bweights[0] = 0;
    EXPECT_DEATH(assert_bweights(g, p_manager, 2), "");
}

#else
TEST(AssertStateTest, AssertsDisabled) {
    // This test will run if asserts are not enabled.
    // It just confirms the test file compiles.
    SUCCEED();
}
#endif

} // namespace HeiProMap
