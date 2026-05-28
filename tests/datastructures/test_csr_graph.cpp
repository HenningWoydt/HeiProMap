#include <gtest/gtest.h>
#include "src/datastructures/csr_graph.h"

namespace HeiProMap {

TEST(CSRGraphTest, BasicInitialization) {
    vertex_t n = 5;
    vertex_t m = 8; // 4 undirected edges
    weight_t g_weight = 10;
    
    CSRGraph g(n, m, g_weight);
    
    EXPECT_EQ(g.n, n);
    EXPECT_EQ(g.m, m);
    EXPECT_EQ(g.g_weight, g_weight);
    EXPECT_FALSE(g.uniform_v_weights);
    EXPECT_FALSE(g.uniform_e_weights);
}

TEST(CSRGraphTest, DegreeCalculation) {
    vertex_t n = 3;
    vertex_t m = 4;
    weight_t g_weight = 3;
    
    CSRGraph g(n, m, g_weight);
    
    // Create a simple triangle graph
    // 0: 1, 2
    // 1: 0, 2
    // 2: 0, 1
    // Actually m=6 for triangle, let's just do a path 0-1-2
    // 0: 1
    // 1: 0, 2
    // 2: 1
    // m = 4
    
    g.neighborhoods[0] = 0;
    g.neighborhoods[1] = 1;
    g.neighborhoods[2] = 3;
    g.neighborhoods[3] = 4;
    
    EXPECT_EQ(g.deg(0), 1);
    EXPECT_EQ(g.deg(1), 2);
    EXPECT_EQ(g.deg(2), 1);
}

} // namespace HeiProMap
