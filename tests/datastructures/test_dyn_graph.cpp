#include <gtest/gtest.h>
#include "src/datastructures/dyn_graph.h"

namespace HeiProMap {

TEST(DynGraphTest, AddVertex) {
    DynGraph g;
    g.add_vertex(0, 10);
    g.add_vertex(1, 20);
    
    EXPECT_EQ(g.n, 2);
    EXPECT_EQ(g.g_weight, 30);
    EXPECT_EQ(g.get_vertex_weight(0), 10);
    EXPECT_EQ(g.get_vertex_weight(1), 20);
    EXPECT_TRUE(g.vertex_exists(0));
    EXPECT_TRUE(g.vertex_exists(1));
    EXPECT_FALSE(g.vertex_exists(2));
}

TEST(DynGraphTest, AddEdge) {
    DynGraph g;
    g.add_edge(0, 1, 5);
    
    EXPECT_EQ(g.n, 2);
    EXPECT_EQ(g.m, 2); // 2 directed edges
    EXPECT_TRUE(g.edge_exists(0, 1));
    EXPECT_TRUE(g.edge_exists(1, 0));
    EXPECT_EQ(g.get_edge_weight(0, 1), 5);
    
    g.add_edge(0, 1, 3);
    EXPECT_EQ(g.get_edge_weight(0, 1), 8);
    EXPECT_EQ(g.m, 2);
}

TEST(DynGraphTest, RemoveEdge) {
    DynGraph g;
    g.add_edge(0, 1, 5);
    g.add_edge(0, 2, 3);
    
    g.remove_edge(0, 1);
    
    EXPECT_FALSE(g.edge_exists(0, 1));
    EXPECT_FALSE(g.edge_exists(1, 0));
    EXPECT_TRUE(g.edge_exists(0, 2));
    EXPECT_EQ(g.m, 2);
}

TEST(DynGraphTest, RemoveVertex) {
    DynGraph g;
    g.add_edge(0, 1);
    g.add_edge(0, 2);
    g.add_edge(1, 2);
    
    g.remove_vertex(0);
    
    EXPECT_FALSE(g.vertex_exists(0));
    EXPECT_TRUE(g.vertex_exists(1));
    EXPECT_TRUE(g.vertex_exists(2));
    EXPECT_FALSE(g.edge_exists(0, 1));
    EXPECT_FALSE(g.edge_exists(0, 2));
    EXPECT_TRUE(g.edge_exists(1, 2));
    EXPECT_EQ(g.m, 2);
}

TEST(DynGraphTest, DirtyTracking) {
    DynGraph g;
    g.add_edge(0, 1);
    g.add_edge(1, 2);
    
    EXPECT_EQ(g.dirty_list.size(), 3);
    
    g.clear_dirty_status();
    EXPECT_EQ(g.dirty_list.size(), 0);
    
    g.add_edge(2, 3);
    EXPECT_EQ(g.dirty_list.size(), 2); // 2 and 3 are dirty
}

} // namespace HeiProMap
