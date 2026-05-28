#include <gtest/gtest.h>
#include "src/datastructures/boundary_vertex_manger.h"

namespace HeiProMap {

TEST(BoundaryVertexManagerTest, BasicFunctionality) {
    BoundaryVertexManager bvm;
    bvm.initialize(10, 2);
    
    bvm.add(0, 0);
    bvm.add(1, 0);
    bvm.add(2, 1);
    
    EXPECT_TRUE(bvm.is_boundary(0));
    EXPECT_TRUE(bvm.is_boundary(1));
    EXPECT_TRUE(bvm.is_boundary(2));
    EXPECT_FALSE(bvm.is_boundary(3));
    
    EXPECT_EQ(bvm.size(0), 2);
    EXPECT_EQ(bvm.size(1), 1);
    
    EXPECT_EQ(bvm.get(0, 0), 0);
    EXPECT_EQ(bvm.get(0, 1), 1);
    EXPECT_EQ(bvm.get(1, 0), 2);
}

} // namespace HeiProMap
