#include <gtest/gtest.h>
#include "src/datastructures/partition_manager.h"

namespace HeiProMap {

TEST(PartitionManagerTest, BasicFunctionality) {
    PartitionManager pm;
    pm.initialize(10, 2, 100); // 10 vertices, 2 blocks, total weight 100
    
    EXPECT_EQ(pm.n, 10);
    EXPECT_EQ(pm.k, 2);
    EXPECT_EQ(pm.get_bweight(0), 100);
    EXPECT_EQ(pm.get_bweight(1), 0);
    
    pm.move_serial(0, 10, 0, 1);
    EXPECT_EQ(pm.get_bweight(0), 90);
    EXPECT_EQ(pm.get_bweight(1), 10);
    EXPECT_EQ(pm[0], 1);
    
    EXPECT_EQ(pm.max_weight(), 90);
    EXPECT_EQ(pm.n_empty_blocks(), 0);
}

} // namespace HeiProMap
