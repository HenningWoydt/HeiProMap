#include <gtest/gtest.h>
#include "src/datastructures/distance_oracle.h"

namespace HeiProMap {

TEST(DistanceOracleTest, BasicHierarchy) {
    DistanceOracle oracle;
    // 2 nodes, 2 sockets per node, 2 cores per socket = 8 cores total
    std::vector<partition_t> hierarchy = {2, 2, 2}; 
    // distances: core=1, socket=10, node=100
    std::vector<weight_t> distances = {1, 10, 100};
    
    oracle.initialize(hierarchy, distances);
    
    EXPECT_EQ(oracle.get_k(), 8);
    
    // Within same socket (cores 0 and 1) -> distance 1
    EXPECT_EQ(oracle.get(0, 1), 1);
    // Different socket, same node (cores 0 and 2) -> distance 10
    EXPECT_EQ(oracle.get(0, 2), 10);
    // Different node (cores 0 and 4) -> distance 100
    EXPECT_EQ(oracle.get(0, 4), 100);
    
    EXPECT_EQ(oracle.get(0, 0), 0);
}

} // namespace HeiProMap
