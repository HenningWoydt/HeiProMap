#include <gtest/gtest.h>
#include "src/utility/small_map.h"

namespace HeiProMap {

TEST(FlatMapTest, BasicFunctionality) {
    FlatMap<int, int> map;
    map[1] = 100;
    map[2] = 200;
    
    EXPECT_EQ(map[1], 100);
    EXPECT_EQ(map[2], 200);
    EXPECT_EQ(map.size(), 2);
    
    map.add(1, 50);
    EXPECT_EQ(map[1], 150);
}

} // namespace HeiProMap
