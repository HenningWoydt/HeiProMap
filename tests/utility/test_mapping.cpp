#include <gtest/gtest.h>
#include "src/utility/mapping.h"

namespace HeiProMap {

TEST(MappingTest, BasicFunctionality) {
    Mapping m;
    m.initialize(5);
    m.set(0, 10);
    m.set(1, 11);
    m.set_coarse_n(2);
    
    EXPECT_EQ(m.get(0), 10);
    EXPECT_EQ(m.get(1), 11);
    EXPECT_EQ(m.get_old_n(), 5);
    EXPECT_EQ(m.get_coarse_n(), 2);
}

} // namespace HeiProMap
