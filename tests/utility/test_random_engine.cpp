#include <gtest/gtest.h>
#include "src/utility/random_engine.h"

namespace HeiProMap {

TEST(RandomEngineTest, Seeding) {
    RandomEngine re1(42);
    RandomEngine re2(42);
    RandomEngine re3(43);

    EXPECT_EQ(re1.get_u64(), re2.get_u64());
    EXPECT_NE(re1.get_u64(), re3.get_u64());
}

TEST(RandomEngineTest, Range) {
    RandomEngine re(0);

    for (int i = 0; i < 100; ++i) {
        float f = re.get_f32(10.0f, 20.0f);
        EXPECT_GE(f, 10.0f);
        EXPECT_LE(f, 20.0f);
    }
    
    for (int i = 0; i < 100; ++i) {
        double d = re.get_f64(100.0, 200.0);
        EXPECT_GE(d, 100.0);
        EXPECT_LE(d, 200.0);
    }
}

} // namespace HeiProMap
