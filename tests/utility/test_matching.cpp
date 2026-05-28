#include <gtest/gtest.h>
#include "src/utility/matching.h"

namespace HeiProMap {

TEST(MatchingTest, Initialization) {
    Matching m;
    m.initialize(10);
    
    EXPECT_EQ(m.get_n(), 10);
    EXPECT_EQ(m.size(), 0);
    for (vertex_t i = 0; i < 10; ++i) {
        EXPECT_FALSE(m.is_matched(i));
        EXPECT_EQ(m.get_partner(i), i);
    }
}

TEST(MatchingTest, AddAndMatch) {
    Matching m;
    m.initialize(10);
    
    m.add(0, 1);
    m.add(2, 3);
    
    EXPECT_EQ(m.size(), 2);
    EXPECT_TRUE(m.is_matched(0));
    EXPECT_TRUE(m.is_matched(1));
    EXPECT_EQ(m.get_partner(0), 1);
    EXPECT_EQ(m.get_partner(1), 0);
    
    EXPECT_FALSE(m.is_matched(4));
}

TEST(MatchingTest, Translation) {
    Matching m;
    m.initialize(6);
    
    m.add(0, 1);
    m.add(3, 5);
    
    m.set_translation();
    
    EXPECT_EQ(m.get_n_coarse_nodes(), 4);
    
    // Matched pairs
    EXPECT_EQ(m.get_n(0), m.get_n(1));
    EXPECT_EQ(m.get_n(3), m.get_n(5));
    
    // Unmatched nodes
    EXPECT_NE(m.get_n(2), m.get_n(0));
    EXPECT_NE(m.get_n(4), m.get_n(0));
    EXPECT_NE(m.get_n(2), m.get_n(4));
}

} // namespace HeiProMap
