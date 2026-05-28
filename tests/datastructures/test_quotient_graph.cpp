#include <gtest/gtest.h>
#include "src/datastructures/quotient_graph.h"

namespace HeiProMap {

TEST(QuotientGraphTest, BasicFunctionality) {
    QuotientGraph qg;
    qg.initialize(4);
    
    qg.add_edge(0, 1, 10);
    qg.add_edge(0, 2, 5);
    qg.add_edge(1, 2, 3);
    
    EXPECT_EQ(qg.get_weight(0, 1), 10);
    EXPECT_EQ(qg.get_weight(1, 0), 10);
    EXPECT_EQ(qg.get_weight(0, 2), 5);
    
    qg.remove_edge(0, 1, 4);
    EXPECT_EQ(qg.get_weight(0, 1), 6);
    
    EXPECT_TRUE(qg.has_edge(0, 1));
    qg.remove_edge(0, 1, 6);
    EXPECT_FALSE(qg.has_edge(0, 1));
}

} // namespace HeiProMap
