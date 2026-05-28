#include <gtest/gtest.h>
#include "src/utility/flow.h"

namespace HeiProMap {

TEST(FlowTest, ResidualFlowNetworkInitialization) {
    ResidualFlowNetwork fn;
    fn.initialize(10);
    
    fn.add_directed_edge(0, 1, 1);
    fn.add_edge_to_source(2, 1);
    fn.add_edge_from_target(3, 1);
    
    fn.finalize();
    
    EXPECT_EQ(fn.get_n(), 10);
    EXPECT_EQ(fn.get_source(), 10);
    EXPECT_EQ(fn.get_target(), 11);
    
    // Check if neighbors were added
    EXPECT_GT(fn.neighbor_count(0), 0);
    EXPECT_GT(fn.neighbor_count(2), 0);
}

// Full testing of SCC and closure algorithms is highly complex
// and would require significant test infrastructure.
// This basic test ensures the component can be built and initialized.

} // namespace HeiProMap
