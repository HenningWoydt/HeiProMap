#include <gtest/gtest.h>
#include <vector>
#include "src/utility/flow.h"
#include "src/utility/translation_table.h"
#include "src/utility/random_engine.h"

namespace HeiProMap {

struct MockGraph {
    std::vector<weight_t> v_weights;
};

class SCCGraphTest : public ::testing::Test {
protected:
    void SetUp() override {
    }
};

TEST_F(SCCGraphTest, IdentificationBasic) {
    // 0 -> 1 -> 2
    ResidualFlowNetwork fn;
    fn.initialize(3);
    fn.add_directed_edge(0, 1, 1);
    fn.add_directed_edge(1, 2, 1);
    fn.finalize();

    MockGraph g;
    g.v_weights = {10, 20, 30};
    
    TranslationTable<vertex_t> tt;
    tt.reserve(3, 3);
    tt.add(0, 0); tt.add(1, 1); tt.add(2, 2);

    SCCGraph scc;
    scc.initialize(fn, g, tt);
    
    // Source(3), Target(4), and 0, 1, 2 are all in separate SCCs because it's a DAG.
    EXPECT_EQ(scc.get_n_scc(), 5);
}

TEST_F(SCCGraphTest, IdentificationCycle) {
    // 0 -> 1 -> 2 -> 0
    ResidualFlowNetwork fn;
    fn.initialize(3);
    fn.add_directed_edge(0, 1, 1);
    fn.add_directed_edge(1, 2, 1);
    fn.add_directed_edge(2, 0, 1);
    fn.finalize();

    MockGraph g;
    g.v_weights = {10, 20, 30};
    
    TranslationTable<vertex_t> tt;
    tt.reserve(3, 3);
    tt.add(0, 0); tt.add(1, 1); tt.add(2, 2);

    SCCGraph scc;
    scc.initialize(fn, g, tt);
    
    // {0, 1, 2} should be one SCC. Source(3) and Target(4) are separate.
    EXPECT_EQ(scc.get_n_scc(), 3);
}

TEST_F(SCCGraphTest, ReductionReachability) {
    // S -> 0 -> 1 -> T
    // 2 is isolated
    ResidualFlowNetwork fn;
    fn.initialize(3);
    fn.add_edge_from_source(0, 1);
    fn.add_directed_edge(0, 1, 1);
    fn.add_edge_to_target(1, 1);
    fn.finalize();

    MockGraph g;
    g.v_weights = {1, 1, 1};
    TranslationTable<vertex_t> tt;
    tt.reserve(3, 3);
    tt.add(0,0); tt.add(1,1); tt.add(2,2);

    SCCGraph scc;
    scc.initialize(fn, g, tt);
    scc.reduce();
}

TEST_F(SCCGraphTest, ClosureFindingSimple) {
    // S -> 0 (weight 10), 1 (weight 20) -> T
    // Active nodes: None (0 is s-successor, 1 is t-predecessor)
    ResidualFlowNetwork fn;
    fn.initialize(2);
    fn.add_edge_from_source(0, 1);
    fn.add_edge_to_target(1, 1);
    fn.finalize();

    MockGraph g;
    g.v_weights = {10, 20};
    TranslationTable<vertex_t> tt;
    tt.reserve(2, 2);
    tt.add(0,0); tt.add(1,1);

    SCCGraph scc;
    scc.initialize(fn, g, tt);
    scc.reduce();

    RandomEngine re(42);
    std::vector<u8> is_left;
    // Lmax = 15. non-region = 0.
    // Fixed source side: {0}. weight 10.
    // Fixed sink side: {1}. weight 20.
    // Total left weight = 0 (non-region) + 10 = 10 <= 15.
    // Total right weight = 0 (non-region) + 20 = 20 > 15.
    bool found = scc.find_best_closure(0, 0, 15, 15, 15.0, 10, re, is_left);
    
    EXPECT_FALSE(found); // Right side (20) exceeds Lmax (15)

    // Increase Lmax to 25
    found = scc.find_best_closure(0, 0, 25, 25, 15.0, 10, re, is_left);
    EXPECT_TRUE(found);
    ASSERT_EQ(is_left.size(), 2);
    EXPECT_EQ(is_left[0], 1); // 0 is on source side
    EXPECT_EQ(is_left[1], 0); // 1 is on sink side
}

TEST_F(SCCGraphTest, ClosureProperty) {
    ResidualFlowNetwork fn;
    fn.initialize(3);
    fn.add_edge_from_source(0, 1);
    fn.add_directed_edge(0, 1, 1);
    fn.add_directed_edge(1, 2, 1);
    fn.add_edge_to_target(2, 1);
    fn.finalize();

    MockGraph g;
    g.v_weights = {1, 1, 1};
    TranslationTable<vertex_t> tt;
    tt.reserve(3, 3);
    tt.add(0,0); tt.add(1,1); tt.add(2,2);

    SCCGraph scc;
    scc.initialize(fn, g, tt);
    scc.reduce();
}

TEST_F(SCCGraphTest, ClosureDirection) {
    // Test if the closure property follows residual edges correctly.
    // S -> 0, 1 -> T.
    // Residual edge: 1 -> 0  (meaning 0 was source-side, 1 was sink-side, and we can push back)
    // If we have 1 -> 0 in residual, then:
    // If 1 is in SOURCE side, 0 MUST be in SOURCE side.
    
    ResidualFlowNetwork fn;
    fn.initialize(2);
    // Nodes 0, 1. Source 2, Sink 3.
    fn.add_directed_edge(1, 0, 1); // Residual edge 1 -> 0
    fn.finalize();

    MockGraph g;
    g.v_weights = {10, 10};
    TranslationTable<vertex_t> tt;
    tt.reserve(2, 2);
    tt.add(0,0); tt.add(1,1);

    SCCGraph scc;
    scc.initialize(fn, g, tt);
    scc.reduce();
    
    RandomEngine re(42);
    std::vector<u8> is_left;
    
    bool found = scc.find_best_closure(0, 0, 11, 11, 10.0, 100, re, is_left);
    EXPECT_TRUE(found);
    EXPECT_EQ(is_left[0], 1); // 0 must be in if 1 is in, or 0 can be in alone.
    // {0} weight 10. Left 10, Right 10. Cost 0.
    // {1} weight 10. NOT A VALID CLOSURE.
}

TEST_F(SCCGraphTest, ClosureDirectionRandomized) {
    ResidualFlowNetwork fn;
    fn.initialize(25); // Nodes 0..24. Source 25, Sink 26.
    fn.add_directed_edge(1, 0, 1); // 1 -> 0 residual
    // nodes 2..24 are isolated and active
    fn.finalize();

    MockGraph g;
    g.v_weights.assign(25, 0);
    g.v_weights[0] = 10;
    g.v_weights[1] = 10;
    // others weight 0
    
    TranslationTable<vertex_t> tt;
    tt.reserve(25, 25);
    for(vertex_t i=0; i<25; ++i) tt.add(i, i);

    SCCGraph scc;
    scc.initialize(fn, g, tt);
    scc.reduce();
    
    RandomEngine re(42);
    std::vector<u8> is_left;
    
    // n_active = 25. Should use randomized version.
    bool found = scc.find_best_closure(0, 0, 11, 11, 10.0, 100, re, is_left);
    EXPECT_TRUE(found);
    EXPECT_EQ(is_left[0], 1); 
}

} // namespace HeiProMap
