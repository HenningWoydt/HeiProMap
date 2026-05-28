#include <gtest/gtest.h>
#include "src/utility/qap.h"
#include "src/datastructures/csr_graph.h"
#include "src/datastructures/partition_manager.h"
#include "src/datastructures/distance_oracle.h"

namespace HeiProMap {

CSRGraph create_qap_test_graph() {
    CSRGraph g;
    g.n = 3;
    g.m = 6;
    g.g_weight = 3;
    g.v_weights.initialize(3, 1);
    g.edges_w.initialize(6);
    g.edges_w[0] = 2; g.edges_w[1] = 3; // 0 -> 1(w:2), 2(w:3)
    g.edges_w[2] = 2; g.edges_w[3] = 4; // 1 -> 0(w:2), 2(w:4)
    g.edges_w[4] = 3; g.edges_w[5] = 4; // 2 -> 0(w:3), 1(w:4)

    g.neighborhoods.initialize(4);
    g.neighborhoods[0] = 0; g.neighborhoods[1] = 2; g.neighborhoods[2] = 4; g.neighborhoods[3] = 6;
    
    g.edges_v.initialize(6);
    g.edges_v[0] = 1; g.edges_v[1] = 2;
    g.edges_v[2] = 0; g.edges_v[3] = 2;
    g.edges_v[4] = 0; g.edges_v[5] = 1;

    return g;
}

TEST(QAPTest, GetQAP) {
    auto g = create_qap_test_graph();
    PartitionManager pm;
    pm.initialize(g.n, 2, g.g_weight);
    pm.partition[0] = 0;
    pm.partition[1] = 1;
    pm.partition[2] = 0;
    
    DistanceOracle oracle;
    oracle.initialize({2}, {10});

    // Edges:
    // 0-1 (w:2), p0-p1, d:10 -> 20
    // 0-2 (w:3), p0-p0, d:0  -> 0
    // 1-2 (w:4), p1-p0, d:10 -> 40
    // Total QAP = (20 + 0 + 40) * 2 = 120
    EXPECT_EQ(get_qap(g, pm, oracle), 120);
}

TEST(QAPTest, GetQAPDelta) {
    auto g = create_qap_test_graph();
    PartitionManager pm;
    pm.initialize(g.n, 2, g.g_weight);
    pm.partition[0] = 0;
    pm.partition[1] = 1;
    pm.partition[2] = 0;
    
    DistanceOracle oracle;
    oracle.initialize({2}, {10});

    // Move vertex 1 from partition 1 to 0
    // old_id=1, new_id=0
    // Neighbors of 1 are 0(p0) and 2(p0)
    // v=0: old_d=d(p0,p1)=10, new_d=d(p0,p0)=0. delta -= (10-0)*w(1,0) = -20
    // v=2: old_d=d(p0,p1)=10, new_d=d(p0,p0)=0. delta -= (10-0)*w(1,2) = -40
    // Total delta is -60. QAP is symmetric, so we double it.
    // wait, get_u_qap_delta is one-sided. So it should be -60.
    weight_t delta = get_u_qap_delta(g, 1, 1, 0, pm, oracle);
    EXPECT_EQ(delta, 60); // The formula is (old - new), so gain is positive.
}

} // namespace HeiProMap
