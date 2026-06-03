#include <gtest/gtest.h>
#include "src/refinement/swap_refinement.h"
#include "src/datastructures/csr_graph.h"
#include "src/datastructures/partition_manager.h"
#include "src/datastructures/distance_oracle.h"
#include "src/datastructures/boundary_vertex_manger.h"
#include "src/datastructures/quotient_graph.h"
#include "src/datastructures/block_conn.h"
#include <vector>

namespace HeiProMap {

TEST(SwapRefinementTest, SimpleSwap) {
    // 4-vertex path: 0-1-2-3
    // Partition: 0: {0, 2}, 1: {1, 3}
    // Edges: (0,1), (1,2), (2,3)
    // All edges are cut!
    
    vertex_t n = 4;
    vertex_t m = 6;
    std::vector<weight_t> v_weights = {1, 1, 1, 1};
    std::vector<size_t> neighborhoods = {0, 1, 3, 5, 6};
    std::vector<vertex_t> edges_v = {1, 0, 2, 1, 3, 2};
    std::vector<weight_t> edges_w = {1, 1, 1, 1, 1, 1};

    CSRGraph g;
    g.initialize(n, m, v_weights.data(), neighborhoods.data(), edges_w.data(), edges_v.data());

    PartitionManager pm;
    pm.initialize(n, 2, 4);
    pm.reset_weights();
    pm.set(0, 1, 0);
    pm.set(2, 1, 0);
    pm.set(1, 1, 1);
    pm.set(3, 1, 1);

    DistanceOracle do_oracle;
    do_oracle.initialize({2}, {1}); // k=2, dist=1

    bv_manager_t bv;
    bv.initialize(n, 2);
    q_graph_t qg;
    qg.initialize(2);
    block_conn_t bc;
    bc.initialize(n, m, 2);

    for (vertex_t u = 0; u < n; ++u) {
        bc.begin_vertex(g, u);
        partition_t u_id = pm[u];
        for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
            vertex_t v = g.edges_v[i];
            partition_t v_id = pm[v];
            bc.add_connection(u, v_id, g.edges_w[i]);
            if (u_id != v_id) {
                bv.add(u, u_id);
                if (u < v) qg.add_edge(u_id, v_id, g.edges_w[i]);
            }
        }
    }

    AlignedArray<weight_t> lmax;
    lmax.initialize(2);
    lmax[0] = 2;
    lmax[1] = 2;

    SwapRefinementConfiguration config("Swap");
    config.enabled = true;
    config.max_iteration = 1;

    SwapRefinement refine;
    refine.initialize(n, m, 2, 1, 42, config);

    weight_t qap_before = get_qap(g, pm, do_oracle);
    refine.refine(g, do_oracle, bv, pm, qg, bc, lmax, true, true);
    weight_t qap_after = get_qap(g, pm, do_oracle);

    EXPECT_LT(qap_after, qap_before);
}

} // namespace HeiProMap
