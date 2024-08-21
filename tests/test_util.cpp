#include "test_util.h"

namespace SPM {

    void graphs_equal(const Graph &g1, const Graph &g2) {
        EXPECT_EQ(g1.get_n(), g2.get_n());
        EXPECT_EQ(g1.get_m(), g2.get_m());
        EXPECT_EQ(g1.get_sum_vertex_weights(), g2.get_sum_vertex_weights());
        EXPECT_EQ(g1.get_sum_edge_weights(), g2.get_sum_edge_weights());

        for (vertex_t u = 0; u < g1.get_n(); ++u) {
            EXPECT_EQ(g1.get_vertex_weight(u), g2.get_vertex_weight(u));
            EXPECT_EQ(g1.get_vertex_state(u), g2.get_vertex_state(u));
            EXPECT_EQ(g1.get_vertex_n_edge(u), g2.get_vertex_n_edge(u));
            EXPECT_EQ(g1[u].size(), g2[u].size());

            for (vertex_t v = u + 1; v < g1.get_n(); ++v) {
                EXPECT_EQ(g1.edge_exists(u, v), g2.edge_exists(u, v));
                EXPECT_EQ(g1.edge_exists(v, u), g2.edge_exists(v, u));
                EXPECT_EQ(g1.edge_exists_2way(u, v), g2.edge_exists_2way(u, v));
                if(g1.edge_exists(u, v)){
                    EXPECT_EQ(g1.get_edge_weight(u, v), g2.get_edge_weight(u, v));
                }
            }
        }
    }
}
