#include "test_utils.h"

namespace HeiProMap {
    void graphs_are_equal(const ISerialGraph& g1, const ISerialGraph& g2) {
        // compare number of vertices
        EXPECT_EQ(g1.get_n(), g2.get_n());

        // compare number of edges
        EXPECT_EQ(g1.get_m(), g2.get_m());

        // compare graph weight
        EXPECT_EQ(g1.get_weight(), g2.get_weight());

        // compare vertex weights
        for (vertex_t u = 0; u < g1.get_n(); u++) {
            EXPECT_EQ(g1.get_weight(u), g2.get_weight(u));
        }

        // compare neighborhood sizes
        for (vertex_t u = 0; u < g1.get_n(); u++) {
            EXPECT_EQ(g1.size(u), g2.size(u));
        }

        // compare neighborhoods
        std::vector<EdgeVW> g1_neighborhood;
        std::vector<EdgeVW> g2_neighborhood;
        for (vertex_t u = 0; u < g1.get_n(); u++) {
            g1_neighborhood.clear();
            g2_neighborhood.clear();
            for (size_t i = 0; i < g1.size(u); i++) {
                vertex_t v = g1.neighbor(u, i);
                weight_t w = g1.get_weight(u, i);
                g1_neighborhood.emplace_back(v, w);
            }
            for (size_t i = 0; i < g2.size(u); i++) {
                vertex_t v = g2.neighbor(u, i);
                weight_t w = g2.get_weight(u, i);
                g2_neighborhood.emplace_back(v, w);
            }
            std::sort(g1_neighborhood.begin(), g1_neighborhood.end());
            std::sort(g2_neighborhood.begin(), g2_neighborhood.end());
            EXPECT_EQ(g1_neighborhood, g2_neighborhood);
        }
    }
}
