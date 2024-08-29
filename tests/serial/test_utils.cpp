#include "test_utils.h"

namespace HeiProMap {

    void graphs_are_equal(IGraph &g1, IGraph &g2){
        EXPECT_EQ(g1.get_n(), g2.get_n());
        EXPECT_EQ(g1.get_weight(), g2.get_weight());

        for(vertex_t u = 0; u < g1.get_n(); ++u){
            EXPECT_EQ(g1.get_weight(u), g2.get_weight(u));
        }

        for(vertex_t u = 0; u < g1.get_n(); ++u){
            for(size_t i = 0; i < g1.size(u); ++i){
                vertex_t v = g1.neighbor(u, i);
                EXPECT_TRUE(g2.edge_exists(u, v));
            }
        }

        for(vertex_t u = 0; u < g2.get_n(); ++u){
            for(size_t i = 0; i < g2.size(u); ++i){
                vertex_t v = g2.neighbor(u, i);
                EXPECT_TRUE(g1.edge_exists(u, v));
            }
        }
    }

    void graphs_are_equal(IGraph &g1, IActiveVertexManager &av_manager1,
                          IGraph &g2, IActiveVertexManager &av_manager2){
        EXPECT_EQ(av_manager1.get_n_active(), av_manager2.get_n_active());
        EXPECT_EQ(g1.get_weight(), g2.get_weight());

        for (av_manager1.reset_iterator(); av_manager1.available(); av_manager1.next()) {
            vertex_t u = av_manager1.get();
            EXPECT_EQ(g1.get_weight(u), g2.get_weight(u));
        }

        for (av_manager1.reset_iterator(); av_manager1.available(); av_manager1.next()) {
            vertex_t u = av_manager1.get();
            for(size_t i = 0; i < g1.size(u); ++i){
                vertex_t v = g1.neighbor(u, i);
                EXPECT_TRUE(g2.edge_exists(u, v));
            }
        }

        for (av_manager2.reset_iterator(); av_manager2.available(); av_manager2.next()) {
            vertex_t u = av_manager2.get();
            for(size_t i = 0; i < g2.size(u); ++i){
                vertex_t v = g2.neighbor(u, i);
                EXPECT_TRUE(g1.edge_exists(u, v));
            }
        }
    }

    void matchings_are_equal(const std::vector<EdgeUV> &match1, const std::vector<EdgeUV> &match2){
        EXPECT_EQ(match1.size(), match2.size());

        for(size_t i = 0; i < match1.size(); ++i){
            EXPECT_EQ(match1[i].u, match2[i].u);
            EXPECT_EQ(match1[i].v, match2[i].v);
        }
    }

}
