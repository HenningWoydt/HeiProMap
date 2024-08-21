#include <vector>
#include <string>

#include <gtest/gtest.h>

#include "../src/utility/definitions.h"
#include "../src/utility/macros.h"
#include "../src/utility/utils.h"
#include "../src/datastructures/graph.h"
#include "../src/coarsening/simple_edge_matcher.h"

#include "test_util.h"

namespace SPM {

    TEST(Coarse, Graph1) {
        std::string graph_in = "../data/mapping/PGPgiantcompo.graph";

        Graph g(graph_in);
        Graph g_copy = g.copy();
        std::vector<s32> marker(g.get_n());
        std::fill(marker.begin(), marker.end(), -1);
        std::vector<vertex_t> partition(g.get_n());

        std::vector<Edge> matches;
        SimpleEdgeMatcher sem;
        sem.match(g, matches, marker, 0);

        // TODO : coarse and uncoarse

        graphs_equal(g, g_copy);
    }

}
