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

    TEST(Copy, Graph1) {
        std::string graph_in = "../data/mapping/PGPgiantcompo.graph";

        Graph g(graph_in);
        Graph g_copy = g.copy();

        graphs_equal(g, g_copy);
    }

}
