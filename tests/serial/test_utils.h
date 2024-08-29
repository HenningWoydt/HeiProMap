#ifndef HEIDELBERGPROCESSMAPPING_TEST_UTILS_H
#define HEIDELBERGPROCESSMAPPING_TEST_UTILS_H

#include <gtest/gtest.h>

#include "../../src/interfaces/IGraph.h"
#include "../../src/interfaces/IActiveVertexManager.h"

namespace HeiProMap {

    void graphs_are_equal(IGraph &g1, IGraph &g2);

    void graphs_are_equal(IGraph &g1, IActiveVertexManager &av_manager1,
                          IGraph &g2, IActiveVertexManager &av_manager2);

    void matchings_are_equal(const std::vector<EdgeUV> &match1, const std::vector<EdgeUV> &match2);

}

#endif //HEIDELBERGPROCESSMAPPING_TEST_UTILS_H
