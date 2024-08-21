#ifndef SERIALPROCESSMAPPING_TEST_UTIL_H
#define SERIALPROCESSMAPPING_TEST_UTIL_H

#include <gtest/gtest.h>

#include "../src/datastructures/graph.h"

namespace SPM {

    void graphs_equal(const Graph &g1, const Graph &g2);

}

#endif //SERIALPROCESSMAPPING_TEST_UTIL_H
