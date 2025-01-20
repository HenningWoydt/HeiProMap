#ifndef HEIDELBERGPROCESSMAPPING_TEST_UTILS_H
#define HEIDELBERGPROCESSMAPPING_TEST_UTILS_H

#include <gtest/gtest.h>

#include "../../src/definitions.h"
#include "../../src/serial/interfaces/ISerialGraph.h"

namespace HeiProMap {
    void graphs_are_equal(const ISerialGraph& g1, const ISerialGraph& g2);
}

#endif //HEIDELBERGPROCESSMAPPING_TEST_UTILS_H
