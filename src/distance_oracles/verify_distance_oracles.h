/*******************************************************************************
 * MIT License
 *
 * This file is part of HeiProMap.
 *
 * Copyright (C) 2025 Henning Woydt <henning.woydt@informatik.uni-heidelberg.de>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#ifndef HEIPROMAP_VERIFY_DISTANCE_ORACLES_H
#define HEIPROMAP_VERIFY_DISTANCE_ORACLES_H

#include <iostream>
#include <vector>

#include "../definitions.h"
#include "matrix_distance_oracle.h"
#include "binary_distance_oracle.h"
#include "division_distance_oracle.h"
#include "stored_division_distance_oracle.h"

namespace HeiProMap {

    inline bool verify_distance_oracles_case(const std::vector<partition_t> &hierarchy,
                                              const std::vector<weight_t> &distance) {
        DistanceOracle matrix_oracle;
        BinaryDistanceOracle binary_oracle;
        DivisionDistanceOracle division_oracle;
        StoredDivisionDistanceOracle stored_division_oracle;

        matrix_oracle.initialize(hierarchy, distance);
        binary_oracle.initialize(hierarchy, distance);
        division_oracle.initialize(hierarchy, distance);
        stored_division_oracle.initialize(hierarchy, distance);

        partition_t k = matrix_oracle.get_k();
        bool ok = true;

        for (partition_t u = 0; u < k && ok; ++u) {
            for (partition_t v = u; v < k && ok; ++v) {
                weight_t d_mat = matrix_oracle.get(u, v);
                weight_t d_bin = binary_oracle.get(u, v);
                weight_t d_div = division_oracle.get(u, v);
                weight_t d_sto = stored_division_oracle.get(u, v);

                if (d_mat != d_bin || d_mat != d_div || d_mat != d_sto) {
                    std::cerr << "Distance oracle mismatch for PEs (" << u << ", " << v << "): "
                              << "matrix=" << d_mat << " binary=" << d_bin
                              << " division=" << d_div << " stored_division=" << d_sto << std::endl;
                    ok = false;
                }

                if (u != v) {
                    partition_t h_mat = matrix_oracle.get_h(u, v);
                    partition_t h_bin = binary_oracle.get_h(u, v);
                    partition_t h_div = division_oracle.get_h(u, v);
                    partition_t h_sto = stored_division_oracle.get_h(u, v);

                    if (h_mat != h_bin || h_mat != h_div || h_mat != h_sto) {
                        std::cerr << "Hierarchy oracle mismatch for PEs (" << u << ", " << v << "): "
                                  << "matrix=" << h_mat << " binary=" << h_bin
                                  << " division=" << h_div << " stored_division=" << h_sto << std::endl;
                        ok = false;
                    }
                }
            }
        }
        return ok;
    }

    inline bool verify_distance_oracles() {
        struct TestCase {
            std::vector<partition_t> hierarchy;
            std::vector<weight_t> distance;
        };

        std::vector<TestCase> cases = {
            {{4, 4},          {1, 10}},             // k=16
            {{2, 3, 5},       {1, 10, 100}},        // k=30
            {{8, 8, 8, 2},    {1, 5, 25, 100}},     // k=1024
        };

        bool ok = true;
        for (const auto &tc : cases) {
            if (!verify_distance_oracles_case(tc.hierarchy, tc.distance)) {
                ok = false;
            }
        }

        if (ok) {
            std::cout << "Distance oracle verification passed for all test cases." << std::endl;
        }
        return ok;
    }

}

#endif //HEIPROMAP_VERIFY_DISTANCE_ORACLES_H
