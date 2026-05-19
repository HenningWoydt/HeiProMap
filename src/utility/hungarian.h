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

#ifndef HEIPROMAP_HUNGARIAN_H
#define HEIPROMAP_HUNGARIAN_H

#include <vector>
#include <limits>
#include <algorithm>
#include "../definitions.h"

namespace Dyn_HeiProMap {
    using namespace HeiProMap;
    /**
     * Simple implementation of the Hungarian algorithm for min-weight perfect matching.
     * cost_matrix: n x n cost matrix. 
     * Returns vector where result[i] is the assignment for the i-th row.
     */
    inline std::vector<int> solve_hungarian(const std::vector<std::vector<weight_t> > &cost_matrix) {
        int n = (int) cost_matrix.size();
        if (n == 0) return {};

        std::vector<weight_t> u(n + 1, 0), v(n + 1, 0);
        std::vector<int> p(n + 1, 0), way(n + 1, 0);

        for (int i = 1; i <= n; ++i) {
            p[0] = i;
            int j0 = 0;
            std::vector<weight_t> minv(n + 1, std::numeric_limits<weight_t>::max());
            std::vector<bool> used(n + 1, false);
            do {
                used[j0] = true;
                int i0 = p[j0], j1 = 0;
                weight_t delta = std::numeric_limits<weight_t>::max();
                for (int j = 1; j <= n; ++j) {
                    if (!used[j]) {
                        weight_t cur = cost_matrix[i0 - 1][j - 1] - u[i0] - v[j];
                        if (cur < minv[j]) {
                            minv[j] = cur;
                            way[j] = j0;
                        }
                        if (minv[j] < delta) {
                            delta = minv[j];
                            j1 = j;
                        }
                    }
                }
                for (int j = 0; j <= n; ++j) {
                    if (used[j]) {
                        u[p[j]] += delta;
                        v[j] -= delta;
                    } else {
                        minv[j] -= delta;
                    }
                }
                j0 = j1;
            } while (p[j0] != 0);
            do {
                int j1 = way[j0];
                p[j0] = p[j1];
                j0 = j1;
            } while (j0);
        }

        std::vector<int> result(n);
        for (int j = 1; j <= n; ++j) {
            result[p[j] - 1] = j - 1;
        }
        return result;
    }
}

#endif // HEIPROMAP_HUNGARIAN_H
