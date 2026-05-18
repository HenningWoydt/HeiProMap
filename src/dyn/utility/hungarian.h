#ifndef DYN_HEIPROMAP_HUNGARIAN_H
#define DYN_HEIPROMAP_HUNGARIAN_H

#include <vector>
#include <limits>
#include <algorithm>
#include "../../definitions.h"

namespace HeiProMap {
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

#endif // DYN_HEIPROMAP_HUNGARIAN_H
