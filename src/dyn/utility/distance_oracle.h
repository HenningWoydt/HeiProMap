#ifndef DYN_HEIPROMAP_DISTANCE_ORACLE_H
#define DYN_HEIPROMAP_DISTANCE_ORACLE_H

#include <vector>
#include "../../definitions.h"

namespace HeiProMap {
    class DistanceOracle {
        std::vector<std::vector<weight_t> > distance_matrix;

    public:
        DistanceOracle() = default;

        DistanceOracle(const std::vector<partition_t> &hierarchy, const std::vector<vertex_t> &distances) {
            u64 n = 1;
            for (auto h: hierarchy) n *= h;

            distance_matrix.assign(n, std::vector<weight_t>(n, 0));

            for (u64 i = 0; i < n; ++i) {
                for (u64 j = 0; j < n; ++j) {
                    if (i == j) {
                        distance_matrix[i][j] = 0;
                    } else {
                        // Decompose
                        std::vector<partition_t> p1 = get_coords(i, hierarchy);
                        std::vector<partition_t> p2 = get_coords(j, hierarchy);

                        for (int l = (int) hierarchy.size() - 1; l >= 0; --l) {
                            if (p1[l] != p2[l]) {
                                distance_matrix[i][j] = (weight_t) distances[l];
                                break;
                            }
                        }
                    }
                }
            }
        }

        weight_t query(partition_t c1, partition_t c2) const {
            if (c1 >= distance_matrix.size() || c2 >= distance_matrix.size()) return 0;
            return distance_matrix[c1][c2];
        }

    private:
        static std::vector<partition_t> get_coords(partition_t c, const std::vector<partition_t> &hierarchy) {
            std::vector<partition_t> coords(hierarchy.size());
            partition_t rest = c;
            for (size_t i = 0; i < hierarchy.size(); ++i) {
                coords[i] = rest % hierarchy[i];
                rest /= hierarchy[i];
            }
            return coords;
        }
    };
}

#endif //DYN_HEIPROMAP_DISTANCE_ORACLE_H
