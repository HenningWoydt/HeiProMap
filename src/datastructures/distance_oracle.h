#ifndef SERIALPROCESSMAPPING_DISTANCE_ORACLE_H
#define SERIALPROCESSMAPPING_DISTANCE_ORACLE_H

#include "../utility/definitions.h"
#include "../utility/macros.h"
#include "../utility/utils.h"

namespace SPM {

    class DistanceOracle{
    private:
        partition_t k = 0;
        std::vector<u64> hierarchy;
        std::vector<u64> distance;

        std::vector<u64> mtx;
        std::vector<u8> h_mtx;

    public:
        DistanceOracle() = default;

        DistanceOracle(u64 k,
                       std::vector<u64> &hierarchy,
                       std::vector<u64> &distance) : k(k), hierarchy(hierarchy), distance(distance) {
            mtx.resize(k*k);
            h_mtx.resize(k*k);

            std::vector<std::vector<u64>> locs(k, std::vector<u64>(hierarchy.size()));
            for(vertex_t u = 0; u < k; ++u){
                determine_loc(u, locs[u]);
            }

            for(vertex_t u_id = 0; u_id < k; ++u_id){
                mtx[u_id * k + u_id] = 0;
                h_mtx[u_id * k + u_id] = 0;
                for(vertex_t v_id = u_id + 1; v_id < k; ++v_id){
                    u64 d = determine_distance(locs[u_id], locs[v_id]);
                    mtx[u_id * k + v_id] = d;
                    mtx[v_id * k + u_id] = d;

                    u8 h = determine_hierarchy(locs[u_id], locs[v_id]);
                    h_mtx[u_id * k + v_id] = h;
                    h_mtx[v_id * k + u_id] = h;
                }
            }
        }

        u64 get(vertex_t u_id, vertex_t v_id) const {
            ASSERT(u_id < k);
            ASSERT(v_id < k);
            ASSERT(u_id * k + v_id < k*k);
            return mtx[u_id * k + v_id];
        }

        u8 get_h(vertex_t u_id, vertex_t v_id) const {
            ASSERT(u_id < k);
            ASSERT(v_id < k);
            ASSERT(u_id * k + v_id < k*k);
            return h_mtx[u_id * k + v_id];
        }

    private:
        void determine_loc(u64 u_id,
                           std::vector<u64> &u_loc) {
            ASSERT(prod<u64>(hierarchy) == k);
            ASSERT(u_id < k);
            ASSERT(!hierarchy.empty());
            ASSERT(u_loc.size() == hierarchy.size());

            u64 r_start = 0;
            u64 r_end = k;

            u64 s = hierarchy.size();
            for (u64 i = 0; i < hierarchy.size(); ++i) {
                u64 n_parts = hierarchy[s - 1 - i];
                u64 add = (r_end - r_start) / n_parts;

                for (u64 j = 0; j < n_parts; ++j) {
                    if (r_start <= u_id && u_id < r_start + add) {
                        // we have found the part of the partition
                        u_loc[s - 1 - i] = j;
                        r_end = r_start + add;
                        break;
                    } else {
                        r_start += add;
                    }
                }
            }
        }

        u64 determine_distance(std::vector<u64> &u_loc,
                               std::vector<u64> &v_loc) {
            ASSERT(prod<u64>(hierarchy) == k);
            ASSERT(!hierarchy.empty());
            ASSERT(hierarchy.size() == distance.size());
            ASSERT(u_loc.size() == hierarchy.size());
            ASSERT(v_loc.size() == hierarchy.size());

            // determine the distance
            u64 s = hierarchy.size();
            for (u64 i = 0; i < hierarchy.size(); ++i) {
                if (u_loc[s - 1 - i] != v_loc[s - 1 - i]) {
                    return distance[s - 1 - i];
                }
            }
            // unreachable
            abort();
        }

        u8 determine_hierarchy(std::vector<u64> &u_loc,
                               std::vector<u64> &v_loc) {
            ASSERT(prod<u64>(hierarchy) == k);
            ASSERT(!hierarchy.empty());
            ASSERT(hierarchy.size() == distance.size());
            ASSERT(u_loc.size() == hierarchy.size());
            ASSERT(v_loc.size() == hierarchy.size());

            // determine the distance
            u64 s = hierarchy.size();
            for (u64 i = 0; i < hierarchy.size(); ++i) {
                if (u_loc[s - 1 - i] != v_loc[s - 1 - i]) {
                    return s - 1 - i;
                }
            }
            // unreachable
            abort();
        }

    };

}

#endif //SERIALPROCESSMAPPING_DISTANCE_ORACLE_H
