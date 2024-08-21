#include "qap.h"

namespace SPM {

    void determine_loc(u64 u_id,
                       u64 k,
                       const std::vector<u64> &hierarchy,
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

    u64 determine_distance(u64 u_id,
                           u64 v_id,
                           u64 k,
                           const std::vector<u64> &hierarchy,
                           const std::vector<u64> &distance,
                           std::vector<u64> &u_loc,
                           std::vector<u64> &v_loc) {
        ASSERT(prod<u64>(hierarchy) == k);
        ASSERT(u_id < k);
        ASSERT(v_id < k);
        ASSERT(!hierarchy.empty());
        ASSERT(hierarchy.size() == distance.size());
        ASSERT(u_loc.size() == hierarchy.size());
        ASSERT(v_loc.size() == hierarchy.size());

        // special case
        if (u_id == v_id) {
            return 0;
        }

        // determine the location of both partitions
        determine_loc(u_id, k, hierarchy, u_loc);
        determine_loc(v_id, k, hierarchy, v_loc);

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

    u64 get_qap(Graph &g,
                PartitionManager &pm,
                std::vector<u64> &hierarchy,
                std::vector<u64> &distance) {
        u64 k = 1;
        for (u64 x: hierarchy) { k *= x; }
        std::vector<u64> u_loc(hierarchy.size());
        std::vector<u64> v_loc(hierarchy.size());

        u64 qap = 0;

        for (ActiveVertexIterator avi(g); avi.not_end(); avi.next()) {
            vertex_t u = avi.get();
            vertex_t u_p_id = pm[u];

            for (EdgeW &e: g[u]) {
                vertex_t v_p_id = pm[e.v];
                u64 d = determine_distance(u_p_id, v_p_id, k, hierarchy, distance, u_loc, v_loc);
                qap += (d * e.w);
            }
        }

        return qap;
    }

    u64 get_qap(Graph &g,
                vertex_t u,
                std::vector<vertex_t> &partition,
                std::vector<u64> &hierarchy,
                u64 k,
                std::vector<u64> &distance,
                std::vector<u64> &u_loc,
                std::vector<u64> &v_loc) {
        u64 qap = 0;
        vertex_t u_p_id = partition[u];

        for (EdgeW &e: g[u]) {
            vertex_t v_p_id = partition[e.v];
            u64 d = determine_distance(u_p_id, v_p_id, k, hierarchy, distance, u_loc, v_loc);
            qap += (d * e.w);
        }

        return qap;
    }

    s64 get_u_qap(Graph &g,
                  vertex_t u,
                  PartitionManager &pm,
                  DistanceOracle &dist_o){
        s64 qap = 0;
        partition_t u_id = pm[u];

        for (EdgeW &e: g[u]) {
            partition_t v_id = pm[e.v];
            u64 d = dist_o.get(u_id, v_id);
            qap += (s64) (d * e.w);
        }

        return qap;
    }

    s64 get_u_qap(const Graph &g,
                  vertex_t u,
                  partition_t new_u_id,
                  const PartitionManager &pm,
                  const DistanceOracle &dist_o){

        s64 qap = 0;
        for (const EdgeW &e: g[u]) {
            partition_t v_id = pm[e.v];
            u64 d = dist_o.get(new_u_id, v_id);
            qap += (s64) (d * e.w);
        }

        return qap;
    }

    void get_u_qap(Graph &g,
                  vertex_t u,
                  const std::vector<partition_t> &moves_to_check,
                  std::vector<u64> &qap_deltas,
                  PartitionManager &pm,
                  DistanceOracle &dist_o){

        qap_deltas.resize(moves_to_check.size());
        std::fill(qap_deltas.begin(), qap_deltas.end(), 0);

        for (EdgeW &e: g[u]) {
            partition_t v_id = pm[e.v];

            for(size_t i = 0; i < moves_to_check.size(); ++i){
                qap_deltas[i] += (s64) (dist_o.get(moves_to_check[i], v_id) * e.w);
            }
        }
    }

    void get_pweights(Graph &g,
                      PartitionManager &pm,
                      std::vector<u64> &pweights) {
        std::fill(pweights.begin(), pweights.end(), 0);
        for (ActiveVertexIterator avi(g); avi.not_end(); avi.next()) {
            vertex_t u = avi.get();
            pweights[pm[u]] += g.get_vertex_weight(u);
        }
    }

    u64 get_lmax(f64 imbalance, u64 k, u64 g_weight) {
        return ceil((1.0 + imbalance) * ((f64) g_weight / (f64) k));
    }

}
