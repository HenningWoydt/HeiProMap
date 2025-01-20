#ifndef HEIDELBERGPROCESSMAPPING_KAFFPA_PARTITIONER_H
#define HEIDELBERGPROCESSMAPPING_KAFFPA_PARTITIONER_H

#include "../../definitions.h"
#include "../../macros.h"
#include "../datastructures/translation_table.h"
#include "../interfaces/ISerialPartitioner.h"
#include "interface/kaHIP_interface.h"

namespace HeiProMap {
    class KaffpaPartitioner final : public ISerialPartitioner {
    public:
        template<typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialPartitionManager>
        void partition(TSerialGraph &g,
                       TSerialActiveVertexManager &av_manager,
                       TSerialPartitionManager &p_manager,
                       const std::vector<partition_t> &hierarchy,
                       const std::vector<weight_t> &distance,
                       const f64 imbalance) {
            // number of vertices and edges
            int n = 0;
            int m = 0;

            // build translation table
            TranslationTable tt;
            vertex_t         translate = 0;
            for (vertex_t    u: av_manager) {
                ASSERT(av_manager.is_active(u));

                tt.add(u, translate);
                translate += 1;

                n += 1;
                m += (int) g.size(u);
            }

            // vertex weights
            int *v_weights = (int *) malloc(n * sizeof(int));
            for (int i = 0; i < n; ++i) { v_weights[i] = (int) g.get_weight(tt.get_o(i)); }

            // pointer to adjacency lists
            int *adj_ptr   = (int *) malloc((n + 1) * sizeof(int));
            int *adj       = (int *) malloc(m * sizeof(int));
            int *e_weights = (int *) malloc(m * sizeof(int));

            // set adj_ptr
            adj_ptr[0] = 0;
            for (int new_u            = 0; new_u < n; ++new_u) {
                vertex_t    old_u      = tt.get_o(new_u);
                int         insert_idx = 0;
                for (size_t i          = 0; i < g.size(old_u); ++i) {
                    vertex_t v                             = g.neighbor(old_u, i);
                    weight_t ew                            = g.get_weight(old_u, i);
                    adj[adj_ptr[new_u] + insert_idx]       = (int) tt.get_n(v);
                    e_weights[adj_ptr[new_u] + insert_idx] = (int) ew;
                    insert_idx += 1;
                }
                adj_ptr[new_u + 1] = adj_ptr[new_u] + insert_idx;
            }
            // imbalance
            double   kaffpa_imbalance = imbalance;

            // hierarchy
            int *kaffpa_hierarchy = (int *) malloc(hierarchy.size() * sizeof(int));
            for (u64 i = 0; i < hierarchy.size(); ++i) { kaffpa_hierarchy[i] = (int) hierarchy[i]; }

            // distance
            int *kaffpa_distance = (int *) malloc(distance.size() * sizeof(int));
            for (u64 i = 0; i < distance.size(); ++i) { kaffpa_distance[i] = (int) distance[i]; }

            // mode
            int kaffpa_map_mode       = MAPMODE_BISECTION; // TODO: Figure out why MAPMODE_MULTISECTION does not work
            int kaffpa_partition_mode = KAFFPA_STRONG;

            // partition result
            int *kaffpa_partition = (int *) malloc(n * sizeof(int));

            int kaffpa_edgecut, kaffpa_qap;

            // execute kaffpa
            process_mapping(&n, v_weights, adj_ptr, e_weights, adj, kaffpa_hierarchy, kaffpa_distance, (int) hierarchy.size(), kaffpa_partition_mode, kaffpa_map_mode, &kaffpa_imbalance, false, 0, &kaffpa_edgecut, &kaffpa_qap, kaffpa_partition);

            // first read partition
            for (int new_u = 0; new_u < n; ++new_u) {
                p_manager.set(tt.get_o(new_u), g.get_weight(tt.get_o(new_u)), kaffpa_partition[new_u]);
            }

            free(v_weights);
            free(adj_ptr);
            free(adj);
            free(e_weights);
            free(kaffpa_hierarchy);
            free(kaffpa_distance);
            free(kaffpa_partition);
        }
    };
}

#endif //HEIDELBERGPROCESSMAPPING_KAFFPA_PARTITIONER_H
