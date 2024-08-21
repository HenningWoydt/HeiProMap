#ifndef SERIALPROCESSMAPPING_KAFFPA_PARTITIONER_H
#define SERIALPROCESSMAPPING_KAFFPA_PARTITIONER_H

#include "../utility/definitions.h"
#include "../utility/macros.h"
#include "../utility/utils.h"
#include "../datastructures/graph.h"
#include "../datastructures/translation_table.h"
#include "../utility/qap.h"

#include "interface/kaHIP_interface.h"

namespace SPM {

    class KaffpaPartitioner {

    private:
        Graph *p_g = nullptr;
        std::vector<u64> hierarchy;
        std::vector<u64> distance;
        u64 k = 0;
        f64 imbalance = 0;

    public:
        KaffpaPartitioner() = default;

        void initialize(Graph *t_g,
                        std::vector<u64> &t_hierarchy,
                        std::vector<u64> &t_distance,
                        u64 t_k,
                        f64 t_imbalance){
            p_g = t_g;
            hierarchy = t_hierarchy;
            distance = t_distance;
            k = t_k;
            imbalance = t_imbalance;
        }

        void partition(PartitionManager &pm,
                       int kaffpa_config) {
            ASSERT(p_g != nullptr);
            Graph &g = *p_g;
            ASSERT(imbalance >= 0.0);

            // build translation table
            TranslationTable tt;
            vertex_t translate = 0;
            for (vertex_t u = 0; u < g.get_n(); ++u) {
                if (g.get_vertex_state(u) == 1) {
                    tt.add(u, translate);
                    translate += 1;
                }
            }

            // number of vertices and edges
            int n = (int) g.get_n_active();
            int m = (int) g.get_m_active();

            // vertex weights
            int *v_weights = (int *) malloc(n * sizeof(int));
            for (int i = 0; i < n; ++i) { v_weights[i] = (int) g.get_vertex_weight(tt.get_o(i)); }

            // pointer to adjacency lists
            int *adj_ptr = (int *) malloc((n + 1) * sizeof(int));
            int *adj = (int *) malloc(m * sizeof(int));
            int *e_weights = (int *) malloc(m * sizeof(int));

            // set adj_ptr
            adj_ptr[0] = 0;
            for (int u = 0; u < n; ++u) {
                int insert_idx = 0;
                for (EdgeW e: g[tt.get_o(u)]) {
                    adj[adj_ptr[u] + insert_idx] = (int) tt.get_n(e.v);
                    e_weights[adj_ptr[u] + insert_idx] = (int) e.w;
                    insert_idx += 1;
                }
                adj_ptr[u + 1] = adj_ptr[u] + insert_idx;
            }
            // imbalance
            double kaffpa_imbalance = imbalance;

            // hierarchy
            int *kaffpa_hierarchy = (int *) malloc(hierarchy.size() * sizeof(int));
            for (u64 i = 0; i < hierarchy.size(); ++i) { kaffpa_hierarchy[i] = (int) hierarchy[i]; }

            // distance
            int *kaffpa_distance = (int *) malloc(distance.size() * sizeof(int));
            for (u64 i = 0; i < distance.size(); ++i) { kaffpa_distance[i] = (int) distance[i]; }

            // mode
            int kaffpa_map_mode = MAPMODE_BISECTION; // TODO: Figure out why MAPMODE_MULTISECTION does not work
            int kaffpa_partition_mode = kaffpa_config;

            // partition result
            int *kaffpa_partition = (int *) malloc(n * sizeof(int));

            int kaffpa_edgecut, kaffpa_qap;

            // execute kaffpa
            process_mapping(&n, v_weights, adj_ptr, e_weights, adj, kaffpa_hierarchy, kaffpa_distance, (int) hierarchy.size(), kaffpa_partition_mode, kaffpa_map_mode, &kaffpa_imbalance, false, 0, &kaffpa_edgecut, &kaffpa_qap, kaffpa_partition);

            // get result
            for (int i = 0; i < n; ++i) {
                pm[tt.get_o(i)] = kaffpa_partition[i];
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

#endif //SERIALPROCESSMAPPING_KAFFPA_PARTITIONER_H
