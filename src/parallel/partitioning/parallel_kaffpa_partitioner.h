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

#ifndef HEIPROMAP_PARALLEL_KAFFPA_PARTITIONER_H
#define HEIPROMAP_PARALLEL_KAFFPA_PARTITIONER_H

#include "../interfaces/IParallelPartitioner.h"
#include "../../definitions.h"
#include "../parallel_definitions_1.h"

#include "interface/kaHIP_interface.h"
#include "../../commons/translation_table.h"

namespace HeiProMap {

    enum ParallelKaffpaPartitionerMode {
        PARALLEL_KAFFPA_PARTITIONER_MODE_UNDEFINED,
        PARALLEL_KAFFPA_PARTITIONER_MODE_STRONG,
        PARALLEL_KAFFPA_PARTITIONER_MODE_ECO,
        PARALLEL_KAFFPA_PARTITIONER_MODE_FAST,
    };

    inline ParallelKaffpaPartitionerMode parallel_string_to_kaffpa_partitioner_mode(const std::string &str) {
        if (str == "UNDEFINED") return PARALLEL_KAFFPA_PARTITIONER_MODE_UNDEFINED;
        if (str == "strong") return PARALLEL_KAFFPA_PARTITIONER_MODE_STRONG;
        if (str == "eco") return PARALLEL_KAFFPA_PARTITIONER_MODE_ECO;
        if (str == "fast") return PARALLEL_KAFFPA_PARTITIONER_MODE_FAST;
        return PARALLEL_KAFFPA_PARTITIONER_MODE_UNDEFINED;
    }

    inline std::string parallel_kaffpa_partitioner_mode_to_string(ParallelKaffpaPartitionerMode mode) {
        switch (mode) {
            case PARALLEL_KAFFPA_PARTITIONER_MODE_UNDEFINED:
                return "UNDEFINED";
            case PARALLEL_KAFFPA_PARTITIONER_MODE_STRONG:
                return "strong";
            case PARALLEL_KAFFPA_PARTITIONER_MODE_ECO:
                return "eco";
            case PARALLEL_KAFFPA_PARTITIONER_MODE_FAST:
                return "fast";
            default:
                return "UNDEFINED";
        }
    }

    enum ParallelKaffpaPartitionerMethod {
        PARALLEL_KAFFPA_PARTITIONER_METHOD_UNDEFINED,
        PARALLEL_KAFFPA_PARTITIONER_METHOD_BISECTION,
        PARALLEL_KAFFPA_PARTITIONER_METHOD_MULTISECTION,
    };

    inline ParallelKaffpaPartitionerMethod parallel_string_to_kaffpa_partitioner_method(const std::string &str) {
        if (str == "UNDEFINED") return PARALLEL_KAFFPA_PARTITIONER_METHOD_UNDEFINED;
        if (str == "bisection") return PARALLEL_KAFFPA_PARTITIONER_METHOD_BISECTION;
        if (str == "multisection") return PARALLEL_KAFFPA_PARTITIONER_METHOD_MULTISECTION;
        return PARALLEL_KAFFPA_PARTITIONER_METHOD_UNDEFINED;
    }

    inline std::string parallel_kaffpa_partitioner_method_to_string(ParallelKaffpaPartitionerMethod mode) {
        switch (mode) {
            case PARALLEL_KAFFPA_PARTITIONER_METHOD_UNDEFINED:
                return "UNDEFINED";
            case PARALLEL_KAFFPA_PARTITIONER_METHOD_BISECTION:
                return "bisection";
            case PARALLEL_KAFFPA_PARTITIONER_METHOD_MULTISECTION:
                return "multisection";
            default:
                return "UNDEFINED";
        }
    }

    class ParallelKaffpaPartitionerConfiguration final : public IParallelPartitionerConfiguration {
    public:
        std::string                     mode_string;
        ParallelKaffpaPartitionerMode   mode; // Which mode to use: strong, eco, fast
        std::string                     method_string;
        ParallelKaffpaPartitionerMethod method; // Which method to use: bisection, multisection
    };

    class KaffpaPartitioner final : public IParallelPartitioner {
    public:
        void partition(IParallelPartitionerConfiguration &i_config,
                       p_graph_t &g,
                       p_av_manager_t &av_manager,
                       p_p_manager_t &p_manager,
                       const std::vector<partition_t> &hierarchy,
                       const std::vector<weight_t> &distance,
                       const f64 imbalance,
                       u64 seed) {
            ParallelKaffpaPartitionerConfiguration &config = dynamic_cast<ParallelKaffpaPartitionerConfiguration&>(i_config);

            // number of vertices and edges
            int n = 0;
            int m = 0;

            // build translation table
            TranslationTable<int> tt;
            tt.reserve(av_manager.get_n_active(), g.get_n());
            vertex_t      translate = 0;
            for (vertex_t u: av_manager) {
                ASSERT(av_manager.is_active(u));

                tt.add(u, translate);
                translate += 1;

                n += 1;
                m += (int) g.size(u);
            }

            // vertex weights
            int      *v_weights = (int *) malloc(n * sizeof(int));
            for (int i          = 0; i < n; ++i) { v_weights[i] = (int) g.get_weight(tt.get_o(i)); }

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
            int      *kaffpa_hierarchy = (int *) malloc(hierarchy.size() * sizeof(int));
            for (u64 i                 = 0; i < hierarchy.size(); ++i) { kaffpa_hierarchy[i] = (int) hierarchy[i]; }

            // distance
            int      *kaffpa_distance = (int *) malloc(distance.size() * sizeof(int));
            for (u64 i                = 0; i < distance.size(); ++i) { kaffpa_distance[i] = (int) distance[i]; }

            // mode
            int kaffpa_map_mode;
            int kaffpa_partition_mode;

            if (config.method == PARALLEL_KAFFPA_PARTITIONER_METHOD_BISECTION) {
                kaffpa_map_mode = MAPMODE_BISECTION;
            } else if (config.method == PARALLEL_KAFFPA_PARTITIONER_METHOD_MULTISECTION) {
                kaffpa_map_mode = MAPMODE_MULTISECTION; // TODO: Figure out why MAPMODE_MULTISECTION does not work
            } else {
                std::cout << "KaFFPa Partitioning method " << parallel_kaffpa_partitioner_method_to_string(config.method) << " with id " << config.method << " not known!" << std::endl;
                exit(EXIT_FAILURE);
            }

            if (config.mode == PARALLEL_KAFFPA_PARTITIONER_MODE_FAST) {
                kaffpa_partition_mode = FAST;
            } else if (config.mode == PARALLEL_KAFFPA_PARTITIONER_MODE_ECO) {
                kaffpa_partition_mode = ECO;
            } else if (config.mode == PARALLEL_KAFFPA_PARTITIONER_MODE_STRONG) {
                kaffpa_partition_mode = STRONG;
            } else {
                std::cout << "KaFFPa Partitioning mode " << parallel_kaffpa_partitioner_mode_to_string(config.mode) << " with id " << config.mode << " not known!" << std::endl;
                exit(EXIT_FAILURE);
            }

            // partition result
            int *kaffpa_partition = (int *) malloc(n * sizeof(int));

            int kaffpa_edgecut, kaffpa_qap;

            // std::string file_path = "error.graph";
            // write_graph(n, m/2, v_weights, adj_ptr, e_weights, adj, file_path);

            // execute kaffpa
            process_mapping(&n, v_weights, adj_ptr, e_weights, adj, kaffpa_hierarchy, kaffpa_distance, (int) hierarchy.size(), kaffpa_partition_mode, kaffpa_map_mode, &kaffpa_imbalance, true, (int) seed, &kaffpa_edgecut, &kaffpa_qap, kaffpa_partition);

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

        void write_graph(int n, int m, int *vwgt, int *xadj, int *adjwgt, int *adjncy, std::string &file_path) {
            std::ofstream file(file_path);

            file << n << " " << m << " 011" << std::endl;
            for (int i = 0; i < n; ++i) {
                file << vwgt[i] << " ";
                for (int j = xadj[i]; j < xadj[i + 1]; ++j) {
                    file << adjncy[j] + 1 << " " << adjwgt[j] << " ";
                }
                file << std::endl;
            }
        }
    };
}

#endif //HEIPROMAP_PARALLEL_KAFFPA_PARTITIONER_H
