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

#ifndef HEIPROMAP_KAFFPA_PARTITIONER_H
#define HEIPROMAP_KAFFPA_PARTITIONER_H

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>

#include <kaHIP_interface.h>

#include "../definitions.h"
#include "../datastructures/csr_graph.h"
#include "../utility/aligned_array.h"
#include "../utility/assert_state.h"


namespace HeiProMap {
    enum KaffpaPartitionMode {
        KAFFPA_PARTITION_STRONG,
        KAFFPA_PARTITION_ECO,
        KAFFPA_PARTITION_FAST
    };

    inline KaffpaPartitionMode string_to_kaffpa_partition_mode(const std::string &str) {
        if (str == "strong") return KAFFPA_PARTITION_STRONG;
        if (str == "eco") return KAFFPA_PARTITION_ECO;
        if (str == "fast") return KAFFPA_PARTITION_FAST;
        return KAFFPA_PARTITION_FAST;
    }

    inline void kaffpa_partition(graph_t &g,
                                 partition_t k,
                                 f64 imb,
                                 KaffpaPartitionMode kaffpa_mode,
                                 u64 seed,
                                 AlignedArray<partition_t> &partition,
                                 u64 kappa,
                                 bool collect_dataset = false,
                                 const std::string &data_dir = "data") {
        ASSERT(assert_graph(g));

        u64 g_hash = 0;
        if (collect_dataset) {
            g_hash = g.hash();
            std::filesystem::create_directories(data_dir + "/graphs");
            std::filesystem::create_directories(data_dir + "/results");
            std::filesystem::create_directories(data_dir + "/partitions");

            std::string graph_path = data_dir + "/graphs/" + std::to_string(g_hash) + ".graph";
            if (!std::filesystem::exists(graph_path)) {
                g.write_graph(graph_path);
            }
        }

        int n = (int) g.n;
        int m = (int) g.m;
        int *vwgt = (int *) malloc(n * sizeof(int));
        int *xadj = (int *) malloc((n + 1) * sizeof(int));
        int *adjcwgt = (int *) malloc(m * sizeof(int));
        int *adjncy = (int *) malloc(m * sizeof(int));
        int nparts = (int) k;
        double imbalance = std::max(imb, 0.0);
        bool suppress_output = true;
        int mode = FAST;
        int edge_cut;
        int best_edge_cut = std::numeric_limits<int>::max();
        int *part = (int *) malloc(n * sizeof(int));

        if (kaffpa_mode == KAFFPA_PARTITION_STRONG) { mode = STRONG; }
        if (kaffpa_mode == KAFFPA_PARTITION_ECO) { mode = ECO; }
        if (kaffpa_mode == KAFFPA_PARTITION_FAST) { mode = FAST; }


        // set vertex weights
        for (int i = 0; i < n; i++) { vwgt[i] = (int) g.v_weights[i]; }

        // set xadj
        for (int i = 0; i < n + 1; i++) { xadj[i] = (int) g.neighborhoods[i]; }

        // set adjncy
        for (int i = 0; i < m; i++) { adjncy[i] = (int) g.edges_v[i]; }

        // set adjcwgt
        for (int i = 0; i < m; i++) { adjcwgt[i] = (int) g.edges_w[i]; }

        auto start_time = get_time_point();

        for (u64 j = 0; j < kappa; j++) {
            kaffpa(&n, vwgt, xadj, adjcwgt, adjncy, &nparts, &imbalance, suppress_output, (int) seed + j, mode, &edge_cut, part);

            if (edge_cut < best_edge_cut) {
                best_edge_cut = edge_cut;
                for (int i = 0; i < n; i++) { partition[i] = part[i]; }
            }
        }

        auto end_time = get_time_point();
        f64 time_ms = get_milli_seconds(start_time, end_time);

        if (collect_dataset) {
            std::string base_name = std::to_string(g_hash) + "_" + std::to_string(k) + "_" + std::to_string(seed);

            // Calculate max block weight
            std::vector<weight_t> block_weights(k, 0);
            for (int i = 0; i < n; i++) {
                block_weights[partition[i]] += g.v_weights[i];
            }
            weight_t max_block_weight = 0;
            for (partition_t i = 0; i < k; i++) {
                if (block_weights[i] > max_block_weight) {
                    max_block_weight = block_weights[i];
                }
            }

            // Save partition
            std::string part_path = data_dir + "/partitions/" + base_name + ".partition";
            std::ofstream p_out(part_path);
            for (int i = 0; i < n; i++) {
                p_out << partition[i] << "\n";
            }

            // Save metadata
            std::string json_path = data_dir + "/results/" + base_name + ".json";
            std::ofstream j_out(json_path);
            j_out << "{\n";
            j_out << "  \"graph_hash\": " << g_hash << ",\n";
            j_out << "  \"k\": " << k << ",\n";
            j_out << "  \"imbalance\": " << imbalance << ",\n";
            j_out << "  \"max_block_weight\": " << max_block_weight << ",\n";
            j_out << "  \"graph_weight\": " << g.g_weight << ",\n";
            j_out << "  \"kaffpa_mode\": " << (int) kaffpa_mode << ",\n";
            j_out << "  \"seed\": " << seed << ",\n";
            j_out << "  \"kappa\": " << kappa << ",\n";
            j_out << "  \"edge_cut\": " << best_edge_cut << ",\n";
            j_out << "  \"time_ms\": " << time_ms << "\n";
            j_out << "}\n";
        }

        free(vwgt);
        free(xadj);
        free(adjcwgt);
        free(adjncy);
        free(part);
    }
}

#endif //HEIPROMAP_KAFFPA_PARTITIONER_H
