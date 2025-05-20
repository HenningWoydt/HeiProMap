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

#ifndef HEIPROMAP_KAFFPA_K_WAY_PARTITIONER_H
#define HEIPROMAP_KAFFPA_K_WAY_PARTITIONER_H

#include <string>

#include "kaHIP_interface.h"

#include "../serial_definitions_1.h"
#include "../serial_definitions_2.h"
#include "../serial_definitions_3.h"
#include "../../commons/random_engine.h"
#include "../../commons/statistic_collector.h"

namespace HeiProMap {
    enum KaffpaKWayPartitionerMode {
        KAFFPA_KWAY_PARTITIONER_MODE_UNDEFINED,
        KAFFPA_KWAY_PARTITIONER_MODE_STRONG,
        KAFFPA_KWAY_PARTITIONER_MODE_ECO,
        KAFFPA_KWAY_PARTITIONER_MODE_FAST,
    };

    inline KaffpaKWayPartitionerMode string_to_kaffpa_kway_partitioner_mode(const std::string &str) {
        if (str == "UNDEFINED") return KAFFPA_KWAY_PARTITIONER_MODE_UNDEFINED;
        if (str == "strong") return KAFFPA_KWAY_PARTITIONER_MODE_STRONG;
        if (str == "eco") return KAFFPA_KWAY_PARTITIONER_MODE_ECO;
        if (str == "fast") return KAFFPA_KWAY_PARTITIONER_MODE_FAST;
        return KAFFPA_KWAY_PARTITIONER_MODE_UNDEFINED;
    }

    inline std::string kaffpa_kway_partitioner_mode_to_string(KaffpaKWayPartitionerMode mode) {
        switch (mode) {
            case KAFFPA_KWAY_PARTITIONER_MODE_UNDEFINED:
                return "UNDEFINED";
            case KAFFPA_KWAY_PARTITIONER_MODE_STRONG:
                return "strong";
            case KAFFPA_KWAY_PARTITIONER_MODE_ECO:
                return "eco";
            case KAFFPA_KWAY_PARTITIONER_MODE_FAST:
                return "fast";
            default:
                return "UNDEFINED";
        }
    }

    class KaffpaKWayPartitionerConfiguration {
    public:
        std::string               mode_string;
        KaffpaKWayPartitionerMode mode; // Which mode to use: strong, eco, fast
    };

    class KaffpaKWayPartitioner {
    private:
        // vars needed for kaffpa
        int *vwgt    = nullptr;
        int *xadj    = nullptr;
        int *adjcwgt = nullptr;
        int *adjncy  = nullptr;
        int *part    = nullptr;

        size_t n_vertices = 0;
        size_t n_edges    = 0;

        TranslationTable<vertex_t> tt;

        f64 time_total  = 0.0;
        f64 time_kaffpa = 0.0;

    public:
        ~KaffpaKWayPartitioner() {
            free(vwgt);
            free(xadj);
            free(adjcwgt);
            free(adjncy);
            free(part);
        }

        /**
         * Partitions the subgraph of g where the vertices have the corresponding id.
         *
         * @param g The graph to partition.
         * @param p_manager The partition manager.
         * @param id The id of the subgraph to partition.
         * @param k The number of partitions.
         * @param lmax The maximum weight of a partition.
         * @param t_random_engine The random engine.
         * @param i_config The configuration.
         * @param t_stat_collect The statistic collector.
         */
        void partition(const graph_t &g,
                       deep_p_manager_t &p_manager,
                       deep_bv_manager_t &bv_manager,
                       deep_q_graph_t &q_graph,
                       partition_t id,
                       partition_t id_increment,
                       partition_t k,
                       weight_t lmax,
                       s32 hierarchy_level,
                       RandomEngine &t_random_engine,
                       const KaffpaKWayPartitionerConfiguration &i_config,
                       StatisticCollector &t_stat_collect) {
            auto sp = std::chrono::high_resolution_clock::now();

            tt.reserve(g.get_n(), g.get_n());

            vertex_t subgraph_n_vertices = 0;
            vertex_t subgraph_n_edges    = 0;
            weight_t subgraph_weight     = 0;
            // step 1: extract the number of vertices and edges, build tt
            forall_gu(g, u)
                {
                    if (p_manager[u] == id) {
                        tt.add(u, subgraph_n_vertices);
                        subgraph_n_vertices += 1;
                        subgraph_weight += g.weight(u);
                        forall_guiv(g, u, i, v)
                            {
                                if (p_manager[v] == id) { subgraph_n_edges += 1; }
                            }
                        endfor
                    }
                }
            endfor

            allocate(subgraph_n_vertices, subgraph_n_edges);
            xadj[0] = 0;
            // step 2: build graph in kaffpa variables
            forall_gu(g, u)
                {
                    if (p_manager[u] == id) {
                        vertex_t new_u  = tt.get_n(u);
                        vwgt[new_u]     = g.weight(u);
                        xadj[new_u + 1] = xadj[new_u];
                        forall_guivw(g, u, i, v, w)
                            {
                                if (p_manager[v] == id) {
                                    vertex_t new_v           = tt.get_n(v);
                                    adjncy[xadj[new_u + 1]]  = (int) new_v;
                                    adjcwgt[xadj[new_u + 1]] = w;
                                    xadj[new_u + 1] += 1;
                                }
                            }
                        endfor
                    }
                }
            endfor

            // step 3: calculate allowed imbalance
            f64 imbalance = ((f64) (lmax * k) / (f64) subgraph_weight) - 1.0;
            if (imbalance <= 0.0) { // Failsafe if subgraph weight is too large // TODO: would be cool if not needed
                imbalance = 0.03;
            }

            // step 4: partition
            int n       = (int) subgraph_n_vertices;
            int n_parts = (int) k;
            int edge_cut;

            auto sp_kaffpa = std::chrono::high_resolution_clock::now();
            kaffpa(&n, vwgt, xadj, adjcwgt, adjncy, &n_parts, &imbalance, true, 0, FAST, &edge_cut, part);
            auto ep_kaffpa = std::chrono::high_resolution_clock::now();

            // step 5: extract partition into p_manager
            for (vertex_t kaffpa_vertex = 0; kaffpa_vertex < subgraph_n_vertices; ++kaffpa_vertex) {
                if (part[kaffpa_vertex] == 0) { continue; }
                partition_t move_id = id + id_increment * part[kaffpa_vertex];

                partition_t vertex_id     = id;
                vertex_t    vertex        = tt.get_o(kaffpa_vertex);
                weight_t    vertex_weight = g.weight(vertex);

                bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                p_manager.move(vertex, vertex_weight, vertex_id, move_id);
            }

            for (partition_t i = 0; i < k; ++i) {
                partition_t move_id = id + id_increment * i;
                p_manager.set_lmax(move_id, lmax);
                p_manager.set_hierarchy_level(move_id, hierarchy_level - 1);
            }

            auto ep = std::chrono::high_resolution_clock::now();

            time_total += get_seconds(sp, ep);
            time_kaffpa += get_seconds(sp_kaffpa, ep_kaffpa);
        }

        void allocate(size_t new_n_vertices, size_t new_n_edges) {
            if (new_n_vertices > n_vertices) {
                n_vertices = new_n_vertices;

                free(vwgt);
                free(xadj);
                free(part);

                vwgt = (int *) malloc(n_vertices * sizeof(int));
                xadj = (int *) malloc((n_vertices + 1) * sizeof(int));
                part = (int *) malloc(n_vertices * sizeof(int));
            }

            if (new_n_edges > n_edges) {
                n_edges = new_n_edges;

                free(adjcwgt);
                free(adjncy);

                adjcwgt = (int *) malloc(n_edges * sizeof(int));
                adjncy  = (int *) malloc(n_edges * sizeof(int));
            }
        }
    };
}

#endif //HEIPROMAP_KAFFPA_K_WAY_PARTITIONER_H
