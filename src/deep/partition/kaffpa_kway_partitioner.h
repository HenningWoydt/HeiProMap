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

#include <atomic>
#include <string>

#include "kaHIP_interface.h"

#include "../../serial/serial_definitions_1.h"
#include "../../serial/serial_definitions_2.h"
#include "../../serial/serial_definitions_3.h"
#include "../../commons/random_engine.h"
#include "../../commons/statistic_collector.h"

namespace HeiProMap {
    enum KaffpaKWayPartitionerMode {
        KAFFPA_KWAY_PARTITIONER_MODE_UNDEFINED,
        KAFFPA_KWAY_PARTITIONER_MODE_STRONG,
        KAFFPA_KWAY_PARTITIONER_MODE_ECO,
        KAFFPA_KWAY_PARTITIONER_MODE_FAST,
    };

    inline KaffpaKWayPartitionerMode string_to_kaffpa_kway_partitioner_mode(const std::string& str) {
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
        std::string mode_string;
        KaffpaKWayPartitionerMode mode; // Which mode to use: strong, eco, fast
    };

    class KaffpaKWayPartitioner {
        vertex_t m_n    = 0;
        partition_t m_k = 0;
        u64 m_threads   = 1;

        struct KaffpaVars {
            // vars needed for kaffpa
            int* vwgt    = nullptr;
            int* xadj    = nullptr;
            int* adjcwgt = nullptr;
            int* adjncy  = nullptr;
            int* part    = nullptr;

            size_t n_vertices = 0;
            size_t n_edges    = 0;

            TranslationTable<vertex_t> tt;
        };

        std::vector<KaffpaVars> m_kaffpa_vars;

        std::atomic<bool> blocks_up_to_date = false;
        std::vector<std::vector<vertex_t>> blocks;

        std::mutex m_mutex;

        RandomEngine* random_engine                      = nullptr;
        const KaffpaKWayPartitionerConfiguration* config = nullptr;

    public:
        void initialize(const vertex_t t_n,
                        const partition_t t_k,
                        const u64 t_threads,
                        RandomEngine& t_random_engine,
                        const KaffpaKWayPartitionerConfiguration& t_config) {
            m_n       = t_n;
            m_k       = t_k;
            m_threads = t_threads;

            m_kaffpa_vars.resize(m_threads);

            random_engine    = &t_random_engine;
            config           = dynamic_cast<const KaffpaKWayPartitionerConfiguration*>(&t_config);
        }

        ~KaffpaKWayPartitioner() {
            for (auto& var : m_kaffpa_vars) {
                free(var.vwgt);
                free(var.xadj);
                free(var.adjcwgt);
                free(var.adjncy);
                free(var.part);
            }
        }

        /**
         * Partitions the subgraph of g where the vertices have the corresponding id.
         *
         * @param g The graph to partition.
         * @param p_manager The partition manager.
         * @param id The id of the subgraph to partition.
         * @param k The number of partitions.
         * @param lmax The maximum weight of a partition.
         */
        void partition(const graph_t& g,
                       deep_p_manager_t& p_manager,
                       deep_bv_manager_t& bv_manager,
                       deep_q_graph_t& q_graph,
                       u64 thread_id,
                       partition_t id,
                       partition_t id_increment,
                       partition_t k,
                       weight_t lmax,
                       partition_t hierarchy_level) {
            KaffpaVars& var = m_kaffpa_vars[thread_id];

            var.tt.reserve(g.get_n(), g.get_n());

            vertex_t subgraph_n_vertices = 0;
            vertex_t subgraph_n_edges    = 0;
            weight_t subgraph_weight     = 0;
            // step 1: extract the number of vertices and edges, build tt
            forall_gu(g, u)
                {
                    if (p_manager[u] == id) {
                        var.tt.add(u, subgraph_n_vertices);
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

            allocate(subgraph_n_vertices, subgraph_n_edges, thread_id);
            var.xadj[0] = 0;
            // step 2: build graph in kaffpa variables
            forall_gu(g, u)
                {
                    if (p_manager[u] == id) {
                        vertex_t new_u      = var.tt.get_n(u);
                        var.vwgt[new_u]     = g.weight(u);
                        var.xadj[new_u + 1] = var.xadj[new_u];
                        forall_guivw(g, u, i, v, w)
                            {
                                if (p_manager[v] == id) {
                                    vertex_t new_v                   = var.tt.get_n(v);
                                    var.adjncy[var.xadj[new_u + 1]]  = (int)new_v;
                                    var.adjcwgt[var.xadj[new_u + 1]] = w;
                                    var.xadj[new_u + 1] += 1;
                                }
                            }
                        endfor
                    }
                }
            endfor

            // step 3: calculate allowed imbalance
            f64 imbalance = ((f64)(lmax * k) / (f64)subgraph_weight) - 1.0;
            if (imbalance <= 0.0) {
                // Failsafe if subgraph weight is too large // TODO: would be cool if not needed
                // std::cout << "failsafe " << imbalance << " " << lmax << " " << k << " " << subgraph_weight << std::endl;
                imbalance = 0.0001;
            }

            // step 4: partition
            int n       = (int)subgraph_n_vertices;
            int n_parts = (int)k;
            int edge_cut;

            kaffpa(&n, var.vwgt, var.xadj, var.adjcwgt, var.adjncy, &n_parts, &imbalance, true, 0, FAST, &edge_cut, var.part);

            m_mutex.lock();
            // step 5: extract partition into p_manager
            for (vertex_t kaffpa_vertex = 0; kaffpa_vertex < subgraph_n_vertices; ++kaffpa_vertex) {
                if (var.part[kaffpa_vertex] == 0) { continue; }
                partition_t move_id = id + id_increment * var.part[kaffpa_vertex];

                partition_t vertex_id  = id;
                vertex_t vertex        = var.tt.get_o(kaffpa_vertex);
                weight_t vertex_weight = g.weight(vertex);

                bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                p_manager.move(vertex, vertex_weight, vertex_id, move_id);
            }
            m_mutex.unlock();

            for (partition_t i = 0; i < k; ++i) {
                partition_t move_id = id + id_increment * i;
                p_manager.set_lmax(move_id, lmax);
                p_manager.set_hierarchy_level(move_id, hierarchy_level - 1);
            }

            blocks_up_to_date = false;
        }

        void allocate(size_t new_n_vertices,
                      size_t new_n_edges,
                      u64 thread_id) {
            KaffpaVars& var = m_kaffpa_vars[thread_id];

            if (new_n_vertices > var.n_vertices) {
                var.n_vertices = new_n_vertices;

                free(var.vwgt);
                free(var.xadj);
                free(var.part);

                var.vwgt = (int*)malloc(var.n_vertices * sizeof(int));
                var.xadj = (int*)malloc((var.n_vertices + 1) * sizeof(int));
                var.part = (int*)malloc(var.n_vertices * sizeof(int));
            }

            if (new_n_edges > var.n_edges) {
                var.n_edges = new_n_edges;

                free(var.adjcwgt);
                free(var.adjncy);

                var.adjcwgt = (int*)malloc(var.n_edges * sizeof(int));
                var.adjncy  = (int*)malloc(var.n_edges * sizeof(int));
            }
        }

        void determine_all_blocks(const graph_t& g,
                                  deep_p_manager_t& p_manager) {
            for (partition_t id = 0; id < std::min((partition_t)blocks.size(), m_k); ++id) {
                blocks[id].clear();
            }
            if (blocks.size() < m_k) {
                blocks.resize(m_k);
            }

            forall_gu(g, u)
                {
                    partition_t u_id = p_manager[u];
                    blocks[u_id].push_back(u);
                }
            endfor
        }
    };
}

#endif //HEIPROMAP_KAFFPA_K_WAY_PARTITIONER_H
