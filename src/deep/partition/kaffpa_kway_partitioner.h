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

#include "../../commons/random_engine.h"
#include "../../src/commons/small_translation_table.h"

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
        std::string mode_string;
        KaffpaKWayPartitionerMode mode; // Which mode to use: strong, eco, fast
    };

    class KaffpaKWayPartitioner {
        vertex_t m_n = 0;
        partition_t m_k = 0;
        u64 m_threads = 1;

        struct KaffpaVars {
            // vars needed for kaffpa
            int *vwgt = nullptr;
            int *xadj = nullptr;
            int *adjcwgt = nullptr;
            int *adjncy = nullptr;
            int *part = nullptr;
            int *best_part = nullptr;

            int edge_cut = 0;
            int best_edge_cut = 0;

            size_t n_vertices = 0;
            size_t n_edges = 0;

            SmallTranslationTable<vertex_t> tt;
        };

        std::vector<KaffpaVars> m_kaffpa_vars;

        std::vector<std::vector<vertex_t> > blocks;

        RandomEngine *random_engine = nullptr;
        const KaffpaKWayPartitionerConfiguration *config = nullptr;

    public:
        void initialize(const vertex_t t_n,
                        const partition_t t_k,
                        const u64 t_threads,
                        RandomEngine &t_random_engine,
                        const KaffpaKWayPartitionerConfiguration &t_config) {
            m_n = t_n;
            m_k = t_k;
            m_threads = t_threads;

            m_kaffpa_vars.resize(m_threads);

            blocks.resize(m_k);

            random_engine = &t_random_engine;
            config = dynamic_cast<const KaffpaKWayPartitionerConfiguration *>(&t_config);
        }

        ~KaffpaKWayPartitioner() {
            for (auto &var: m_kaffpa_vars) {
                free(var.vwgt);
                free(var.xadj);
                free(var.adjcwgt);
                free(var.adjncy);
                free(var.part);
                free(var.best_part);
            }
        }

        void initial_partition(const deep_graph_t &g,
                               deep_p_manager_t &p_manager,
                               partition_t k,
                               partition_t k_spacing,
                               weight_t lmax,
                               partition_t hierarchy_level,
                               size_t kappa = 10) {

#pragma omp parallel for num_threads(m_threads) schedule(static)
            for (u64 thread_id = 0; thread_id < m_threads; thread_id++) {
                KaffpaVars &var = m_kaffpa_vars[thread_id];

                // step 1: extract the number of vertices and edges
                vertex_t subgraph_n_vertices = g.get_n();
                vertex_t subgraph_n_edges = g.get_m();
                weight_t subgraph_weight = g.weight();

                allocate(subgraph_n_vertices, subgraph_n_edges, thread_id);
                var.xadj[0] = 0;
                // step 2: build graph in kaffpa variables
                forall_gu(g, u)
                    {
                        var.vwgt[u] = (int) g.weight(u);
                        var.xadj[u + 1] = var.xadj[u];
                        forall_guivw(g, u, i, v, w) {
                                var.adjncy[var.xadj[u + 1]] = (int) v;
                                var.adjcwgt[var.xadj[u + 1]] = (int) w;
                                var.xadj[u + 1] += 1;
                            }
                        endfor
                    }
                endfor

                // step 3: calculate allowed imbalance
                f64 imbalance = ((f64) (lmax * k) / (f64) subgraph_weight) - 1.0;
                if (imbalance <= 0.0) {
                    // Failsafe if subgraph weight is too large // TODO: would be cool if not needed
                    imbalance = 0.0001;
                }

                // step 4: partition
                int n = (int) subgraph_n_vertices;
                int n_parts = (int) k;

                var.best_edge_cut = std::numeric_limits<int>::max();

                for (size_t i = 0; i < std::max((size_t) 1, kappa); ++i) {
                    kaffpa(&n, var.vwgt, var.xadj, var.adjcwgt, var.adjncy, &n_parts, &imbalance, true, 0, FAST, &var.edge_cut, var.part);
                    if (var.edge_cut < var.best_edge_cut) {
                        std::swap(var.edge_cut, var.best_edge_cut);
                        std::swap(var.part, var.best_part);
                    }
                }
            }

            // choose the best across all partitioning
            for (u64 thread_id = 1; thread_id < m_threads; thread_id++) {
                if (m_kaffpa_vars[thread_id].best_edge_cut < m_kaffpa_vars[0].best_edge_cut) {
                    std::swap(m_kaffpa_vars[0].best_edge_cut, m_kaffpa_vars[thread_id].best_edge_cut);
                    std::swap(m_kaffpa_vars[0].best_part, m_kaffpa_vars[thread_id].best_part);
                }
            }

            // step 5: extract partition into p_manager
            forall_gu(g, u)
                {
                    if (m_kaffpa_vars[0].best_part[u] == 0) { continue; }

                    partition_t old_id = 0;
                    partition_t new_id = m_kaffpa_vars[0].best_part[u] * k_spacing;

                    p_manager.move(u, g.weight(u), old_id, new_id);
                }
            endfor

            for (partition_t i = 0; i < k; ++i) {
                partition_t move_id = k_spacing * i;
                p_manager.set_lmax(move_id, lmax);
                p_manager.set_hierarchy_level(move_id, hierarchy_level - 1);
            }
        }

        void intermediate_partition(const deep_graph_t &g,
                                    deep_p_manager_t &p_manager,
                                    u64 thread_id,
                                    partition_t block_id,
                                    partition_t k_spacing,
                                    partition_t block_k,
                                    weight_t block_lmax,
                                    partition_t block_level,
                                    size_t kappa = 1) {
            if (blocks[block_id].empty()) {
                for (partition_t i = 0; i < block_k; ++i) {
                    partition_t move_id = block_id + k_spacing * i;
                    p_manager.set_lmax(move_id, block_lmax);
                    p_manager.set_hierarchy_level(move_id, block_level - 1);
                }
            }

            KaffpaVars &var = m_kaffpa_vars[thread_id];

            var.tt.clear();

            vertex_t subgraph_n_vertices = 0;
            vertex_t subgraph_n_edges = 0;
            weight_t subgraph_weight = 0;
            // step 1: extract the number of vertices and edges, build tt
            for (vertex_t u: blocks[block_id]) {
                var.tt.add(u, subgraph_n_vertices);
                subgraph_n_vertices += 1;
                subgraph_weight += g.weight(u);
                forall_guiv(g, u, i, v) {
                        if (p_manager[v] == block_id) { subgraph_n_edges += 1; }
                    }
                endfor
            }

            allocate(subgraph_n_vertices, subgraph_n_edges, thread_id);
            var.xadj[0] = 0;
            // step 2: build graph in kaffpa variables
            for (vertex_t u: blocks[block_id]) {
                vertex_t new_u = var.tt.get_n(u);
                var.vwgt[new_u] = (int) g.weight(u);
                var.xadj[new_u + 1] = var.xadj[new_u];
                forall_guivw(g, u, i, v, w) {
                        if (p_manager[v] == block_id) {
                            vertex_t new_v = var.tt.get_n(v);
                            var.adjncy[var.xadj[new_u + 1]] = (int) new_v;
                            var.adjcwgt[var.xadj[new_u + 1]] = (int) w;
                            var.xadj[new_u + 1] += 1;
                        }
                    }
                endfor
            }

            // step 3: calculate allowed imbalance
            f64 imbalance = ((f64) (block_lmax * block_k) / (f64) subgraph_weight) - 1.0;
            if (imbalance <= 0.0) {
                // Failsafe if subgraph weight is too large // TODO: would be cool if not needed
                imbalance = 0.0001;
            }

            // step 4: partition
            int n = (int) subgraph_n_vertices;
            int n_parts = (int) block_k;

            var.best_edge_cut = std::numeric_limits<int>::max();
            for (size_t i = 0; i < kappa; ++i) {
                kaffpa(&n, var.vwgt, var.xadj, var.adjcwgt, var.adjncy, &n_parts, &imbalance, true, 0, FAST, &var.edge_cut, var.part);
                if (var.edge_cut < var.best_edge_cut) {
                    std::swap(var.edge_cut, var.best_edge_cut);
                    std::swap(var.part, var.best_part);
                }
            }

            // step 5: extract partition into p_manager
            for (vertex_t kaffpa_vertex = 0; kaffpa_vertex < subgraph_n_vertices; ++kaffpa_vertex) {
                if (var.best_part[kaffpa_vertex] == 0) { continue; }
                partition_t move_id = block_id + k_spacing * var.best_part[kaffpa_vertex];

                partition_t vertex_id = block_id;
                vertex_t vertex = var.tt.get_o(kaffpa_vertex);
                weight_t vertex_weight = g.weight(vertex);

                p_manager.move(vertex, vertex_weight, vertex_id, move_id);
            }

            for (partition_t i = 0; i < block_k; ++i) {
                partition_t move_id = block_id + k_spacing * i;
                p_manager.set_lmax(move_id, block_lmax);
                p_manager.set_hierarchy_level(move_id, block_level - 1);
            }
        }

        void allocate(size_t new_n_vertices,
                      size_t new_n_edges,
                      u64 thread_id) {
            KaffpaVars &var = m_kaffpa_vars[thread_id];

            if (new_n_vertices > var.n_vertices) {
                var.n_vertices = new_n_vertices;

                free(var.vwgt);
                free(var.xadj);
                free(var.part);
                free(var.best_part);

                var.vwgt = (int *) malloc(var.n_vertices * sizeof(int));
                var.xadj = (int *) malloc((var.n_vertices + 1) * sizeof(int));
                var.part = (int *) malloc(var.n_vertices * sizeof(int));
                var.best_part = (int *) malloc(var.n_vertices * sizeof(int));
            }

            if (new_n_edges > var.n_edges) {
                var.n_edges = new_n_edges;

                free(var.adjcwgt);
                free(var.adjncy);

                var.adjcwgt = (int *) malloc(var.n_edges * sizeof(int));
                var.adjncy = (int *) malloc(var.n_edges * sizeof(int));
            }
        }

        void determine_all_blocks(const deep_graph_t &g,
                                  const deep_p_manager_t &p_manager) {
            for (partition_t id = 0; id < m_k; ++id) {
                blocks[id].clear();
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
