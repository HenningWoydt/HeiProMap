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

#ifndef HEIPROMAP_ASSERT_STATE_H
#define HEIPROMAP_ASSERT_STATE_H

#include <algorithm>
#include <map>
#include <utility>
#include <vector>
#include <omp.h>

#include "../definitions.h"
#include "utils.h"
#include "../datastructures/csr_graph.h"
#include "../distance_oracles/distance_oracle.h"
#include "../datastructures/partition_manager.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/quotient_graph.h"
#include "../datastructures/large_quotient_graph.h"
#include "../datastructures/block_conn.h"

namespace HeiProMap {
    inline bool assert_csr_structure(const graph_t &g, const u64 num_threads = 1) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_csr_structure");

        ASSERT(g.neighborhoods[0] == 0);
        ASSERT(g.neighborhoods[g.n] == g.m);
        #pragma omp parallel for num_threads(num_threads) schedule(static)
        for (vertex_t u = 0; u < g.n; ++u) {
            ASSERT(g.neighborhoods[u] <= g.neighborhoods[u + 1]);
        }
        return true;
    }

    inline bool assert_no_self_loops(const graph_t &g, const u64 num_threads = 1) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_no_self_loops");

        #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 1024)
        for (vertex_t u = 0; u < g.n; ++u) {
            for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                for (size_t j = i + 1; j < g.neighborhoods[u + 1]; ++j) {
                    ASSERT(g.edges_v[i] != g.edges_v[j]);
                }
            }
        }
        return true;
    }

    inline bool assert_no_double_edges(const graph_t &g, const u64 num_threads = 1) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_no_double_edges");

        #pragma omp parallel num_threads(num_threads)
        {
            std::vector<vertex_t> manual;
            #pragma omp for schedule(dynamic, 1024)
            for (vertex_t u = 0; u < g.n; ++u) {
                manual.clear();
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    manual.push_back(g.edges_v[i]);
                }
                std::sort(manual.begin(), manual.end());
                ASSERT(no_duplicates_sorted(manual));
            }
        }
        return true;
    }

    inline bool assert_correct_partition_size(const graph_t &g,
                                              const p_manager_t &p_manager,
                                              const partition_t k,
                                              const u64 num_threads = 1) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_correct_partition_size");

        std::vector<std::vector<size_t>> thread_sizes(num_threads, std::vector<size_t>(k, 0));

        #pragma omp parallel num_threads(num_threads)
        {
            const u64 tid = omp_get_thread_num();
            auto &local_sizes = thread_sizes[tid];
            #pragma omp for schedule(static)
            for (vertex_t u = 0; u < g.n; ++u) {
                local_sizes[p_manager[u]] += 1;
            }
        }

        std::vector<size_t> sizes(k, 0);
        #pragma omp parallel for num_threads(num_threads) schedule(static)
        for (partition_t id = 0; id < k; ++id) {
            for (u64 t = 0; t < num_threads; ++t) {
                sizes[id] += thread_sizes[t][id];
            }
            ASSERT(sizes[id] == p_manager.size(id));
        }

        return true;
    }

    inline bool assert_bweights(const graph_t &g,
                                const p_manager_t &p_manager,
                                const partition_t k,
                                const u64 num_threads = 1) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_bweights");

        std::vector<std::vector<weight_t>> thread_weights(num_threads, std::vector<weight_t>(k, 0));

        #pragma omp parallel num_threads(num_threads)
        {
            const u64 tid = omp_get_thread_num();
            auto &local_weights = thread_weights[tid];
            #pragma omp for schedule(static)
            for (vertex_t u = 0; u < g.n; ++u) {
                const partition_t u_id = p_manager[u];
                local_weights[u_id] += g.uniform_v_weights ? 1 : g.v_weights[u];
            }
        }

        #pragma omp parallel for num_threads(num_threads) schedule(static)
        for (partition_t id = 0; id < k; ++id) {
            weight_t w = 0;
            for (u64 t = 0; t < num_threads; ++t) {
                w += thread_weights[t][id];
            }
            ASSERT(w == p_manager.get_bweight(id));
        }

        return true;
    }

    inline bool assert_bweights(const graph_t &g,
                                const p_manager_t &p_manager,
                                const TranslationTable<vertex_t> &tt,
                                const u64 offset,
                                const partition_t k,
                                const u64 num_threads = 1) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_bweights_partial");

        #pragma omp parallel for num_threads(num_threads) schedule(static)
        for (vertex_t u = 0; u < g.n; ++u) {
            partition_t u_id = p_manager[tt.get_o(u)];
            ASSERT(u_id >= offset && u_id < offset + k);
        }
        return true;
    }

    inline bool assert_state_partial(const graph_t &g,
                                     const p_manager_t &p_manager,
                                     const TranslationTable<vertex_t> &tt,
                                     const u64 offset,
                                     const partition_t k,
                                     const u64 num_threads = 1) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_state_partial");

        ASSERT(assert_csr_structure(g, num_threads));
        ASSERT(assert_no_self_loops(g, num_threads));
        ASSERT(assert_no_double_edges(g, num_threads));
        ASSERT(assert_bweights(g, p_manager, tt, offset, k, num_threads));

        return true;
    }

    inline bool assert_correct_vertices_boundary(const graph_t &g,
                                                 const p_manager_t &p_manager,
                                                 const bv_manager_t &bv_manager,
                                                 const u64 num_threads = 1) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_correct_vertices_boundary");

        std::vector<std::vector<vertex_t>> thread_manual(num_threads);
        #pragma omp parallel num_threads(num_threads)
        {
            const u64 tid = omp_get_thread_num();
            auto &local_manual = thread_manual[tid];
            #pragma omp for schedule(dynamic, 1024)
            for (vertex_t u = 0; u < g.n; ++u) {
                const partition_t u_id = p_manager[u];
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    const partition_t v_id = p_manager[g.edges_v[i]];
                    if (u_id != v_id) {
                        local_manual.push_back(u);
                        break;
                    }
                }
            }
        }

        std::vector<vertex_t> manual;
        for (u64 t = 0; t < num_threads; ++t) {
            manual.insert(manual.end(), thread_manual[t].begin(), thread_manual[t].end());
        }

        std::vector<std::vector<vertex_t>> thread_auto(num_threads);
        #pragma omp parallel num_threads(num_threads)
        {
            const u64 tid = omp_get_thread_num();
            auto &local_auto = thread_auto[tid];
            #pragma omp for schedule(dynamic, 64)
            for (partition_t id = 0; id < bv_manager.get_k(); ++id) {
                for (size_t i = 0; i < bv_manager.size(id); ++i) {
                    local_auto.push_back(bv_manager.get(id, i));
                }
            }
        }

        std::vector<vertex_t> automatic;
        for (u64 t = 0; t < num_threads; ++t) {
            automatic.insert(automatic.end(), thread_auto[t].begin(), thread_auto[t].end());
        }

        std::sort(manual.begin(), manual.end());
        std::sort(automatic.begin(), automatic.end());
        ASSERT(no_duplicates_sorted(manual));
        ASSERT(no_duplicates_sorted(automatic));
        ASSERT(manual == automatic);
        return manual == automatic;
    }

    inline bool assert_correct_vertices_boundary_per_block(const graph_t &g,
                                                           const p_manager_t &p_manager,
                                                           const bv_manager_t &bv_manager,
                                                           const partition_t k,
                                                           const u64 num_threads = 1) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_correct_vertices_boundary_per_block");

        std::vector<std::vector<vertex_t>> manual(k);
        std::vector<std::vector<std::vector<vertex_t>>> thread_manual(num_threads, std::vector<std::vector<vertex_t>>(k));

        #pragma omp parallel num_threads(num_threads)
        {
            const u64 tid = omp_get_thread_num();
            auto &local_manual = thread_manual[tid];
            #pragma omp for schedule(dynamic, 1024)
            for (vertex_t u = 0; u < g.n; ++u) {
                const partition_t u_id = p_manager[u];
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    const partition_t v_id = p_manager[g.edges_v[i]];
                    if (u_id != v_id) {
                        local_manual[u_id].push_back(u);
                        break;
                    }
                }
            }
        }

        #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 64)
        for (partition_t id = 0; id < k; ++id) {
            std::vector<vertex_t> block_manual;
            for (u64 t = 0; t < num_threads; ++t) {
                block_manual.insert(block_manual.end(), thread_manual[t][id].begin(), thread_manual[t][id].end());
            }

            std::vector<vertex_t> automatic;
            automatic.reserve(bv_manager.size(id));
            for (size_t i = 0; i < bv_manager.size(id); ++i) {
                automatic.push_back(bv_manager.get(id, i));
            }

            std::sort(block_manual.begin(), block_manual.end());
            std::sort(automatic.begin(), automatic.end());
            ASSERT(no_duplicates_sorted(block_manual));
            ASSERT(no_duplicates_sorted(automatic));
            ASSERT(block_manual == automatic);
        }

        return true;
    }

    inline bool assert_correct_boundary(const graph_t &g,
                                        const p_manager_t &p_manager,
                                        const bv_manager_t &bv_manager,
                                        const partition_t k,
                                        const u64 num_threads = 1) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_correct_boundary");

        ASSERT(assert_correct_vertices_boundary(g, p_manager, bv_manager, num_threads));
        ASSERT(assert_correct_vertices_boundary_per_block(g, p_manager, bv_manager, k, num_threads));
        return true;
    }

    template<typename QGraphT>
    inline bool assert_correct_quotient_graph(const graph_t &g,
                                              const p_manager_t &p_manager,
                                              const bv_manager_t &bv_manager,
                                              const QGraphT &q_graph,
                                              const partition_t k,
                                              const u64 num_threads = 1) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_correct_quotient_graph");

        QGraphT manual_q;
        manual_q.compute_from_scratch(g, p_manager, bv_manager, num_threads);

        #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 64)
        for (partition_t id = 0; id < k; ++id) {
            q_graph.for_each_neighbor(id, [&](const partition_t v, const weight_t w) {
                ASSERT(w == manual_q.get_weight(id, v));
            });
            manual_q.for_each_neighbor(id, [&](const partition_t v, const weight_t w) {
                ASSERT(w == q_graph.get_weight(id, v));
            });
        }
        return true;
    }

    template<typename QGraphT>
    inline bool assert_correct_quotient_graph(const graph_t &g,
                                              const p_manager_t &p_manager,
                                              const QGraphT &q_graph,
                                              const partition_t k,
                                              const u64 num_threads = 1) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_correct_quotient_graph");

        std::vector<std::map<std::pair<partition_t, partition_t>, weight_t>> thread_maps(num_threads);

        #pragma omp parallel num_threads(num_threads)
        {
            const u64 tid = omp_get_thread_num();
            auto &local_map = thread_maps[tid];
            #pragma omp for schedule(dynamic, 1024)
            for (vertex_t u = 0; u < g.n; ++u) {
                const partition_t u_id = p_manager[u];
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    const vertex_t v = g.edges_v[i];
                    if (u < v) {
                        const partition_t v_id = p_manager[v];
                        if (u_id != v_id) {
                            const partition_t lo = std::min(u_id, v_id);
                            const partition_t hi = std::max(u_id, v_id);
                            local_map[{lo, hi}] += g.edges_w[i];
                        }
                    }
                }
            }
        }

        std::map<std::pair<partition_t, partition_t>, weight_t> manual;
        for (u64 t = 0; t < num_threads; ++t) {
            for (const auto &[pair, w] : thread_maps[t]) {
                manual[pair] += w;
            }
        }

        for (const auto &[pair, w] : manual) {
            ASSERT(w == q_graph.get_weight(pair.first, pair.second));
        }
        return true;
    }

    inline bool assert_correct_block_conn(const graph_t &g,
                                          const p_manager_t &p_manager,
                                          const block_conn_t &block_conn,
                                          [[maybe_unused]] const partition_t k,
                                          const u64 num_threads = 1) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_correct_block_conn");

        #pragma omp parallel num_threads(num_threads)
        {
            std::vector<std::pair<partition_t, weight_t>> manual;
            std::vector<std::pair<partition_t, weight_t>> automatic;

            #pragma omp for schedule(dynamic, 1024)
            for (vertex_t u = 0; u < g.n; ++u) {
                manual.clear();
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    const vertex_t v = g.edges_v[i];
                    const weight_t w = g.edges_w[i];
                    const partition_t v_id = p_manager[v];
                    manual.emplace_back(v_id, w);
                }

                std::sort(manual.begin(), manual.end(), [](const auto &a, const auto &b) {
                    return a.first < b.first;
                });

                size_t write_idx = 0;
                for (size_t i = 0; i < manual.size(); ++i) {
                    if (write_idx > 0 && manual[write_idx - 1].first == manual[i].first) {
                        manual[write_idx - 1].second += manual[i].second;
                    } else {
                        manual[write_idx++] = manual[i];
                    }
                }
                manual.resize(write_idx);

                automatic.clear();
                for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) {
                    automatic.emplace_back(block_conn.get_id(i), block_conn.get_w(i));
                }

                std::sort(automatic.begin(), automatic.end(), [](const auto &a, const auto &b) {
                    return a.first < b.first;
                });

                ASSERT(manual.size() == automatic.size());
                for (size_t i = 0; i < manual.size(); ++i) {
                    ASSERT(manual[i].first == automatic[i].first);
                    ASSERT(manual[i].second == automatic[i].second);
                }
            }
        }

        return true;
    }

    inline bool assert_graph(const graph_t &g, const u64 num_threads = 1) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_graph");

        ASSERT(assert_csr_structure(g, num_threads));
        ASSERT(assert_no_self_loops(g, num_threads));
        ASSERT(assert_no_double_edges(g, num_threads));

        return true;
    }

    inline bool assert_state_pre_partitioning(const graph_t &g,
                                              const p_manager_t &p_manager,
                                              const partition_t k,
                                              const u64 num_threads = 1) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_state_pre_partitioning");

        ASSERT(assert_csr_structure(g, num_threads));
        ASSERT(assert_no_self_loops(g, num_threads));
        ASSERT(assert_no_double_edges(g, num_threads));
        ASSERT(assert_correct_partition_size(g, p_manager, k, num_threads));

        return true;
    }

    template<typename QGraphT>
    inline bool assert_state_after_partitioning(const graph_t &g,
                                                const p_manager_t &p_manager,
                                                const bv_manager_t &bv_manager,
                                                const QGraphT &q_graph,
                                                const block_conn_t &block_conn,
                                                const partition_t k,
                                                const u64 num_threads = 1) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_state_after_partitioning");

        ASSERT(assert_csr_structure(g, num_threads));
        ASSERT(assert_no_self_loops(g, num_threads));
        ASSERT(assert_no_double_edges(g, num_threads));
        ASSERT(assert_correct_partition_size(g, p_manager, k, num_threads));
        ASSERT(assert_bweights(g, p_manager, k, num_threads));
        ASSERT(assert_correct_vertices_boundary(g, p_manager, bv_manager, num_threads));
        ASSERT(assert_correct_vertices_boundary_per_block(g, p_manager, bv_manager, k, num_threads));
        ASSERT(assert_correct_quotient_graph(g, p_manager, bv_manager, q_graph, k, num_threads));
        ASSERT(assert_correct_block_conn(g, p_manager, block_conn, k, num_threads));

        return true;
    }

    inline bool assert_state_after_partitioning(const graph_t &g,
                                                const p_manager_t &p_manager,
                                                const partition_t k,
                                                const u64 num_threads = 1) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_state_after_partitioning");

        ASSERT(assert_csr_structure(g, num_threads));
        ASSERT(assert_no_self_loops(g, num_threads));
        ASSERT(assert_no_double_edges(g, num_threads));
        ASSERT(assert_correct_partition_size(g, p_manager, k, num_threads));
        ASSERT(assert_bweights(g, p_manager, k, num_threads));

        return true;
    }
}

#endif //HEIPROMAP_ASSERT_STATE_H
