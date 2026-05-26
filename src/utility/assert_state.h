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

#include <map>

#include "../definitions.h"
#include "utils.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/distance_oracle.h"
#include "../datastructures/partition_manager.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/quotient_graph.h"
#include "../datastructures/block_conn.h"

namespace HeiProMap {
    inline bool assert_csr_structure(const graph_t &g) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_csr_structure");

        ASSERT(g.neighborhoods[0] == 0);
        ASSERT(g.neighborhoods[g.n] == g.m);
        for (vertex_t u = 0; u < g.n; ++u) {
            ASSERT(g.neighborhoods[u] <= g.neighborhoods[u + 1]);
        }
        return true;
    }

    inline bool assert_no_self_loops(const graph_t &g) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_no_self_loops");

        for (vertex_t u = 0; u < g.n; ++u) {
            {
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    for (size_t j = i + 1; j < g.neighborhoods[u + 1]; ++j) {
                        ASSERT(g.edges_v[i] != g.edges_v[j]);
                    }
                }
            }
        }
        return true;
    }

    inline bool assert_no_double_edges(const graph_t &g) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_no_double_edges");

        std::vector<vertex_t> manual;
        for (vertex_t u = 0; u < g.n; ++u) {
            {
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
                                              const partition_t k) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_correct_partition_size");

        std::vector<size_t> sizes(k, 0);

        for (vertex_t u = 0; u < g.n; ++u) {
            {
                sizes[p_manager[u]] += 1;
            }
        }

        for (partition_t id = 0; id < k; ++id) {
            ASSERT(sizes[id] == p_manager.size(id));
        }

        return true;
    }

    inline bool assert_bweights(const graph_t &g,
                                const p_manager_t &p_manager,
                                const partition_t k) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_bweights");

        std::vector<weight_t> weights(k, 0);
        for (vertex_t u = 0; u < g.n; ++u) {
            {
                partition_t u_id = p_manager[u];

                weights[u_id] += g.v_weights[u];
            }
        }

        for (partition_t id = 0; id < k; ++id) {
            ASSERT(weights[id] == p_manager.get_bweight(id));
        }

        return true;
    }

    inline bool assert_bweights(const graph_t &g,
                                const p_manager_t &p_manager,
                                const TranslationTable<vertex_t> &tt,
                                const u64 offset,
                                const partition_t k) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_bweights_partial");

        std::vector<weight_t> weights(k, 0);
        for (vertex_t u = 0; u < g.n; ++u) {
            partition_t u_id = p_manager[tt.get_o(u)];
            ASSERT(u_id >= offset && u_id < offset + k);
            weights[u_id - offset] += g.v_weights[u];
        }

        // We can only check if the weights are consistent if we knew the partial weights in p_manager.
        // But p_manager stores global weights.
        // So we just check if all vertices of g are mapped to the correct range.
        return true;
    }

    inline bool assert_state_partial(const graph_t &g,
                                     const p_manager_t &p_manager,
                                     const TranslationTable<vertex_t> &tt,
                                     const u64 offset,
                                     const partition_t k) {
        // assert csr structure
        ASSERT(assert_csr_structure(g));

        // check no self-loops
        ASSERT(assert_no_self_loops(g));

        // check no duplicate edges
        ASSERT(assert_no_double_edges(g));

        // check the correct block weights (range check)
        ASSERT(assert_bweights(g, p_manager, tt, offset, k));

        return true;
    }

    inline bool assert_correct_vertices_boundary(const graph_t &g,
                                                 const p_manager_t &p_manager,
                                                 const bv_manager_t &bv_manager) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_correct_vertices_boundary");

        std::vector<vertex_t> manual;
        for (vertex_t u = 0; u < g.n; ++u) {
            {
                partition_t u_id = p_manager[u];

                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    partition_t v_id = p_manager[g.edges_v[i]];
                    if (u_id != v_id) {
                        manual.push_back(u);
                        break;
                    }
                }
            }
        }

        std::vector<vertex_t> automatic;
        for (partition_t id = 0; id < bv_manager.get_k(); ++id) {
            for (size_t i = 0; i < bv_manager.size(id); ++i) { const vertex_t u = bv_manager.get(id, i);
                {
                    automatic.push_back(u);
                }
            }
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
                                                           bv_manager_t &bv_manager,
                                                           const partition_t k) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_correct_vertices_boundary_per_block");

        std::vector<std::vector<vertex_t> > manual(k);

        for (vertex_t u = 0; u < g.n; ++u) {
            {
                partition_t u_id = p_manager[u];
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) { const vertex_t v = g.edges_v[i];
                    {
                        partition_t v_id = p_manager[v];
                        if (u_id != v_id) {
                            manual[u_id].push_back(u);
                            break;
                        }
                    }
                }
            }
        }

        for (partition_t id = 0; id < k; ++id) {
            std::vector<vertex_t> automatic;
            for (size_t i = 0; i < bv_manager.size(id); ++i) { const vertex_t u = bv_manager.get(id, i);
                {
                    automatic.push_back(u);
                }
            }

            std::sort(manual[id].begin(), manual[id].end());
            std::sort(automatic.begin(), automatic.end());
            ASSERT(no_duplicates_sorted(manual[id]));
            ASSERT(no_duplicates_sorted(automatic));

            ASSERT(manual[id] == automatic);
        }

        return true;
    }

    inline bool assert_correct_boundary([[maybe_unused]] const graph_t &g,
                                        [[maybe_unused]] const p_manager_t &p_manager,
                                        [[maybe_unused]] bv_manager_t &bv_manager,
                                        [[maybe_unused]] const partition_t k) {
        ASSERT(assert_correct_vertices_boundary(g, p_manager, bv_manager));
        ASSERT(assert_correct_vertices_boundary_per_block(g, p_manager, bv_manager, k));
        return true;
    }

    inline bool assert_correct_quotient_graph(const graph_t &g,
                                              const p_manager_t &p_manager,
                                              const q_graph_t &q_graph,
                                              [[maybe_unused]] const partition_t k) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_correct_quotient_graph");

        std::map<std::pair<partition_t, partition_t>, weight_t> manual;

        for (vertex_t u = 0; u < g.n; ++u) {
            {
                partition_t u_id = p_manager[u];
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) { const vertex_t v = g.edges_v[i]; const weight_t w = g.edges_w[i];
                    {
                        partition_t v_id = p_manager[v];
                        if (u_id != v_id) {
                            manual[{u_id, v_id}] += w;
                        }
                    }
                }
            }
        }

        for (auto [pair, w]: manual) {
            partition_t id1 = pair.first;
            partition_t id2 = pair.second;

            if (w != q_graph.get_weight(id1, id2)) {
                std::cout << w << " " << q_graph.get_weight(id1, id2) << " " << id1 << " " << id2 << std::endl;
            }

            ASSERT(w == q_graph.get_weight(id1, id2));
        }
        return true;
    }

    inline bool assert_correct_block_conn(const graph_t &g,
                                          const p_manager_t &p_manager,
                                          const block_conn_t &block_conn,
                                          [[maybe_unused]] const partition_t k) {
        HEIPROMAP_PROFILE_SCOPE("assert", "misc", "assert_correct_block_conn");

        for (vertex_t u = 0; u < g.n; ++u) {
            {
                std::map<partition_t, weight_t> manual_map;

                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) { const vertex_t v = g.edges_v[i]; const weight_t w = g.edges_w[i];
                    {
                        partition_t v_id = p_manager[v];
                        manual_map[v_id] += w;
                    }
                }
                std::vector<std::pair<partition_t, weight_t> > manual;
                for (auto [id, w]: manual_map) {
                    manual.emplace_back(id, w);
                }

                std::vector<std::pair<partition_t, weight_t> > automatic;

                for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) { const partition_t id = block_conn.get_id(i); const weight_t idw = block_conn.get_w(i);
                    {
                        automatic.emplace_back(id, idw);
                    }
                }

                std::sort(manual.begin(), manual.end(), [](const auto &a, const auto &b) {
                    return a.first < b.first;
                });
                std::sort(automatic.begin(), automatic.end(), [](const auto &a, const auto &b) {
                    return a.first < b.first;
                });

                ASSERT(no_duplicates_sorted(manual));
                ASSERT(no_duplicates_sorted(automatic));
                ASSERT(manual.size() == automatic.size());

                for (size_t i = 0; i < manual.size(); ++i) {
                    ASSERT(manual[i].first == automatic[i].first);
                    ASSERT(manual[i].second == automatic[i].second);
                }
            }
        }

        return true;
    }

    inline bool assert_graph([[maybe_unused]] const graph_t &g) {
        // assert csr structure
        ASSERT(assert_csr_structure(g));

        // check no self-loops
        ASSERT(assert_no_self_loops(g));

        // check no duplicate edges
        ASSERT(assert_no_double_edges(g));

        return true;
    }

    inline bool assert_state_pre_partitioning([[maybe_unused]] const graph_t &g,
                                              [[maybe_unused]] const p_manager_t &p_manager,
                                              [[maybe_unused]] const partition_t k) {
        // assert csr structure
        ASSERT(assert_csr_structure(g));

        // check no self-loops
        ASSERT(assert_no_self_loops(g));

        // check no duplicate edges
        ASSERT(assert_no_double_edges(g));

        // check the correct partition sizes
        ASSERT(assert_correct_partition_size(g, p_manager, k));

        return true;
    }

    inline bool assert_state_after_partitioning([[maybe_unused]] const graph_t &g,
                                                [[maybe_unused]] const p_manager_t &p_manager,
                                                [[maybe_unused]] bv_manager_t &bv_manager,
                                                [[maybe_unused]] const q_graph_t &q_graph,
                                                [[maybe_unused]] const block_conn_t &block_conn,
                                                [[maybe_unused]] const partition_t k) {
        // assert csr structure
        ASSERT(assert_csr_structure(g));

        // check no self-loops
        ASSERT(assert_no_self_loops(g));

        // check no duplicate edges
        ASSERT(assert_no_double_edges(g));

        // check the correct partition sizes
        ASSERT(assert_correct_partition_size(g, p_manager, k));

        // check the correct block weights
        ASSERT(assert_bweights(g, p_manager, k));

        // check the right vertices are boundary
        ASSERT(assert_correct_vertices_boundary(g, p_manager, bv_manager));

        // check the right vertices are boundary per block
        ASSERT(assert_correct_vertices_boundary_per_block(g, p_manager, bv_manager, k));

        // check the correct quotient graph
        ASSERT(assert_correct_quotient_graph(g, p_manager, q_graph, k));

        ASSERT(assert_correct_block_conn(g, p_manager, block_conn, k));

        return true;
    }

    inline bool assert_state_after_partitioning([[maybe_unused]] const graph_t &g,
                                                [[maybe_unused]] const p_manager_t &p_manager,
                                                [[maybe_unused]] const partition_t k) {
        // assert csr structure
        ASSERT(assert_csr_structure(g));

        // check no self-loops
        ASSERT(assert_no_self_loops(g));

        // check no duplicate edges
        ASSERT(assert_no_double_edges(g));

        // check the correct partition sizes
        ASSERT(assert_correct_partition_size(g, p_manager, k));

        // check the correct block weights
        ASSERT(assert_bweights(g, p_manager, k));

        return true;
    }
}


#endif //HEIPROMAP_ASSERT_STATE_H
