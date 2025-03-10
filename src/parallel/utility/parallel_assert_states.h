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

#ifndef HEIPROMAP_PARALLEL_ASSERT_STATES_H
#define HEIPROMAP_PARALLEL_ASSERT_STATES_H

#include "../parallel_definitions_1.h"
#include "../parallel_definitions_2.h"

namespace HeiProMap {

    bool assert_n_active_vertices(p_graph_t &g,
                                  p_av_manager_t &av_manager,
                                  u64 threads) {
        vertex_t      count = 0;
        for (vertex_t u     = 0; u < g.get_n(); ++u) {
            count += av_manager.is_active(u);
        }

        ASSERT(count == av_manager.get_n_active());
        return count == av_manager.get_n_active();
    }

    bool assert_correct_vertices_active(p_graph_t &g,
                                        p_av_manager_t &av_manager,
                                        u64 threads) {
        std::vector<vertex_t> manual;
        for (vertex_t         u = 0; u < g.get_n(); ++u) {
            if (av_manager.is_active(u)) { manual.push_back(u); }
        }

        std::vector<vertex_t> automatic;
        for (vertex_t         v: av_manager) {
            automatic.push_back(v);
        }

        std::sort(manual.begin(), manual.end());
        std::sort(automatic.begin(), automatic.end());
        ASSERT(no_duplicates_sorted(manual));
        ASSERT(no_duplicates_sorted(automatic));

        ASSERT(manual == automatic);
        return manual == automatic;
    }

    bool assert_active_edges(p_graph_t &g,
                             p_av_manager_t &av_manager,
                             u64 threads) {
        for (vertex_t u: av_manager) {
            for (size_t i = 0; i < g.size(u); ++i) {
                const vertex_t v = g.neighbor(u, i);
                ASSERT(av_manager.is_active(v));
            }
        }

        return true;
    }

    bool assert_no_self_loops(p_graph_t &g,
                              p_av_manager_t &av_manager,
                              u64 threads) {
        for (vertex_t u: av_manager) {
            for (size_t i = 0; i < g.size(u); ++i) {
                for (size_t j = i + 1; j < g.size(u); ++j) {
                    ASSERT(g.neighbor(u, i) != g.neighbor(u, j));
                }
            }
        }
        return true;
    }

    bool assert_no_double_edges(p_graph_t &g,
                                p_av_manager_t &av_manager,
                                u64 threads) {
        std::vector<vertex_t> manual;
        for (vertex_t         u: av_manager) {
            manual.clear();
            for (size_t i = 0; i < g.size(u); ++i) {
                manual.push_back(g.neighbor(u, i));
            }
            std::sort(manual.begin(), manual.end());
            ASSERT(no_duplicates_sorted(manual));
        }
        return true;
    }

    bool assert_bweights(p_graph_t &g,
                         p_av_manager_t &av_manager,
                         p_p_manager_t &p_manager,
                         partition_t k,
                         u64 threads) {
        std::vector<weight_t> weights(k, 0);
        for (vertex_t         u: av_manager) {
            partition_t u_id = p_manager[u];

            weights[u_id] += g.get_weight(u);
        }

        for (partition_t id = 0; id < k; ++id) {
            ASSERT(weights[id] == p_manager.get_bweight(id));
        }

        return true;
    }

    bool assert_n_boundary_vertices(p_graph_t &g,
                                    p_av_manager_t &av_manager,
                                    p_p_manager_t &p_manager,
                                    p_bv_manager_t &bv_manager,
                                    u64 threads) {
        vertex_t      count = 0;
        for (vertex_t u: av_manager) {
            partition_t u_id = p_manager[u];

            for (size_t i = 0; i < g.size(u); ++i) {
                partition_t v_id = p_manager[g.neighbor(u, i)];

                if (u_id != v_id) {
                    count += 1;
                    break;
                }
            }
        }

        std::cout << count << " " << bv_manager.get_n_boundary() << " " << av_manager.get_n_active() << std::endl;

        ASSERT(count == bv_manager.get_n_boundary());
        return count == bv_manager.get_n_boundary();
    }

    bool assert_correct_vertices_boundary(p_graph_t &g,
                                          p_av_manager_t &av_manager,
                                          p_p_manager_t &p_manager,
                                          p_bv_manager_t &bv_manager,
                                          u64 threads) {
        std::vector<vertex_t> manual;
        for (vertex_t         u: av_manager) {
            partition_t u_id = p_manager[u];

            for (size_t i = 0; i < g.size(u); ++i) {
                partition_t v_id = p_manager[g.neighbor(u, i)];
                if (u_id != v_id) {
                    manual.push_back(u);
                    break;
                }
            }
        }

        std::vector<vertex_t> automatic;
        for (size_t           i = 0; i < bv_manager.get_n_boundary(); ++i) {
            automatic.push_back(bv_manager.get(i));
        }

        std::sort(manual.begin(), manual.end());
        std::sort(automatic.begin(), automatic.end());
        ASSERT(no_duplicates_sorted(manual));
        ASSERT(no_duplicates_sorted(automatic));

        ASSERT(manual == automatic);
        return manual == automatic;
    }

    bool assert_correct_vertices_boundary_per_block(p_graph_t &g,
                                                    p_av_manager_t &av_manager,
                                                    p_p_manager_t &p_manager,
                                                    p_bv_manager_t &bv_manager,
                                                    partition_t k,
                                                    u64 threads) {
        for (partition_t id = 0; id < k; ++id) {
            std::vector<vertex_t> manual;

            for (vertex_t u: av_manager) {
                partition_t u_id = p_manager[u];
                if (u_id == id) {
                    for (size_t i = 0; i < g.size(u); ++i) {
                        partition_t v_id = p_manager[g.neighbor(u, i)];
                        if (u_id != v_id) {
                            manual.push_back(u);
                            break;
                        }
                    }
                }
            }

            std::vector<vertex_t> automatic;
            for (size_t           i = 0; i < bv_manager.get_n_boundary(id); ++i) {
                automatic.push_back(bv_manager.get(id, i));
            }

            std::sort(manual.begin(), manual.end());
            std::sort(automatic.begin(), automatic.end());
            ASSERT(no_duplicates_sorted(manual));
            ASSERT(no_duplicates_sorted(automatic));

            ASSERT(manual == automatic);
        }
        return true;
    }

    bool assert_correct_quotient_graph(p_graph_t &g,
                                       p_av_manager_t &av_manager,
                                       p_p_manager_t &p_manager,
                                       p_q_graph_t &q_graph,
                                       partition_t k,
                                       u64 threads) {
        std::vector<weight_t> manual(k * k, 0);

        for (vertex_t u: av_manager) {
            partition_t u_id = p_manager[u];
            for (size_t i    = 0; i < g.size(u); ++i) {
                vertex_t v = g.neighbor(u, i);
                weight_t w = g.get_weight(u, i);

                partition_t v_id = p_manager[v];

                manual[u_id * k + v_id] += w;
                manual[v_id * k + u_id] += w;
            }
        }

        for (partition_t id_1 = 0; id_1 < k; ++id_1) {
            for (partition_t id_2 = id_1 + 1; id_2 < k; ++id_2) {
                ASSERT(manual[id_1 * k + id_2] == q_graph.get_weight(id_1, id_2));
            }
        }
        return true;
    }

    bool assert_state_pre_partitioning(p_graph_t &g, p_av_manager_t &av_manager, u64 threads) {
        // check the right number of vertices active
        ASSERT(assert_n_active_vertices(g, av_manager, threads));

        // check the right vertices are active
        ASSERT(assert_correct_vertices_active(g, av_manager, threads));

        // check only active edges
        ASSERT(assert_active_edges(g, av_manager, threads));

        // check no self-loops
        ASSERT(assert_no_self_loops(g, av_manager, threads));

        // check no duplicate edges
        ASSERT(assert_no_double_edges(g, av_manager, threads));

        return true;
    }

    bool assert_state_after_partitioning(p_graph_t &g,
                                         p_av_manager_t &av_manager,
                                         p_p_manager_t &p_manager,
                                         p_bv_manager_t &bv_manager,
                                         p_q_graph_t &q_graph,
                                         partition_t k,
                                         u64 threads) {
        // check the right number of vertices active
        ASSERT(assert_n_active_vertices(g, av_manager, threads));

        // check the right vertices are active
        ASSERT(assert_correct_vertices_active(g, av_manager, threads));

        // check only active edges
        ASSERT(assert_active_edges(g, av_manager, threads));

        // check no self-loops
        ASSERT(assert_no_self_loops(g, av_manager, threads));

        // check no duplicate edges
        ASSERT(assert_no_double_edges(g, av_manager, threads));

        // check the correct block weights
        ASSERT(assert_bweights(g, av_manager, p_manager, k, threads));

        // check the correct number of partitions
        ASSERT(assert_n_boundary_vertices(g, av_manager, p_manager, bv_manager, threads));

        // check the right vertices are boundary
        ASSERT(assert_correct_vertices_boundary(g, av_manager, p_manager, bv_manager, threads));

        // check the right vertices are boundary per block
        ASSERT(assert_correct_vertices_boundary_per_block(g, av_manager, p_manager, bv_manager, k, threads));

        // check the correct quotient graph
        ASSERT(assert_correct_quotient_graph(g, av_manager, p_manager, q_graph, k, threads));

        return true;
    }

}

#endif //HEIPROMAP_PARALLEL_ASSERT_STATES_H
