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

#include "../../definitions.h"
#include "../../macros.h"
#include "../../commons/utils.h"
#include "../../parallel/parallel_definitions_1.h"

namespace HeiProMap {
    bool assert_no_self_loops(graph_t& g) {
        forall_gu(g, u)
            {
                for (size_t i = 0; i < g.size(u); ++i) {
                    for (size_t j = i + 1; j < g.size(u); ++j) {
                        ASSERT(g.neighbor(u, i) != g.neighbor(u, j));
                    }
                }
            }
        endfor
        return true;
    }

    bool assert_no_double_edges(graph_t& g) {
        std::vector<vertex_t> manual;
        forall_gu(g, u)
            {
                manual.clear();
                for (size_t i = 0; i < g.size(u); ++i) {
                    manual.push_back(g.neighbor(u, i));
                }
                std::sort(manual.begin(), manual.end());
                ASSERT(no_duplicates_sorted(manual));
            }
        endfor
        return true;
    }

    bool assert_bweights(graph_t& g,
                         p_manager_t& p_manager,
                         partition_t k) {
        std::vector<weight_t> weights(k, 0);
        forall_gu(g, u)
            {
                partition_t u_id = p_manager[u];

                weights[u_id] += g.weight(u);
            }
        endfor

        for (partition_t id = 0; id < k; ++id) {
            ASSERT(weights[id] == p_manager.get_bweight(id));
        }

        return true;
    }

    bool assert_n_boundary_vertices(graph_t& g,
                                    p_manager_t& p_manager,
                                    bv_manager_t& bv_manager) {
        vertex_t count = 0;
        forall_gu(g, u)
            {
                partition_t u_id = p_manager[u];

                for (size_t i = 0; i < g.size(u); ++i) {
                    partition_t v_id = p_manager[g.neighbor(u, i)];

                    if (u_id != v_id) {
                        count += 1;
                        break;
                    }
                }
            }
        endfor

        ASSERT(count == bv_manager.size());
        return count == bv_manager.size();
    }

    bool assert_correct_vertices_boundary(graph_t& g,
                                          p_manager_t& p_manager,
                                          bv_manager_t& bv_manager) {
        std::vector<vertex_t> manual;
        forall_gu(g, u)
            {
                partition_t u_id = p_manager[u];

                for (size_t i = 0; i < g.size(u); ++i) {
                    partition_t v_id = p_manager[g.neighbor(u, i)];
                    if (u_id != v_id) {
                        manual.push_back(u);
                        break;
                    }
                }
            }
        endfor

        std::vector<vertex_t> automatic;
        forall_bv_iu(bv_manager, i, u)
            {
                automatic.push_back(u);
            }
        endfor

        std::sort(manual.begin(), manual.end());
        std::sort(automatic.begin(), automatic.end());
        ASSERT(no_duplicates_sorted(manual));
        ASSERT(no_duplicates_sorted(automatic));

        ASSERT(manual == automatic);
        return manual == automatic;
    }

    bool assert_correct_vertices_boundary_per_block(graph_t& g,
                                                    p_manager_t& p_manager,
                                                    bv_manager_t& bv_manager,
                                                    partition_t k) {
        for (partition_t id = 0; id < k; ++id) {
            std::vector<vertex_t> manual;

            forall_gu(g, u)
                {
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
            endfor

            std::vector<vertex_t> automatic;
            forall_bv_id_iu(bv_manager, id, i, u)
                {
                    automatic.push_back(u);
                }
            endfor

            std::sort(manual.begin(), manual.end());
            std::sort(automatic.begin(), automatic.end());
            ASSERT(no_duplicates_sorted(manual));
            ASSERT(no_duplicates_sorted(automatic));

            ASSERT(manual == automatic);
        }
        return true;
    }

    bool assert_correct_quotient_graph(graph_t& g,
                                       p_manager_t& p_manager,
                                       q_graph_t& q_graph,
                                       partition_t k) {
        std::vector<weight_t> manual(k * k, 0);

        forall_gu(g, u)
            {
                partition_t u_id = p_manager[u];
                forall_guivw(g, u, i, v, w)
                    {
                        partition_t v_id = p_manager[v];

                        manual[u_id * k + v_id] += w;
                        manual[v_id * k + u_id] += w;
                    }
                endfor
            }
        endfor

        for (partition_t id_1 = 0; id_1 < k; ++id_1) {
            for (partition_t id_2 = id_1 + 1; id_2 < k; ++id_2) {
                ASSERT(manual[id_1 * k + id_2] == q_graph.get_weight(id_1, id_2));
            }
        }
        return true;
    }

    bool assert_state_pre_partitioning(graph_t& g) {
        // check no self-loops
        ASSERT(assert_no_self_loops(g));

        // check no duplicate edges
        ASSERT(assert_no_double_edges(g));

        return true;
    }

    bool assert_state_after_partitioning(graph_t& g,
                                         p_manager_t& p_manager,
                                         bv_manager_t& bv_manager,
                                         q_graph_t& q_graph,
                                         partition_t k) {
        // check no self-loops
        ASSERT(assert_no_self_loops(g));

        // check no duplicate edges
        ASSERT(assert_no_double_edges(g));

        // check the correct block weights
        ASSERT(assert_bweights(g,  p_manager, k));

        // check the correct number of partitions
        ASSERT(assert_n_boundary_vertices(g,  p_manager, bv_manager));

        // check the right vertices are boundary
        ASSERT(assert_correct_vertices_boundary(g,  p_manager, bv_manager));

        // check the right vertices are boundary per block
        ASSERT(assert_correct_vertices_boundary_per_block(g,  p_manager, bv_manager, k));

        // check the correct quotient graph
        ASSERT(assert_correct_quotient_graph(g, p_manager, q_graph, k));

        return true;
    }
}


#endif //HEIPROMAP_ASSERT_STATE_H
