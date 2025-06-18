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

#ifndef HEIPROMAP_DEEP_ASSERT_STATE_H
#define HEIPROMAP_DEEP_ASSERT_STATE_H

#include <map>

#include "../../../commons/definitions.h"
#include "../../../commons/utils.h"
#include "../../../serial/serial_definitions_1.h"
#include "../../../serial/serial_definitions_2.h"
#include "../../../serial/serial_definitions_3.h"

namespace HeiProMap {
    bool deep_assert_no_self_loops(const graph_t &g) {
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

    bool deep_assert_no_double_edges(const graph_t &g) {
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

    bool deep_assert_correct_partition_size(const graph_t &g,
                                            const deep_p_manager_t &p_manager,
                                            const partition_t k) {
        std::vector<size_t> sizes(k, 0);

        forall_gu(g, u)
            {
                sizes[p_manager[u]] += 1;
            }
        endfor

        for (partition_t id = 0; id < k; ++id) {
            ASSERT(sizes[id] == p_manager.size(id));
        }

        return true;
    }

    bool deep_assert_bweights(const graph_t &g,
                              const deep_p_manager_t &p_manager,
                              const partition_t k) {
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

    bool deep_assert_n_boundary_vertices(const graph_t &g,
                                         const deep_p_manager_t &p_manager,
                                         const deep_bv_manager_t &bv_manager) {
        vertex_t count = 0;
        forall_gu(g, u)
            {
                partition_t u_id = p_manager[u];

                forall_guiv(g, u, i, v)
                    {
                        partition_t v_id = p_manager[v];

                        if (u_id != v_id) {
                            ASSERT(bv_manager.is_boundary(u));
                            ASSERT(bv_manager.is_boundary(v));
                            count += 1;
                            break;
                        }
                    }
                endfor
            }
        endfor

        ASSERT(count == bv_manager.size());
        return count == bv_manager.size();
    }

    bool deep_assert_correct_vertices_boundary(const graph_t &g,
                                               const deep_p_manager_t &p_manager,
                                               const deep_bv_manager_t &bv_manager) {
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

    bool deep_assert_correct_vertices_boundary_per_block(const graph_t &g,
                                                         const deep_p_manager_t &p_manager,
                                                         deep_bv_manager_t &bv_manager,
                                                         const partition_t k) {
        std::vector<std::vector<vertex_t>> manual(k);

        forall_gu(g, u)
            {
                partition_t u_id = p_manager[u];
                forall_guiv(g, u, i, v)
                    {
                        partition_t v_id = p_manager[v];
                        if (u_id != v_id) {
                            manual[u_id].push_back(u);
                            break;
                        }
                    }
                endfor
            }
        endfor

        for (partition_t id = 0; id < k; ++id) {
            std::vector<vertex_t> automatic;
            forall_bv_id_iu(bv_manager, id, i, u)
                {
                    automatic.push_back(u);
                }
            endfor

            std::sort(manual[id].begin(), manual[id].end());
            std::sort(automatic.begin(), automatic.end());
            ASSERT(no_duplicates_sorted(manual[id]));
            ASSERT(no_duplicates_sorted(automatic));

            ASSERT(manual[id] == automatic);
        }

        return true;
    }

    bool deep_assert_correct_boundary(const graph_t &g,
                                      const deep_p_manager_t &p_manager,
                                      deep_bv_manager_t &bv_manager,
                                      const partition_t k) {
        deep_assert_n_boundary_vertices(g, p_manager, bv_manager);
        deep_assert_correct_vertices_boundary(g, p_manager, bv_manager);
        deep_assert_correct_vertices_boundary_per_block(g, p_manager, bv_manager, k);
        return true;
    }

    bool deep_assert_correct_quotient_graph(const graph_t &g,
                                            const deep_p_manager_t &p_manager,
                                            const deep_q_graph_t &q_graph,
                                            const partition_t k) {
        std::map<std::pair<partition_t, partition_t>, weight_t> manual;

        forall_gu(g, u)
            {
                partition_t u_id = p_manager[u];
                forall_guivw(g, u, i, v, w)
                    {
                        partition_t v_id = p_manager[v];
                        if (u_id != v_id) {
                            manual[{u_id, v_id}] += w;
                        }
                    }
                endfor
            }
        endfor

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

    bool deep_assert_state_pre_partitioning(const graph_t &g,
                                            const deep_p_manager_t &p_manager,
                                            const partition_t k) {
        // check no self-loops
        deep_assert_no_self_loops(g);

        // check no duplicate edges
        deep_assert_no_double_edges(g);

        // check the correct partition sizes
        deep_assert_correct_partition_size(g, p_manager, k);

        return true;
    }

    bool deep_assert_state_after_partitioning(const graph_t &g,
                                              const deep_p_manager_t &p_manager,
                                              deep_bv_manager_t &bv_manager,
                                              const deep_q_graph_t &q_graph,
                                              const partition_t k) {
        // check no self-loops
        deep_assert_no_self_loops(g);

        // check no duplicate edges
        deep_assert_no_double_edges(g);

        // check the correct partition sizes
        deep_assert_correct_partition_size(g, p_manager, k);

        // check the correct block weights
        deep_assert_bweights(g, p_manager, k);

        // check the correct number of partitions
        deep_assert_n_boundary_vertices(g, p_manager, bv_manager);

        // check the right vertices are boundary
        deep_assert_correct_vertices_boundary(g, p_manager, bv_manager);

        // check the right vertices are boundary per block
        deep_assert_correct_vertices_boundary_per_block(g, p_manager, bv_manager, k);

        // check the correct quotient graph
        deep_assert_correct_quotient_graph(g, p_manager, q_graph, k);

        return true;
    }

    bool deep_assert_state_after_partitioning(const graph_t &g,
                                              const deep_p_manager_t &p_manager,
                                              const partition_t k) {
        // check no self-loops
        deep_assert_no_self_loops(g);

        // check no duplicate edges
        deep_assert_no_double_edges(g);

        // check the correct partition sizes
        deep_assert_correct_partition_size(g, p_manager, k);

        // check the correct block weights
        deep_assert_bweights(g, p_manager, k);

        return true;
    }
}

#endif //HEIPROMAP_DEEP_ASSERT_STATE_H
