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
#include "../definitions_1.h"
#include "../definitions_2.h"
#include "../definitions_3.h"

namespace HeiProMap {
    bool assert_no_self_loops(const graph_t &g) {
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

    bool assert_no_double_edges(const graph_t &g) {
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

    bool assert_correct_partition_size(const graph_t &g,
                                       const p_manager_t &p_manager,
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

    bool assert_bweights(const graph_t &g,
                         const p_manager_t &p_manager,
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

    bool assert_correct_vertices_boundary(const graph_t &g,
                                          const p_manager_t &p_manager,
                                          const bv_manager_t &bv_manager) {
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
        for (partition_t id = 0; id < bv_manager.get_k(); ++id) {
            forall_bv_id_iu(bv_manager, id, i, u)
            {
                automatic.push_back(u);
            }
            endfor
        }

        std::sort(manual.begin(), manual.end());
        std::sort(automatic.begin(), automatic.end());
        ASSERT(no_duplicates_sorted(manual));
        ASSERT(no_duplicates_sorted(automatic));

        ASSERT(manual == automatic);
        return manual == automatic;
    }

    bool assert_correct_vertices_boundary_per_block(const graph_t &g,
                                                    const p_manager_t &p_manager,
                                                    bv_manager_t &bv_manager,
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

    bool assert_correct_boundary([[maybe_unused]] const graph_t &g,
                                 [[maybe_unused]] const p_manager_t &p_manager,
                                 [[maybe_unused]] bv_manager_t &bv_manager,
                                 [[maybe_unused]] const partition_t k) {
        ASSERT(assert_correct_vertices_boundary(g, p_manager, bv_manager));
        ASSERT(assert_correct_vertices_boundary_per_block(g, p_manager, bv_manager, k));
        return true;
    }

    bool assert_correct_quotient_graph(const graph_t &g,
                                       const p_manager_t &p_manager,
                                       const q_graph_t &q_graph,
                                       [[maybe_unused]] const partition_t k) {
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

    bool assert_state_pre_partitioning([[maybe_unused]] const graph_t &g,
                                       [[maybe_unused]] const p_manager_t &p_manager,
                                       [[maybe_unused]] const partition_t k) {
        // check no self-loops
        ASSERT(assert_no_self_loops(g));

        // check no duplicate edges
        ASSERT(assert_no_double_edges(g));

        // check the correct partition sizes
        ASSERT(assert_correct_partition_size(g, p_manager, k));

        return true;
    }

    bool assert_state_after_partitioning([[maybe_unused]] const graph_t &g,
                                         [[maybe_unused]] const p_manager_t &p_manager,
                                         [[maybe_unused]] bv_manager_t &bv_manager,
                                         [[maybe_unused]] const q_graph_t &q_graph,
                                         [[maybe_unused]] const partition_t k) {
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

        return true;
    }

    bool assert_state_after_partitioning([[maybe_unused]] const graph_t &g,
                                         [[maybe_unused]] const p_manager_t &p_manager,
                                         [[maybe_unused]] const partition_t k) {
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
