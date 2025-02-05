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

#include "utils.h"
#include "../../definitions.h"
#include "../../macros.h"

namespace HeiProMap {
    template <typename TGraph, typename TActiveVertexManager>
    bool assert_n_active_vertices(TGraph& g,
                                  TActiveVertexManager& av_manager) {
        vertex_t count = 0;
        for (vertex_t u = 0; u < g.get_n(); ++u) {
            count += av_manager.is_active(u);
        }

        ASSERT(count == av_manager.get_n_active());
        return count == av_manager.get_n_active();
    }

    template <typename TGraph, typename TActiveVertexManager>
    bool assert_correct_vertices_active(TGraph& g,
                                        TActiveVertexManager& av_manager) {
        std::vector<vertex_t> manual;
        for (vertex_t u = 0; u < g.get_n(); ++u) {
            if (av_manager.is_active(u)) { manual.push_back(u); }
        }

        std::vector<vertex_t> automatic;
        for (vertex_t v : av_manager) {
            automatic.push_back(v);
        }

        std::sort(manual.begin(), manual.end());
        std::sort(automatic.begin(), automatic.end());
        ASSERT(no_duplicates_sorted(manual));
        ASSERT(no_duplicates_sorted(automatic));

        ASSERT(manual == automatic);
        return manual == automatic;
    }

    template <typename TGraph, typename TActiveVertexManager>
    bool assert_active_edges(TGraph& g,
                             TActiveVertexManager& av_manager) {
        for (vertex_t u : av_manager) {
            for (size_t i = 0; i < g.size(u); ++i) {
                const vertex_t v = g.neighbor(u, i);
                ASSERT(av_manager.is_active(v));
            }
        }

        return true;
    }

    template <typename TGraph, typename TActiveVertexManager>
    bool assert_no_self_loops(TGraph& g,
                              TActiveVertexManager& av_manager) {
        for (vertex_t u : av_manager) {
            for (size_t i = 0; i < g.size(u); ++i) {
                for (size_t j = i + 1; j < g.size(u); ++j) {
                    ASSERT(g.neighbor(u, i) != g.neighbor(u, j));
                }
            }
        }
        return true;
    }

    template <typename TGraph, typename TActiveVertexManager>
    bool assert_no_double_edges(TGraph& g,
                                TActiveVertexManager& av_manager) {
        std::vector<vertex_t> manual;
        for (vertex_t u : av_manager) {
            manual.clear();
            for (size_t i = 0; i < g.size(u); ++i) {
                manual.push_back(g.neighbor(u, i));
            }
            std::sort(manual.begin(), manual.end());
            ASSERT(no_duplicates_sorted(manual));
        }
        return true;
    }

    template <typename TGraph, typename TActiveVertexManager, typename TPartitionManager>
    bool assert_bweights(TGraph& g,
                         TActiveVertexManager& av_manager,
                         TPartitionManager& p_manager,
                         partition_t k) {
        std::vector<weight_t> weights(k, 0);
        for (vertex_t u : av_manager) {
            partition_t u_id = p_manager[u];

            weights[u_id] += g.get_weight(u);
        }

        for (partition_t id = 0; id < k; ++id) {
            ASSERT(weights[id] == p_manager.get_bweight(id));
        }

        return true;
    }

    template <typename TGraph, typename TActiveVertexManager, typename TPartitionManager, typename TBoundaryVertexManager>
    bool assert_n_boundary_vertices(TGraph& g,
                                    TActiveVertexManager& av_manager,
                                    TPartitionManager& p_manager,
                                    TBoundaryVertexManager& bv_manager) {
        vertex_t count = 0;
        for (vertex_t u : av_manager) {
            partition_t u_id = p_manager[u];

            for (size_t i = 0; i < g.size(u); ++i) {
                partition_t v_id = p_manager[g.neighbor(u, i)];

                if (u_id != v_id) {
                    count += 1;
                    break;
                }
            }
        }

        ASSERT(count == bv_manager.get_n_boundary());
        return count == bv_manager.get_n_boundary();
    }

    template <typename TGraph, typename TActiveVertexManager, typename TPartitionManager, typename TBoundaryVertexManager>
    bool assert_correct_vertices_boundary(TGraph& g,
                                          TActiveVertexManager& av_manager,
                                          TPartitionManager& p_manager,
                                          TBoundaryVertexManager& bv_manager) {
        std::vector<vertex_t> manual;
        for (vertex_t u : av_manager) {
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
        for (vertex_t u : bv_manager) {
            automatic.push_back(u);
        }

        std::sort(manual.begin(), manual.end());
        std::sort(automatic.begin(), automatic.end());
        ASSERT(no_duplicates_sorted(manual));
        ASSERT(no_duplicates_sorted(automatic));

        ASSERT(manual == automatic);
        return manual == automatic;
    }

    template <typename TGraph, typename TActiveVertexManager, typename TPartitionManager, typename TBoundaryVertexManager>
    bool assert_correct_vertices_boundary_per_block(TGraph& g,
                                                    TActiveVertexManager& av_manager,
                                                    TPartitionManager& p_manager,
                                                    TBoundaryVertexManager& bv_manager,
                                                    partition_t k) {
        for (partition_t id = 0; id < k; ++id) {
            std::vector<vertex_t> manual;

            for (vertex_t u : av_manager) {
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
            for (vertex_t u : bv_manager[id]) {
                automatic.push_back(u);
            }

            std::sort(manual.begin(), manual.end());
            std::sort(automatic.begin(), automatic.end());
            ASSERT(no_duplicates_sorted(manual));
            ASSERT(no_duplicates_sorted(automatic));

            ASSERT(manual == automatic);
        }
        return true;
    }

    template <typename TGraph, typename TActiveVertexManager>
    bool assert_state_pre_partitioning([[maybe_unused]] TGraph& g,
                                       [[maybe_unused]] TActiveVertexManager& av_manager) {
        // check the right number of vertices active
        ASSERT(assert_n_active_vertices(g, av_manager));

        // check the right vertices are active
        ASSERT(assert_correct_vertices_active(g, av_manager));

        // check only active edges
        ASSERT(assert_active_edges(g, av_manager));

        // check no self-loops
        ASSERT(assert_no_self_loops(g, av_manager));

        // check no duplicate edges
        ASSERT(assert_no_double_edges(g, av_manager));

        return true;
    }

    template <typename TGraph, typename TActiveVertexManager, typename TPartitionManager, typename TBoundaryVertexManager>
    bool assert_state_after_partitioning([[maybe_unused]] TGraph& g,
                                         [[maybe_unused]] TActiveVertexManager& av_manager,
                                         [[maybe_unused]] TPartitionManager& p_manager,
                                         [[maybe_unused]] TBoundaryVertexManager& bv_manager,
                                         [[maybe_unused]] partition_t k) {
        // check the right number of vertices active
        ASSERT(assert_n_active_vertices(g, av_manager));

        // check the right vertices are active
        ASSERT(assert_correct_vertices_active(g, av_manager));

        // check only active edges
        ASSERT(assert_active_edges(g, av_manager));

        // check no self-loops
        ASSERT(assert_no_self_loops(g, av_manager));

        // check no duplicate edges
        ASSERT(assert_no_double_edges(g, av_manager));

        // check the correct block weights
        ASSERT(assert_bweights(g, av_manager, p_manager, k));

        // check the correct number of partitions
        ASSERT(assert_n_boundary_vertices(g, av_manager, p_manager, bv_manager));

        // check the right vertices are boundary
        ASSERT(assert_correct_vertices_boundary(g, av_manager, p_manager, bv_manager));

        // check the right vertices are boundary per block
        ASSERT(assert_correct_vertices_boundary_per_block(g, av_manager, p_manager, bv_manager, k));

        return true;
    }
}


#endif //HEIPROMAP_ASSERT_STATE_H
