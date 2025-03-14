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

#ifndef HEIPROMAP_PARALLEL_LABEL_PROPAGATION_REFINEMENT_H
#define HEIPROMAP_PARALLEL_LABEL_PROPAGATION_REFINEMENT_H

#include <omp.h>

#include "../interfaces/IParallelRefiner.h"

namespace HeiProMap {

    template<typename TParallelGraph,
             typename TParallelActiveVertexManager,
             typename TParallelBoundaryVertexManager,
             typename TParallelPartitionManager,
             typename TParallelDistanceOracle,
             typename TParallelQuotientGraph>
    class ParallelLabelPropagationRefinement : public IParallelRefiner<TParallelGraph, TParallelActiveVertexManager, TParallelBoundaryVertexManager, TParallelPartitionManager, TParallelDistanceOracle, TParallelQuotientGraph> {
    private:
        TParallelGraph *m_p_g = nullptr;
        TParallelActiveVertexManager *m_p_av_manager = nullptr;
        TParallelBoundaryVertexManager *m_p_bv_manager = nullptr;
        TParallelPartitionManager *m_p_p_manger = nullptr;
        TParallelDistanceOracle *m_p_d_oracle = nullptr;
        TParallelQuotientGraph *m_p_qgraph = nullptr;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;
        weight_t m_lmax = 0;
        u64 m_n_threads = 1;

        std::vector<s32> m_used;
        std::vector<u8> m_neighbor_changed;
        s32 m_mark = -1;

    public:
        ParallelLabelPropagationRefinement() = default;

        void initialize(TParallelGraph *t_p_g,
                        TParallelActiveVertexManager *t_p_av_manager,
                        TParallelBoundaryVertexManager *t_p_bv_manager,
                        TParallelPartitionManager *t_p_p_manger,
                        TParallelDistanceOracle *t_p_d_oracle,
                        TParallelQuotientGraph *t_p_qgraph,
                        std::vector<partition_t> &t_hierarchy,
                        std::vector<weight_t> &t_distance,
                        weight_t t_lmax,
                        u64 t_n_threads) final {
            m_p_g = t_p_g;
            m_p_av_manager = t_p_av_manager;
            m_p_bv_manager = t_p_bv_manager;
            m_p_p_manger = t_p_p_manger;
            m_p_d_oracle = t_p_d_oracle;
            m_p_qgraph = t_p_qgraph;
            m_hierarchy = t_hierarchy;
            m_distance = t_distance;
            m_lmax = t_lmax;
            m_n_threads = t_n_threads;

            m_used.resize(m_p_g->get_n(), -1);
            m_neighbor_changed.resize(m_p_g->get_n());
        }

        void refine() final {
            refine_pm_aware();
            return;

            ASSERT(m_p_g != nullptr);
            ASSERT(m_p_av_manager != nullptr);
            ASSERT(m_p_bv_manager != nullptr);
            ASSERT(m_p_p_manger != nullptr);
            ASSERT(m_p_d_oracle != nullptr);
            ASSERT(m_p_qgraph != nullptr);

            std::vector<vertex_t> best_us(m_n_threads);
            std::vector<partition_t> best_ids(m_n_threads);
            std::vector<s64> best_qap_deltas(m_n_threads);

            u64 max_iterations = 100;
            bool move_occurred = true;

            for (u64 iteration = 0; iteration < max_iterations && move_occurred; ++iteration) {
                move_occurred = false;
                std::fill(best_qap_deltas.begin(), best_qap_deltas.end(), -std::numeric_limits<s64>::max());

#pragma omp parallel default(none) shared(best_us, best_ids, best_qap_deltas, std::cerr) num_threads(m_n_threads)
                {
                    u64 thread_id = omp_get_thread_num();
                    vertex_t local_best_u;
                    partition_t local_best_id;
                    s64 local_best_qap_delta = -std::numeric_limits<s64>::max();

                    for (size_t i = thread_id; i < m_p_bv_manager->get_n_boundary(); i += m_n_threads) {
                        vertex_t u = m_p_bv_manager->get_vertex(i);
                        ASSERT(m_p_av_manager->is_active(u));
                        ASSERT(m_p_bv_manager->is_boundary(u));

                        // this vertex can be moved
                        partition_t u_id = (*m_p_p_manger)[u];
                        weight_t u_weight = m_p_g->get_weight(u);
                        s64 u_qap = get_u_qap(*m_p_g, u, *m_p_p_manger, *m_p_d_oracle);

                        for (size_t j = 0; j < m_p_g->size(u); ++j) {
                            vertex_t v = m_p_g->neighbor(u, j);
                            partition_t v_id = (*m_p_p_manger)[v];
                            weight_t v_id_weight = m_p_p_manger->get_bweight(v_id);

                            if (u_id != v_id && v_id_weight + u_weight <= m_lmax) {
                                // check if v_id is a better option
                                s64 qap_delta = u_qap - get_u_qap(*m_p_g, u, v_id, *m_p_p_manger, *m_p_d_oracle);
                                if (qap_delta > local_best_qap_delta && qap_delta > 0) {
                                    local_best_u = u;
                                    local_best_id = v_id;
                                    local_best_qap_delta = qap_delta;
                                }
                            }
                        }

                        best_us[thread_id] = local_best_u;
                        best_ids[thread_id] = local_best_id;
                        best_qap_deltas[thread_id] = local_best_qap_delta;
                    }
                }

                size_t max_idx = argmax(best_qap_deltas);

                if (best_qap_deltas[max_idx] > 0) {
                    // change the partition of the best vertex
                    vertex_t u = best_us[max_idx];
                    partition_t old_id = (*m_p_p_manger)[u];
                    partition_t new_id = best_ids[max_idx];
                    m_p_p_manger->move(u, old_id, new_id);
                    m_p_bv_manager->move(u, old_id, new_id);
                    move_occurred = true;
                }
            }
        }

        void refine_pm_aware() {
            ASSERT(m_p_g != nullptr);
            ASSERT(m_p_av_manager != nullptr);
            ASSERT(m_p_bv_manager != nullptr);
            ASSERT(m_p_p_manger != nullptr);
            ASSERT(m_p_d_oracle != nullptr);
            ASSERT(m_p_qgraph != nullptr);

            std::vector<vertex_t> best_us(m_n_threads);
            std::vector<partition_t> best_ids(m_n_threads);
            std::vector<s64> best_qap_deltas(m_n_threads);

            u64 max_global_iterations = 10;
            bool global_move_occurred = true;

            std::vector<u64> local_max_iterations = {1, 3, 5};
            bool local_move_occurred = true;

            for (u64 global_iteration = 0; global_iteration < max_global_iterations && global_move_occurred; ++global_iteration) {
                global_move_occurred = false;

                for(partition_t h_i = 0; h_i < local_max_iterations.size(); ++h_i) {
                    partition_t h = local_max_iterations.size() - 1 - h_i;
                    for (u64 local_iteration = 0; local_iteration < local_max_iterations[h] && local_move_occurred; ++local_iteration) {
                        m_mark += 1;
                        local_move_occurred = false;
                        std::fill(best_qap_deltas.begin(), best_qap_deltas.end(), -std::numeric_limits<s64>::max());

#pragma omp parallel default(none) firstprivate(h) shared(best_us, best_ids, best_qap_deltas, std::cerr) num_threads(m_n_threads)
                        {
                            u64 thread_id = omp_get_thread_num();
                            vertex_t local_best_u = 0;
                            partition_t local_best_id = 0;
                            s64 local_best_qap_delta = -std::numeric_limits<s64>::max();

                            for (size_t i = thread_id; i < m_p_bv_manager->get_n_boundary(); i += m_n_threads) {
                                vertex_t u = m_p_bv_manager->get_vertex(i);
                                ASSERT(m_p_av_manager->is_active(u));
                                ASSERT(m_p_bv_manager->is_boundary(u));

                                // this vertex can be moved
                                partition_t u_id = (*m_p_p_manger)[u];
                                weight_t u_weight = m_p_g->get_weight(u);
                                s64 u_qap = get_u_qap(*m_p_g, u, *m_p_p_manger, *m_p_d_oracle);

                                for (size_t j = 0; j < m_p_g->size(u); ++j) {
                                    vertex_t v = m_p_g->neighbor(u, j);
                                    partition_t v_id = (*m_p_p_manger)[v];
                                    weight_t v_id_weight = m_p_p_manger->get_bweight(v_id);

                                    if (u_id != v_id && v_id_weight + u_weight <= m_lmax && m_p_d_oracle->get_h(u_id, v_id) == h) {
                                        // check if v_id is a better option
                                        s64 qap_delta = u_qap - get_u_qap(*m_p_g, u, v_id, *m_p_p_manger, *m_p_d_oracle);
                                        if (qap_delta > local_best_qap_delta && qap_delta > 0) {
                                            local_best_u = u;
                                            local_best_id = v_id;
                                            local_best_qap_delta = qap_delta;
                                        }
                                    }
                                }

                                best_us[thread_id] = local_best_u;
                                best_ids[thread_id] = local_best_id;
                                best_qap_deltas[thread_id] = local_best_qap_delta;
                            }
                        }

                        size_t max_idx = argmax(best_qap_deltas);

                        if (best_qap_deltas[max_idx] > 0) {
                            // change the partition of the best vertex
                            vertex_t u = best_us[max_idx];
                            partition_t old_id = (*m_p_p_manger)[u];
                            partition_t new_id = best_ids[max_idx];
                            m_p_p_manger->move(u, old_id, new_id);
                            m_p_bv_manager->move(u, old_id, new_id);
                            m_p_qgraph->move(u, old_id, new_id);
                            local_move_occurred = true;
                        }
                    }
                    global_move_occurred |= local_move_occurred;
                }
            }
        }

    };
}

#endif //HEIPROMAP_PARALLEL_LABEL_PROPAGATION_REFINEMENT_H
