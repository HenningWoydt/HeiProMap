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

#ifndef HEIPROMAP_PARALLEL_QUOTIENT_GRAPH_REFINEMENT_H
#define HEIPROMAP_PARALLEL_QUOTIENT_GRAPH_REFINEMENT_H

#include <queue>
#include <omp.h>
#include <atomic>

#include "../interfaces/IParallelRefiner.h"
#include "../../serial/datastructures/indexed_max_heap.h"

namespace HeiProMap {

    template<typename TParallelGraph,
             typename TParallelActiveVertexManager,
             typename TParallelBoundaryVertexManager,
             typename TParallelPartitionManager,
             typename TParallelDistanceOracle,
             typename TParallelQuotientGraph>
    class ParallelQuotientGraphRefinement : public IParallelRefiner<TParallelGraph, TParallelActiveVertexManager, TParallelBoundaryVertexManager, TParallelPartitionManager, TParallelDistanceOracle, TParallelQuotientGraph> {
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
        std::atomic<s32> m_mark = -1;

    public:
        ParallelQuotientGraphRefinement() = default;

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
            ASSERT(m_p_g != nullptr);
            ASSERT(m_p_av_manager != nullptr);
            ASSERT(m_p_bv_manager != nullptr);
            ASSERT(m_p_p_manger != nullptr);
            ASSERT(m_p_d_oracle != nullptr);
            ASSERT(m_p_qgraph != nullptr);

            std::vector<std::vector<s32>> threads_vertex_moved(m_n_threads, std::vector<s32>(m_p_g->get_n(), -1));
            std::vector<std::vector<vertex_t>> threads_moves(m_n_threads);
            std::vector<IndexedMaxHeap<Swap>> threads_max_hep(m_n_threads, IndexedMaxHeap<Swap>(m_p_g->get_n()));
            m_mark = 0;

            std::vector<std::vector<QGraphUV>> matching_hierarchy;

            size_t max_global_iterations = 3;
            size_t max_iterations = 1;

            for (size_t global_iteration = 0; global_iteration < max_global_iterations; ++global_iteration) {
                matching_hierarchy.clear();
                m_p_qgraph->get_matching_hierarchy(matching_hierarchy);

                while(matching_hierarchy.back().size() <= m_n_threads){
                    matching_hierarchy.pop_back();
                }

                for (std::vector<QGraphUV> &matchings: matching_hierarchy) {
#pragma omp parallel for default(none) firstprivate(global_iteration, max_iterations) shared(matchings, threads_vertex_moved, threads_moves, threads_max_hep, std::cout, std::cerr) num_threads(m_n_threads)
                    for (QGraphUV uv: matchings) {
                        u64 thread_id = omp_get_thread_num();

                        partition_t b1_id = uv.u;
                        partition_t b2_id = uv.v;

                        // do fm refinement on blocks b1_id and b2_id
                        std::vector<s32> &vertex_moved = threads_vertex_moved[thread_id];
                        std::vector<vertex_t> &moves = threads_moves[thread_id];
                        IndexedMaxHeap<Swap> &max_heap = threads_max_hep[thread_id];

                        for (size_t iteration = 0; iteration < max_iterations; ++iteration) {
                            s32 mark = m_mark.fetch_add(1);

                            moves.clear();
                            max_heap.clear();

                            // initialize all swaps and their deltas from b1
                            for (m_p_bv_manager->reset_iterator(b1_id); m_p_bv_manager->available(b1_id); m_p_bv_manager->next(b1_id)) {
                                vertex_t u = m_p_bv_manager->get(b1_id);
                                ASSERT((*m_p_p_manger)[u] == b1_id);

                                s64 qap_delta = get_u_qap_delta(*m_p_g, u, b1_id, b2_id, *m_p_p_manger, *m_p_d_oracle);
                                Swap s(u, qap_delta);
                                max_heap.push(u, s);
                            }

                            // initialize all swaps and their deltas from b2
                            for (m_p_bv_manager->reset_iterator(b2_id); m_p_bv_manager->available(b2_id); m_p_bv_manager->next(b2_id)) {
                                vertex_t v = m_p_bv_manager->get(b2_id);
                                ASSERT((*m_p_p_manger)[v] == b2_id);

                                s64 qap_delta = get_u_qap_delta(*m_p_g, v, b2_id, b1_id, *m_p_p_manger, *m_p_d_oracle);
                                Swap s(v, qap_delta);
                                max_heap.push(v, s);
                            }

                            size_t best_move_idx = std::numeric_limits<size_t>::max();
                            s64 best_cum_qap_delta = 0;
                            s64 curr_cum_qap_delta = 0;
                            s64 n_neg_deltas = 0;
                            // pop all elements until the max heap is empty
                            while (!max_heap.empty()) {
                                vertex_t u = max_heap.top().u;
                                s64 qap_delta = max_heap.top().qap_delta;
                                max_heap.pop();

                                partition_t u_old_id = (*m_p_p_manger)[u];
                                partition_t u_new_id = (*m_p_p_manger)[u] == b1_id ? b2_id : b1_id;

                                ASSERT(u_old_id == b1_id || u_old_id == b2_id);
                                ASSERT(u_new_id == b1_id || u_new_id == b2_id);

                                if (!m_p_p_manger->is_boundary(u)) {
                                // if (!m_p_bv_manager->is_boundary(u)) {
                                    // vertex could not be boundary anymore
                                    continue;
                                }

                                if(qap_delta < 0){
                                    n_neg_deltas += 1;
                                    if(n_neg_deltas >= 4){
                                        break;
                                    }
                                } else{
                                    n_neg_deltas = 0;
                                }

                                // make the move
                                m_p_p_manger->move(u, u_old_id, u_new_id);
                                // m_p_bv_manager->move(u, u_old_id, u_new_id);

                                moves.push_back(u);
                                curr_cum_qap_delta += qap_delta;
                                if (curr_cum_qap_delta > best_cum_qap_delta && m_p_p_manger->get_bweight(b1_id) <= m_lmax && m_p_p_manger->get_bweight(b2_id) <= m_lmax) {
                                    // new best and not overloaded
                                    best_cum_qap_delta = curr_cum_qap_delta;
                                    best_move_idx = moves.size();
                                }
                                vertex_moved[u] = mark;

                                // update the neighbors of u if they are in b1 or b2, are boundary vertices and have not been moved
                                for (size_t i = 0; i < m_p_g->size(u); ++i) {
                                    vertex_t v = m_p_g->neighbor(u, i);
                                    partition_t v_old_id = (*m_p_p_manger)[v];
                                    if (vertex_moved[v] == mark || (v_old_id != b1_id && v_old_id != b2_id) || !m_p_p_manger->is_boundary(v)) { //!m_p_bv_manager->is_boundary(v)) {
                                        // vertex cannot be moved
                                        continue;
                                    }

                                    // update the qap delta
                                    partition_t v_new_id = (*m_p_p_manger)[v] == b1_id ? b2_id : b1_id;
                                    s64 v_qap_delta = get_u_qap_delta(*m_p_g, v, v_old_id, v_new_id, *m_p_p_manger, *m_p_d_oracle);

                                    // insert/update the qap delta
                                    max_heap.push_update(v, {v, v_qap_delta});
                                }
                            }

                            size_t n_reverts = best_move_idx == std::numeric_limits<size_t>::max() ? moves.size() : moves.size() - best_move_idx;
                            // revert to best state
                            /*
                            for (size_t idx = 0; idx < n_reverts; ++idx) {
                                vertex_t v = moves[moves.size() - 1 - idx];
                                partition_t v_old_id = (*m_p_p_manger)[v];
                                partition_t v_new_id = (*m_p_p_manger)[v] == b1_id ? b2_id : b1_id;

                                ASSERT(v_old_id == b1_id || v_old_id == b2_id);
                                ASSERT(v_new_id == b1_id || v_new_id == b2_id);

                                // revert the move
                                m_p_p_manger->move(v, v_old_id, v_new_id);
                                // m_p_bv_manager->move(v, v_old_id, v_new_id);
                            }
                             */

                            // make boundary vertices move to best vertex
                            for(size_t idx = 0; idx < moves.size(); ++idx){
                                vertex_t v = moves[moves.size() - 1 - idx];
                                partition_t v_old_id = (*m_p_p_manger)[v];
                                partition_t v_new_id = (*m_p_p_manger)[v] == b1_id ? b2_id : b1_id;

                                ASSERT(v_old_id == b1_id || v_old_id == b2_id);
                                ASSERT(v_new_id == b1_id || v_new_id == b2_id);

                                // revert the move
                                m_p_p_manger->move(v, v_old_id, v_new_id);
                                // m_p_bv_manager->move(v, v_old_id, v_new_id);
                            }
                            for(size_t idx = 0; idx < moves.size() - n_reverts; ++idx){
                                vertex_t v = moves[idx];
                                partition_t v_old_id = (*m_p_p_manger)[v];
                                partition_t v_new_id = (*m_p_p_manger)[v] == b1_id ? b2_id : b1_id;

                                ASSERT(v_old_id == b1_id || v_old_id == b2_id);
                                ASSERT(v_new_id == b1_id || v_new_id == b2_id);

                                // make the move
                                m_p_p_manger->move(v, v_old_id, v_new_id);
                                m_p_bv_manager->move(v, v_old_id, v_new_id);
                                m_p_qgraph->move(v, v_old_id, v_new_id);
                            }
                        }
                    }
                }
            }
        }

        void refine_inaccurate() {
            ASSERT(m_p_g != nullptr);
            ASSERT(m_p_av_manager != nullptr);
            ASSERT(m_p_bv_manager != nullptr);
            ASSERT(m_p_p_manger != nullptr);
            ASSERT(m_p_d_oracle != nullptr);
            ASSERT(m_p_qgraph != nullptr);

            std::vector<std::vector<s32>> threads_vertex_moved(m_n_threads, std::vector<s32>(m_p_g->get_n(), -1));
            std::vector<std::vector<vertex_t>> threads_moves(m_n_threads);
            std::vector<IndexedMaxHeap<Swap>> threads_max_hep(m_n_threads, IndexedMaxHeap<Swap>(m_p_g->get_n()));
            m_mark = 0;

            std::vector<std::vector<QGraphUV>> matching_hierarchy;

            size_t max_global_iterations = 10;
            size_t max_iterations = 1;

            for (size_t global_iteration = 0; global_iteration < max_global_iterations; ++global_iteration) {
                matching_hierarchy.clear();
                m_p_qgraph->get_matching_hierarchy(matching_hierarchy);

                for (std::vector<QGraphUV> &matchings: matching_hierarchy) {
#pragma omp parallel for default(none) firstprivate(global_iteration, max_iterations) shared(matchings, threads_vertex_moved, threads_moves, threads_max_hep, std::cout, std::cerr) num_threads(m_n_threads)
                    for (QGraphUV uv: matchings) {
                        u64 thread_id = omp_get_thread_num();

                        partition_t b1_id = uv.u;
                        partition_t b2_id = uv.v;

                        // do fm refinement on blocks b1_id and b2_id
                        std::vector<s32> &vertex_moved = threads_vertex_moved[thread_id];
                        std::vector<vertex_t> &moves = threads_moves[thread_id];
                        IndexedMaxHeap<Swap> &max_heap = threads_max_hep[thread_id];

                        for (size_t iteration = 0; iteration < max_iterations; ++iteration) {
                            s32 mark = m_mark.fetch_add(1);

                            moves.clear();
                            max_heap.clear();

                            // initialize all swaps and their deltas from b1
                            for (m_p_bv_manager->reset_iterator(b1_id); m_p_bv_manager->available(b1_id); m_p_bv_manager->next(b1_id)) {
                                vertex_t u = m_p_bv_manager->get(b1_id);
                                ASSERT((*m_p_p_manger)[u] == b1_id);

                                s64 qap_delta = get_u_qap_delta(*m_p_g, u, b1_id, b2_id, *m_p_p_manger, *m_p_d_oracle);
                                Swap s(u, qap_delta);
                                max_heap.push(u, s);
                            }

                            // initialize all swaps and their deltas from b2
                            for (m_p_bv_manager->reset_iterator(b2_id); m_p_bv_manager->available(b2_id); m_p_bv_manager->next(b2_id)) {
                                vertex_t v = m_p_bv_manager->get(b2_id);
                                ASSERT((*m_p_p_manger)[v] == b2_id);

                                s64 qap_delta = get_u_qap_delta(*m_p_g, v, b2_id, b1_id, *m_p_p_manger, *m_p_d_oracle);
                                Swap s(v, qap_delta);
                                max_heap.push(v, s);
                            }

                            size_t best_move_idx = std::numeric_limits<size_t>::max();
                            s64 best_cum_qap_delta = 0;
                            s64 curr_cum_qap_delta = 0;
                            // pop all elements until the max heap is empty
                            while (!max_heap.empty()) {
                                vertex_t u = max_heap.top().u;
                                s64 qap_delta = max_heap.top().qap_delta;
                                max_heap.pop();

                                if (!m_p_p_manger->is_boundary(u)) {
                                    // vertex could not be boundary anymore
                                    continue;
                                }

                                partition_t u_old_id = (*m_p_p_manger)[u];
                                partition_t u_new_id = (*m_p_p_manger)[u] == b1_id ? b2_id : b1_id;
                                ASSERT(u_old_id == b1_id || u_old_id == b2_id);
                                ASSERT(u_new_id == b1_id || u_new_id == b2_id);

                                s64 true_qap_delta = get_u_qap_delta(*m_p_g, u, u_old_id, u_new_id, *m_p_p_manger, *m_p_d_oracle);
                                if(qap_delta > true_qap_delta){
                                    // reinsert into heap, since true delta was smaller
                                    max_heap.push(u, {u, true_qap_delta});
                                    continue;
                                }
                                qap_delta = true_qap_delta;

                                // make the move
                                m_p_p_manger->move(u, u_old_id, u_new_id);
                                // m_p_bv_manager->move(u, u_old_id, u_new_id);

                                moves.push_back(u);
                                curr_cum_qap_delta += qap_delta;
                                if (curr_cum_qap_delta > best_cum_qap_delta && m_p_p_manger->get_bweight(b1_id) <= m_lmax && m_p_p_manger->get_bweight(b2_id) <= m_lmax) {
                                    // new best and not overloaded
                                    best_cum_qap_delta = curr_cum_qap_delta;
                                    best_move_idx = moves.size();
                                }
                                vertex_moved[u] = mark;

                                // update the neighbors of u if they are in b1 or b2, are boundary vertices and have not been moved
                                for (size_t i = 0; i < m_p_g->size(u); ++i) {
                                    vertex_t v = m_p_g->neighbor(u, i);
                                    partition_t v_old_id = (*m_p_p_manger)[v];
                                    if (vertex_moved[v] == mark || (v_old_id != b1_id && v_old_id != b2_id) || !m_p_p_manger->is_boundary(v) || max_heap.entry_exists(v)) {
                                        // vertex cannot be moved
                                        continue;
                                    }

                                    // update the qap delta
                                    partition_t v_new_id = (*m_p_p_manger)[v] == b1_id ? b2_id : b1_id;
                                    s64 v_qap_delta = get_u_qap_delta(*m_p_g, v, v_old_id, v_new_id, *m_p_p_manger, *m_p_d_oracle);
                                    max_heap.push(v, {v, v_qap_delta});
                                }
                            }

                            size_t n_reverts = best_move_idx == std::numeric_limits<size_t>::max() ? moves.size() : moves.size() - best_move_idx;
                            // revert to best state
                            /*
                            for (size_t idx = 0; idx < n_reverts; ++idx) {
                                vertex_t v = moves[moves.size() - 1 - idx];
                                partition_t v_old_id = (*m_p_p_manger)[v];
                                partition_t v_new_id = (*m_p_p_manger)[v] == b1_id ? b2_id : b1_id;

                                ASSERT(v_old_id == b1_id || v_old_id == b2_id);
                                ASSERT(v_new_id == b1_id || v_new_id == b2_id);

                                // revert the move
                                m_p_p_manger->move(v, v_old_id, v_new_id);
                                // m_p_bv_manager->move(v, v_old_id, v_new_id);
                            }
                             */

                            // make boundary vertices move to best vertex
                            for(size_t idx = 0; idx < moves.size(); ++idx){
                                vertex_t v = moves[moves.size() - 1 - idx];
                                partition_t v_old_id = (*m_p_p_manger)[v];
                                partition_t v_new_id = (*m_p_p_manger)[v] == b1_id ? b2_id : b1_id;

                                ASSERT(v_old_id == b1_id || v_old_id == b2_id);
                                ASSERT(v_new_id == b1_id || v_new_id == b2_id);

                                // revert the move
                                m_p_p_manger->move(v, v_old_id, v_new_id);
                                // m_p_bv_manager->move(v, v_old_id, v_new_id);
                            }
                            for(size_t idx = 0; idx < moves.size() - n_reverts; ++idx){
                                vertex_t v = moves[idx];
                                partition_t v_old_id = (*m_p_p_manger)[v];
                                partition_t v_new_id = (*m_p_p_manger)[v] == b1_id ? b2_id : b1_id;

                                ASSERT(v_old_id == b1_id || v_old_id == b2_id);
                                ASSERT(v_new_id == b1_id || v_new_id == b2_id);

                                // make the move
                                m_p_p_manger->move(v, v_old_id, v_new_id);
                                m_p_bv_manager->move(v, v_old_id, v_new_id);
                                m_p_qgraph->move(v, v_old_id, v_new_id);
                            }
                        }
                    }
                    // break;
                }
            }
        }
    };
}

#endif //HEIPROMAP_PARALLEL_QUOTIENT_GRAPH_REFINEMENT_H
