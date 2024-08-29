#ifndef HEIDELBERGPROCESSMAPPING_PARALLEL_QUOTIENT_GRAPH_REFINEMENT_H
#define HEIDELBERGPROCESSMAPPING_PARALLEL_QUOTIENT_GRAPH_REFINEMENT_H

#include <queue>
#include <omp.h>
#include <atomic>

#include "../interfaces/IParallelRefiner.h"
#include "../../serial/datastructures/indexed_max_heap.h"

namespace HeiProMap {

    class ParallelQuotientGraphRefinement : public IParallelRefiner {
    private:
        IParallelGraph *m_p_g = nullptr;
        IParallelActiveVertexManager *m_p_av_manager = nullptr;
        IParallelBoundaryVertexManager *m_p_bv_manager = nullptr;
        IParallelPartitionManager *m_p_p_manger = nullptr;
        IParallelDistanceOracle *m_p_d_oracle = nullptr;
        IParallelQuotientGraph *m_p_qgraph = nullptr;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;
        weight_t m_lmax = 0;
        u64 m_n_threads = 1;

        std::vector<s32> m_used;
        std::vector<u8> m_neighbor_changed;
        s32 m_mark = -1;

    public:
        ParallelQuotientGraphRefinement() = default;

        void initialize(IParallelGraph *t_p_g,
                        IParallelActiveVertexManager *t_p_av_manager,
                        IParallelBoundaryVertexManager *t_p_bv_manager,
                        IParallelPartitionManager *t_p_p_manger,
                        IParallelDistanceOracle *t_p_d_oracle,
                        IParallelQuotientGraph *t_p_qgraph,
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

            std::vector<std::vector<QGraphUV>> matching_hierarchy;
            m_p_qgraph->get_matching_hierarchy(matching_hierarchy);

            size_t max_global_iterations = 10;
            for (size_t global_iteration = 0; global_iteration < max_global_iterations; ++global_iteration) {
                for (std::vector<QGraphUV> &matchings: matching_hierarchy) {
#pragma omp parallel for default(none) shared(matchings, std::cerr) num_threads(m_n_threads)
                    for (QGraphUV uv: matchings) {
                        partition_t b1_id = uv.u;
                        partition_t b2_id = uv.v;

                        // do fm refinement on blocks b1_id and b2_id
                        std::vector<u8> vertex_moved(m_p_g->get_n());
                        std::vector<vertex_t> moves;
                        IndexedMaxHeap<Swap> max_heap(m_p_g->get_n());

                        size_t max_iterations = 1;
                        for (size_t iteration = 0; iteration < max_iterations; ++iteration) {
                            std::fill(vertex_moved.begin(), vertex_moved.end(), 0);
                            moves.clear();
                            max_heap.clear();

                            // initialize all swaps and their deltas from b1
                            for (m_p_bv_manager->reset_iterator(b1_id); m_p_bv_manager->available(b1_id); m_p_bv_manager->next(b1_id)) {
                                vertex_t u = m_p_bv_manager->get(b1_id);
                                ASSERT((*m_p_p_manger)[u] == b1_id);

                                s64 qap_delta = get_u_qap(*m_p_g, u, *m_p_p_manger, *m_p_d_oracle) - get_u_qap(*m_p_g, u, b2_id, *m_p_p_manger, *m_p_d_oracle);
                                Swap s(u, qap_delta);
                                max_heap.push(u, s);
                            }

                            // initialize all swaps and their deltas from b2
                            for (m_p_bv_manager->reset_iterator(b2_id); m_p_bv_manager->available(b2_id); m_p_bv_manager->next(b2_id)) {
                                vertex_t v = m_p_bv_manager->get(b2_id);
                                ASSERT((*m_p_p_manger)[v] == b2_id);

                                s64 qap_delta = get_u_qap(*m_p_g, v, *m_p_p_manger, *m_p_d_oracle) - get_u_qap(*m_p_g, v, b1_id, *m_p_p_manger, *m_p_d_oracle);
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

                                if (!m_p_bv_manager->is_boundary(u)) {
                                    // vertex could not be boundary anymore
                                    continue;
                                }

                                if(qap_delta < 0){
                                    n_neg_deltas += 1;
                                    if(n_neg_deltas >= 10){
                                        break;
                                    }
                                } else{
                                    n_neg_deltas = 0;
                                }

                                // make the move
                                m_p_p_manger->move(u, u_old_id, u_new_id);
                                m_p_bv_manager->move(u, u_old_id, u_new_id);

                                moves.push_back(u);
                                curr_cum_qap_delta += qap_delta;
                                if (curr_cum_qap_delta > best_cum_qap_delta && m_p_p_manger->get_bweight(b1_id) <= m_lmax && m_p_p_manger->get_bweight(b2_id) <= m_lmax) {
                                    // new best and not overloaded
                                    best_cum_qap_delta = curr_cum_qap_delta;
                                    best_move_idx = moves.size();
                                }
                                vertex_moved[u] = 1;

                                // update the neighbors of u if they are in b1 or b2, are boundary vertices and have not been moved
                                for (size_t i = 0; i < m_p_g->size(u); ++i) {
                                    vertex_t v = m_p_g->neighbor(u, i);
                                    partition_t v_old_id = (*m_p_p_manger)[v];
                                    if (vertex_moved[v] == 1 || (v_old_id != b1_id && v_old_id != b2_id) || !m_p_bv_manager->is_boundary(v)) {
                                        // vertex cannot be moved
                                        continue;
                                    }

                                    // update the qap delta
                                    partition_t v_new_id = (*m_p_p_manger)[v] == b1_id ? b2_id : b1_id;
                                    s64 v_qap_delta = get_u_qap(*m_p_g, v, *m_p_p_manger, *m_p_d_oracle) - get_u_qap(*m_p_g, v, v_new_id, *m_p_p_manger, *m_p_d_oracle);

                                    // insert/update the qap delta
                                    max_heap.push_update(v, {v, v_qap_delta});
                                }
                            }

                            // revert to best state
                            size_t n_reverts = best_move_idx == std::numeric_limits<size_t>::max() ? moves.size() : moves.size() - best_move_idx;
                            for (size_t idx = 0; idx < n_reverts; ++idx) {
                                vertex_t v = moves[moves.size() - 1 - idx];
                                partition_t v_old_id = (*m_p_p_manger)[v];
                                partition_t v_new_id = (*m_p_p_manger)[v] == b1_id ? b2_id : b1_id;

                                ASSERT(v_old_id == b1_id || v_old_id == b2_id);
                                ASSERT(v_new_id == b1_id || v_new_id == b2_id);

                                // revert the move
                                m_p_p_manger->move(v, v_old_id, v_new_id);
                                m_p_bv_manager->move(v, v_old_id, v_new_id);
                            }
                        }
                    }
                }
            }
        }
    };
}

#endif //HEIDELBERGPROCESSMAPPING_PARALLEL_QUOTIENT_GRAPH_REFINEMENT_H
