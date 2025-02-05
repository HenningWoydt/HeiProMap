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

#ifndef HEIPROMAP_PARALLEL_KAFFPA_PARTITIONER_H
#define HEIPROMAP_PARALLEL_KAFFPA_PARTITIONER_H

#include "../interfaces/IParallelPartitioner.h"
#include "../../serial/datastructures/translation_table.h"

#include "interface/kaHIP_interface.h"

namespace HeiProMap {

    template<typename TParallelGraph, typename TParallelActiveVertexManager, typename TParallelBoundaryVertexManager, typename TParallelPartitionManager>
    class ParallelKaffpaPartitioner : public IParallelPartitioner<TParallelGraph, TParallelActiveVertexManager, TParallelBoundaryVertexManager, TParallelPartitionManager> {
    private:
        TParallelGraph *m_p_g = nullptr;
        TParallelActiveVertexManager *m_p_av_manager = nullptr;
        TParallelBoundaryVertexManager *m_p_bv_manager = nullptr;
        TParallelPartitionManager *m_p_p_manager = nullptr;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;
        f64 m_imbalance = 0;

        u64 m_n_threads = 1;

    public:
        // initialization
        void initialize(TParallelGraph *t_p_g,
                        TParallelActiveVertexManager *t_p_av_manager,
                        TParallelBoundaryVertexManager *t_p_bv_manager,
                        TParallelPartitionManager *t_p_p_manager,
                        std::vector<partition_t> &t_hierarchy,
                        std::vector<weight_t> &t_distance,
                        f64 t_imbalance,
                        u64 n_threads) final {
            ASSERT(t_p_g != nullptr);
            ASSERT(t_p_av_manager != nullptr);
            ASSERT(t_p_bv_manager != nullptr);
            ASSERT(t_p_p_manager != nullptr);
            ASSERT(t_imbalance >= 0.0);

            m_p_g = t_p_g;
            m_p_av_manager = t_p_av_manager;
            m_p_bv_manager = t_p_bv_manager;
            m_p_p_manager = t_p_p_manager;
            m_hierarchy = t_hierarchy;
            m_distance = t_distance;
            m_imbalance = t_imbalance;

            m_n_threads = n_threads;
        }

        // partition
        void partition() final {
            ASSERT(m_p_g != nullptr);
            ASSERT(m_p_av_manager != nullptr);
            ASSERT(m_p_bv_manager != nullptr);
            ASSERT(m_p_p_manager != nullptr);

            // number of vertices and edges
            int n = 0;
            int m = 0;

            // build translation table
            TranslationTable tt;
            vertex_t translate = 0;
            for (m_p_av_manager->reset_iterator(); m_p_av_manager->available(); m_p_av_manager->next()) {
                vertex_t u = m_p_av_manager->get();
                ASSERT(m_p_av_manager->is_active(u));

                tt.add(u, translate);
                translate += 1;

                n += 1;
                m += (int) m_p_g->size(u);
            }

            // vertex weights
            int *v_weights = (int *) malloc(n * sizeof(int));
            for (int i = 0; i < n; ++i) { v_weights[i] = (int) m_p_g->get_weight(tt.get_o(i)); }

            // pointer to adjacency lists
            int *adj_ptr = (int *) malloc((n + 1) * sizeof(int));
            int *adj = (int *) malloc(m * sizeof(int));
            int *e_weights = (int *) malloc(m * sizeof(int));

            // set adj_ptr
            adj_ptr[0] = 0;
            for (int new_u = 0; new_u < n; ++new_u) {
                vertex_t old_u = tt.get_o(new_u);
                int insert_idx = 0;
                for (size_t i = 0; i < m_p_g->size(old_u); ++i) {
                    vertex_t v = m_p_g->neighbor(old_u, i);
                    weight_t ew = m_p_g->get_weight(old_u, i);
                    adj[adj_ptr[new_u] + insert_idx] = (int) tt.get_n(v);
                    e_weights[adj_ptr[new_u] + insert_idx] = (int) ew;
                    insert_idx += 1;
                }
                adj_ptr[new_u + 1] = adj_ptr[new_u] + insert_idx;
            }
            // imbalance
            double kaffpa_imbalance = m_imbalance;

            // hierarchy
            int *kaffpa_hierarchy = (int *) malloc(m_hierarchy.size() * sizeof(int));
            for (u64 i = 0; i < m_hierarchy.size(); ++i) { kaffpa_hierarchy[i] = (int) m_hierarchy[i]; }

            // distance
            int *kaffpa_distance = (int *) malloc(m_distance.size() * sizeof(int));
            for (u64 i = 0; i < m_distance.size(); ++i) { kaffpa_distance[i] = (int) m_distance[i]; }

            // mode
            int kaffpa_map_mode = MAPMODE_BISECTION; // TODO: Figure out why MAPMODE_MULTISECTION does not work
            int kaffpa_partition_mode = KAFFPA_FAST;

            // partition result
            int *kaffpa_partition = (int *) malloc(n * sizeof(int));

            int kaffpa_edgecut, kaffpa_qap;

            // execute kaffpa
            process_mapping(&n, v_weights, adj_ptr, e_weights, adj, kaffpa_hierarchy, kaffpa_distance, (int) m_hierarchy.size(), kaffpa_partition_mode, kaffpa_map_mode, &kaffpa_imbalance, false, 0, &kaffpa_edgecut, &kaffpa_qap, kaffpa_partition);

            // first read partition
            for (int new_u = 0; new_u < n; ++new_u) {
                m_p_p_manager->set(tt.get_o(new_u), kaffpa_partition[new_u]);
            }
            // then initialize boundary vertices
            for (int new_u = 0; new_u < n; ++new_u) {
                m_p_bv_manager->insert(tt.get_o(new_u), kaffpa_partition[new_u]);
            }

            free(v_weights);
            free(adj_ptr);
            free(adj);
            free(e_weights);
            free(kaffpa_hierarchy);
            free(kaffpa_distance);
            free(kaffpa_partition);
        }
    };


}

#endif //HEIPROMAP_PARALLEL_KAFFPA_PARTITIONER_H
