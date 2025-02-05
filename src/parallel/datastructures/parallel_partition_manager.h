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

#ifndef HEIPROMAP_PARALLEL_PARTITION_MANAGER_H
#define HEIPROMAP_PARALLEL_PARTITION_MANAGER_H

#include <atomic>

#include "../../definitions.h"
#include "../../macros.h"
#include "../../serial/utility/utils.h"
#include "../../serial/utility/qap.h"
#include "../../interfaces/IPartitionManager.h"
#include "../interfaces/IParallelActiveVertexManager.h"
#include "../interfaces/IParallelPartitionManager.h"

namespace HeiProMap {

    template<typename TParallelGraph, typename TParallelActiveVertexManager>
    class ParallelPartitionManager : public IParallelPartitionManager<TParallelGraph, TParallelActiveVertexManager> {
    private:
        TParallelGraph *m_p_g = nullptr;
        TParallelActiveVertexManager *m_p_av_manager = nullptr;
        partition_t m_k = 0;

        u64 m_n_threads = 1;

        // actual partition
        std::vector<partition_t> partition;

        // partition weights
        std::atomic<weight_t> *m_p_bweights = nullptr;

    public:

        ~ParallelPartitionManager() { free(m_p_bweights); }

        void initialize(TParallelGraph *t_p_g,
                        TParallelActiveVertexManager *t_p_av_manager,
                        partition_t t_k,
                        u64 t_n_threads) final {
            ASSERT(t_p_g != nullptr);
            ASSERT(t_p_av_manager != nullptr);

            m_p_g = t_p_g;
            m_p_av_manager = t_p_av_manager;
            m_k = t_k;

            m_n_threads = t_n_threads;

            // actual partition
            partition.resize(m_p_g->get_n());

            // partition weights
            m_p_bweights = (std::atomic<weight_t> *) malloc(m_k * sizeof(std::atomic<weight_t>));
            for(size_t i = 0; i < m_k; ++i){ m_p_bweights[i] = 0; }
        }

        // read
        const partition_t &operator[](vertex_t u) const final { return partition[u]; }

        // write
        void set(vertex_t u, partition_t id) final {
            ASSERT(m_p_g != nullptr);
            m_p_bweights[id] += m_p_g->get_weight(u);
            partition[u] = id;
        }

        bool is_boundary(vertex_t u) final {
            partition_t u_id = partition[u];
            for(size_t i = 0; i < m_p_g->size(u); ++i){
                if(u_id != partition[m_p_g->neighbor(u, i)]){
                    return true;
                }
            }
            return false;
        }

        void move(vertex_t u, partition_t old_id, partition_t new_id) final {
            ASSERT(m_p_g != nullptr);
            m_p_bweights[old_id] -= m_p_g->get_weight(u);
            m_p_bweights[new_id] += m_p_g->get_weight(u);
            partition[u] = new_id;
        }

        // weights
        weight_t get_bweight(partition_t id) const final { return m_p_bweights[id]; }

        std::vector<weight_t> get_bweights() const final {
            std::vector<weight_t> temp(m_k);
            for (size_t i = 0; i < m_k; ++i) { temp[i] = m_p_bweights[i]; }
            return temp;
        }

        // uncoarsing
        void uncontract(vertex_t u, vertex_t v) final { partition[v] = partition[u]; }
    };

}

#endif //HEIPROMAP_PARALLEL_PARTITION_MANAGER_H
