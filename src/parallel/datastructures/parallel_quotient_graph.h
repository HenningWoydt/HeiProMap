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

#ifndef HEIPROMAP_PARALLEL_QUOTIENT_GRAPH_H
#define HEIPROMAP_PARALLEL_QUOTIENT_GRAPH_H

#include "../parallel_definitions_1.h"
#include "../interfaces/IParallelQuotientGraph.h"

namespace HeiProMap {
    class ParallelQuotientGraph final : public IParallelQuotientGraph {
        partition_t k = 0;

        std::vector<weight_t> m_adj_mtx;

    public:
        ParallelQuotientGraph() = default;

        void initialize(partition_t t_k) override {
            k = t_k;

            m_adj_mtx.resize(k * k, 0);
        }

        void add_edge(partition_t u, partition_t v, weight_t w) override {
            partition_t min = std::min(u, v);
            partition_t max = std::max(u, v);
            m_adj_mtx[min * k + max] += w;
        }

        void remove_edge(partition_t u, partition_t v, weight_t w) override {
            partition_t min = std::min(u, v);
            partition_t max = std::max(u, v);
            m_adj_mtx[min * k + max] -= w;
        }

        bool has_edge(partition_t u, partition_t v) override {
            partition_t min = std::min(u, v);
            partition_t max = std::max(u, v);
            return m_adj_mtx[min * k + max] > 0;
        }

        weight_t get_weight(partition_t u, partition_t v) override {
            partition_t min = std::min(u, v);
            partition_t max = std::max(u, v);
            return m_adj_mtx[min * k + max];
        }

        void move(p_graph_t &g, p_p_manager_t &p_manager, vertex_t u, partition_t old_id, partition_t new_id) {
            ASSERT(new_id < k);
            ASSERT(new_id != old_id);

            for (size_t i = 0; i < g.size(u); ++i) {
                vertex_t    v    = g.neighbor(u, i);
                weight_t    w    = g.get_weight(u, i);
                partition_t v_id = p_manager[v];

                // remove old edge, if existed
                if (old_id != v_id) {
                    remove_edge(old_id, v_id, w * 2);
                }
                // add new edge, if has to exist
                if (new_id != v_id) {
                    add_edge(new_id, v_id, w * 2);
                }
            }
        }
    };
}


#endif //HEIPROMAP_PARALLEL_QUOTIENT_GRAPH_H
