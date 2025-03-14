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

#ifndef HEIPROMAP_QUOTIENT_GRAPH_H
#define HEIPROMAP_QUOTIENT_GRAPH_H

#include "../../commons/definitions.h"
#include "../interfaces/ISerialQuotientGraph.h"

namespace HeiProMap {
    class QuotientGraph final : public ISerialQuotientGraph {
        partition_t m_k = 0;

        weight_t* m_adj_mtx = nullptr;

    public:
        QuotientGraph() = default;
        ~QuotientGraph() override { free(m_adj_mtx); }

        void initialize(const partition_t t_k) override {
            partition_t k_k_64 = round_up_64(t_k * t_k);
            m_k                = t_k;

            m_adj_mtx = (weight_t*)aligned_alloc(64, k_k_64 * sizeof(weight_t));
            std::fill_n(m_adj_mtx, k_k_64, 0);
        }

        void add_edge(const partition_t u_id, const partition_t v_id, const weight_t w) override {
            partition_t min = std::min(u_id, v_id);
            partition_t max = std::max(u_id, v_id);
            m_adj_mtx[min * m_k + max] += w;
        }

        void remove_edge(const partition_t u_id, const partition_t v_id, const weight_t w) override {
            partition_t min = std::min(u_id, v_id);
            partition_t max = std::max(u_id, v_id);
            m_adj_mtx[min * m_k + max] -= w;
        }

        bool has_edge(const partition_t u_id, const partition_t v_id) const override {
            partition_t min = std::min(u_id, v_id);
            partition_t max = std::max(u_id, v_id);
            return m_adj_mtx[min * m_k + max] > 0;
        }

        weight_t get_weight(const partition_t u_id, const partition_t v_id) const override {
            partition_t min = std::min(u_id, v_id);
            partition_t max = std::max(u_id, v_id);
            return m_adj_mtx[min * m_k + max];
        }

        void move(const graph_t& g,
                  const p_manager_t& p_manager,
                  const vertex_t u,
                  const partition_t old_id,
                  const partition_t new_id) override {
            ASSERT(new_id < m_k);
            ASSERT(new_id != old_id);

            forall_guivw(g, u, i, v, w)
                {
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
            endfor
        }
    };
}

#endif //HEIPROMAP_QUOTIENT_GRAPH_H
