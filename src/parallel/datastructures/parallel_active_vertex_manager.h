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

#ifndef HEIPROMAP_PARALLEL_ACTIVE_VERTEX_MANAGER_H
#define HEIPROMAP_PARALLEL_ACTIVE_VERTEX_MANAGER_H

#include "../../definitions.h"
#include "../../serial/utility/utils.h"
#include "../../macros.h"
#include "../../interfaces/IGraph.h"
#include "../../interfaces/IActiveVertexManager.h"
#include "../interfaces/IParallelActiveVertexManager.h"

namespace HeiProMap {

    template<typename TParallelGraph>
    class ParallelActiveVertexManager : public IParallelActiveVertexManager<TParallelGraph> {
    private:
        TParallelGraph *m_p_g = nullptr;

        u64 m_n_threads = 1;

        // active states
        std::vector<bool> m_states;
        std::vector<vertex_t> m_vertices;
        std::atomic<vertex_t> m_n_active = 0;

        // iterator
        size_t idx = 0;

    public:
        // initialize
        void initialize(TParallelGraph *t_p_g,
                        u64 n_threads) final {
            ASSERT(t_p_g != nullptr);

            m_p_g = t_p_g;

            m_n_threads = n_threads;

            m_states.resize(m_p_g->get_n(), true);
            m_vertices.resize(m_p_g->get_n());
            std::iota(m_vertices.begin(), m_vertices.end(), 0);
            m_n_active = m_p_g->get_n();
            idx = 0;
        }

        // active vertex manipulation
        vertex_t get_n_active() const final { return m_n_active; }

        void activate_vertex(vertex_t u) final {
            if (!m_states[u]) {
                m_states[u] = true;
                m_vertices.push_back(u);
                m_n_active += 1;
            }
        }

        void disable_vertex(vertex_t u) final {
            if (m_states[u]) {
                m_states[u] = false;
                m_n_active -= 1;
            }
        }

        bool is_active(vertex_t u) const final { return m_states[u]; }

        bool is_disabled(vertex_t u) const final { return !m_states[u]; }

        bool get_state(vertex_t u) const final { return m_states[u]; }

        // coarsing and uncoarsing
        void contract([[maybe_unused]] vertex_t u, vertex_t v) final { disable_vertex(v); }

        void uncontract([[maybe_unused]] vertex_t u, vertex_t v) final { activate_vertex(v); }

        // iteration
        void reset_iterator() final { idx = 0; }

        vertex_t get() final { return m_vertices[idx]; }

        void next() final { idx += 1; }

        bool available() final {
            while (idx < m_vertices.size() && is_disabled(m_vertices[idx])) {
                m_vertices[idx] = m_vertices.back();
                m_vertices.pop_back();
            }

            return idx < m_vertices.size();
        }
    };

}

#endif //HEIPROMAP_PARALLEL_ACTIVE_VERTEX_MANAGER_H
