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

#ifndef HEIPROMAP_ACTIVE_VERTEX_MANAGER_H
#define HEIPROMAP_ACTIVE_VERTEX_MANAGER_H

#include <numeric>

#include "../../definitions.h"
#include "../interfaces/ISerialActiveVertexManager.h"

namespace HeiProMap {
    class ActiveVertexManager final : public ISerialActiveVertexManager {
        // active states
        std::vector<bool> m_states;
        std::vector<vertex_t> m_vertices;
        vertex_t m_n_active = 0;

    public:
        // initialize
        void initialize(const size_t n) override {
            m_states.resize(n, true);
            m_vertices.resize(n);
            std::iota(m_vertices.begin(), m_vertices.end(), 0);
            m_n_active = n;
        }

        // active vertex manipulation
        vertex_t get_n_active() const override { return m_n_active; }

        void activate_vertex(const vertex_t u) override {
            if (!m_states[u]) {
                m_states[u] = true;
                m_vertices.push_back(u);
                m_n_active += 1;
            }
        }

        void disable_vertex(const vertex_t u) override {
            if (m_states[u]) {
                m_states[u] = false;
                m_n_active -= 1;
            }
        }

        bool is_active(const vertex_t u) const override { return m_states[u]; }
        bool is_disabled(const vertex_t u) const override { return !m_states[u]; }
        bool get_state(const vertex_t u) const override { return m_states[u]; }

        void contract(const std::vector<EdgeUV>& matches) override {
            for (const auto [u, v] : matches) {
                disable_vertex(v);
            }
        }

        void uncontract(const std::vector<EdgeUV>& matches) override {
            for (const auto [u, v] : matches) {
                activate_vertex(v);
            }
        }

        // iteration
        // void reset_iterator() override { idx = 0; }
        // vertex_t get() const override { return m_vertices[idx]; }
        // void next() override { idx += 1; }

        /*
        bool available() override {
            while (idx < m_vertices.size() && is_disabled(m_vertices[idx])) {
                m_vertices[idx] = m_vertices.back();
                m_vertices.pop_back();
            }

            return idx < m_vertices.size();
        }
        */

        class Iterator {
            std::vector<vertex_t>& m_vertices;
            std::vector<bool>& m_states;
            size_t m_idx = 0;

        public:
            // Constructor
            Iterator(std::vector<vertex_t>& vertices,
                     std::vector<bool>& states) : m_vertices(vertices), m_states(states) {
                // Skip disabled vertices during initialization
                advance_to_next_valid();
            }

            // Dereference operator
            vertex_t operator*() const {
                return m_vertices[m_idx];
            }

            // Pre-increment operator
            Iterator& operator++() {
                // Remove the current inactive vertex from the vector
                if (!m_states[m_vertices[m_idx]]) {
                    m_vertices[m_idx] = m_vertices.back();
                    m_vertices.pop_back();
                } else {
                    ++m_idx;
                }

                // Advance to the next valid vertex
                advance_to_next_valid();
                return *this;
            }

            bool operator!=(const Iterator& other) const {
                return m_idx != other.m_vertices.size();
            }

        private:
            // Helper to advance the iterator to the next active vertex
            void advance_to_next_valid() const {
                while (m_idx < m_vertices.size() && !m_states[m_vertices[m_idx]]) {
                    // Pop inactive vertices out of the container
                    m_vertices[m_idx] = m_vertices.back();
                    m_vertices.pop_back();
                }
            }
        };

        Iterator begin() { return {m_vertices, m_states}; }
        Iterator end() { return {m_vertices, m_states}; }
    };
}

#endif //HEIPROMAP_ACTIVE_VERTEX_MANAGER_H
