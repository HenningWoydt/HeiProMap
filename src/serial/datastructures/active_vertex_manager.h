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
        std::vector<u8>       m_states;
        std::vector<vertex_t> m_vertices;
        vertex_t              m_n_active = 0;

    public:
        // initialize
        void initialize(const size_t n) override {
            m_states.resize(n, 1);
            m_vertices.resize(n);
            std::iota(m_vertices.begin(), m_vertices.end(), 0);
            m_n_active = n;
        }

        // active vertex manipulation
        vertex_t get_n_active() const override { return m_n_active; }

        bool is_active(const vertex_t u) const override { return m_states[u] == 1; }

        bool is_disabled(const vertex_t u) const override { return m_states[u] == 0; }

        bool get_state(const vertex_t u) const override { return m_states[u]; }

        void contract(const EdgeUV *matches, size_t &matches_size) override {
            matches = ASSUME_ALIGNED(EdgeUV*, matches, 64);

            m_n_active -= matches_size;
            for (size_t i = 0; i < matches_size; ++i) { m_states[matches[i].v] = 0; }
        }

        void uncontract(const EdgeUV *matches, size_t &matches_size) override {
            matches = ASSUME_ALIGNED(EdgeUV*, matches, 64);

            m_n_active += matches_size;
            for (size_t i = 0; i < matches_size; ++i) {
                vertex_t v = matches[i].v;
                m_states[v] = 1;
                m_vertices.push_back(v);
            }
        }

        class Iterator {
            std::vector<vertex_t> &m_vertices;
            std::vector<u8>       &m_states;
            size_t                m_idx = 0;

        public:
            // Constructor
            Iterator(std::vector<vertex_t> &vertices,
                     std::vector<u8> &states) : m_vertices(vertices), m_states(states) {
                // Skip disabled vertices during initialization
                advance_to_next_valid();
            }

            // Dereference operator
            vertex_t operator*() const {
                return m_vertices[m_idx];
            }

            // Pre-increment operator
            Iterator &operator++() {
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

            bool operator!=(const Iterator &other) const {
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
