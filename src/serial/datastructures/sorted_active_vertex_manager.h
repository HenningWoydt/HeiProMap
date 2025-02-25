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

#ifndef HEIPROMAP_SORTED_ACTIVE_VERTEX_MANAGER_H
#define HEIPROMAP_SORTED_ACTIVE_VERTEX_MANAGER_H

#include <numeric>
#include <algorithm>

#include "../../definitions.h"
#include "../interfaces/ISerialActiveVertexManager.h"

namespace HeiProMap {
    class SortedActiveVertexManager final : public ISerialActiveVertexManager {
        u8       *m_states       = nullptr;
        vertex_t *m_vertices     = nullptr;
        size_t   m_vertices_size = 0;
        vertex_t m_n_active      = 0;

        vertex_t            *m_vertices_temp = nullptr;
        std::vector<size_t> temp_points;

    public:
        SortedActiveVertexManager() = default;

        ~SortedActiveVertexManager() override {
            free(m_states);
            free(m_vertices);

            free(m_vertices_temp);
        }

        // initialize
        void initialize(const size_t t_n) override {
            size_t t_n_64 = round_up_64(t_n);

            m_states = (u8 *) aligned_alloc(64, t_n_64 * sizeof(u8));
            std::fill(m_states, m_states + t_n, (u8) 1);

            m_vertices      = (vertex_t *) aligned_alloc(64, t_n_64 * sizeof(vertex_t));
            m_vertices_size = t_n;
            std::iota(m_vertices, m_vertices + t_n, 0);

            m_n_active = t_n;

            m_vertices_temp = (vertex_t *) aligned_alloc(64, t_n_64 * sizeof(vertex_t));
            temp_points.push_back(t_n);
        }

        // active vertex manipulation
        vertex_t get_n_active() const override { return m_n_active; }

        bool is_active(const vertex_t u) const override { return m_states[u] == 1; }

        bool is_disabled(const vertex_t u) const override { return m_states[u] == 0; }

        bool get_state(const vertex_t u) const override { return m_states[u]; }

        void contract(const EdgeUV *matches, size_t &matches_size) override {
            matches = ASSUME_ALIGNED(EdgeUV*, matches, 64);

            size_t temp_idx = 0;
            m_n_active -= matches_size;
            for (size_t i = 0; i < matches_size; ++i) { m_states[matches[i].v] = 0; }

            size_t      write_idx = 0;
            for (size_t read_idx  = 0; read_idx < m_vertices_size; ++read_idx) {
                if (!is_disabled(m_vertices[read_idx])) {
                    m_vertices[write_idx] = m_vertices[read_idx];
                    ++write_idx;
                } else {
                    m_vertices_temp[temp_idx] = m_vertices[read_idx];
                    ++temp_idx;
                }
            }
            m_vertices_size = write_idx;

            std::copy(m_vertices_temp, m_vertices_temp + temp_idx, m_vertices + write_idx);
            temp_points.push_back(write_idx);
        }

        void uncontract(const EdgeUV *matches, size_t &matches_size) override {
            matches = ASSUME_ALIGNED(EdgeUV*, matches, 64);

            m_n_active += matches_size;
            m_vertices_size = temp_points[temp_points.size() - 2];
            for (size_t i   = 0; i < matches_size; ++i) { m_states[matches[i].v] = 1; }
            std::inplace_merge(m_vertices, m_vertices + temp_points[temp_points.size() - 1], m_vertices + temp_points[temp_points.size() - 2]); // merge
            temp_points.pop_back();
        }

        class Iterator {
            vertex_t *m_vertices = nullptr;
            size_t   m_idx;

        public:
            // Constructor
            explicit Iterator(vertex_t *vertices, size_t idx) {
                m_vertices = ASSUME_ALIGNED(vertex_t *, vertices, 64);
                m_idx      = idx;
            }

            // Dereference operator
            vertex_t operator*() const {
                return m_vertices[m_idx];
            }

            // Pre-increment operator
            Iterator &operator++() {
                m_idx++;
                return *this;
            }

            bool operator!=(const Iterator &other) const { return m_idx != other.m_idx; }
        };

        Iterator begin() { return Iterator(m_vertices, 0); }

        Iterator end() { return Iterator(m_vertices, m_vertices_size); }
    };
}

#endif //HEIPROMAP_SORTED_ACTIVE_VERTEX_MANAGER_H
