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

#include <algorithm>
#include <numeric>

#include "../../definitions.h"
#include "../../macros.h"
#include "../../commons/utils.h"
#include "../interfaces/IParallelActiveVertexManager.h"

namespace HeiProMap {
    class ParallelActiveVertexManager final : public IParallelActiveVertexManager {
        u8* m_states           = nullptr;
        vertex_t* m_vertices   = nullptr;
        size_t m_vertices_size = 0;
        vertex_t m_n_active    = 0;

        vertex_t* m_vertices_temp = nullptr;
        std::vector<size_t> temp_points;

    public:
        ~ParallelActiveVertexManager() override {
            free(m_states);
            free(m_vertices);
            free(m_vertices_temp);
        }

        void initialize(const size_t t_n) override {
            const size_t t_n_64 = round_up_64(t_n);

            m_states = (u8*)aligned_alloc(64, t_n_64 * sizeof(u8));
            std::fill_n(m_states, t_n, (u8)1);

            m_vertices      = (vertex_t*)aligned_alloc(64, t_n_64 * sizeof(vertex_t));
            m_vertices_size = t_n;
            std::iota(m_vertices, m_vertices + t_n, 0);

            m_n_active = t_n;

            m_vertices_temp = (vertex_t*)aligned_alloc(64, t_n_64 * sizeof(vertex_t));
            temp_points.push_back(t_n);
        }

        size_t size() const override { return m_n_active; }
        vertex_t get(const size_t i) const override { return m_vertices[i]; }
        bool is_active(const vertex_t u) const override { return m_states[u] == 1; }
        bool is_disabled(const vertex_t u) const override { return m_states[u] == 0; }
        bool get_state(const vertex_t u) const override { return m_states[u]; }

        void contract(const EdgeUV* matches, const size_t& matches_size) override {
            matches = ASSUME_ALIGNED(EdgeUV*, matches, 64);

            size_t temp_idx = 0;
            m_n_active -= matches_size;
            for (size_t i = 0; i < matches_size; ++i) { m_states[matches[i].v] = 0; }

            size_t write_idx = 0;
            for (size_t read_idx = 0; read_idx < m_vertices_size; ++read_idx) {
                if (!is_disabled(m_vertices[read_idx])) {
                    m_vertices[write_idx] = m_vertices[read_idx];
                    ++write_idx;
                } else {
                    m_vertices_temp[temp_idx] = m_vertices[read_idx];
                    ++temp_idx;
                }
            }
            m_vertices_size = write_idx;

            std::copy_n(m_vertices_temp, temp_idx, m_vertices + write_idx);
            temp_points.push_back(write_idx);
        }

        void uncontract(const EdgeUV* matches, const size_t& matches_size) override {
            matches = ASSUME_ALIGNED(EdgeUV*, matches, 64);

            m_n_active += matches_size;
            m_vertices_size = temp_points[temp_points.size() - 2];
            for (size_t i = 0; i < matches_size; ++i) { m_states[matches[i].v] = 1; }
            std::inplace_merge(m_vertices, m_vertices + temp_points[temp_points.size() - 1], m_vertices + temp_points[temp_points.size() - 2]); // merge
            temp_points.pop_back();
        }
    };
}

#endif //HEIPROMAP_PARALLEL_ACTIVE_VERTEX_MANAGER_H
