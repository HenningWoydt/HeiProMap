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

#ifndef HEIPROMAP_MATCHING_H
#define HEIPROMAP_MATCHING_H

#include <atomic>
#include <numeric>

#include "definitions.h"
#include "translation_table.h"
#include "utils.h"

namespace HeiProMap {
    class Matching {
    private:
        vertex_t m_n = 0;

        // AlignedArray<EdgeUV> matches;           // O(n)
        std::atomic<size_t> matches_size = 0;

        AlignedArray<vertex_t> partner;         // O(n)

        AlignedArray<vertex_t> o_to_n;          // O(n)

    public:
        Matching() = default;

        void initialize(vertex_t n) {
            vertex_t n_64 = round_up_64(n);
            m_n = n;

            // matches.initialize((n_64 / 2));
            matches_size = 0;

            partner.initialize(m_n);
            std::iota(partner.get_ptr(), partner.get_ptr() + n_64, 0);

            o_to_n.initialize(m_n);
            std::iota(o_to_n.get_ptr(), o_to_n.get_ptr() + n_64, 0);
        }

        // Move constructor
        Matching(Matching &&other) noexcept {
            m_n = other.m_n;
            // std::swap(matches, other.matches);
            size_t temp1 = matches_size;
            size_t temp2 = other.matches_size;
            matches_size = temp2;
            other.matches_size = temp1;
            std::swap(partner, other.partner);
            std::swap(o_to_n, other.o_to_n);
        }

        // Optionally disable copying.
        Matching(const Matching &) = delete;

        Matching &operator=(const Matching &) = delete;

        ~Matching() = default;

        void add(vertex_t u, vertex_t v) {
            // ASSERT(matches_size < (m_n / 2));
            ASSERT(u != v);
            matches_size.fetch_add(1);
            // matches[matches_size.fetch_add(1)] = {u, v};
            partner[u] = v;
            partner[v] = u;
        }

        // EdgeUV operator[](size_t i) const {
        // ASSERT(i < matches_size);
        // return matches[i];
        // }

        size_t size() const { return matches_size; }

        vertex_t get_n() const { return m_n; }

        void clear() {
            matches_size = 0;
        }

        vertex_t get_partner(vertex_t u) const { return partner[u]; }

        void set_translation() {
            vertex_t new_u = 0;
            for (vertex_t u = 0; u < m_n; ++u) {
                vertex_t v = partner[u];
                if (u == v || u < v) {
                    o_to_n[u] = new_u;
                    o_to_n[v] = new_u; // v gets same new ID
                    new_u += 1;
                }
            }
            ASSERT(new_u == get_n_coarse_nodes());
        }

        vertex_t get_n_coarse_nodes() const { return m_n - matches_size; }

        vertex_t get_n(vertex_t o) const { return o_to_n[o]; }
    };
}

#endif //HEIPROMAP_MATCHING_H
