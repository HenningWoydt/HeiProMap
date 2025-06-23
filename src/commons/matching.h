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

        AlignedArray<EdgeUV> matches;
        std::atomic<size_t> matches_size = 0;

        AlignedArray<vertex_t> partner;

        TranslationTable<vertex_t> tt;

    public:
        Matching() = default;

        void initialize(vertex_t n) {
            vertex_t n_64 = round_up_64(n);
            m_n           = n;

            matches.initialize((n_64 / 2));
            matches_size = 0;

            partner.initialize(m_n);
            std::iota(partner.get_ptr(), partner.get_ptr() + n_64, 0);

            tt.reserve(n, n);
        }

        // Move constructor
        Matching(Matching&& other) noexcept {
            m_n = other.m_n;
            std::swap(matches, other.matches);
            size_t temp1= matches_size;
            size_t temp2= other.matches_size;
            matches_size = temp2;
            other.matches_size = temp1;
            std::swap(partner, other.partner);
            tt.swap(other.tt);
        }

        // Optionally disable copying.
        Matching(const Matching&) = delete;

        Matching& operator=(const Matching&) = delete;

        ~Matching() = default;

        void add(vertex_t u, vertex_t v) {
            ASSERT(matches_size < (m_n / 2));
            ASSERT(u != v);
            ASSERT(partner[u] == u);
            ASSERT(partner[v] == v);
            matches[matches_size.fetch_add(1)] = {u, v};
            partner[u]              = v;
            partner[v]              = u;
        }

        EdgeUV operator[](size_t i) const {
            ASSERT(i < matches_size);
            return matches[i];
        }

        size_t size() const { return matches_size; }

        vertex_t get_n() const { return m_n; }

        void clear() {
            matches_size = 0;
        }

        bool is_matched(vertex_t u) { return u != partner[u]; }

        vertex_t get_partner(vertex_t u) const { return partner[u]; }

        void set_translation() {
            vertex_t new_u = 0;
            for (vertex_t old_u = 0; old_u < m_n; ++old_u) {
                if (old_u == partner[old_u] || old_u < partner[old_u]) {
                    tt.add(old_u, new_u);
                    new_u += 1;
                }
            }
            ASSERT(new_u == get_n_coarse_nodes());
        }

        vertex_t get_n_coarse_nodes() const { return m_n - matches_size; }

        vertex_t get_o(vertex_t n) const { return tt.get_o(n); }

        vertex_t get_n(vertex_t o) const { return tt.get_n(o); }
    };
}

#endif //HEIPROMAP_MATCHING_H
