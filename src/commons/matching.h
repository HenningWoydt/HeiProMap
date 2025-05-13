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

#include <numeric>

#include "definitions.h"
#include "translation_table.h"
#include "utils.h"

namespace HeiProMap {
    class Matching {
    private:
        vertex_t m_n = 0;

        EdgeUV* matches     = nullptr;
        size_t matches_size = 0;

        vertex_t* partner = nullptr;

        TranslationTable<vertex_t> tt;

    public:
        Matching() = default;

        void initialize(vertex_t n) {
            vertex_t n_64 = round_up_64(n);
            m_n           = n;

            matches      = (EdgeUV*)aligned_alloc(64, (n_64 / 2) * sizeof(EdgeUV));
            matches_size = 0;

            partner = (vertex_t*)aligned_alloc(64, n_64 * sizeof(vertex_t));
            std::iota(partner, partner + n_64, 0);

            tt.reserve(n, n);
        }

        // Move constructor
        Matching(Matching&& other) noexcept {
            m_n          = other.m_n;
            matches      = other.matches;
            matches_size = other.matches_size;
            partner      = other.partner;
            std::swap(tt, other.tt);

            other.m_n          = 0;
            other.matches      = nullptr;
            other.matches_size = 0;
            other.partner      = nullptr;
        }

        // Optionally disable copying.
        Matching(const Matching&) = delete;

        Matching& operator=(const Matching&) = delete;

        ~Matching() {
            free(matches);
            free(partner);
        }

        void add(vertex_t u, vertex_t v) {
            ASSERT(matches_size < (m_n / 2));
            matches[matches_size++] = {u, v};
            partner[u]  = v;
            partner[v]  = u;
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
