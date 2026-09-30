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

#ifndef HEIPROMAP_BIT_VECTOR_ARRAY_H
#define HEIPROMAP_BIT_VECTOR_ARRAY_H

#include <cstddef>
#include <cstring>
#include <bit>
#include <limits>
#include "../definitions.h"
#include "../utility/aligned_array.h"

namespace HeiProMap {
    /**
     * Stores L bits for N elements in a contiguous memory layout.
     * Each element has ceil(L / 64) 64-bit words.
     */
    class BitVectorArray {
        size_t m_n = 0;
        size_t m_num_bits = 0;
        size_t m_words_per_element = 0;
        AlignedArray<u64> m_data;

    public:
        BitVectorArray() = default;

        void initialize(const size_t n, const size_t num_bits) {
            m_n = n;
            m_num_bits = num_bits;
            m_words_per_element = (num_bits + 63) / 64;
            if (m_words_per_element == 0) {
                m_words_per_element = 1;
            }
            m_data.initialize(m_n * m_words_per_element, 0);
        }

        void clear() {
            if (m_data.size() > 0) {
                m_data.fill(0);
            }
        }

        inline size_t num_elements() const { return m_n; }
        inline size_t num_bits() const { return m_num_bits; }
        inline size_t words_per_element() const { return m_words_per_element; }

        inline void set_bit(const size_t element_idx, const size_t bit_idx) {
            ASSERT(element_idx < m_n);
            ASSERT(bit_idx < m_num_bits);
            const size_t word_offset = element_idx * m_words_per_element + (bit_idx / 64);
            m_data[word_offset] |= (1ULL << (bit_idx % 64));
        }

        inline bool test_bit(const size_t element_idx, const size_t bit_idx) const {
            ASSERT(element_idx < m_n);
            ASSERT(bit_idx < m_num_bits);
            const size_t word_offset = element_idx * m_words_per_element + (bit_idx / 64);
            return (m_data[word_offset] & (1ULL << (bit_idx % 64))) != 0;
        }

        inline const u64 *get_words(const size_t element_idx) const {
            ASSERT(element_idx < m_n);
            return &m_data[element_idx * m_words_per_element];
        }

        inline u64 *get_words(const size_t element_idx) {
            ASSERT(element_idx < m_n);
            return &m_data[element_idx * m_words_per_element];
        }

        inline void or_into(const size_t element_idx, u64 *dst) const {
            ASSERT(element_idx < m_n);
            const u64 *src = &m_data[element_idx * m_words_per_element];
            for (size_t w = 0; w < m_words_per_element; ++w) {
                dst[w] |= src[w];
            }
        }

        inline void or_from(const size_t element_idx, const u64 *src) {
            ASSERT(element_idx < m_n);
            u64 *dst = &m_data[element_idx * m_words_per_element];
            for (size_t w = 0; w < m_words_per_element; ++w) {
                dst[w] |= src[w];
            }
        }

        /**
         * Finds the first bit (0 .. m_num_bits - 1) that is 0 in the combined mask.
         * Returns m_num_bits if all bits are 1.
         */
        inline size_t find_first_zero_bit(const u64 *mask) const {
            for (size_t w = 0; w < m_words_per_element; ++w) {
                const u64 inverted = ~mask[w];
                if (inverted != 0) {
                    const size_t bit = w * 64 + __builtin_ctzll(inverted);
                    return (bit < m_num_bits) ? bit : m_num_bits;
                }
            }
            return m_num_bits;
        }
    };
}

#endif //HEIPROMAP_BIT_VECTOR_ARRAY_H
