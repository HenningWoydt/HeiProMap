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

#ifndef HEIPROMAP_ALIGNED_ARRAY_H
#define HEIPROMAP_ALIGNED_ARRAY_H

#include "macros.h"
#include "utils.h"

namespace HeiProMap {
    template<typename T>
    class AlignedArray {
        T *m_ptr = nullptr;
        size_t m_n = 0;

        static_assert(std::is_trivially_destructible<T>::value,
                      "AlignedArray requires trivially destructible types");

    public:
        AlignedArray() = default;

        void initialize(const size_t n) {
            size_t size = round_up_64(n);

            if (size > m_n) {
                m_n = size;
                free(m_ptr);
                m_ptr = (T *) aligned_alloc(64, size * sizeof(T));
            }
        }

        void initialize(const size_t n, const T fill_value) {
            size_t size = round_up_64(n);

            if (size > m_n) {
                m_n = size;
                free(m_ptr);
                m_ptr = (T *) aligned_alloc(64, size * sizeof(T));
            }
            std::fill_n(m_ptr, size, fill_value);
        }

        void free_memory() {
            free(m_ptr);
            m_ptr = nullptr;
            m_n = 0;
        }

        ~AlignedArray() { free(m_ptr); }

        // Copy constructor
        AlignedArray(const AlignedArray &other) {
            m_n = other.m_n;
            if (m_n > 0) {
                m_ptr = (T *) aligned_alloc(64, m_n * sizeof(T));
                std::memcpy(m_ptr, other.m_ptr, m_n * sizeof(T));
            }
        }

        // Copy assignment
        AlignedArray &operator=(const AlignedArray &other) {
            if (this != &other) {
                if (m_n != other.m_n) {
                    free(m_ptr);
                    m_n = other.m_n;
                    m_ptr = m_n > 0 ? (T *) aligned_alloc(64, m_n * sizeof(T)) : nullptr;
                }
                if (m_n > 0) {
                    std::memcpy(m_ptr, other.m_ptr, m_n * sizeof(T));
                }
            }
            return *this;
        }


        AlignedArray(AlignedArray &&other) noexcept {
            free(m_ptr);
            m_n = 0;

            m_ptr = other.m_ptr;
            m_n = other.m_n;

            other.m_ptr = nullptr;
            other.m_n = 0;
        }

        AlignedArray &operator=(AlignedArray &&other) noexcept {
            if (this != &other) {
                free(m_ptr);
                m_n = 0;

                m_ptr = other.m_ptr;
                m_n = other.m_n;

                other.m_ptr = nullptr;
                other.m_n = 0;
            }
            return *this;
        }

        T &operator[](size_t index) { return m_ptr[index]; }

        const T &operator[](size_t index) const { return m_ptr[index]; }

        T *get_ptr() { return m_ptr; }
    };

    template<typename T>
    void swap(AlignedArray<T> &a, AlignedArray<T> &b) noexcept {
        using std::swap;
        swap(a.m_ptr, b.m_ptr);
        swap(a.m_n, b.m_n);
    }
}

#endif //HEIPROMAP_ALIGNED_ARRAY_H
