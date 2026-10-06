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

#include <algorithm>
#include <cstring>
#include <vector>
#include "macros.h"
#include "utils.h"
#include "memory_pool.h"

namespace HeiProMap {
    template<typename T>
    class AlignedArray {
        T *m_ptr = nullptr;
        size_t m_n = 0;
        MemoryPool *m_pool = nullptr;
        bool m_own_memory = true;

        static_assert(std::is_trivially_destructible<T>::value,
                      "AlignedArray requires trivially destructible types");

        void do_alloc(size_t size) {
            m_n = size;
            if (m_pool) {
                m_ptr = (T *) m_pool->borrow(size * sizeof(T));
                m_own_memory = false;
            } else {
                m_ptr = (T *) aligned_alloc(64, size * sizeof(T));
                m_own_memory = true;
            }
        }

        void do_free() {
            if (m_ptr) {
                if (m_own_memory) {
                    free(m_ptr);
                } else {
                    m_pool->give_back(m_ptr, m_n * sizeof(T));
                }
                m_ptr = nullptr;
                m_n = 0;
            }
        }

    public:
        AlignedArray() = default;

        explicit AlignedArray(MemoryPool *pool) : m_pool(pool), m_own_memory(true) {}

        void set_pool(MemoryPool *pool) { m_pool = pool; }
        MemoryPool *get_pool() const { return m_pool; }

        void initialize(const size_t n) {
            size_t size = round_up_64(n);

            if (size > m_n) {
                do_free();
                do_alloc(size);
            }
        }

        void initialize(const size_t n, const T fill_value) {
            size_t size = round_up_64(n);

            if (size > m_n) {
                do_free();
                do_alloc(size);
            }
            std::fill_n(m_ptr, size, fill_value);
        }

        void free_memory() {
            do_free();
        }

        ~AlignedArray() { do_free(); }

        AlignedArray(const AlignedArray &other) {
            m_pool = other.m_pool;
            m_n = other.m_n;
            if (m_n > 0) {
                do_alloc(m_n);
                std::memcpy(m_ptr, other.m_ptr, m_n * sizeof(T));
            }
        }

        AlignedArray &operator=(const AlignedArray &other) {
            if (this != &other) {
                if (m_n != other.m_n) {
                    do_free();
                    m_pool = other.m_pool;
                    m_n = 0;
                    if (other.m_n > 0) {
                        do_alloc(other.m_n);
                    }
                }
                if (m_n > 0) {
                    std::memcpy(m_ptr, other.m_ptr, m_n * sizeof(T));
                }
            }
            return *this;
        }

        AlignedArray(AlignedArray &&other) noexcept {
            m_ptr = other.m_ptr;
            m_n = other.m_n;
            m_pool = other.m_pool;
            m_own_memory = other.m_own_memory;

            other.m_ptr = nullptr;
            other.m_n = 0;
        }

        AlignedArray &operator=(AlignedArray &&other) noexcept {
            if (this != &other) {
                do_free();
                m_ptr = other.m_ptr;
                m_n = other.m_n;
                m_pool = other.m_pool;
                m_own_memory = other.m_own_memory;

                other.m_ptr = nullptr;
                other.m_n = 0;
            }
            return *this;
        }

        T &operator[](size_t index) {
            T * HEIPROMAP_RESTRICT ptr = HEIPROMAP_ASSUME_ALIGNED(m_ptr, 64);
            return ptr[index];
        }

        const T &operator[](size_t index) const {
            const T * HEIPROMAP_RESTRICT ptr = HEIPROMAP_ASSUME_ALIGNED(m_ptr, 64);
            return ptr[index];
        }

        T *get_ptr() {
            return HEIPROMAP_ASSUME_ALIGNED(m_ptr, 64);
        }

        const T *get_ptr() const {
            return HEIPROMAP_ASSUME_ALIGNED(m_ptr, 64);
        }

        std::vector<T> get_vector() const {
            std::vector<T> vec(m_n);
            if (m_n > 0) {
                std::memcpy(vec.data(), m_ptr, m_n * sizeof(T));
            }
            return vec;
        }

        size_t size() const { return m_n; }
        size_t heap_bytes() const { return m_own_memory ? m_n * sizeof(T) : 0; }

        T *begin() { return get_ptr(); }
        const T *begin() const { return get_ptr(); }
        T *end() { return get_ptr() + m_n; }
        const T *end() const { return get_ptr() + m_n; }

        void fill(const T value) {
            if (m_n > 0) {
                if constexpr (sizeof(T) == 1) {
                    std::memset(m_ptr, static_cast<int>(value), m_n);
                } else {
                    if (value == 0) {
                        std::memset(m_ptr, 0, m_n * sizeof(T));
                    } else {
                        std::fill_n(get_ptr(), m_n, value);
                    }
                }
            }
        }

        void swap(AlignedArray &other) noexcept {
            std::swap(m_ptr, other.m_ptr);
            std::swap(m_n, other.m_n);
            std::swap(m_pool, other.m_pool);
            std::swap(m_own_memory, other.m_own_memory);
        }
    };

    template<typename T>
    void swap(AlignedArray<T> &a, AlignedArray<T> &b) noexcept {
        a.swap(b);
    }
}

#endif //HEIPROMAP_ALIGNED_ARRAY_H
