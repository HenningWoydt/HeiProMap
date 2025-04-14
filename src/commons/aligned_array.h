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

#include "utils.h"

namespace HeiProMap {
    template <typename T>
    class AlignedArray {
        T* ptr = nullptr;

    public:
        AlignedArray() = default;

        explicit AlignedArray(const size_t n) {
            size_t n_64 = round_up_64(n);
            ptr         = (T*)aligned_alloc(64, n_64 * sizeof(T));
        }

        explicit AlignedArray(const size_t n, const T fill_value) {
            size_t n_64 = round_up_64(n);
            ptr         = (T*)aligned_alloc(64, n_64 * sizeof(T));
            std::fill_n(ptr, n_64, fill_value);
        }

        ~AlignedArray() { free(ptr); }

        AlignedArray(const AlignedArray&) = delete;
        AlignedArray& operator=(const AlignedArray&) = delete;
        // Move constructor
        AlignedArray(AlignedArray&& other) noexcept {
            ptr = other.ptr;
            other.ptr = nullptr;
        }

        // Move assignment
        AlignedArray& operator=(AlignedArray&& other) noexcept {
            if (this != &other) {
                free(ptr);
                ptr = other.ptr;
                other.ptr = nullptr;
            }
            return *this;
        }

        T& operator[](size_t index) { return ptr[index]; }
        const T& operator[](size_t index) const { return ptr[index]; }
    };
}

#endif //HEIPROMAP_ALIGNED_ARRAY_H
