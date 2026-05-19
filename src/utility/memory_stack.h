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

#ifndef HEIPROMAP_MEMORY_STACK_H
#define HEIPROMAP_MEMORY_STACK_H

#include <cstdint>
#include <cstdlib>

namespace HeiProMap {
    inline uint64_t align64(uint64_t n) { return (n + 63) & ~uint64_t(63); }

    class MemoryStack {
        static constexpr size_t ALIGNMENT = 64;
        char *buffer = nullptr;
        size_t capacity = 0;
        size_t offset = 0;

    public:
        ~MemoryStack() { std::free(buffer); }

        MemoryStack() = default;

        MemoryStack(const MemoryStack &) = delete;

        MemoryStack &operator=(const MemoryStack &) = delete;

        MemoryStack(MemoryStack &&other) noexcept: buffer(other.buffer), capacity(other.capacity), offset(other.offset) {
            other.buffer = nullptr;
            other.capacity = 0;
            other.offset = 0;
        }

        MemoryStack &operator=(MemoryStack &&other) noexcept {
            if (this != &other) {
                std::free(buffer);
                buffer = other.buffer;
                capacity = other.capacity;
                offset = other.offset;
                other.buffer = nullptr;
                other.capacity = 0;
                other.offset = 0;
            }
            return *this;
        }

        void ensure(size_t total_bytes) {
            total_bytes = align64(total_bytes);
            if (total_bytes > capacity) {
                std::free(buffer);
                buffer = static_cast<char *>(std::aligned_alloc(ALIGNMENT, total_bytes));
                capacity = total_bytes;
            }
        }

        void *get_memory(size_t size) {
            size = align64(size);
            void *ptr = buffer + offset;
            offset += size;
            return ptr;
        }

        void clear() { offset = 0; }
    };
}

#endif // HEIPROMAP_MEMORY_STACK_H
