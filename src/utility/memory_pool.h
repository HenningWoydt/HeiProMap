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

#ifndef HEIPROMAP_MEMORY_POOL_H
#define HEIPROMAP_MEMORY_POOL_H

#include <cstdlib>
#include <cstddef>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <mutex>

namespace HeiProMap {

    class MemoryPool {
        struct FreeBlock {
            void *ptr;
            size_t size_bytes;
        };

        std::vector<FreeBlock> m_free_list;
        std::vector<void *> m_all_allocations;
        size_t m_total_allocated = 0;
        size_t m_total_in_use = 0;
        size_t m_peak_in_use = 0;
        std::mutex m_mutex;

    public:
        MemoryPool() = default;

        ~MemoryPool() {
            for (void *ptr : m_all_allocations) {
                free(ptr);
            }
        }

        MemoryPool(const MemoryPool &) = delete;
        MemoryPool &operator=(const MemoryPool &) = delete;

        void *borrow(size_t size_bytes) {
            size_bytes = round_up_64(size_bytes);
            if (size_bytes == 0) return nullptr;

            std::lock_guard<std::mutex> lock(m_mutex);

            // Find best-fit block in free list
            size_t best_idx = m_free_list.size();
            size_t best_size = std::numeric_limits<size_t>::max();
            for (size_t i = 0; i < m_free_list.size(); ++i) {
                if (m_free_list[i].size_bytes >= size_bytes && m_free_list[i].size_bytes < best_size) {
                    best_idx = i;
                    best_size = m_free_list[i].size_bytes;
                    if (best_size == size_bytes) break;
                }
            }

            if (best_idx < m_free_list.size()) {
                void *ptr = m_free_list[best_idx].ptr;
                size_t block_size = m_free_list[best_idx].size_bytes;
                m_free_list[best_idx] = m_free_list.back();
                m_free_list.pop_back();
                m_total_in_use += block_size;
                m_peak_in_use = std::max(m_peak_in_use, m_total_in_use);
                return ptr;
            }

            // No suitable block found, allocate new
            void *ptr = aligned_alloc(64, size_bytes);
            m_all_allocations.push_back(ptr);
            m_total_allocated += size_bytes;
            m_total_in_use += size_bytes;
            m_peak_in_use = std::max(m_peak_in_use, m_total_in_use);
            return ptr;
        }

        void give_back(void *ptr, size_t size_bytes) {
            size_bytes = round_up_64(size_bytes);
            if (ptr == nullptr || size_bytes == 0) return;

            std::lock_guard<std::mutex> lock(m_mutex);
            m_free_list.push_back({ptr, size_bytes});
            m_total_in_use -= size_bytes;
        }

        size_t total_allocated() const { return m_total_allocated; }
        size_t total_in_use() const { return m_total_in_use; }
        size_t total_free() const { return m_total_allocated - m_total_in_use; }
        size_t peak_in_use() const { return m_peak_in_use; }
        size_t free_list_size() const { return m_free_list.size(); }

        void print_stats() const {
            auto mb = [](size_t bytes) { return (double) bytes / (1024.0 * 1024.0); };
            std::cout << "------- Memory Pool -------" << std::endl;
            std::cout << std::fixed << std::setprecision(1);
            std::cout << "Total allocated       : " << mb(m_total_allocated) << " MB" << std::endl;
            std::cout << "Currently in use      : " << mb(m_total_in_use) << " MB" << std::endl;
            std::cout << "Free (reusable)       : " << mb(m_total_allocated - m_total_in_use) << " MB" << std::endl;
            std::cout << "Peak in use           : " << mb(m_peak_in_use) << " MB" << std::endl;
            std::cout << "Free list entries     : " << m_free_list.size() << std::endl;
        }

    private:
        static size_t round_up_64(size_t n) {
            return (n + 63) & ~size_t(63);
        }
    };

}

#endif //HEIPROMAP_MEMORY_POOL_H
