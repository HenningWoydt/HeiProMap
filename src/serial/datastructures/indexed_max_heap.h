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

#ifndef HEIPROMAP_INDEXED_MAX_HEAP_H
#define HEIPROMAP_INDEXED_MAX_HEAP_H

#include "../../commons/definitions.h"
#include "../../commons/macros.h"

namespace HeiProMap {
    /**
     * One entry in the IndexedMaxHeap.
     *
     * @tparam T The value.
     */
    template <typename T>
    class IndexedMaxHeapEntry {
    public:
        size_t key = 0;
        T val{};

        IndexedMaxHeapEntry() = default;

        IndexedMaxHeapEntry(const size_t t_key, const T t_val) {
            key = t_key;
            val = t_val;
        }
    };

    /**
     * A datastructure that acts like a heap and additionally has O(1) access
     * to all heap elements. This comes at a cost of additional O(n) memory, and
     * each element must be identifiable by a key in the range [0, n-1].
     *
     * @tparam T The value.
     */
    template <typename T>
    class IndexedMaxHeap {
    private:
        size_t m_n                     = 0;
        size_t m_heap_size             = 0;
        IndexedMaxHeapEntry<T>* m_heap = nullptr;
        size_t* m_indices              = nullptr; // Mapping of keys to heap indices

        u64 m_iteration          = 0;
        u64* m_iteration_counter = nullptr;

    public:
        IndexedMaxHeap() = default;

        ~IndexedMaxHeap() {
            free(m_heap);
            free(m_indices);
            free(m_iteration_counter);
        }

        void initialize(const size_t t_n) {
            size_t m_n_64 = round_up_64(t_n);

            m_n         = t_n;
            m_heap_size = 0;
            m_heap      = (IndexedMaxHeapEntry<T>*)aligned_alloc(64, m_n_64 * sizeof(IndexedMaxHeapEntry<T>));
            m_indices   = (size_t*)aligned_alloc(64, m_n_64 * sizeof(size_t));

            m_iteration          = 0;
            m_iteration_counter = (u64*)aligned_alloc(64, m_n_64 * sizeof(u64));
            std::fill_n(m_iteration_counter, m_n_64, 0);
        }

        void push(const size_t key, const T t) {
            ASSERT(!entry_exists(key));
            m_indices[key]           = m_heap_size;
            m_iteration_counter[key] = m_iteration;
            m_heap[m_heap_size]      = {key, t};
            m_heap_size += 1;
            bubble_up(m_heap_size - 1);
        }

        void update(const size_t key, const T t) {
            ASSERT(entry_exists(key));
            m_heap[m_indices[key]].val = t;
            bubble_up(m_indices[key]);
            bubble_down(m_indices[key]);
        }

        void increment(const size_t key, const T t) {
            ASSERT(entry_exists(key));
            m_heap[m_indices[key]].val += t;
            bubble_up(m_indices[key]);
            bubble_down(m_indices[key]);
        }

        void push_update(const size_t key, const T t) {
            if (entry_exists(key)) {
                update(key, t);
            } else {
                push(key, t);
            }
        }

        void push_increment(const size_t key, const T t) {
            if (entry_exists(key)) {
                increment(key, t);
            } else {
                push(key, t);
            }
        }

        bool entry_exists(const size_t key) const {
            ASSERT(key < m_n);
            return m_iteration_counter[key] == m_iteration && m_indices[key] != HEAP_TOMBSTONE;
        }

        T get(const size_t key) {
            ASSERT(key < m_n);
            return m_heap[m_indices[key]].val;
        }

        void pop() {
            ASSERT(!empty());

            size_t last_index        = m_heap_size - 1;
            m_indices[m_heap[0].key] = HEAP_TOMBSTONE;
            if (last_index > 0) {
                m_heap[0]                = m_heap[last_index];
                m_indices[m_heap[0].key] = 0;
                m_heap_size -= 1;
                bubble_down(0);
            } else {
                m_heap_size -= 1;
            }
        }

        size_t top_key() {
            ASSERT(!empty());
            return m_heap[0].key;
        }

        T& top() {
            ASSERT(!empty());
            return m_heap[0].val;
        }

        bool empty() const {
            return m_heap_size == 0;
        }

        size_t size() const {
            return m_heap_size;
        }

        void clear() {
            m_heap_size = 0;
            m_iteration += 1;
        }

    private:
        // Bubbles up the element at the given index to restore the heap property
        void bubble_up(size_t index) {
            while (index > 0) {
                size_t parent_index = (index - 1) / 2;
                if (m_heap[index].val <= m_heap[parent_index].val) break;

                swap(index, parent_index);
                index = parent_index;
            }
        }

        // Bubbles down the element at the given index to restore the heap property
        void bubble_down(size_t index) {
            size_t last_index = m_heap_size - 1;

            while (true) {
                size_t left_child_index  = 2 * index + 1;
                size_t right_child_index = 2 * index + 2;
                size_t largest_index     = index;

                if (left_child_index <= last_index && m_heap[left_child_index].val > m_heap[largest_index].val) {
                    largest_index = left_child_index;
                }
                if (right_child_index <= last_index && m_heap[right_child_index].val > m_heap[largest_index].val) {
                    largest_index = right_child_index;
                }
                if (largest_index == index) break;

                swap(index, largest_index);
                index = largest_index;
            }
        }

        // Swaps two elements in the heap and updates the indices
        void swap(size_t i, size_t j) {
            std::swap(m_heap[i], m_heap[j]);
            m_indices[m_heap[i].key] = i;
            m_indices[m_heap[j].key] = j;
        }
    };
}

#endif //HEIPROMAP_INDEXED_MAX_HEAP_H
