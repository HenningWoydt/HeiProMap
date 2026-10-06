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

#include "../definitions.h"
#include "macros.h"
#include "aligned_array.h"
#include <vector>

namespace HeiProMap {
    /**
     * One entry in the IndexedMaxHeap.
     *
     * @tparam T The value (payload).
     */
    template<typename T>
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
     * A unified indexed max-heap that allows O(1) existence checks and O(log n) updates.
     * Replaces both IndexedMaxHeap and IndexedUpdateHeap.
     *
     * @tparam T The payload type. Must support standard comparison operators if used for sorting.
     */
    template<typename T>
    class IndexedMaxHeap {
    private:
        size_t m_n = 0;
        size_t m_heap_size = 0;
        AlignedArray<IndexedMaxHeapEntry<T> > m_heap;
        AlignedArray<size_t> m_indices;

        u64 m_iteration = 1;
        AlignedArray<u64> m_iteration_counter;

    public:
        IndexedMaxHeap() = default;

        size_t heap_bytes() const {
            return m_heap.heap_bytes() + m_indices.heap_bytes() + m_iteration_counter.heap_bytes();
        }

        void set_pool(MemoryPool *pool) {
            m_heap.set_pool(pool);
            m_indices.set_pool(pool);
            m_iteration_counter.set_pool(pool);
        }

        ~IndexedMaxHeap() = default;

        void initialize(const size_t t_n) {
            m_n = t_n;
            m_heap_size = 0;
            m_heap.initialize(m_n);
            m_indices.initialize(m_n);
            m_iteration = 1;
            m_iteration_counter.initialize(m_n, 0);
        }

        void ensure_capacity(const size_t t_n) {
            if (t_n <= m_n) {
                clear();
                return;
            }
            initialize(t_n);
        }

        void grow(const size_t t_n) {
            if (t_n <= m_n) return;
            AlignedArray<IndexedMaxHeapEntry<T>> new_heap(m_heap.get_pool());
            AlignedArray<size_t> new_indices(m_indices.get_pool());
            AlignedArray<u64> new_iteration_counter(m_iteration_counter.get_pool());
            new_heap.initialize(t_n);
            new_indices.initialize(t_n);
            new_iteration_counter.initialize(t_n, 0);
            std::memcpy(new_heap.get_ptr(), m_heap.get_ptr(), m_heap_size * sizeof(IndexedMaxHeapEntry<T>));
            for (size_t i = 0; i < m_heap_size; ++i) {
                size_t key = m_heap[i].key;
                new_indices[key] = m_indices[key];
                new_iteration_counter[key] = m_iteration_counter[key];
            }
            m_heap = std::move(new_heap);
            m_indices = std::move(new_indices);
            m_iteration_counter = std::move(new_iteration_counter);
            m_n = t_n;
        }

        size_t capacity() const { return m_n; }

        void push(const size_t key, const T val) {
            ASSERT(!entry_exists(key));
            m_indices[key] = m_heap_size;
            m_iteration_counter[key] = m_iteration;
            m_heap[m_heap_size] = {key, val};
            m_heap_size += 1;
            bubble_up(m_heap_size - 1);
        }

        void push_many_heapify(const std::vector<std::pair<size_t, T> > &entries) {
            push_many_heapify(entries.data(), entries.size());
        }

        void push_many_heapify(const std::pair<size_t, T> *entries, const size_t count) {
            for (size_t i = 0; i < count; ++i) {
                const size_t key = entries[i].first;
                const T &val = entries[i].second;
                ASSERT(!entry_exists(key));
                m_indices[key] = m_heap_size;
                m_iteration_counter[key] = m_iteration;
                m_heap[m_heap_size] = {key, val};
                ++m_heap_size;
            }
            if (m_heap_size > 1) {
                for (size_t i = (m_heap_size - 2) / 2 + 1; i > 0; --i) {
                    bubble_down(i - 1);
                }
            }
        }

        void update(const size_t key, const T val) {
            ASSERT(entry_exists(key));
            size_t idx = m_indices[key];
            m_heap[idx].val = val;
            bubble_up(idx);
            bubble_down(idx);
        }

        void increment(const size_t key, const T val) {
            ASSERT(entry_exists(key));
            size_t idx = m_indices[key];
            m_heap[idx].val += val;
            bubble_up(idx);
            bubble_down(idx);
        }

        void push_update(const size_t key, const T val) {
            if (entry_exists(key)) {
                update(key, val);
            } else {
                push(key, val);
            }
        }

        void push_increment(const size_t key, const T val) {
            if (entry_exists(key)) {
                increment(key, val);
            } else {
                push(key, val);
            }
        }

        bool entry_exists(const size_t key) const {
            ASSERT(key < m_n);
            return m_iteration_counter[key] == m_iteration && m_indices[key] != HEAP_TOMBSTONE;
        }

        void remove(const size_t key) {
            ASSERT(entry_exists(key));
            size_t idx = m_indices[key];
            m_indices[key] = HEAP_TOMBSTONE;
            size_t last_idx = m_heap_size - 1;
            if (idx != last_idx) {
                m_heap[idx] = m_heap[last_idx];
                m_indices[m_heap[idx].key] = idx;
                m_heap_size -= 1;
                bubble_up(idx);
                bubble_down(idx);
            } else {
                m_heap_size -= 1;
            }
        }

        void pop() {
            ASSERT(!empty());
            m_indices[m_heap[0].key] = HEAP_TOMBSTONE;
            size_t last_idx = m_heap_size - 1;
            if (last_idx > 0) {
                m_heap[0] = m_heap[last_idx];
                m_indices[m_heap[0].key] = 0;
                m_heap_size -= 1;
                bubble_down(0);
            } else {
                m_heap_size -= 1;
            }
        }

        T get(const size_t key) const {
            ASSERT(entry_exists(key));
            return m_heap[m_indices[key]].val;
        }

        size_t top_key() const {
            ASSERT(!empty());
            return m_heap[0].key;
        }

        const T &top() const {
            ASSERT(!empty());
            return m_heap[0].val;
        }

        T &top() {
            ASSERT(!empty());
            return m_heap[0].val;
        }

        bool empty() const { return m_heap_size == 0; }
        size_t size() const { return m_heap_size; }

        void clear() {
            m_heap_size = 0;
            m_iteration += 1;
            if (m_iteration == 0) {
                // handle overflow
                m_iteration_counter.initialize(m_n, 0);
                m_iteration = 1;
            }
        }

    private:
        void bubble_up(size_t index) {
            while (index > 0) {
                size_t parent_index = (index - 1) / 2;
                if (m_heap[index].val <= m_heap[parent_index].val) break;
                swap_nodes(index, parent_index);
                index = parent_index;
            }
        }

        void bubble_down(size_t index) {
            size_t last_index = m_heap_size - 1;
            while (true) {
                size_t left = 2 * index + 1;
                size_t right = 2 * index + 2;
                size_t largest = index;
                if (left <= last_index && m_heap[left].val > m_heap[largest].val) largest = left;
                if (right <= last_index && m_heap[right].val > m_heap[largest].val) largest = right;
                if (largest == index) break;
                swap_nodes(index, largest);
                index = largest;
            }
        }

        void swap_nodes(size_t i, size_t j) {
            std::swap(m_heap[i], m_heap[j]);
            m_indices[m_heap[i].key] = i;
            m_indices[m_heap[j].key] = j;
        }
    };
}

#endif //HEIPROMAP_INDEXED_MAX_HEAP_H
