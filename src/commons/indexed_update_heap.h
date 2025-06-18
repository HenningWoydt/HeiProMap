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

#ifndef HEIPROMAP_INDEXED_UPDATE_HEAP_H
#define HEIPROMAP_INDEXED_UPDATE_HEAP_H

#include "definitions.h"
#include "macros.h"
#include "utils.h"

namespace HeiProMap {
    class IndexedUpdateHeapEntry {
    public:
        vertex_t    u;
        partition_t id;
        s64         qap_delta;

        bool operator>(const IndexedUpdateHeapEntry &m) const { return qap_delta > m.qap_delta; }

        bool operator>=(const IndexedUpdateHeapEntry &m) const { return qap_delta >= m.qap_delta; }

        bool operator<(const IndexedUpdateHeapEntry &m) const { return qap_delta < m.qap_delta; }

        bool operator<=(const IndexedUpdateHeapEntry &m) const { return qap_delta <= m.qap_delta; }
    };

    class IndexedUpdateHeap {
        size_t                               m_n         = 0;
        size_t                               m_heap_size = 0;
        AlignedArray<IndexedUpdateHeapEntry> m_heap;
        AlignedArray<size_t>                 m_indices;

        u64               m_iteration = 0;
        AlignedArray<u64> m_iteration_counter;

    public:
        IndexedUpdateHeap() = default;

        ~IndexedUpdateHeap() = default;

        void initialize(const size_t t_n) {
            m_n         = t_n;
            m_heap_size = 0;
            m_heap.initialize(m_n);
            m_indices.initialize(m_n);

            m_iteration = 1;
            m_iteration_counter.initialize(m_n, 0);
        }

        void push(const vertex_t u, const partition_t move_id, const s64 qap_delta) {
            ASSERT(!entry_exists(u));
            m_indices[u]           = m_heap_size;
            m_iteration_counter[u] = m_iteration;
            m_heap[m_heap_size]    = {u, move_id, qap_delta};
            m_heap_size += 1;
            bubble_up(m_heap_size - 1);
        }

        void update(const vertex_t u, const partition_t move_id, const s64 qap_delta) {
            ASSERT(entry_exists(u));
            m_heap[m_indices[u]].id        = move_id;
            m_heap[m_indices[u]].qap_delta = qap_delta;
            bubble_up(m_indices[u]);
            bubble_down(m_indices[u]);
        }

        void push_update(const vertex_t u, const partition_t move_id, const s64 qap_delta) {
            if (entry_exists(u)) {
                update(u, move_id, qap_delta);
            } else {
                push(u, move_id, qap_delta);
            }
        }

        bool entry_exists(const size_t u) const {
            ASSERT(u < m_n);
            return m_iteration_counter[u] == m_iteration && m_indices[u] != HEAP_TOMBSTONE;
        }

        vertex_t top_u() const { return m_heap[0].u; }

        partition_t top_id() const { return m_heap[0].id; }

        s64 top_qap_delta() const { return m_heap[0].qap_delta; }

        bool empty() const { return m_heap_size == 0; }

        size_t size() const { return m_heap_size; }

        void pop() {
            ASSERT(!empty());

            size_t last_index = m_heap_size - 1;
            m_indices[m_heap[0].u] = HEAP_TOMBSTONE;
            if (last_index > 0) {
                m_heap[0]              = m_heap[last_index];
                m_indices[m_heap[0].u] = 0;
                m_heap_size -= 1;
                bubble_down(0);
            } else {
                m_heap_size -= 1;
            }
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
                if (m_heap[index] <= m_heap[parent_index]) break;

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

                if (left_child_index <= last_index && m_heap[left_child_index] > m_heap[largest_index]) {
                    largest_index = left_child_index;
                }
                if (right_child_index <= last_index && m_heap[right_child_index] > m_heap[largest_index]) {
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
            m_indices[m_heap[i].u] = i;
            m_indices[m_heap[j].u] = j;
        }
    };
}

#endif //HEIPROMAP_INDEXED_UPDATE_HEAP_H
