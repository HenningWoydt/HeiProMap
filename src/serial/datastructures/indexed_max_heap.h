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

#include <vector>

#include "../../definitions.h"
#include "../../macros.h"

namespace HeiProMap {
    /**
     * One entry in the IndexedMaxHeap.
     *
     * @tparam T The value.
     */
    template <typename T>
    class IndexedMaxHeapEntry {
    public:
        size_t key;
        T val;

        IndexedMaxHeapEntry() = default;

        IndexedMaxHeapEntry(size_t t_key, T t_val) {
            key = t_key;
            val = t_val;
        }
    };

    /**
     * A datastructure that acts like a heap and additionally has O(1) access
     * to all heap elements. This comes at a cost of additional O(n) memory and
     * each element must be identifiable by a key in the range [0, n-1].
     *
     * @tparam T The value.
     */
    template <typename T>
    class IndexedMaxHeap {
    private:
        size_t n = 0;
        std::vector<IndexedMaxHeapEntry<T>> heap; // The heap array
        std::vector<size_t> indices; // Mapping of keys to heap indices

    public:
        IndexedMaxHeap() = default;

        explicit IndexedMaxHeap(size_t t_n) {
            n = t_n;
            heap.reserve(n);
            indices.resize(n, HEAP_TOMBSTONE);
        }

        void push(size_t key, T t) {
            ASSERT(!entry_exists(key));
            indices[key] = heap.size();
            heap.emplace_back(key, t);
            bubble_up(heap.size() - 1);
        }

        void update(size_t key, T t) {
            ASSERT(entry_exists(key));
            heap[indices[key]].val = t;
            bubble_up(indices[key]);
            bubble_down(indices[key]);
        }

        void increment(size_t key, T t) {
            ASSERT(entry_exists(key));
            heap[indices[key]].val += t;
            bubble_up(indices[key]);
            bubble_down(indices[key]);
        }

        void push_update(size_t key, T t) {
            if (entry_exists(key)) {
                update(key, t);
            } else {
                push(key, t);
            }
        }

        void push_increment(size_t key, T t) {
            if (entry_exists(key)) {
                increment(key, t);
            } else {
                push(key, t);
            }
        }

        bool entry_exists(size_t key) {
            ASSERT(key < n);
            return indices[key] != HEAP_TOMBSTONE;
        }

        T get(size_t key) {
            ASSERT(key < n);
            return heap[indices[key]].val;
        }

        void pop() {
            ASSERT(!empty());

            size_t last_index    = heap.size() - 1;
            indices[heap[0].key] = HEAP_TOMBSTONE;
            if (last_index > 0) {
                heap[0]              = heap[last_index];
                indices[heap[0].key] = 0;
                heap.pop_back();
                bubble_down(0);
            } else {
                heap.pop_back();
            }
        }

        size_t top_key() {
            ASSERT(!heap.empty());
            return heap[0].key;
        }

        T& top() {
            ASSERT(!heap.empty());
            return heap[0].val;
        }

        bool empty() const {
            return heap.empty();
        }

        size_t size() const {
            return heap.size();
        }

        void clear() {
            heap.clear();
            std::fill(indices.begin(), indices.end(), HEAP_TOMBSTONE);
        }

    private:
        // Bubbles up the element at the given index to restore the heap property
        void bubble_up(size_t index) {
            while (index > 0) {
                size_t parent_index = (index - 1) / 2;
                if (heap[index].val <= heap[parent_index].val) break;

                swap(index, parent_index);
                index = parent_index;
            }
        }

        // Bubbles down the element at the given index to restore the heap property
        void bubble_down(size_t index) {
            size_t last_index = heap.size() - 1;

            while (true) {
                size_t left_child_index  = 2 * index + 1;
                size_t right_child_index = 2 * index + 2;
                size_t largest_index     = index;

                if (left_child_index <= last_index && heap[left_child_index].val > heap[largest_index].val) {
                    largest_index = left_child_index;
                }
                if (right_child_index <= last_index && heap[right_child_index].val > heap[largest_index].val) {
                    largest_index = right_child_index;
                }
                if (largest_index == index) break;

                swap(index, largest_index);
                index = largest_index;
            }
        }

        // Swaps two elements in the heap and updates the indices
        void swap(size_t i, size_t j) {
            std::swap(heap[i], heap[j]);
            indices[heap[i].key] = i;
            indices[heap[j].key] = j;
        }
    };
}

#endif //HEIPROMAP_INDEXED_MAX_HEAP_H
