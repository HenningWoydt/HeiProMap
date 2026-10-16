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

#ifndef HEIPROMAP_ACTIVE_BLOCK_SCHEDULING_H
#define HEIPROMAP_ACTIVE_BLOCK_SCHEDULING_H

#include "../utility/aligned_array.h"

namespace HeiProMap {
    class ActiveBlockScheduling {
    public:
        partition_t m_k = 0;
        // active block scheduling
        AlignedArray<u8> active_this_round;
        AlignedArray<u8> active_next_round;

        size_t heap_bytes() const {
            return active_this_round.heap_bytes() + active_next_round.heap_bytes();
        }

        void set_pool(MemoryPool *pool) {
            active_this_round.set_pool(pool);
            active_next_round.set_pool(pool);
        }

        void initialize(const partition_t t_k, MemoryPool *pool = nullptr) {
            if (pool) set_pool(pool);
            m_k = t_k;
            active_this_round.initialize(t_k);
            active_next_round.initialize(t_k);
        }

        void reset(const partition_t t_k) {
            active_this_round.initialize(t_k, 1);
            active_next_round.initialize(t_k, 0);
        }

        void next_round() {
            std::swap(active_this_round, active_next_round);
            active_next_round.initialize(m_k, 0);
        }

        void activate(partition_t id1, partition_t id2) {
            active_next_round[id1] = 1;
            active_next_round[id2] = 1;
        }
    };
}

#endif //HEIPROMAP_ACTIVE_BLOCK_SCHEDULING_H
