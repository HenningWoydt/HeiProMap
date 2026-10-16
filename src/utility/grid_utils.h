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

#ifndef HEIPROMAP_GRID_UTILS_H
#define HEIPROMAP_GRID_UTILS_H

#include <cstddef>

#include "../definitions.h"

namespace HeiProMap {
    // Writes the (up to four) grid neighbors of cell b on a grid of width W with k cells (cell id = row * W + col)
    // to nbrs, in the order up, down, left, right, and returns how many there are.
    inline size_t grid_neighbors(const partition_t b, const size_t W, const partition_t k, partition_t *nbrs) {
        size_t cnt = 0;
        const size_t c = b % W;
        if (b >= static_cast<partition_t>(W)) nbrs[cnt++] = b - static_cast<partition_t>(W);
        if (b + static_cast<partition_t>(W) < k) nbrs[cnt++] = b + static_cast<partition_t>(W);
        if (c > 0) nbrs[cnt++] = b - 1;
        if (c + 1 < W && b + 1 < k) nbrs[cnt++] = b + 1;
        return cnt;
    }
}

#endif //HEIPROMAP_GRID_UTILS_H
