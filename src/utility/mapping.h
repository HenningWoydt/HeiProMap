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

#ifndef HEIPROMAP_MAPPING_H
#define HEIPROMAP_MAPPING_H

#include <numeric>

#include "aligned_array.h"
#include "../definitions.h"

namespace HeiProMap {
    class Mapping {
        vertex_t old_n = 0;
        vertex_t coarse_n = 0;
        AlignedArray<vertex_t> map;

    public:
        explicit Mapping() = default;

        void initialize(vertex_t n) {
            old_n = n;
            map.initialize(n);
        }

        void set(vertex_t u, vertex_t v) { map[u] = v; }
        vertex_t get(vertex_t u) const { return map[u]; }
        void set_coarse_n(vertex_t n) { coarse_n = n; }
        vertex_t get_coarse_n() const { return coarse_n; }
        vertex_t get_old_n() const { return old_n; }
    };
}

#endif //HEIPROMAP_MAPPING_H
