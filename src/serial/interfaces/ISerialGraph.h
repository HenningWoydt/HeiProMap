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

#ifndef HEIPROMAP_ISERIALGRAPH_H
#define HEIPROMAP_ISERIALGRAPH_H

#include "../../definitions.h"

namespace HeiProMap {
    class ISerialGraph {
    public:
        virtual ~ISerialGraph() = default;
        virtual vertex_t get_n() const = 0;
        virtual vertex_t get_m() const = 0;
        virtual weight_t get_weight() const = 0;
        virtual weight_t get_weight(vertex_t u) const = 0;
        virtual size_t size(vertex_t u) const = 0;
        virtual vertex_t neighbor(vertex_t u, size_t idx) const = 0;
        virtual weight_t get_weight(vertex_t u, size_t idx) const = 0;
        virtual bool edge_exists(vertex_t u, vertex_t v) const = 0;
    };
}

#endif //HEIPROMAP_ISERIALGRAPH_H
