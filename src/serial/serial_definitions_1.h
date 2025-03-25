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

#ifndef HEIPROMAP_SERIAL_DEFINITIONS_1_H
#define HEIPROMAP_SERIAL_DEFINITIONS_1_H

#include "datastructures/distance_oracle.h"
#include "datastructures/graph_split.h"

namespace HeiProMap {
    typedef GraphSplit graph_t;
    typedef DistanceOracle d_oracle_t;

    // Macro to iterate over the neighborhood of vertex u of a graph
#define forall_guivw(g, u, i, v, w) for (size_t i = 0; i < g.size(u); ++i) { const vertex_t v = g.neighbor(u, i); const weight_t w = g.weight(u, i);
#define forall_guiv(g, u, i, v)  for (size_t i = 0; i < g.size(u); ++i) { const vertex_t v = g.neighbor(u, i);
#define forall_gu(g, u) for (vertex_t u = 0; u < g.get_n(); ++u) {

#define endfor }
}

#endif //HEIPROMAP_SERIAL_DEFINITIONS_1_H
