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

#ifndef HEIPROMAP_SERIAL_DEFINITIONS_3_H
#define HEIPROMAP_SERIAL_DEFINITIONS_3_H

#include "datastructures/deep_boundary_vertex_manager.h"
#include "datastructures/deep_quotient_graph.h"
#include "datastructures/boundary_vertex_manger.h"
#include "datastructures/quotient_graph.h"

namespace HeiProMap {
#if USE_DEEP_DATASTRUCTURES
    typedef DeepBoundaryVertexManager bv_manager_t;
    typedef DeepQuotientGraph         q_graph_t;
#else
    typedef BoundaryVertexManager bv_manager_t;
    typedef QuotientGraph q_graph_t;
#endif

    // Macro to iterate over all boundary vertices
#define forall_bv_iu(bv_manager, i, u) for (size_t i = 0; i < bv_manager.size(); ++i) { const vertex_t u = bv_manager.get(i);
#define forall_bv_id_iu(bv_manager, id, i, u) for (size_t i = 0; i < bv_manager.size(id); ++i) { const vertex_t u = bv_manager.get(id, i);
}

#endif //HEIPROMAP_SERIAL_DEFINITIONS_3_H
