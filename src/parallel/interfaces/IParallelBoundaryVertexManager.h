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

#ifndef HEIPROMAP_IPARALLELBOUNDARYVERTEXMANAGER_H
#define HEIPROMAP_IPARALLELBOUNDARYVERTEXMANAGER_H

#include "../parallel_definitions_1.h"
#include "../parallel_definitions_2.h"
#include "../../definitions.h"

namespace HeiProMap {

    class IParallelBoundaryVertexManager {
    public:
        virtual ~IParallelBoundaryVertexManager() = default;
        virtual void initialize(vertex_t t_n, partition_t t_k) = 0;
        virtual size_t size() = 0;
        virtual size_t size(partition_t id) = 0;
        virtual vertex_t get(size_t i) = 0;
        virtual vertex_t get(partition_t id, size_t i) = 0;
        virtual bool is_boundary(vertex_t u) const = 0;
        virtual void add(vertex_t u, partition_t id) = 0;
        virtual void move(const p_graph_t& g, const p_p_manager_t& p_manager, vertex_t u, partition_t old_id, partition_t new_id) = 0;
        virtual void uncontract(const p_graph_t& g, const p_av_manager_t& av_manager, const p_p_manager_t& p_manager) = 0;
    };

}

#endif //HEIPROMAP_IPARALLELBOUNDARYVERTEXMANAGER_H
