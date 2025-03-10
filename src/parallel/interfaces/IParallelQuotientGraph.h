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

#ifndef HEIPROMAP_IPARALLELQUOTIENTGRAPH_H
#define HEIPROMAP_IPARALLELQUOTIENTGRAPH_H

#include "../parallel_definitions_1.h"

namespace HeiProMap {
    class IParallelQuotientGraph {
    public:
        virtual ~IParallelQuotientGraph() = default;
        virtual void initialize(partition_t k) = 0;
        virtual void add_edge(partition_t u, partition_t v, weight_t w) = 0;
        virtual void remove_edge(partition_t u, partition_t v, weight_t w) = 0;
        virtual bool has_edge(partition_t u, partition_t v) = 0;
        virtual weight_t get_weight(partition_t u, partition_t v) = 0;
        virtual void move(p_graph_t& g, p_p_manager_t& p_manager, vertex_t u, partition_t old_id, partition_t new_id) = 0;
    };
}

#endif //HEIPROMAP_IPARALLELQUOTIENTGRAPH_H
