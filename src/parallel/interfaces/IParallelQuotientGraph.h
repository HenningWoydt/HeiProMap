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

#include "../../interfaces/IQuotientGraph.h"

namespace HeiProMap {

    template<typename TParallelGraph,
             typename TParallelActiveVertexManager,
             typename TParallelPartitionManager,
             typename TParallelDistanceOracle>
    class IParallelQuotientGraph : public IQuotientGraph {
        static_assert(std::is_base_of<IParallelGraph, TParallelGraph>::value, "TParallelGraph must inherit from IParallelGraph");
        static_assert(std::is_base_of<IParallelPartitionManager<TParallelGraph, TParallelActiveVertexManager>, TParallelPartitionManager>::value, "TParallelPartitionManager must inherit from IParallelPartitionManager");
        static_assert(std::is_base_of<IParallelDistanceOracle, TParallelDistanceOracle>::value, "TParallelDistanceOracle must inherit from IParallelDistanceOracle");
    public:
        virtual void initialize(TParallelGraph *t_p_g,
                                TParallelPartitionManager *t_p_p_manager,
                                TParallelDistanceOracle *t_p_d_oracle,
                                partition_t k,
                                u64 n_threads) = 0;

        virtual void move(vertex_t u, partition_t old_id, partition_t new_id) = 0;
    };

}

#endif //HEIPROMAP_IPARALLELQUOTIENTGRAPH_H
