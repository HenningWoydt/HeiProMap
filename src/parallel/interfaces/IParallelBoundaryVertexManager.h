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

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../../interfaces/IBoundaryVertexManager.h"
#include "IParallelGraph.h"
#include "IParallelActiveVertexManager.h"
#include "IParallelPartitionManager.h"

namespace HeiProMap {

    template<typename TParallelGraph, typename TParallelActiveVertexManager, typename TParallelPartitionManager>
    class IParallelBoundaryVertexManager : public IBoundaryVertexManager {
        static_assert(std::is_base_of<IParallelGraph, TParallelGraph>::value, "TParallelGraph must inherit from IParallelGraph");
        static_assert(std::is_base_of<IParallelActiveVertexManager<TParallelGraph>, TParallelActiveVertexManager>::value, "TParallelActiveVertexManager must inherit from IParallelActiveVertexManager");
        static_assert(std::is_base_of<IParallelPartitionManager<TParallelGraph, TParallelActiveVertexManager>, TParallelPartitionManager>::value, "TParallelActiveVertexManager must inherit from IParallelActiveVertexManager");
    public:
        // initialization
        virtual void initialize(TParallelGraph *t_p_g,
                                TParallelActiveVertexManager *t_p_av_manager,
                                TParallelPartitionManager *t_p_p_manger,
                                partition_t t_k,
                                u64 n_threads) = 0;

        virtual vertex_t get_vertex(size_t idx) const = 0;
    };

}

#endif //HEIPROMAP_IPARALLELBOUNDARYVERTEXMANAGER_H
