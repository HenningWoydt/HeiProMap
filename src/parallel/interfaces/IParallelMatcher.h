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

#ifndef HEIPROMAP_IPARALLELMATCHER_H
#define HEIPROMAP_IPARALLELMATCHER_H

#include "../parallel_definitions_1.h"
#include "../../definitions.h"

namespace HeiProMap {
    class IParallelMatcherConfiguration {
    public:
        virtual ~IParallelMatcherConfiguration() = default;
    };

    class IParallelMatcher {
    public:
        virtual ~IParallelMatcher() = default;
        virtual void initialize(vertex_t t_n, vertex_t t_m, partition_t t_k, weight_t t_l_max, u64 t_seed) = 0;
        virtual void match(size_t level, IParallelMatcherConfiguration& config, p_graph_t& g, p_av_manager_t& av_manager, EdgeUV* matches, size_t& matches_size) = 0;
    };
}

#endif //HEIPROMAP_IPARALLELMATCHER_H
