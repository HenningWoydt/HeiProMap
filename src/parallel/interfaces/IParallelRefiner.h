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

#ifndef HEIPROMAP_IPARALLELREFINER_H
#define HEIPROMAP_IPARALLELREFINER_H

#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "IParallelActiveVertexManager.h"
#include "IParallelBoundaryVertexManager.h"
#include "IParallelDistanceOracle.h"
#include "IParallelQuotientGraph.h"
#include "../datastructures/parallel_static_csr_graph.h"
#include "../datastructures/parallel_distance_oracle.h"

namespace HeiProMap {

    class IParallelRefiner {
    public:
        virtual ~IParallelRefiner() = default;

        // initialization
        virtual void initialize(vertex_t t_n,
                                vertex_t t_m,
                                partition_t t_k,
                                weight_t t_lmax,
                                const std::vector<partition_t> &t_hierarchy,
                                const std::vector<weight_t> &t_distance,
                                u64 t_seed) = 0;
    };

}

#endif //HEIPROMAP_IPARALLELREFINER_H
