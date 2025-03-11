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

#ifndef HEIPROMAP_ISERIALPARTITIONER_H
#define HEIPROMAP_ISERIALPARTITIONER_H

#include "../serial_definitions_1.h"
#include "../../definitions.h"

namespace HeiProMap {
    class ISerialPartitionerConfiguration {
    public:
        virtual ~ISerialPartitionerConfiguration() = default;
    };

    class ISerialPartitioner {
    public:
        virtual ~ISerialPartitioner() = default;

        virtual void partition(const graph_t& g,
                               p_manager_t& p_manager,
                               const std::vector<partition_t>& hierarchy,
                               const std::vector<weight_t>& distance,
                               f64 imbalance,
                               RandomEngine& t_random_engine,
                               const ISerialPartitionerConfiguration& i_config,
                               StatisticCollector& t_stat_collect) = 0;
    };
}

#endif //HEIPROMAP_ISERIALPARTITIONER_H
