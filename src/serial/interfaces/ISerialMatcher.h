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

#ifndef HEIPROMAP_ISERIALMATCHER_H
#define HEIPROMAP_ISERIALMATCHER_H

#include "../serial_definitions_1.h"
#include "../../commons/definitions.h"
#include "../../commons/random_engine.h"
#include "../../commons/statistic_collector.h"

namespace HeiProMap {
    class ISerialMatcherConfiguration {
    public:
        virtual ~ISerialMatcherConfiguration() = default;
    };

    class ISerialMatcher {
    public:
        virtual ~ISerialMatcher() = default;

        virtual void initialize(vertex_t t_n,
                                vertex_t t_m,
                                partition_t t_k,
                                weight_t t_l_max,
                                RandomEngine& t_random_engine,
                                const ISerialMatcherConfiguration& i_config,
                                StatisticCollector& t_stat_collect) = 0;

        virtual void match(size_t level,
                           const graph_t& g,
                           p_manager_t& p_manager,
                           Matching &matching) = 0;

        virtual JSONString get_stats() = 0;
    };
}

#endif //HEIPROMAP_ISERIALMATCHER_H
