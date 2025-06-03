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

#ifndef HEIPROMAP_ISERIALREFINER_H
#define HEIPROMAP_ISERIALREFINER_H

#include <vector>

#include "../serial_definitions_1.h"
#include "../serial_definitions_2.h"
#include "../serial_definitions_3.h"
#include "../../commons/definitions.h"

namespace HeiProMap {
    class ISerialRefinerConfiguration {
    public:
        explicit ISerialRefinerConfiguration(const std::string& t_name) { name = t_name; }
        virtual ~ISerialRefinerConfiguration() = default;

        std::string name;
        bool enabled = false;
    };

    class ISerialRefiner {
    public:
        virtual ~ISerialRefiner() = default;
        virtual void initialize(vertex_t t_n,
                                vertex_t t_m,
                                partition_t t_k,
                                f64 t_imbalance,
                                weight_t t_lmax,
                                const std::vector<partition_t>& t_hierarchy,
                                const std::vector<weight_t>& t_distance,
                                RandomEngine& t_random_engine,
                                const ISerialRefinerConfiguration& i_config,
                                StatisticCollector& t_stat_collect) = 0;

        virtual void refine(u64 level,
                            u64 max_level,
                            graph_t& g,
                            d_oracle_t& d_oracle,
                            bv_manager_t& bv_manager,
                            p_manager_t& p_manager,
                            q_graph_t& q_graph) = 0;

        virtual JSONString get_stats() = 0;
    };
}

#endif //HEIPROMAP_ISERIALREFINER_H
