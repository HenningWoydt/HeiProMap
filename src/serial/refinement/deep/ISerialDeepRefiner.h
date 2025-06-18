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

#ifndef HEIPROMAP_ISERIALDEEPREFINER_H
#define HEIPROMAP_ISERIALDEEPREFINER_H

#include <vector>

#include "../../serial_definitions_1.h"
#include "../../serial_definitions_2.h"
#include "../../serial_definitions_3.h"
#include "../../../commons/definitions.h"
#include "../../../commons/random_engine.h"
#include "../../../commons/statistic_collector.h"

namespace HeiProMap {
    class ISerialDeepRefinerConfiguration {
    public:
        explicit ISerialDeepRefinerConfiguration(const std::string &t_name) { name = t_name; }

        virtual ~ISerialDeepRefinerConfiguration() = default;

        std::string name;
        bool        enabled = false;
    };

    class ISerialDeepRefiner {
    public:
        virtual ~ISerialDeepRefiner() = default;

        virtual void initialize(const vertex_t t_n,
                                const vertex_t t_m,
                                const partition_t t_k,
                                const f64 t_imbalance,
                                const u64 t_threads,
                                const std::vector<partition_t> &t_hierarchy,
                                const std::vector<weight_t> &t_distance,
                                RandomEngine &t_random_engine,
                                const ISerialDeepRefinerConfiguration &i_config) = 0;

        virtual void refine(u64 level,
                            u64 max_level,
                            const graph_t &g,
                            deep_d_oracle_t &d_oracle,
                            deep_bv_manager_t &bv_manager,
                            deep_p_manager_t &p_manager,
                            deep_q_graph_t &q_graph) = 0;
    };
}

#endif //HEIPROMAP_ISERIALDEEPREFINER_H
