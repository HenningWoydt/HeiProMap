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

#ifndef HEIPROMAP_DEEP_REBALANCER_H
#define HEIPROMAP_DEEP_REBALANCER_H

#include "../../serial_definitions_1.h"
#include "../../serial_definitions_2.h"
#include "../../serial_definitions_3.h"

namespace HeiProMap {
    class DeepRebalancer {
    public:
        static void rebalance(const graph_t& g,
                       deep_p_manager_t& p_manager,
                       deep_bv_manager_t& bv_manager,
                       deep_q_graph_t& q_graph,
                       partition_t k) {
            for (partition_t id = 0; id < k; ++id) {
                if (p_manager.get_hierarchy_level(id) != k && p_manager.get_bweight(id) > p_manager.get_lmax(id)) {
                    // std::cout << id << " is unbalanced with " << p_manager.get_bweight(id) << " of " << p_manager.get_lmax(id) << std::endl;



                }
            }
        }


    };
}

#endif //HEIPROMAP_DEEP_REBALANCER_H
