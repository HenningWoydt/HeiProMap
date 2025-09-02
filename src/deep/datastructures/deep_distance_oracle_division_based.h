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

#ifndef HEIPROMAP_DEEP_DISTANCE_ORACLE_DIVISION_BASED_H
#define HEIPROMAP_DEEP_DISTANCE_ORACLE_DIVISION_BASED_H

#include "../../commons/definitions.h"
#include "../../commons/macros.h"
#include "../../commons/utils.h"
#include "../../commons/aligned_array.h"

namespace HeiProMap {
    class DeepDistanceOracleDivisionBased {
        std::vector<partition_t> m_hierarchy;           // O(l)
        std::vector<weight_t> m_distance;               // O(l)
        size_t m_l = 0;
        partition_t m_k = 0;

        std::vector<partition_t> divisor;               // O(l)

    public:
        void initialize(const std::vector<partition_t> &t_hierarchy,
                        const std::vector<weight_t> &t_distance,
                        const u64 t_threads) {
            m_hierarchy = t_hierarchy;
            m_distance = t_distance;
            m_l = m_hierarchy.size();
            m_k = prod<partition_t>(m_hierarchy);

            divisor.resize(m_l);
            partition_t p = 1;
            for (size_t i = 0; i < m_l; ++i) {
                p *= m_hierarchy[m_l - 1 - i];
                divisor[m_l - 1 - i] = m_k / p;
            }
        }

        weight_t get(const partition_t u_id,
                     const partition_t v_id) const {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);

            if (u_id == v_id) { return 0; }

            for (size_t i = 0; i < m_l; ++i) {
                if ((u_id / divisor[m_l - 1 - i]) != (v_id / divisor[m_l - 1 - i])) {
                    return m_distance[m_l - 1 - i];
                }
            }
            return 0;
        }

        partition_t get_h(const partition_t u_id,
                          const partition_t v_id) const {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);

            if (u_id == v_id) { return 0; }

            for (size_t i = 0; i < m_l; ++i) {
                if ((u_id / divisor[m_l - 1 - i]) != (v_id / divisor[m_l - 1 - i])) {
                    return m_l - 1 - i;
                }
            }

            return 0;
        }

        bool last_level_pair(const partition_t u_id,
                             const partition_t v_id) const {
            return (u_id / m_hierarchy[0]) == (v_id / m_hierarchy[0]);
        }
    };
}

#endif //HEIPROMAP_DEEP_DISTANCE_ORACLE_DIVISION_BASED_H
