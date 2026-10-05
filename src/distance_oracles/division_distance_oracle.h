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

#ifndef HEIPROMAP_DIVISION_DISTANCE_ORACLE_H
#define HEIPROMAP_DIVISION_DISTANCE_ORACLE_H

#include "../definitions.h"
#include "../utility/macros.h"
#include "../utility/utils.h"
#include "../utility/profiler.h"

namespace HeiProMap {

    // O(ℓ) time per query, O(ℓ) space.
    // Computes distance by dividing PE IDs through the hierarchy.
    class DivisionDistanceOracle {
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;
        partition_t m_k = 0;
        size_t m_l = 0;

        // h[i] = k / prod(hierarchy[0..i])
        // Dividing a PE ID by h[i] gives its group at level i (from top).
        std::vector<partition_t> m_divisors;

    public:
        size_t heap_bytes() const {
            return m_hierarchy.capacity() * sizeof(partition_t)
                 + m_distance.capacity() * sizeof(weight_t)
                 + m_divisors.capacity() * sizeof(partition_t);
        }

        void initialize(const std::vector<partition_t> &t_hierarchy,
                        const std::vector<weight_t> &t_distance) {
            HEIPROMAP_PROFILE_SCOPE("misc", "DivisionDistanceOracle", "initialize");

            m_hierarchy = t_hierarchy;
            m_distance = t_distance;
            m_l = m_hierarchy.size();
            m_k = prod<partition_t>(m_hierarchy);

            // Build divisors from top level down.
            // Level ordering: hierarchy[l-1] is the top level (coarsest).
            // m_divisors[i] corresponds to hierarchy level (l-1-i), i.e. from top to bottom.
            // m_divisors[i] = k / prod(hierarchy[l-1] * hierarchy[l-2] * ... * hierarchy[l-1-i])
            m_divisors.resize(m_l);
            partition_t product = 1;
            for (size_t i = 0; i < m_l; ++i) {
                product *= m_hierarchy[m_l - 1 - i];
                m_divisors[i] = m_k / product;
            }
        }

        weight_t get(partition_t u_id, partition_t v_id) const {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);
            if (u_id == v_id) return 0;

            for (size_t i = 0; i < m_l; ++i) {
                if (u_id / m_divisors[i] != v_id / m_divisors[i]) {
                    return m_distance[m_l - 1 - i];
                }
            }
            return 0;
        }

        partition_t get_h(partition_t u_id, partition_t v_id) const {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);
            if (u_id == v_id) return 0;

            for (size_t i = 0; i < m_l; ++i) {
                if (u_id / m_divisors[i] != v_id / m_divisors[i]) {
                    return m_l - 1 - i;
                }
            }
            return 0;
        }

        partition_t get_k() const { return m_k; }

        bool last_level_pair(partition_t u_id, partition_t v_id) const {
            return (u_id / m_hierarchy[0]) == (v_id / m_hierarchy[0]);
        }
    };

}

#endif //HEIPROMAP_DIVISION_DISTANCE_ORACLE_H
