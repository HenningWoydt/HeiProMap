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

#ifndef HEIPROMAP_STORED_DIVISION_DISTANCE_ORACLE_H
#define HEIPROMAP_STORED_DIVISION_DISTANCE_ORACLE_H

#include "../definitions.h"
#include "../utility/macros.h"
#include "../utility/utils.h"
#include "../utility/profiler.h"

namespace HeiProMap {

    // O(ℓ) time per query, O(kℓ) space.
    // Precomputes division results for every PE to avoid divisions at query time.
    class StoredDivisionDistanceOracle {
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;
        partition_t m_k = 0;
        size_t m_l = 0;

        // m_groups[id * m_l + i] = id / divisor[i], the group of PE id at level i (from top).
        std::vector<partition_t> m_groups;

    public:
        size_t heap_bytes() const {
            return m_hierarchy.capacity() * sizeof(partition_t)
                 + m_distance.capacity() * sizeof(weight_t)
                 + m_groups.capacity() * sizeof(partition_t);
        }

        void initialize(const std::vector<partition_t> &t_hierarchy,
                        const std::vector<weight_t> &t_distance) {
            HEIPROMAP_PROFILE_SCOPE("misc", "StoredDivisionDistanceOracle", "initialize");

            m_hierarchy = t_hierarchy;
            m_distance = t_distance;
            m_l = m_hierarchy.size();
            m_k = prod<partition_t>(m_hierarchy);

            // Build divisors (same as DivisionDistanceOracle)
            std::vector<partition_t> divisors(m_l);
            partition_t product = 1;
            for (size_t i = 0; i < m_l; ++i) {
                product *= m_hierarchy[m_l - 1 - i];
                divisors[i] = m_k / product;
            }

            // Precompute all division results
            m_groups.resize(static_cast<size_t>(m_k) * m_l);
            for (partition_t id = 0; id < m_k; ++id) {
                for (size_t i = 0; i < m_l; ++i) {
                    m_groups[static_cast<size_t>(id) * m_l + i] = id / divisors[i];
                }
            }
        }

        weight_t get(partition_t u_id, partition_t v_id) const {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);
            if (u_id == v_id) return 0;

            const partition_t *u_groups = m_groups.data() + static_cast<size_t>(u_id) * m_l;
            const partition_t *v_groups = m_groups.data() + static_cast<size_t>(v_id) * m_l;

            for (size_t i = 0; i < m_l; ++i) {
                if (u_groups[i] != v_groups[i]) {
                    return m_distance[m_l - 1 - i];
                }
            }
            return 0;
        }

        partition_t get_h(partition_t u_id, partition_t v_id) const {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);
            if (u_id == v_id) return 0;

            const partition_t *u_groups = m_groups.data() + static_cast<size_t>(u_id) * m_l;
            const partition_t *v_groups = m_groups.data() + static_cast<size_t>(v_id) * m_l;

            for (size_t i = 0; i < m_l; ++i) {
                if (u_groups[i] != v_groups[i]) {
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

#endif //HEIPROMAP_STORED_DIVISION_DISTANCE_ORACLE_H
