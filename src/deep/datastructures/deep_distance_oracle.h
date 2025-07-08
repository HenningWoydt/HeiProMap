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

#ifndef HEIPROMAP_DEEP_DISTANCE_ORACLE_H
#define HEIPROMAP_DEEP_DISTANCE_ORACLE_H


#include "../../commons/definitions.h"
#include "../../commons/macros.h"
#include "../../commons/utils.h"

namespace HeiProMap {
    class DeepDistanceOracle {
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t>    m_distance;
        partition_t              m_k = 0;

        u64 m_threads = 1;

        std::vector<AlignedArray<partition_t>> locs;

    public:
        void initialize(const std::vector<partition_t> &t_hierarchy,
                        const std::vector<weight_t> &t_distance,
                        const u64 t_threads) {
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;
            m_k         = prod<partition_t>(m_hierarchy);

            m_threads = t_threads;

            AlignedArray<partition_t> u_loc;
            AlignedArray<partition_t> v_loc;
            u_loc.initialize(m_hierarchy.size());
            v_loc.initialize(m_hierarchy.size());

            locs.resize(m_k);
            for (partition_t i = 0; i < m_k; ++i) {
                locs[i].initialize(m_hierarchy.size());
                determine_loc(i, locs[i]);
            }
        }

        weight_t get(const partition_t u_id,
                     const partition_t v_id) const {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);

            if (u_id == v_id) { return 0; }

            return determine_distance(locs[u_id], locs[v_id]);
        }

        partition_t get_h(const partition_t u_id,
                          const partition_t v_id) const {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);

            if (u_id == v_id) { return 0; }

            return determine_hierarchy(locs[u_id], locs[v_id]);
        }

        bool last_level_pair(const partition_t u_id,
                             const partition_t v_id) const {
            return (u_id / m_hierarchy[0]) == (v_id / m_hierarchy[0]);
        }

    private:
        void determine_loc(partition_t u_id,
                           AlignedArray<partition_t> &loc) const {
            ASSERT(prod<u64>(m_hierarchy) == m_k);
            ASSERT(u_id < m_k);
            ASSERT(!m_hierarchy.empty());

            u64 r_start = 0;
            u64 r_end   = m_k;

            u64      s = m_hierarchy.size();
            for (u64 i = 0; i < m_hierarchy.size(); ++i) {
                u64 n_parts = m_hierarchy[s - 1 - i];
                u64 add     = (r_end - r_start) / n_parts;

                for (u64 j = 0; j < n_parts; ++j) {
                    if (r_start <= u_id && u_id < r_start + add) {
                        // we have found the part of the partition
                        loc[s - 1 - i] = j;
                        r_end = r_start + add;
                        break;
                    } else {
                        r_start += add;
                    }
                }
            }
        }

        weight_t determine_distance(const AlignedArray<partition_t> &loc_1,
                                    const AlignedArray<partition_t> &loc_2) const {
            ASSERT(prod<u64>(m_hierarchy) == m_k);
            ASSERT(!m_hierarchy.empty());
            ASSERT(m_hierarchy.size() == m_distance.size());

            // determine the distance
            u64      s = m_hierarchy.size();
            for (u64 i = 0; i < m_hierarchy.size(); ++i) {
                if (loc_1[s - 1 - i] != loc_2[s - 1 - i]) {
                    return m_distance[s - 1 - i];
                }
            }
            // unreachable
            abort();
        }

        partition_t determine_hierarchy(const AlignedArray<partition_t> &loc_1,
                                        const AlignedArray<partition_t> &loc_2) const {
            ASSERT(prod<u64>(m_hierarchy) == m_k);
            ASSERT(!m_hierarchy.empty());
            ASSERT(m_hierarchy.size() == m_distance.size());

            // determine the distance
            u64      s = m_hierarchy.size();
            for (u64 i = 0; i < m_hierarchy.size(); ++i) {
                if (loc_1[s - 1 - i] != loc_2[s - 1 - i]) {
                    return s - 1 - i;
                }
            }
            // unreachable
            abort();
        }
    };
}

#endif //HEIPROMAP_DEEP_DISTANCE_ORACLE_H
