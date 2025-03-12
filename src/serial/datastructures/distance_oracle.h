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

#ifndef HEIPROMAP_DISTANCE_ORACLE_H
#define HEIPROMAP_DISTANCE_ORACLE_H

#include "../../definitions.h"
#include "../../macros.h"
#include "../../commons/utils.h"
#include "../interfaces/ISerialDistanceOracle.h"

namespace HeiProMap {
    class DistanceOracle final : public ISerialDistanceOracle {
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;
        partition_t m_k = 0;

        weight_t* m_mtx = nullptr;
        partition_t* m_h_mtx = nullptr;

    public:
        DistanceOracle() = default;
        ~DistanceOracle() override {
            free(m_mtx);
            free(m_h_mtx);
        }

        void initialize(const std::vector<partition_t>& t_hierarchy,
                        const std::vector<weight_t>& t_distance) override {
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;
            m_k         = prod<partition_t>(m_hierarchy);

            size_t m_k_m_k_64 = round_up_64(m_k * m_k);
            m_mtx = (weight_t*) aligned_alloc(64, m_k_m_k_64 * sizeof(weight_t));
            m_h_mtx = (partition_t*) aligned_alloc(64, m_k_m_k_64 * sizeof(partition_t));

            std::vector<std::vector<partition_t>> locs(m_k, std::vector<partition_t>(m_hierarchy.size()));
            for (partition_t id = 0; id < m_k; ++id) {
                determine_loc(id, locs[id]);
            }

            for (partition_t u_id = 0; u_id < m_k; ++u_id) {
                m_mtx[u_id * m_k + u_id]   = 0;
                m_h_mtx[u_id * m_k + u_id] = 0;
                for (partition_t v_id = u_id + 1; v_id < m_k; ++v_id) {
                    weight_t d               = determine_distance(locs[u_id], locs[v_id]);
                    m_mtx[u_id * m_k + v_id] = d;
                    m_mtx[v_id * m_k + u_id] = d;

                    partition_t h              = determine_hierarchy(locs[u_id], locs[v_id]);
                    m_h_mtx[u_id * m_k + v_id] = h;
                    m_h_mtx[v_id * m_k + u_id] = h;
                }
            }
        }

        weight_t get(partition_t u_id, partition_t v_id) const {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);
            ASSERT(u_id * m_k + v_id < m_k * m_k);
            weight_t* m_mtx_copy = ASSUME_ALIGNED(weight_t*, m_mtx, 64);
            return m_mtx_copy[u_id * m_k + v_id];
        }

        partition_t get_h(partition_t u_id, partition_t v_id) const {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);
            ASSERT(u_id * m_k + v_id < m_k * m_k);
            partition_t* m_h_mtx_copy = ASSUME_ALIGNED(partition_t*, m_h_mtx, 64);
            return m_h_mtx_copy[u_id * m_k + v_id];
        }

    private:
        void determine_loc(partition_t u_id,
                           std::vector<partition_t>& u_loc) const {
            ASSERT(prod<u64>(m_hierarchy) == m_k);
            ASSERT(u_id < m_k);
            ASSERT(!m_hierarchy.empty());
            ASSERT(u_loc.size() == m_hierarchy.size());

            u64 r_start = 0;
            u64 r_end   = m_k;

            u64 s = m_hierarchy.size();
            for (u64 i = 0; i < m_hierarchy.size(); ++i) {
                u64 n_parts = m_hierarchy[s - 1 - i];
                u64 add     = (r_end - r_start) / n_parts;

                for (u64 j = 0; j < n_parts; ++j) {
                    if (r_start <= u_id && u_id < r_start + add) {
                        // we have found the part of the partition
                        u_loc[s - 1 - i] = j;
                        r_end            = r_start + add;
                        break;
                    } else {
                        r_start += add;
                    }
                }
            }
        }

        weight_t determine_distance(const std::vector<partition_t>& u_loc,
                                    const std::vector<partition_t>& v_loc) const {
            ASSERT(prod<u64>(m_hierarchy) == m_k);
            ASSERT(!m_hierarchy.empty());
            ASSERT(m_hierarchy.size() == m_distance.size());
            ASSERT(u_loc.size() == m_hierarchy.size());
            ASSERT(v_loc.size() == m_hierarchy.size());

            // determine the distance
            u64 s = m_hierarchy.size();
            for (u64 i = 0; i < m_hierarchy.size(); ++i) {
                if (u_loc[s - 1 - i] != v_loc[s - 1 - i]) {
                    return m_distance[s - 1 - i];
                }
            }
            // unreachable
            abort();
        }

        partition_t determine_hierarchy(const std::vector<partition_t>& u_loc,
                                        const std::vector<partition_t>& v_loc) const {
            ASSERT(prod<u64>(m_hierarchy) == m_k);
            ASSERT(!m_hierarchy.empty());
            ASSERT(m_hierarchy.size() == m_distance.size());
            ASSERT(u_loc.size() == m_hierarchy.size());
            ASSERT(v_loc.size() == m_hierarchy.size());

            // determine the distance
            u64 s = m_hierarchy.size();
            for (u64 i = 0; i < m_hierarchy.size(); ++i) {
                if (u_loc[s - 1 - i] != v_loc[s - 1 - i]) {
                    return s - 1 - i;
                }
            }
            // unreachable
            abort();
        }
    };
}

#endif //HEIPROMAP_DISTANCE_ORACLE_H
