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

#ifndef HEIPROMAP_PARALLEL_DISTANCE_ORACLE_H
#define HEIPROMAP_PARALLEL_DISTANCE_ORACLE_H

#include "../../definitions.h"
#include "../../macros.h"
#include "../../serial/utility/utils.h"
#include "../../interfaces/IDistanceOracle.h"
#include "../interfaces/IParallelDistanceOracle.h"

namespace HeiProMap {

    class ParallelDistanceOracle : public IParallelDistanceOracle {
    private:
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;
        partition_t m_k = 0;

        u64 m_n_threads = 1;

        std::vector<weight_t> m_mtx;
        std::vector<partition_t> m_h_mtx;

    public:

        void initialize(std::vector<partition_t> &t_hierarchy,
                        std::vector<weight_t> &t_distance,
                        u64 n_threads) final {
            m_hierarchy = t_hierarchy;
            m_distance = t_distance;
            m_k = prod<partition_t>(m_hierarchy);

            m_n_threads = n_threads;

            m_mtx.resize(m_k * m_k);
            m_h_mtx.resize(m_k * m_k);

            std::vector<std::vector<partition_t>> locs(m_k, std::vector<partition_t>(m_hierarchy.size()));
// #pragma omp parallel for default(none) shared(m_mtx, m_h_mtx, locs) num_threads(n_threads)
            for (partition_t u_id = 0; u_id < m_k; ++u_id) {
                determine_loc(u_id, locs[u_id]);
                m_mtx[u_id * m_k + u_id] = 0;
                m_h_mtx[u_id * m_k + u_id] = 0;
            }

#pragma omp parallel for default(none) shared(m_mtx, m_h_mtx, m_k, locs) num_threads(m_n_threads)
            for (partition_t u_id = 0; u_id < m_k; ++u_id) {
                for (partition_t v_id = u_id + 1; v_id < m_k; ++v_id) {
                    weight_t d = determine_distance(locs[u_id], locs[v_id]);
                    m_mtx[u_id * m_k + v_id] = d;
                    m_mtx[v_id * m_k + u_id] = d;

                    partition_t h = determine_hierarchy(locs[u_id], locs[v_id]);
                    m_h_mtx[u_id * m_k + v_id] = h;
                    m_h_mtx[v_id * m_k + u_id] = h;
                }
            }
        }

        weight_t get(partition_t u_id, partition_t v_id) const final {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);
            ASSERT(u_id * m_k + v_id < m_k * m_k);
            return m_mtx[u_id * m_k + v_id];
        }

        partition_t get_h(partition_t u_id, partition_t v_id) const final {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);
            ASSERT(u_id * m_k + v_id < m_k * m_k);
            return m_h_mtx[u_id * m_k + v_id];
        }

    private:
        void determine_loc(partition_t u_id,
                           std::vector<partition_t> &u_loc) {
            ASSERT(prod<u64>(m_hierarchy) == m_k);
            ASSERT(u_id < m_k);
            ASSERT(!m_hierarchy.empty());
            ASSERT(u_loc.size() == m_hierarchy.size());

            u64 r_start = 0;
            u64 r_end = m_k;

            u64 s = m_hierarchy.size();
            for (u64 i = 0; i < m_hierarchy.size(); ++i) {
                u64 n_parts = m_hierarchy[s - 1 - i];
                u64 add = (r_end - r_start) / n_parts;

                for (u64 j = 0; j < n_parts; ++j) {
                    if (r_start <= u_id && u_id < r_start + add) {
                        // we have found the part of the partition
                        u_loc[s - 1 - i] = j;
                        r_end = r_start + add;
                        break;
                    } else {
                        r_start += add;
                    }
                }
            }
        }

        weight_t determine_distance(std::vector<partition_t> &u_loc,
                                    std::vector<partition_t> &v_loc) {
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

        partition_t determine_hierarchy(std::vector<partition_t> &u_loc,
                                        std::vector<partition_t> &v_loc) {
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

#endif //HEIPROMAP_PARALLEL_DISTANCE_ORACLE_H
