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

#ifndef HEIPROMAP_DEEP_DISTANCE_ORACLE_BINARY_H
#define HEIPROMAP_DEEP_DISTANCE_ORACLE_BINARY_H


#include "../../commons/definitions.h"
#include "../../commons/macros.h"
#include "../../commons/utils.h"

namespace HeiProMap {
    class DeepDistanceOracleBinaryBased {
        std::vector<partition_t> m_hierarchy; // O(l)
        std::vector<weight_t> m_distance; // O(l)
        size_t m_l = 0;
        partition_t m_k = 0;

        std::vector<u64> identifier; // O(k)
        std::vector<weight_t> dist_lookup; // O(64)
        std::vector<partition_t> hierarchy_lookup; // O(64)

    public:
        void initialize(const std::vector<partition_t> &t_hierarchy,
                        const std::vector<weight_t> &t_distance,
                        const u64 t_threads) {
            m_hierarchy = t_hierarchy;
            m_distance = t_distance;
            m_l = m_hierarchy.size();
            m_k = prod<partition_t>(m_hierarchy);

            std::vector<u64> bit_sizes;
            bit_sizes.reserve(m_hierarchy.size());

            u64 total_bits = 0;
            for (partition_t size: m_hierarchy) {
                const u64 w = bitsNeeded64(size - 1);
                bit_sizes.push_back(w);
                total_bits += w;
            }

            if (total_bits > 64) {
                std::cout << "Hierarchy too large! k = " << m_k << std::endl;
                print(m_hierarchy);
                std::exit(EXIT_FAILURE);
            }

            // --- MSB-first packing plan ---
            // We place the packed region in the top total_bits of the 64-bit word.
            // This makes clz(xor) equal to "number of equal leading bits" before first difference.
            const u64 msb_start = 64 - total_bits;

            identifier.resize(m_k);
            std::vector<partition_t> loc;
            loc.resize(m_hierarchy.size());

            for (partition_t id = 0; id < m_k; ++id) {
                determine_loc(id, loc);
                u64 ident = loc[m_hierarchy.size() - 1];

                for (size_t i = 1; i < m_hierarchy.size(); ++i) {
                    ident = ident << bit_sizes[m_hierarchy.size() - 1 - i];
                    ident |= loc[m_hierarchy.size() - 1 - i];
                }

                // print(m_hierarchy);
                // print(bit_sizes);
                // print(loc);
                // std::cout << id << std::endl;
                // std::cout << toBinary(ident) << std::endl;
                identifier[id] = ident; // high bits [0..msb_start-1] are zero
            }

            dist_lookup.assign(64, m_distance.back());
            hierarchy_lookup.assign(64, m_hierarchy.size() - 1);
            size_t idx = 0;
            for (size_t i = 0; i < m_hierarchy.size(); ++i) {
                for (size_t j = 0; j < bit_sizes[i]; ++j) {
                    dist_lookup[63 - idx] = m_distance[i];
                    hierarchy_lookup[63 - idx] = i;
                    idx += 1;
                }
            }
        }

        weight_t get(const partition_t u_id, const partition_t v_id) const {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);
            if (u_id == v_id) return 0;

            const u64 x = identifier[u_id] ^ identifier[v_id];
            const int lz = __builtin_clzll(x); // 0..64, x!=0 here

            // std::cout << u_id << " " << v_id << std::endl;
            // std::cout << toBinary(identifier[u_id]) << std::endl;
            // std::cout << toBinary(identifier[v_id]) << std::endl;
            // std::cout << toBinary(x) << std::endl;
            // std::cout << lz << " " << dist_lookup[lz] << std::endl;

            return dist_lookup[lz];
        }

        partition_t get_h(const partition_t u_id, const partition_t v_id) const {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);
            if (u_id == v_id) return 0;

            const u64 x = identifier[u_id] ^ identifier[v_id];
            const int lz = __builtin_clzll(x);
            return hierarchy_lookup[static_cast<size_t>(lz)];
        }

        bool last_level_pair(const partition_t u_id,
                             const partition_t v_id) const {
            return (u_id / m_hierarchy[0]) == (v_id / m_hierarchy[0]);
        }

    private:
        void determine_loc(partition_t u_id,
                           std::vector<partition_t> &loc) const {
            ASSERT(prod<u64>(m_hierarchy) == m_k);
            ASSERT(u_id < m_k);
            ASSERT(!m_hierarchy.empty());

            u64 r_start = 0;
            u64 r_end = m_k;

            u64 s = m_hierarchy.size();
            for (u64 i = 0; i < m_hierarchy.size(); ++i) {
                u64 n_parts = m_hierarchy[s - 1 - i];
                u64 add = (r_end - r_start) / n_parts;

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
    };
}


#endif //HEIPROMAP_DEEP_DISTANCE_ORACLE_BINARY_H
