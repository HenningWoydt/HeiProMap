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

#ifndef HEIPROMAP_PARTITION_MANAGER_H
#define HEIPROMAP_PARTITION_MANAGER_H

#include "../../commons/definitions.h"

namespace HeiProMap {
    class PartitionManager {
        vertex_t    m_n  = 0;
        partition_t m_k  = 0;
        weight_t    lmax = 0;

        AlignedArray<partition_t> partition;
        AlignedArray<partition_t> partition_temp;
        AlignedArray<partition_t> bweights;
        AlignedArray<size_t>      n_vertices;

    public:
        PartitionManager() = default;

        ~PartitionManager() = default;

        void initialize(const vertex_t t_n,
                        const partition_t t_k,
                        const weight_t t_lmax) {
            m_n  = t_n;
            m_k  = t_k;
            lmax = t_lmax;

            partition.initialize(m_n, 0);
            partition_temp.initialize(m_n);
            bweights.initialize(m_k, 0);
            n_vertices.initialize(m_k, 0);
            n_vertices[0] = t_n;
        }

        // read
        const partition_t &operator[](const vertex_t u) const { return partition[u]; }

        // write
        void set(const vertex_t u,
                 const weight_t w,
                 const partition_t id) {
            n_vertices[id] += 1;
            bweights[id] += w;
            partition[u] = id;
        }

        void move(const vertex_t u,
                  const weight_t w,
                  const partition_t old_id,
                  const partition_t new_id) {
            n_vertices[old_id] -= 1;
            n_vertices[new_id] += 1;
            bweights[old_id] -= w;
            bweights[new_id] += w;
            partition[u] = new_id;
        }

        weight_t get_bweight(const partition_t id) const { return bweights[id]; }

        size_t size(const partition_t id) const { return n_vertices[id]; }

        std::vector<weight_t> get_bweights() const {
            std::vector<weight_t> weights(m_k);
            for (size_t           i = 0; i < m_k; ++i) {
                weights[i] = bweights[i];
            }
            return weights;
        }

        void contract(const Matching &matching) {
            for (size_t i = 0; i < matching.size(); ++i) {
                auto     [u, v]         = matching[i];
                vertex_t smaller_vertex = std::min(u, v);
                vertex_t new_u          = matching.get_n(smaller_vertex);
                partition_temp[new_u] = partition[smaller_vertex];
                n_vertices[partition[smaller_vertex]] -= 1;
            }
            std::swap(partition, partition_temp);
        }

        void uncontract(const Matching &matching) {
            for (vertex_t new_u = 0; new_u < matching.get_n_coarse_nodes(); ++new_u) {
                vertex_t old_u         = matching.get_o(new_u);
                vertex_t old_u_partner = matching.get_partner(old_u);
                partition_temp[old_u]         = partition[new_u];
                partition_temp[old_u_partner] = partition[new_u];
                if (old_u != old_u_partner) {
                    n_vertices[partition[new_u]] += 1;
                }
            }
            std::swap(partition, partition_temp);
        }

        bool is_overloaded() {
            for (size_t i = 0; i < m_k; ++i) {
                if (bweights[i] > lmax) { return true; }
            }
            return false;
        }

        void reset_weights() {
            bweights.initialize(m_k, 0);
            n_vertices.initialize(m_k, 0);
        }

        void set_lmax(const partition_t id, const weight_t new_lmax) {}

        weight_t get_lmax(const partition_t id) { return 0; }

        void set_hierarchy_level(const partition_t id, const partition_t level) {}

        partition_t get_hierarchy_level(const partition_t id) { return 0; }
    };
}

#endif //HEIPROMAP_PARTITION_MANAGER_H
