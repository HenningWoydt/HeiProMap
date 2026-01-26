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

#include "../definitions.h"
#include "../utility/aligned_array.h"

namespace HeiProMap {
    class PartitionManager {
        vertex_t m_n = 0;
        partition_t m_k = 0;
        weight_t lmax = 0;

        AlignedArray<partition_t> partition;
        AlignedArray<partition_t> partition_temp;
        AlignedArray<weight_t> bweights;
        AlignedArray<size_t> n_vertices;

    public:
        void initialize(const vertex_t t_n,
                        const partition_t t_k,
                        const weight_t t_lmax) {
            ScopedTimer _t("io", "PartitionManager", "initialize");

            m_n = t_n;
            m_k = t_k;
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
            for (size_t i = 0; i < m_k; ++i) {
                weights[i] = bweights[i];
            }
            return weights;
        }

        weight_t max_weight() const {
            weight_t m = 0;
            for (size_t i = 0; i < m_k; ++i) {
                m = std::max(m, bweights[i]);
            }
            return m;
        }

        partition_t n_empty_blocks() const {
            partition_t n = 0;
            for (size_t i = 0; i < m_k; ++i) {
                n += bweights[i] == 0;
            }
            return n;
        }

        partition_t n_oload_blocks() const {
            partition_t n = 0;
            for (size_t i = 0; i < m_k; ++i) {
                n += bweights[i] > lmax;
            }
            return n;
        }

        weight_t sum_oload_weight() const {
            weight_t w = 0;
            for (size_t i = 0; i < m_k; ++i) {
                w += std::max((weight_t) 0 , bweights[i] - lmax);
            }
            return w;
        }

        void contract(const Mapping &mapping) {
            ScopedTimer _t("contraction", "PartitionManager", "contract");

            for (vertex_t u = 0; u < mapping.get_old_n(); ++u) {
                vertex_t map_u = mapping.get(u);
                partition_temp[map_u] = partition[u];
            }
            std::swap(partition, partition_temp);

            n_vertices.initialize(m_k, 0);
            for (vertex_t u = 0; u < mapping.get_coarse_n(); ++u) {
                n_vertices[partition[u]] += 1;
            }
        }

        void uncontract(const Mapping &mapping) {
            ScopedTimer _t("uncontraction", "PartitionManager", "uncontract");
            for (vertex_t u = 0; u < mapping.get_old_n(); ++u) {
                vertex_t map_u = mapping.get(u);
                partition_temp[u] = partition[map_u];
            }
            std::swap(partition, partition_temp);

            n_vertices.initialize(m_k, 0);
            for (vertex_t u = 0; u < mapping.get_old_n(); ++u) {
                n_vertices[partition[u]] += 1;
            }
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
    };
}

#endif //HEIPROMAP_PARTITION_MANAGER_H
