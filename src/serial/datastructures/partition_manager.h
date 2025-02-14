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

#include "../../definitions.h"
#include "../interfaces/ISerialPartitionManager.h"

namespace HeiProMap {
    class PartitionManager final : public ISerialPartitionManager {
        vertex_t m_n    = 0;
        partition_t m_k = 0;
        weight_t lmax   = 0;

        std::vector<partition_t> partition;
        std::vector<weight_t> bweights;

    public:
        void initialize(const vertex_t t_n,
                        const partition_t t_k,
                        const weight_t t_lmax) override {
            m_n = t_n;
            m_k = t_k;
            lmax = t_lmax;

            partition.resize(m_n);
            bweights.resize(m_k, 0);
        }

        // read
        const partition_t& operator[](const vertex_t u) const override { return partition[u]; }

        // write
        void set(const vertex_t u, const weight_t w, const partition_t id) override {
            bweights[id] += w;
            partition[u] = id;
        }

        void move(const vertex_t u, const weight_t w, const partition_t old_id, const partition_t new_id) override {
            bweights[old_id] -= w;
            bweights[new_id] += w;
            partition[u] = new_id;
        }

        weight_t get_bweight(const partition_t id) const override { return bweights[id]; }
        std::vector<weight_t> get_bweights() const override { return bweights; }
        void uncontract(const vertex_t u, const vertex_t v) override { partition[v] = partition[u]; }

        void uncontract(const std::vector<EdgeUV>& matches) override {
            for (const auto [u, v] : matches) {
                partition[v] = partition[u];
            }
        }

        bool is_overloaded() override {
            for (const weight_t w : bweights) {
                if (w > lmax) {
                    return true;
                }
            }
            return false;
        }
    };
}

#endif //HEIPROMAP_PARTITION_MANAGER_H
