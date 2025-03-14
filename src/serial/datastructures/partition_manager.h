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
#include "../interfaces/ISerialPartitionManager.h"

namespace HeiProMap {
    class PartitionManager final : public ISerialPartitionManager {
        vertex_t m_n    = 0;
        partition_t m_k = 0;
        weight_t lmax   = 0;

        partition_t* partition      = nullptr;
        partition_t* partition_temp = nullptr;
        weight_t* bweights          = nullptr;

    public:
        PartitionManager() = default;

        ~PartitionManager() override {
            free(partition);
            free(partition_temp);
            free(bweights);
        }

        void initialize(const vertex_t t_n,
                        const partition_t t_k,
                        const weight_t t_lmax) override {
            vertex_t t_n_64 = round_up_64(t_n);
            vertex_t t_k_64 = round_up_64(t_k);

            m_n  = t_n;
            m_k  = t_k;
            lmax = t_lmax;

            partition      = (partition_t*)aligned_alloc(64, t_n_64 * sizeof(partition_t));
            partition_temp = (partition_t*)aligned_alloc(64, t_n_64 * sizeof(partition_t));
            bweights       = (weight_t*)aligned_alloc(64, t_k_64 * sizeof(weight_t));
            std::fill_n(bweights, t_k_64, 0);
        }

        // read
        const partition_t& operator[](const vertex_t u) const override { return partition[u]; }

        // write
        void set(const vertex_t u,
                 const weight_t w,
                 const partition_t id) override {
            bweights[id] += w;
            partition[u] = id;
        }

        void move(const vertex_t u,
                  const weight_t w,
                  const partition_t old_id,
                  const partition_t new_id) override {
            bweights[old_id] -= w;
            bweights[new_id] += w;
            partition[u] = new_id;
        }

        weight_t get_bweight(const partition_t id) const override { return bweights[id]; }


        std::vector<weight_t> get_bweights() const override {
            std::vector<weight_t> weights(m_k);
            for (size_t i = 0; i < m_k; ++i) {
                weights[i] = bweights[i];
            }
            return weights;
        }

        void uncontract(const Matching& matching) override {
            for (vertex_t new_u = 0; new_u < matching.get_n_coarse_nodes(); ++new_u) {
                vertex_t old_u           = matching.get_o(new_u);
                vertex_t old_u_partner   = matching.get_partner(old_u);
                partition_temp[old_u]         = partition[new_u];
                partition_temp[old_u_partner] = partition[new_u];
            }
            std::swap(partition, partition_temp);
        }

        bool is_overloaded() override {
            for (size_t i = 0; i < m_k; ++i) {
                if (bweights[i] > lmax) { return true; }
            }
            return false;
        }
    };
}

#endif //HEIPROMAP_PARTITION_MANAGER_H
