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

#ifndef HEIPROMAP_PARALLEL_PARTITION_MANAGER_H
#define HEIPROMAP_PARALLEL_PARTITION_MANAGER_H

#include "../../definitions.h"
#include "../../macros.h"
#include "../interfaces/IParallelPartitionManager.h"

namespace HeiProMap {

    class ParallelPartitionManager final : public IParallelPartitionManager {
        vertex_t    m_n  = 0;
        partition_t m_k  = 0;
        weight_t    lmax = 0;

        partition_t *partition = nullptr;
        weight_t    *bweights  = nullptr;
        size_t      *n_nodes   = nullptr;

    public:
        ParallelPartitionManager() = default;

        ~ParallelPartitionManager() override {
            free(partition);
            free(bweights);
            free(n_nodes);
        }

        void initialize(const vertex_t t_n,
                        const partition_t t_k,
                        const weight_t t_lmax) override {
            vertex_t t_n_64 = round_up_64(t_n);
            vertex_t t_k_64 = round_up_64(t_k);

            m_n  = t_n;
            m_k  = t_k;
            lmax = t_lmax;

            partition = (partition_t *) aligned_alloc(64, t_n_64 * sizeof(partition_t));
            bweights  = (weight_t *) aligned_alloc(64, t_k_64 * sizeof(weight_t));
            n_nodes   = (size_t *) aligned_alloc(64, t_k_64 * sizeof(size_t));
            std::fill_n(bweights, t_k_64, 0);
            std::fill_n(n_nodes, t_k_64, 0);
        }

        // read
        const partition_t &operator[](const vertex_t u) const override { return partition[u]; }

        // write
        void set(const vertex_t u, const weight_t w, const partition_t id) override {
            bweights[id] += w;
            partition[u] = id;
            n_nodes[id] += 1;
        }

        void move(const vertex_t u, const weight_t w, const partition_t old_id, const partition_t new_id) override {
            bweights[old_id] -= w;
            bweights[new_id] += w;
            partition[u] = new_id;
            n_nodes[old_id] -= 1;
            n_nodes[new_id] += 1;
        }

        weight_t get_bweight(const partition_t id) const override { return bweights[id]; }

        vertex_t get_n_nodes(const partition_t id) const { return n_nodes[id]; }

        std::vector<weight_t> get_bweights() const override {
            std::vector<weight_t> weights(m_k);
            for (size_t           i = 0; i < m_k; ++i) {
                weights[i] = bweights[i];
            }
            return weights;
        }

        void uncontract(const EdgeUV *matches, size_t &matches_size) override {
            matches = ASSUME_ALIGNED(EdgeUV*, matches, 64);

            for (size_t i = 0; i < matches_size; ++i) {
                partition[matches[i].v] = partition[matches[i].u];
                n_nodes[partition[matches[i].u]] += 1;
            }
        }

        bool is_overloaded() override {
            for (size_t i = 0; i < m_k; ++i) {
                if (bweights[i] > lmax) { return true; }
            }
            return false;
        }
    };

}

#endif //HEIPROMAP_PARALLEL_PARTITION_MANAGER_H
