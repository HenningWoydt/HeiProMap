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

#ifndef HEIPROMAP_IPARALLELPARTITIONMANAGER_H
#define HEIPROMAP_IPARALLELPARTITIONMANAGER_H

#include <vector>

#include "../../definitions.h"

namespace HeiProMap {

    class IParallelPartitionManager {
    public:
        virtual ~IParallelPartitionManager() = default;
        virtual void initialize(vertex_t n, partition_t k, weight_t t_lmax) = 0;
        virtual const partition_t& operator[](vertex_t u) const = 0;
        virtual void set(vertex_t u, weight_t w, partition_t id) = 0;
        virtual void move(vertex_t u, weight_t w, partition_t old_id, partition_t new_id) = 0;
        virtual weight_t get_bweight(partition_t id) const = 0;
        virtual std::vector<weight_t> get_bweights() const = 0;
        virtual void uncontract(const EdgeUV* matches, size_t &matches_size) = 0;
        virtual bool is_overloaded() = 0;
    };

}

#endif //HEIPROMAP_IPARALLELPARTITIONMANAGER_H
