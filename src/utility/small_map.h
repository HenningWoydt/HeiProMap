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

#ifndef HEIPROMAP_SMALL_MAP_H
#define HEIPROMAP_SMALL_MAP_H

#include <vector>
#include <utility>
#include <cstddef>

namespace HeiProMap {
    template<typename Key, typename Val>
    class FlatMap {
    public:
        using value_type = std::pair<Key, Val>;
        using iterator = typename std::vector<value_type>::iterator;
        using const_iterator = typename std::vector<value_type>::const_iterator;

        FlatMap() = default;

        explicit FlatMap(size_t capacity) { reserve(capacity); }

        void reserve(size_t n) { data_.reserve(n); }

        void clear() noexcept { data_.clear(); }

        bool empty() const noexcept { return data_.empty(); }

        size_t size() const noexcept { return data_.size(); }

        // --- lookup and insert ----------------------------------------------------
        Val &operator[](const Key &k) {
            // linear find — works best for small sets
            for (auto &kv: data_) {
                if (kv.first == k) return kv.second;
            }
            data_.emplace_back(k, Val{});
            return data_.back().second;
        }

        // --- iteration ------------------------------------------------------------
        iterator begin() noexcept { return data_.begin(); }
        iterator end() noexcept { return data_.end(); }
        const_iterator begin() const noexcept { return data_.begin(); }
        const_iterator end() const noexcept { return data_.end(); }

    private:
        std::vector<value_type> data_;
    };
}

#endif //HEIPROMAP_SMALL_MAP_H
