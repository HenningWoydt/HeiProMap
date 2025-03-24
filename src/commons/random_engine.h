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

#ifndef HEIPROMAP_RANDOMENGINE_H
#define HEIPROMAP_RANDOMENGINE_H

#include <random>

#include "../commons/definitions.h"

namespace HeiProMap {
    class RandomEngine {
    public:
        std::uniform_int_distribution<int> dis_int;
        std::uniform_real_distribution<float> dis_float;
        std::mt19937 gen;

        RandomEngine() = default;

        explicit RandomEngine(const u64 t_seed) {
            gen.seed(t_seed);
            dis_float = std::uniform_real_distribution<f32>(0.0f, 1.0f);
            dis_int = std::uniform_int_distribution<int>(std::numeric_limits<int>::min(), std::numeric_limits<int>::max());
        }

        f32 get_f32() { return dis_float(gen); }

        int get_int() { return dis_int(gen); }
    };
}

#endif //HEIPROMAP_RANDOMENGINE_H
