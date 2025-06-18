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

#ifndef HEIPROMAP_JSON_UTILS_H
#define HEIPROMAP_JSON_UTILS_H

#include <map>
#include <string>
#include <vector>

#include "definitions.h"

namespace HeiProMap {
    struct JSONString {
        std::string s;
    };

#define to_JSON_MACRO(x) (std::string("\"") + (#x) + "\" : " + to_JSON_value(x) + ",\n")

    std::string to_JSON_value(u8 x);

    std::string to_JSON_value(u16 x);

    std::string to_JSON_value(u32 x);

    std::string to_JSON_value(u64 x);

    std::string to_JSON_value(s8 x);

    std::string to_JSON_value(s16 x);

    std::string to_JSON_value(s32 x);

    std::string to_JSON_value(s64 x);

    std::string to_JSON_value(f32 x);

    std::string to_JSON_value(f64 x);

    std::string to_JSON_value(const std::string &s);

    std::string to_JSON_value(const JSONString &s);

    template<typename T1, typename T2>
    std::string to_JSON_value(const std::map<T1, T2> &m) {
        if (m.empty()) {
            return "{}";
        }
        std::vector<T1> arg1;
        std::vector<T2> arg2;
        for (auto const &x: m) {
            arg1.push_back(x.first);
            arg2.push_back(x.second);
        }

        if (arg1.size() == 1) {
            return "{" + to_JSON_value(arg1[0]) + " : " + to_JSON_value(arg2[0]) + "}";
        }
        std::string s = "{";
        for (size_t i = 0; i < arg1.size() - 1; ++i) {
            s += to_JSON_value(arg1[i]) + " : " + to_JSON_value(arg2[i]) + ", ";
        }
        s += to_JSON_value(arg1[arg1.size() - 1]) + " : " + to_JSON_value(arg2[arg1.size() - 1]) + "}";

        return s;
    }

    template<typename T>
    std::string to_JSON_value(const std::vector<T> &vec) {
        if (vec.empty()) {
            return "[]";
        }
        if (vec.size() == 1) {
            return "[" + to_JSON_value(vec[0]) + "]";
        }
        std::string s = "[";
        for (size_t i = 0; i < vec.size() - 1; ++i) {
            s += to_JSON_value(vec[i]) + ", ";
        }
        s += to_JSON_value(vec.back()) + "]";

        return s;
    }

    template<typename T>
    std::string to_JSON_value(const std::pair<T, T> &p) {
        std::string s = to_JSON_value(p.first);
        s += " : " + to_JSON_value(p.second);
        return s;
    }
}

#endif //HEIPROMAP_JSON_UTILS_H
