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

#ifndef HEIPROMAP_UTILS_H
#define HEIPROMAP_UTILS_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

#include "../../definitions.h"

namespace HeiProMap {
    /**
     * Splits a string into multiple sub-strings. The specified character will
     * serve as the delimiter and will not be present in any string.
     *
     * @param str The string.
     * @param c The character.
     * @return Vector of sub-strings.
     */
    std::vector<std::string> split(const std::string &str,
                                   char c);


    /**
     * Converts a string into the specified datatype. Conversion is done via
     * string stream and ">>" operator.
     *
     * @tparam T The desired datatype.
     * @param str The string.
     * @return The converted string.
     */
    template<typename T>
    T convert_to(const std::string &str) {
        T                  result;
        std::istringstream iss(str);
        iss >> result;
        return result;
    }

    /**
     * Converts the vector of strings into a vector of T's.
     *
     * @tparam T Type of conversion.
     * @param vec The vector.
     * @return Vector of transformed T's.
     */
    template<typename T>
    std::vector<T> convert(const std::vector<std::string> &vec) {
        std::vector<T> v;

        for (auto &s: vec) {
            v.push_back(convert_to<T>(s));
        }

        return v;
    }

    /**
     * Converts the vector of strings into a vector of T's.
     *
     * @tparam T Type of conversion.
     * @param vec The vector.
     * @return Vector of transformed T's.
     */
    template<typename T>
    std::vector<T> convert(const std::vector<std::string> &&vec) {
        std::vector<T> v;

        for (auto &s: vec) {
            v.push_back(convert_to<T>(s));
        }

        return v;
    }

    void line_to_ints(const std::string &line, std::vector<u64> &ints);

    /**
     * Multiplies all elements in the vector.
     *
     * @tparam T1 Resulting type.
     * @tparam T2 Type of vector elements.
     * @param vec The vector.
     * @return The product.
     */
    template<typename T1, typename T2>
    T1 prod(const std::vector<T2> &vec) {
        T1 p = (T1) 1;

        for (auto &x: vec) {
            p *= (T1) x;
        }

        return p;
    }

    /**
     * Sums all the elements in the vector.
     *
     * @tparam T1 Resulting type.
     * @tparam T2 Type of vector elements.
     * @param vec The vector.
     * @return The sum.
     */
    template<typename T1, typename T2>
    T1 sum(const std::vector<T2> &vec) {
        T1 s = (T1) 0;

        for (auto &x: vec) {
            s += (T1) x;
        }

        return s;
    }

    /**
     * Determines the maximum in the vector.
     *
     * @tparam T Resulting type.
     * @param vec The vector.
     * @return The maximum.
     */
    template<typename T>
    T max(const std::vector<T> &vec) {
        T m = vec[0];

        for (auto &x: vec) {
            m = std::max(m, x);
        }

        return m;
    }

    template<typename T>
    T min(const std::vector<T> &vec) {
        T m = vec[0];

        for (auto &x: vec) {
            m = std::min(m, x);
        }

        return m;
    }

    template<typename T>
    T avg(const std::vector<T> &vec) {
        T m = 0;

        for (auto &x: vec) {
            m += x;
        }

        return (f64) m / (f64) vec.size();
    }

    template<typename T>
    size_t argmin(const std::vector<T> &vec) {
        size_t idx = 0;

        for (size_t i = 1; i < vec.size(); ++i) {
            if (vec[i] < vec[idx]) {
                idx = i;
            }
        }

        return idx;
    }

    template<typename T>
    size_t argmax(const std::vector<T> &vec) {
        size_t idx = 0;

        for (size_t i = 1; i < vec.size(); ++i) {
            if (vec[i] > vec[idx]) {
                idx = i;
            }
        }

        return idx;
    }

    /**
     * Determines if an element exists in the vector. The "==" operator is used
     * to determine if elements are equal.
     *
     * @tparam T Type of element and vector elements.
     * @param vec The vector.
     * @param x The element to find.
     * @return Vector of transformed T's.
     */
    template<typename T>
    bool exists(const std::vector<T> &vec, const T &x) {
        return std::find(vec.begin(), vec.end(), x) != vec.end();
    }

    template<typename T>
    void print(const std::vector<T> &vec) {
        std::cout << "[";
        for (size_t i = 0; i < vec.size(); ++i) {
            std::cout << vec[i];
            if (i != vec.size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]" << std::endl;
    }

    /**
     * Determines if a duplicate exists. The "==" operator is used to determine
     * if elements are equal.
     *
     * @tparam T Type of element and vector elements.
     * @param vec The vector.
     * @return True if no duplicate exists, false else.
     */
    template<typename T>
    bool no_duplicates(const std::vector<T> &vec) {
        for (u64 i = 0; i < vec.size(); ++i) {
            for (u64 j = i + 1; j < vec.size(); ++j) {
                if (vec[i] == vec[j]) {
                    return false;
                }
            }
        }
        return true;
    }

    template<typename T>
    bool no_duplicates_sorted(const std::vector<T> &vec) {
        if (vec.empty()) {
            return true;
        }
        for (u64 i = 0; i < vec.size() - 1; ++i) {
            if (vec[i] == vec[i + 1]) {
                return false;
            }
        }
        return true;
    }

    /**
     * Checks if a file exists at the specified path.
     *
     * @param path The file path.
     * @return True if the file exists, false else.
     */
    bool file_exists(const std::string &path);

    std::string read_file(const std::string &path);

    template<typename T>
    std::string to_string(const std::vector<T> &vec) {
        std::string s;
        if (vec.empty()) {
            s = "[]";
            return s;
        }
        s = "[";
        for (size_t i = 0; i < vec.size() - 1; ++i) {
            s += std::to_string(vec[i]) + ", ";
        }
        s += std::to_string(vec.back()) + "]";

        return s;
    }

    std::string to_string(const std::vector<EdgeVW> &vec);

    template<typename T>
    std::string concat(const std::vector<T> &vec) {
        std::string s;
        if (vec.empty()) {
            s = "[]";
            return s;
        }
        s = "[";
        for (size_t i = 0; i < vec.size() - 1; ++i) {
            s += vec[i] + ", ";
        }
        s += vec.back() + "]";

        return s;
    }

    void busyFunction(float duration);

    // Function to trim leading and trailing spaces
    std::string trim(std::string str);

    bool startsWith(const std::string &s, const std::string &start);

    bool endsWith(const std::string &s, const std::string &end);

    void locked_print(std::mutex &lock, const std::string &s);

    f64 get_seconds(std::chrono::high_resolution_clock::time_point sp, std::chrono::high_resolution_clock::time_point ep);

    void counting_sort(std::vector<EdgeUVW> &edges, std::vector<EdgeUVW> &edges_help, std::vector<weight_t> &help, weight_t min_w, weight_t max_w);

    size_t own_lower_bound_guaranteed(const std::vector<EdgeVW> &edges, vertex_t v);

    size_t own_lower_bound_not_guaranteed(const std::vector<EdgeVW> &edges, vertex_t v);

    void write_partition(const std::vector<partition_t> &partition, const std::string &file_path);

    void read_partition(const std::string &file_path, std::vector<partition_t> &p);

    void move_while(const char *arr, size_t &i, const char &x, size_t size);

    void move_while_not(const char *arr, size_t &i, const char &x, size_t size);

    void str_to_ints(const std::string &str,
                     std::vector<u64> &ints);
}

#endif //HEIPROMAP_UTILS_H
