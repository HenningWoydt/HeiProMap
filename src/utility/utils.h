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
#include <cmath>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstring>
#include <filesystem>

#include "../definitions.h"

namespace HeiProMap {
    /**
     * Splits a string into multiple sub-strings. The specified character will
     * serve as the delimiter and will not be present in any string.
     *
     * @param str The string.
     * @param c The character.
     * @return Vector of sub-strings.
     */
    inline std::vector<std::string> split(const std::string &str,
                                   char c) {
        std::vector<std::string> splits;
        std::istringstream iss(str);
        std::string token;
        while (std::getline(iss, token, c)) {
            splits.push_back(token);
        }
        return splits;
    }


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
        T result;
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

    inline void line_to_ints(const std::string &line, std::vector<u64> &ints) {
        ints.resize(line.size());
        u64 idx = 0;
        u64 curr_number = 0;
        for (char c: line) {
            if (c == ' ') {
                ints[idx] = curr_number;
                idx += curr_number != 0;
                curr_number = 0;
            } else {
                curr_number = curr_number * 10 + (c - '0');
            }
        }
        ints[idx] = curr_number;
        idx += curr_number != 0;
        ints.resize(idx);
    }

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
    f64 avg(const std::vector<T> &vec) {
        T m = 0;

        for (auto &x: vec) {
            m += x;
        }

        return static_cast<f64>(m) / static_cast<f64>(vec.size());
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
    bool exists(T *const &arr, const size_t size, const T &x) {
        for (size_t i = 0; i < size; ++i) {
            if (arr[i] == x) {
                return true;
            }
        }
        return false;
    }

    template<typename T>
    void print(const std::vector<T> &vec) {
        std::cout << "[";
        for (size_t i = 0; i < vec.size(); ++i) {
            std::cout << +vec[i];
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

    inline bool file_exists(const std::string &path) {
        std::ifstream f(path.c_str());
        return f.good();
    }

    inline std::string read_file(const std::string &path) {
        std::ifstream t(path);
        std::stringstream buffer;
        buffer << t.rdbuf();
        return buffer.str();
    }

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

    inline std::string to_string(const std::vector<EdgeVW> &vec) {
        std::string s;
        if (vec.empty()) {
            s = "[]";
            return s;
        }
        s = "[";
        for (size_t i = 0; i < vec.size() - 1; ++i) {
            s += "(" + std::to_string(vec[i].v) + ", " + std::to_string(vec[i].w) + "), ";
        }
        s += "(" + std::to_string(vec.back().v) + ", " + std::to_string(vec.back().w) + ")]";
        return s;
    }

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

    inline void busyFunction(float duration) {
        auto start = std::chrono::high_resolution_clock::now();
        auto end = start;
        volatile float uselessResult = 0.0f;

        // Continue running until the specified duration has passed
        while (std::chrono::duration<float>(end - start).count() < duration) {
            // Perform some meaningless calculations
            for (int i = 0; i < 1000; ++i) {
                uselessResult += std::sqrt(static_cast<float>(i)) * std::sqrt(static_cast<float>(i + 1));
                uselessResult -= std::sqrt(static_cast<float>(i + 1)) * std::sqrt(static_cast<float>(i));
            }
            // Update the end time
            end = std::chrono::high_resolution_clock::now();
        }
    }

    // Function to trim leading and trailing spaces
    inline std::string trim(std::string str) {
        if (str.empty()) {
            return str;
        }

        // Find the first non-space character
        size_t start = 0;
        while (start < str.size() && std::isspace(str[start])) {
            ++start;
        }

        // Find the last non-space character
        size_t end = str.size() - 1;
        while (end > start && std::isspace(str[end])) {
            --end;
        }

        // Resize and move the string to contain only the trimmed part
        return str.substr(start, end - start + 1);
    }

    inline bool startsWith(const std::string &s, const std::string &start) {
        if (s.size() < start.size()) return false;
        return s.compare(0, start.size(), start) == 0;
    }

    inline bool endsWith(const std::string &s, const std::string &end) {
        if (s.size() < end.size()) return false;
        return s.compare(s.size() - end.size(), end.size(), end) == 0;
    }

    inline void locked_print(std::mutex &lock, const std::string &s) {
        lock.lock();
        std::cout << s << std::endl;
        lock.unlock();
    }

    inline f64 get_seconds(std::chrono::high_resolution_clock::time_point sp, std::chrono::high_resolution_clock::time_point ep) {
        return (f64) std::chrono::duration_cast<std::chrono::nanoseconds>(ep - sp).count() / 1e9;
    }

    inline f64 get_seconds(std::chrono::steady_clock::time_point sp, std::chrono::steady_clock::time_point ep) {
        return (f64) std::chrono::duration_cast<std::chrono::nanoseconds>(ep - sp).count() / 1e9;
    }

    inline void counting_sort(std::vector<EdgeUVW> &edges, std::vector<EdgeUVW> &edges_help, std::vector<weight_t> &help, weight_t min_w, weight_t max_w) {
        if (edges.empty()) return;

        // Step 2: Create the count array
        help.resize(max_w - min_w + 1);
        std::fill(help.begin(), help.end(), 0);

        // Step 3: Count the occurrences of each weight
        for (const auto &edge: edges) {
            help[edge.w - min_w]++;
        }

        // Step 4: Accumulate the counts
        for (size_t i = 1; i < help.size(); ++i) {
            help[i] += help[i - 1];
        }

        // Step 5: Build the output array
        edges_help.resize(edges.size());
        for (int i = edges.size() - 1; i >= 0; --i) {
            const EdgeUVW &edge = edges[i];
            edges_help[--help[edge.w - min_w]] = edge;
        }

        // Step 6: Copy the sorted array back to the original array
        edges.swap(edges_help);
        std::reverse(edges.begin(), edges.end());
    }

    inline size_t own_lower_bound_guaranteed(const std::vector<EdgeVW> &edges, vertex_t v) {
        if (edges[0].v == v) { return 0; }
        if (edges[1].v == v) { return 1; }
        if (edges[2].v == v) { return 2; }
        if (edges[3].v == v) { return 3; }
        if (edges[4].v == v) { return 4; }
        if (edges[5].v == v) { return 5; }
        if (edges[6].v == v) { return 6; }
        if (edges[7].v == v) { return 7; }

        size_t left = 8;
        size_t right = edges.size();

        while (left < right) {
            size_t mid = left + (right - left) / 2;

            if (edges[mid].v < v) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }

        // At this point, left should be the position of the vertex v.
        return left;
    }

    inline size_t own_lower_bound_not_guaranteed(const std::vector<EdgeVW> &edges, vertex_t v) {
        size_t left = 0;
        size_t right = edges.size();

        while (left < right) {
            size_t mid = left + (right - left) / 2;

            if (edges[mid].v < v) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }

        // After the loop, left is the position where `v` would go if it exists.
        // It could be the first element that is not less than `v`.
        return left;
    }

    inline void write_partition(const std::vector<partition_t> &partition, const std::string &file_path) {
        std::ofstream out(file_path, std::ios::binary); // Open file in binary mode for faster writing
        if (!out.is_open()) return; // Ensure the file is open before proceeding

        // Write each partition element directly to the file, separating them by newlines
        for (const auto &i: partition) {
            out << i << '\n'; // Write each element followed by a newline
        }

        out.close();
    }

    inline void read_partition(const std::string &file_path, std::vector<partition_t> &p) {
        partition_t x;
        std::fstream in(file_path);
        vertex_t u = 0;
        while (in >> x) {
            p[u++] = x;
        }
    }

    inline void move_while(const char *arr, size_t &i, const char &x, size_t size) {
        while (i < size && arr[i] == x) {
            ++i;
        }
    }

    inline void move_while_not(const char *arr, size_t &i, const char &x, size_t size) {
        while (i < size - 4) {
            if (arr[i] == x) { return; }
            if (arr[i + 1] == x) {
                i += 1;
                return;
            }
            if (arr[i + 2] == x) {
                i += 2;
                return;
            }
            if (arr[i + 3] == x) {
                i += 3;
                return;
            }
            i += 4;
        }

        while (i < size && arr[i] != x) {
            ++i;
        }
    }

    inline void str_to_ints(const std::string &str,
                     std::vector<u64> &ints) {
        ints.resize(str.size());

        u64 idx = 0;
        u64 curr_number = 0;

        for (const char c: str) {
            if (c == ' ') {
                ints[idx] = curr_number;
                idx += curr_number != 0;
                curr_number = 0;
            } else {
                curr_number = curr_number * 10 + (c - '0');
            }
        }

        ints[idx] = curr_number;
        idx += curr_number != 0;
        ints.resize(idx);
    }

    inline size_t round_up_64(std::size_t x) {
        return (x + 63) & ~static_cast<std::size_t>(63);
    }

    template<typename T>
    std::vector<T> diff(const std::vector<T> &vec1, const std::vector<T> &vec2) {
        // Make local copies so we can sort them.
        std::vector<T> sorted1 = vec1;
        std::vector<T> sorted2 = vec2;

        std::sort(sorted1.begin(), sorted1.end());
        std::sort(sorted2.begin(), sorted2.end());

        std::vector<T> difference;

        // Compute the symmetric difference: elements in sorted1 or sorted2 but not in both.
        std::set_symmetric_difference(
            sorted1.begin(), sorted1.end(),
            sorted2.begin(), sorted2.end(),
            std::back_inserter(difference)
        );

        return difference;
    }

    inline size_t get_memory_usage_kb() {
        std::ifstream status_file("/proc/self/status");
        std::string line;

        while (std::getline(status_file, line)) {
            if (line.rfind("VmRSS:", 0) == 0) {
                // Resident Set Size
                size_t kb;
                sscanf(line.c_str(), "VmRSS: %zu kB", &kb);
                return kb;
            }
        }
        return 0; // fallback
    }

    inline double get_memory_usage_gb() {
        return (double) get_memory_usage_kb() / 1024.0 / 1024.0;
    }

    inline constexpr u64 bitsNeeded64(u64 x) {
        if (x == 0) return 1;
        return 64u - static_cast<unsigned>(__builtin_clzll(x));
    }

    // Generic for any unsigned type
    template<typename T>
    std::string toBinary(T value) {
        static_assert(std::is_unsigned<T>::value, "Use unsigned types for bit printing");

        constexpr size_t bits = std::numeric_limits<T>::digits; // number of value bits
        std::string out;
        out.reserve(bits);

        for (size_t i = 0; i < bits; ++i) {
            // Check bit from MSB (bits-1) down to LSB (0)
            size_t shift = bits - 1 - i;
            out.push_back((value & (T{1} << shift)) ? '1' : '0');
        }
        return out;
    }

    inline auto get_time_point() {
        return std::chrono::high_resolution_clock::now();
    }

    inline f64 get_milli_seconds(std::chrono::high_resolution_clock::time_point sp, std::chrono::high_resolution_clock::time_point ep) {
        return (f64) std::chrono::duration_cast<std::chrono::nanoseconds>(ep - sp).count() / 1e6;
    }

    inline std::size_t floor_log2(std::size_t x) noexcept {
        if (x == 0) return 0;
#if defined(__GNUC__) || defined(__clang__)
        return static_cast<std::size_t>(8 * sizeof(unsigned long long) - 1 - __builtin_clzll(static_cast<unsigned long long>(x)));
#elif defined(_MSC_VER)
        unsigned long index;
#   if defined(_WIN64)
        _BitScanReverse64(&index, x);
#   else
        _BitScanReverse(&index, static_cast<unsigned long>(x));
#   endif
        return static_cast<std::size_t>(index);
#else
        // Fallback portable loop
        std::size_t res = 0;
        while ((std::size_t(1) << (res + 1)) <= x) ++res;
        return res;
#endif
    }

    // Suggested shape of your helper
    struct MMap {
        char *data = nullptr;
        size_t size = 0;
        int fd = -1; // keep fd so you can close it later
    };

    inline MMap mmap_file_ro(const std::string &path) {
        MMap mm;

        int fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
        if (fd < 0) {
            perror("open");
            std::exit(EXIT_FAILURE);
        }

        struct stat st{};
        if (fstat(fd, &st) != 0) {
            perror("fstat");
            std::exit(EXIT_FAILURE);
        }
        size_t size = static_cast<size_t>(st.st_size);

#ifdef __linux__
        // 1) Tell the kernel we’ll read sequentially (before mmap)
        (void) posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif

        void *addr = ::mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0);
        if (addr == MAP_FAILED) {
            perror("mmap");
            std::exit(EXIT_FAILURE);
        }

#ifdef __linux__
        // 2) Hint that we’ll need these pages, sequentially (right after mmap)
        (void) madvise(addr, size, MADV_SEQUENTIAL | MADV_WILLNEED);
#endif

        mm.data = static_cast<char *>(addr);
        mm.size = size;
        mm.fd = fd; // store; close in your munmap_file(...)
        return mm;
    }

    inline void munmap_file(const MMap &mm) {
        if (mm.data && mm.size) ::munmap(mm.data, mm.size);
        if (mm.fd >= 0) ::close(mm.fd);
    }

    template<class T, class URBG>
    inline void fast_shuffle_unchecked(T *first, T *last, URBG &gen) {
        std::size_t n = last - first;
        std::size_t swaps = n / 8;

        std::uniform_int_distribution<std::size_t> dist(0, n - 1);

        for (std::size_t i = 0; i < swaps; ++i)
            std::swap(first[dist(gen)], first[dist(gen)]);
    }

    static inline uint64_t splitmix64(uint64_t x) {
        x += 0x9E3779B97F4A7C15ull;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
        return x ^ (x >> 31);
    }

    // Hash a vertex id (optionally salted by map_u)
    static inline uint64_t hash_vertex(vertex_t v, vertex_t salt = 0) {
        uint64_t x = static_cast<uint64_t>(v) ^ (static_cast<uint64_t>(salt) * 0x9E3779B97F4A7C15ull);
        return splitmix64(x);
    }
}

#endif //HEIPROMAP_UTILS_H
