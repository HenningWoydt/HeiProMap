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

#include "utils.h"


namespace HeiProMap {
    std::vector<std::string> split(const std::string &str,
                                   char c) {
        std::vector<std::string> splits;

        std::istringstream iss(str);
        std::string token;

        while (std::getline(iss, token, c)) {
            splits.push_back(token);
        }

        return splits;
    }

    bool file_exists(const std::string &path) {
        std::ifstream f(path.c_str());
        return f.good();
    }

    std::string read_file(const std::string &path) {
        std::ifstream t(path);
        std::stringstream buffer;
        buffer << t.rdbuf();
        return buffer.str();
    }

    void line_to_ints(const std::string &line, std::vector<u64> &ints) {
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

    void busyFunction(float duration) {
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

    std::string trim(std::string str) {
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

    bool startsWith(const std::string &s, const std::string &start) {
        if (s.size() < start.size()) return false;
        return s.compare(0, start.size(), start) == 0;
    }

    bool endsWith(const std::string &s, const std::string &end) {
        if (s.size() < end.size()) return false;
        return s.compare(s.size() - end.size(), end.size(), end) == 0;
    }

    void locked_print(std::mutex &lock, const std::string &s) {
        lock.lock();
        std::cout << s << std::endl;
        lock.unlock();
    }

    f64 get_seconds(std::chrono::high_resolution_clock::time_point sp, std::chrono::high_resolution_clock::time_point ep) {
        return (f64) std::chrono::duration_cast<std::chrono::nanoseconds>(ep - sp).count() / 1e9;
    }

    f64 get_seconds(std::chrono::steady_clock::time_point sp, std::chrono::steady_clock::time_point ep) {
        return (f64) std::chrono::duration_cast<std::chrono::nanoseconds>(ep - sp).count() / 1e9;
    }

    std::string to_string(const std::vector<EdgeVW> &vec) {
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

    void counting_sort(std::vector<EdgeUVW> &edges, std::vector<EdgeUVW> &edges_help, std::vector<weight_t> &help, weight_t min_w, weight_t max_w) {
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

    size_t own_lower_bound_guaranteed(const std::vector<EdgeVW> &edges, vertex_t v) {
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

    size_t own_lower_bound_not_guaranteed(const std::vector<EdgeVW> &edges, vertex_t v) {
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

    void write_partition(const std::vector<partition_t> &partition, const std::string &file_path) {
        std::ofstream out(file_path, std::ios::binary); // Open file in binary mode for faster writing
        if (!out.is_open()) return; // Ensure the file is open before proceeding

        // Write each partition element directly to the file, separating them by newlines
        for (const auto &i: partition) {
            out << i << '\n'; // Write each element followed by a newline
        }

        out.close();
    }

    void read_partition(const std::string &file_path, std::vector<partition_t> &p) {
        partition_t x;
        std::fstream in(file_path);
        vertex_t u = 0;
        while (in >> x) {
            p[u++] = x;
        }
    }

    void move_while(const char *arr, size_t &i, const char &x, size_t size) {
        while (i < size && arr[i] == x) {
            ++i;
        }
    }

    void move_while_not(const char *arr, size_t &i, const char &x, size_t size) {
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

    void str_to_ints(const std::string &str,
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

    size_t round_up_64(std::size_t x) {
        return (x + 63) & ~static_cast<std::size_t>(63);
    }
}
