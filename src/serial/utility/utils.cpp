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

    // Function to trim leading and trailing spaces
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

    void locked_print(std::mutex &lock, std::string s) {
        lock.lock();
        std::cout << s << std::endl;
        lock.unlock();
    }

    f64 get_seconds(std::chrono::high_resolution_clock::time_point sp, std::chrono::high_resolution_clock::time_point ep) {
        return (f64) std::chrono::duration_cast<std::chrono::nanoseconds>(ep - sp).count() / 1e9;;
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

    size_t own_lower_bound_guaranteed(const SmallVector<EdgeVW> &edges, vertex_t v) {
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

    size_t own_lower_bound_not_guaranteed(const SmallVector<EdgeVW> &edges, vertex_t v) {
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
        std::ofstream out(file_path, std::ios::binary);  // Open file in binary mode for faster writing
        if (!out.is_open()) return;                      // Ensure the file is open before proceeding

        // Write each partition element directly to the file, separating them by newlines
        for (const auto &i : partition) {
            out << i << '\n';                            // Write each element followed by a newline
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

    void move_while_not_avx2(const char* arr, size_t& i, char x, size_t size) {
        const size_t simd_width = 32;  // AVX2 processes 32 bytes at a time
        __m256i x_vec = _mm256_set1_epi8(x);  // Set all bytes in the vector to x

        /*
        __m256i block0, block1, block2, block3, block4, block5, block6, block7;
        __m256i cmp0, cmp1, cmp2, cmp3, cmp4, cmp5, cmp6, cmp7;
        int mask0, mask1,mask2,mask3,mask4,mask5,mask6,mask7;

        // Process the array in chunks of simd_width
        for (; i + simd_width*8 <= size; i += simd_width*8) {
            block0 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(arr + i + 0*simd_width));
            block1 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(arr + i + 1*simd_width));
            block2 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(arr + i + 2*simd_width));
            block3 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(arr + i + 3*simd_width));
            block4 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(arr + i + 4*simd_width));
            block5 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(arr + i + 5*simd_width));
            block6 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(arr + i + 6*simd_width));
            block7 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(arr + i + 7*simd_width));

            cmp0 = _mm256_cmpeq_epi8(block0, x_vec);
            cmp1 = _mm256_cmpeq_epi8(block1, x_vec);
            cmp2 = _mm256_cmpeq_epi8(block2, x_vec);
            cmp3 = _mm256_cmpeq_epi8(block3, x_vec);
            cmp4 = _mm256_cmpeq_epi8(block4, x_vec);
            cmp5 = _mm256_cmpeq_epi8(block5, x_vec);
            cmp6 = _mm256_cmpeq_epi8(block6, x_vec);
            cmp7 = _mm256_cmpeq_epi8(block7, x_vec);

            mask0 = _mm256_movemask_epi8(cmp0);
            mask1 = _mm256_movemask_epi8(cmp1);
            mask2 = _mm256_movemask_epi8(cmp2);
            mask3 = _mm256_movemask_epi8(cmp3);
            mask4 = _mm256_movemask_epi8(cmp4);
            mask5 = _mm256_movemask_epi8(cmp5);
            mask6 = _mm256_movemask_epi8(cmp6);
            mask7 = _mm256_movemask_epi8(cmp7);

            if (mask0 != 0) {i += 0*simd_width + __builtin_ctz(mask0); return; }
            if (mask1 != 0) {i += 1*simd_width + __builtin_ctz(mask1); return; }
            if (mask2 != 0) {i += 2*simd_width + __builtin_ctz(mask2); return; }
            if (mask3 != 0) {i += 3*simd_width + __builtin_ctz(mask3); return; }
            if (mask4 != 0) {i += 4*simd_width + __builtin_ctz(mask4); return; }
            if (mask5 != 0) {i += 5*simd_width + __builtin_ctz(mask5); return; }
            if (mask6 != 0) {i += 6*simd_width + __builtin_ctz(mask6); return; }
            if (mask7 != 0) {i += 7*simd_width + __builtin_ctz(mask7); return; }
        }
         */

        // Process the array in chunks of simd_width
        for (; i + simd_width <= size; i += simd_width) {
            __m256i block = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(arr + i));
            __m256i cmp = _mm256_cmpeq_epi8(block, x_vec);
            int mask = _mm256_movemask_epi8(cmp);
            if (mask != 0) {
                i += __builtin_ctz(mask);
                return;
            }
        }

        for (; i < size; ++i) {
            if (arr[i] == x) {
                return;
            }
        }
    }

    void move_while(const char *arr, size_t &i, const char &x, size_t size){
        while(i < size && arr[i] == x){
            ++i;
        }
    }

    void move_while_not(const char *arr, size_t &i, const char &x, size_t size){
        move_while_not_avx2(arr, i, x, size);
        return;

        while(i < size - 4){
            if(arr[i] == x){ return; }
            if(arr[i+1] == x){ i += 1; return; }
            if(arr[i+2] == x){ i += 2; return; }
            if(arr[i+3] == x){ i += 3; return; }
            i += 4;
        }

        while(i < size && arr[i] != x){
            ++i;
        }
    }

}
