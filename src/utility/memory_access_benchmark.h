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

#ifndef HEIPROMAP_MEMORY_ACCESS_BENCHMARK_H
#define HEIPROMAP_MEMORY_ACCESS_BENCHMARK_H

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <numeric>
#include <random>
#include <vector>

#include <omp.h>

namespace HeiProMap {
    // Measures parallel scaling for three memory access patterns:
    //   1. Sequential: streaming sum over a contiguous array (bandwidth-bound)
    //   2. Random: sum via randomly permuted indices (latency-bound)
    //   3. Graph-like: random jump to a short neighborhood (2-8 elements),
    //      then sequential read within it (mimics CSR graph traversal)
    //
    // Compares 1 thread vs 2 threads to show hardware-level parallel efficiency
    // without memory bandwidth saturation effects at higher thread counts.
    static void run_memory_access_benchmark() {
        constexpr size_t N = 500'000'000;
        constexpr int ROUNDS = 10;

        std::vector<int32_t> data(N);
        std::vector<int32_t> indices(N);

        std::mt19937 gen(42);
        for (size_t i = 0; i < N; ++i) data[i] = static_cast<int32_t>(gen() % 1000);
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), gen);

        // Graph-like setup: random chunks of size 2-8 simulating neighborhoods
        std::vector<size_t> nb_start;
        std::vector<size_t> nb_end;
        {
            size_t pos = 0;
            while (pos < N) {
                size_t deg = 2 + (gen() % 7);
                size_t end = std::min(pos + deg, N);
                nb_start.push_back(pos);
                nb_end.push_back(end);
                pos = end;
            }
        }
        const size_t n_nodes = nb_start.size();

        std::vector<size_t> node_order(n_nodes);
        std::iota(node_order.begin(), node_order.end(), 0);
        std::shuffle(node_order.begin(), node_order.end(), gen);

        printf("%-8s | %-18s %-10s | %-18s %-10s | %-18s %-10s\n", "Threads",
               "Sequential (ms)", "Speedup", "Random (ms)", "Speedup", "Graph-like (ms)", "Speedup");
        printf("---------+-------------------------------+-------------------------------+------------------------------\n");

        double seq_time = 0, rand_time = 0, graph_time = 0;

        for (int t = 1; t <= 2; t++) {
            double total_ms1 = 0, total_ms2 = 0, total_ms3 = 0;

            for (int r = 0; r < ROUNDS; ++r) {
                // Test 1: Sequential access
                volatile int64_t sum1 = 0;
                auto t0 = std::chrono::high_resolution_clock::now();
                int64_t local_sum1 = 0;
                #pragma omp parallel for num_threads(t) reduction(+:local_sum1)
                for (size_t i = 0; i < N; ++i) {
                    local_sum1 += data[i];
                }
                sum1 = local_sum1;
                auto t1 = std::chrono::high_resolution_clock::now();
                total_ms1 += std::chrono::duration<double, std::milli>(t1 - t0).count();

                // Test 2: Random access
                volatile int64_t sum2 = 0;
                auto t2 = std::chrono::high_resolution_clock::now();
                int64_t local_sum2 = 0;
                #pragma omp parallel for num_threads(t) reduction(+:local_sum2)
                for (size_t i = 0; i < N; ++i) {
                    local_sum2 += data[indices[i]];
                }
                sum2 = local_sum2;
                auto t3 = std::chrono::high_resolution_clock::now();
                total_ms2 += std::chrono::duration<double, std::milli>(t3 - t2).count();

                // Test 3: Graph-like access
                volatile int64_t sum3 = 0;
                auto t4 = std::chrono::high_resolution_clock::now();
                int64_t local_sum3 = 0;
                #pragma omp parallel for num_threads(t) reduction(+:local_sum3)
                for (size_t i = 0; i < n_nodes; ++i) {
                    const size_t node = node_order[i];
                    for (size_t j = nb_start[node]; j < nb_end[node]; ++j) {
                        local_sum3 += data[j];
                    }
                }
                sum3 = local_sum3;
                auto t5 = std::chrono::high_resolution_clock::now();
                total_ms3 += std::chrono::duration<double, std::milli>(t5 - t4).count();
            }

            double ms1 = total_ms1 / ROUNDS;
            double ms2 = total_ms2 / ROUNDS;
            double ms3 = total_ms3 / ROUNDS;

            if (t == 1) { seq_time = ms1; rand_time = ms2; graph_time = ms3; }

            printf("%-8d | %-18.2f %-10.2fx | %-18.2f %-10.2fx | %-18.2f %-10.2fx\n",
                   t, ms1, seq_time / ms1, ms2, rand_time / ms2, ms3, graph_time / ms3);
        }
        printf("\n");
    }
}

#endif //HEIPROMAP_MEMORY_ACCESS_BENCHMARK_H
