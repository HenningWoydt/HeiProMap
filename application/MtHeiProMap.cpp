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

#include <iostream>
#include <chrono>
#include <vector>
#include <string>
#include <cstring>
#include "../src/parallel/utility/parallel_algorithm_configuration.h"
#include "../src/parallel/datastructures/parallel_solver.h"

using namespace HeiProMap;

int main(const int argc, char *argv[]) {
    if (argc == 1) {
        const auto sp = std::chrono::high_resolution_clock::now();
        {
            std::vector<std::pair<std::string, std::string>> input = {
                    // {"--graph", "../data/mapping/2cubes_sphere.mtx.graph"},
                    // {"--mapping", "../data/out/partition/2cubes_sphere.txt"},
                    // {"--statistics", "../data/out/statistics/2cubes_sphere.JSON"},
                    {"--graph",                                                      "../data/dimacs10_delaunay/delaunay_n22.graph"},
                    {"--mapping",                                                    "../data/out/partition/delaunay_n22.txt"},
                    {"--statistics",                                                 "../data/out/statistics/delaunay_n22.JSON"},
                    // {"--graph",                                                      "../data/dimacs10_delaunay/delaunay_n24.graph"},
                    // {"--mapping",                                                    "../data/out/partition/delaunay_n24.txt"},
                    // {"--statistics",                                                 "../data/out/statistics/delaunay_n24.JSON"},
                    // {"--graph", "../data/training/598a.graph"},
                    // {"--mapping", "../data/out/partition/598a.txt"},
                    // {"--statistics", "../data/out/statistics/598a.JSON"},
                    // {"--graph",                                                      "../data/training/rgg_n26.graph"},
                    // {"--mapping",                                                    "../data/out/partition/rgg_n26.txt"},
                    // {"--statistics",                                                 "../data/out/statistics/rgg_n26.JSON"},
                    // {"--graph",                         "../data/training/G3_circuit.graph"},
                    // {"--mapping",                       "../data/out/partition/G3_circuit.txt"},
                    // {"--statistics",                    "../data/out/statistics/G3_circuit.JSON"},
                    // {"--graph", "../data/dimacs10_random/rgg_n_2_15_s0.graph"},
                    // {"--mapping", "../data/out/partition/rgg_n_2_15_s0.txt"},
                    // {"--statistics", "../data/out/statistics/rgg_n_2_15_s0.JSON"},
                    {"--hierarchy",                                                  "4:8:6"},
                    {"--distance",                                                   "1:10:100"},
                    {"--imbalance",                                                  "0.03"},
                    {"--threads",                                                    "16"},
                    {"--config",                                                     "experimental"},
                    {"--seed",                                                       "0"},

                    // coarsening
                    {"--parallel-coarsening-algorithm",                              "heavy-matching"},

                    // coarsening - heavy configuration
                    {"--parallel-coarsening-algorithm-heavy-matching-pendant-first", "1"},

                    // Partitioning
                    {"--parallel-partitioning-algorithm",                            "kaffpa"},
                    // {"--partitioning-algorithm", "multisection"},

                    // partitioning - kaffpa configuration
                    {"--parallel-partitioning-algorithm-method",                     "bisection"},
                    {"--parallel-partitioning-algorithm-mode",                       "fast"},
            };

            std::vector<std::string> args = {"HeiProMap"};
            for (const auto &[key, val]: input) {
                args.push_back(key);
                args.push_back(val);
            }

            // Step 3: Prepare argc and argv.
            int argc = args.size();
            if (argc < 0) {
                std::cerr << "Error: Invalid argc size" << std::endl;
                exit(EXIT_FAILURE);
            }

            // Allocate an array of char* for argv.
            char **argv = new char *[argc];

            for (int i = 0; i < argc; ++i) {
                // Allocate enough space for the string plus the null terminator.
                argv[i] = new char[args[i].size() + 1];
                std::strcpy(argv[i], args[i].c_str());
            }

            ParallelAlgorithmConfiguration ac(argc, argv);

            ParallelSolver parallel_solver(ac);
            parallel_solver.solve();

            for (int i = 0; i < argc; ++i) { delete[] argv[i]; }
            delete[] argv;
        }
        const auto ep = std::chrono::high_resolution_clock::now();
        std::cout << "Total time: " << get_seconds(sp, ep) << std::endl;
    } else {
        ParallelAlgorithmConfiguration ac(argc, argv);

        ParallelSolver parallel_solver(ac);
        parallel_solver.solve();
    }

    return 0;
}