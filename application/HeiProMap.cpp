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

#include <cstring>
#include <iostream>

#include "../src/definitions.h"
#include "../src/macros.h"
#include "../src/serial/datastructures/solver.h"
#include "../src/serial/utility/AlgorithmConfiguration.h"
#include "../src/serial/utility/utils.h"

using namespace HeiProMap;

// TODO: KaFFPa never allows a node v to participate in a contraction if the weight of v exceeds 1.5n / 20k

int main(const int argc, char* argv[]) {
    if (argc == 1) {
        const auto sp = std::chrono::high_resolution_clock::now();
        {
            std::vector<std::pair<std::string, std::string>> input = {
                    // {"--graph", "../data/mapping/2cubes_sphere.mtx.graph"},
                    // {"--mapping", "../data/out/partition/2cubes_sphere.txt"},
                    // {"--statistics", "../data/out/statistics/2cubes_sphere.JSON"},
                    {"--graph", "../data/dimacs10_delaunay/delaunay_n21.graph"},
                    {"--mapping", "../data/out/partition/delaunay_n21.txt"},
                    {"--statistics", "../data/out/statistics/delaunay_n21.JSON"},
                    // {"--graph", "../data/dimacs10_random/rgg_n_2_15_s0.graph"},
                    // {"--mapping", "../data/out/partition/rgg_n_2_15_s0.txt"},
                    // {"--statistics", "../data/out/statistics/rgg_n_2_15_s0.JSON"},
                    {"--hierarchy", "4:8:6"},
                    {"--distance", "1:10:100"},
                    {"--imbalance", "0.03"},
                    {"--config", "Faraj20-fastest"},
                    {"--seed", "0"},

                    // coarsening
                    // {"--coarsening-algorithm", "greedy-matching"},
                    {"--coarsening-algorithm", "global-paths"},

                    // coarsening - greedy configuration
                    {"--coarsening-algorithm-greedy-matching-pendant-first", "0"},
                    {"--coarsening-algorithm-greedy-matching-no-overload", "1"},

                    // Partitioning
                    //{"--partitioning-algorithm",                             "kaffpa-multisection"},
                    // {"--partitioning-algorithm", "multisection"},

                    // {"--partitioning-algorithm-multisection-mode", "fast"},

                    // Rebalancing
                    {"--rebalancing-algorithm", "simple"},

                    // Refinement
                    // {"--refinement-lable-propagation-faraj20-enable", "0"},

                    // {"--refinement-quotient-graph-faraj20-enable", "0"},
                    // {"--refinement-quotient-graph-faraj20-max-iterations", "5"},

                    // {"--refinement-k-way-fm-faraj20-enable", "0"},
                    // {"--refinement-k-way-fm-faraj20-max-iterations", "5"},

                    // {"--refinement-multi-try-fm-faraj20-enable", "0"},
                    // {"--refinement-multi-try-fm-faraj20-max-iterations", "0"},

                    // {"--refinement-hierarchy-aware-cycles-enable", "0"},
                };

            std::vector<std::string> args = {"HeiProMap"};
            for (const auto& [key, val] : input) {
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
            char** argv = new char*[argc];

            for (int i = 0; i < argc; ++i) {
                // Allocate enough space for the string plus the null terminator.
                argv[i] = new char[args[i].size() + 1];
                std::strcpy(argv[i], args[i].c_str());
            }

            AlgorithmConfiguration ac(argc, argv);

            Solver solver(ac);
            solver.solve();

            for (int i = 0; i < argc; ++i) { delete[] argv[i]; }
            delete[] argv;
        }
        const auto ep = std::chrono::high_resolution_clock::now();
        std::cout << "Total time: " << get_seconds(sp, ep) << std::endl;
    } else {
        AlgorithmConfiguration ac(argc, argv);

        Solver solver(ac);
        solver.solve();
    }

    return 0;
}
