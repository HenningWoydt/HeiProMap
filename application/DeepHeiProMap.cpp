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

#include "../../src/commons/definitions.h"
#include "../src/commons/utils.h"
#include "../src/deep/datastructures/deep_solver.h"
#include "../src/serial/utility/algorithm_configuration.h"

using namespace HeiProMap;

int main(const int argc, char* argv[]) {
    if (argc == 1) {
        {
            std::vector<std::pair<std::string, std::string>> input = {
                    // {"--graph", "../data/mapping/2cubes_sphere.mtx.graph"},
                    // {"--mapping", "../data/out/partition/2cubes_sphere.txt"},
                    // {"--statistics", "../data/out/statistics/2cubes_sphere.JSON"},
                    // {"--graph", "../data/dimacs10_matrix/af_shell10.graph"},
                    // {"--mapping", "../data/out/partition/af_shell10.txt"},
                    // {"--statistics", "../data/out/statistics/af_shell10.JSON"},
                    // {"--graph", "../data/training/PGPgiantcompo.graph"},
                    // {"--mapping", "../data/out/partition/PGPgiantcompo.txt"},
                    // {"--statistics", "../data/out/statistics/PGPgiantcompo.JSON"},
                    {"--graph", "../data/dimacs10_delaunay/delaunay_n22.graph"}, // To Beat 2715456 in 220.94 s
                    {"--mapping", "../data/out/partition/delaunay_n22.txt"},
                    {"--statistics", "../data/out/statistics/delaunay_n22.JSON"},
                    // {"--graph", "../data/mapping/del26.graph"},
                    // {"--mapping", "../data/out/partition/del26.txt"},
                    // {"--statistics", "../data/out/statistics/del26.JSON"},
                    // {"--graph", "../data/training/598a.graph"},
                    // {"--mapping", "../data/out/partition/598a.txt"},
                    // {"--statistics", "../data/out/statistics/598a.JSON"},
                    // {"--graph", "../data/training/rgg_n26.graph"},
                    // {"--mapping", "../data/out/partition/rgg_n26.txt"},
                    // {"--statistics", "../data/out/statistics/rgg_n26.JSON"},
                    // {"--graph", "../data/training/G3_circuit.graph"},
                    // {"--mapping", "../data/out/partition/G3_circuit.txt"},
                    // {"--statistics", "../data/out/statistics/G3_circuit.JSON"},
                    // {"--graph", "../data/dimacs10_random/rgg_n_2_15_s0.graph"}, // To beat 207196 in 0.29
                    // {"--mapping", "../data/out/partition/rgg_n_2_15_s0.txt"},
                    // {"--statistics", "../data/out/statistics/rgg_n_2_15_s0.JSON"},
                    {"--hierarchy", "32:16:8:8"},
                    {"--distance", "1:10:50:100"},
                    // {"--hierarchy", "8:8:16"},
                    // {"--distance", "1:10:100"},
                    {"--imbalance", "0.03"},
                    {"--config", "experimental"},
                    {"--threads", "4"},
                    {"--seed", "0"},

                    // coarsening
                    // {"--coarsening-algorithm", "heavy-matching"},
                    // {"--coarsening-algorithm", "greedy-matching"},
                    {"--coarsening-algorithm", "global-paths"},
                    // {"--coarsening-algorithm", "random-matching"},

                    // coarsening - greedy configuration
                    {"--coarsening-algorithm-greedy-matching-pendant-first", "0"},

                    // coarsening - greedy configuration
                    {"--coarsening-algorithm-heavy-matching-pendant-first", "0"},

                    // Partitioning
                    //{"--partitioning-algorithm", "kaffpa-multisection"},
                    // {"--partitioning-algorithm", "multisection"},

                    // {"--partitioning-algorithm-multisection-mode", "fast"},

                    // Rebalancing
                    {"--rebalancing-algorithm", "simple"},

                    // Refinement Label Propagation Faraj20
                    // {"--refinement-lable-propagation-faraj20-enable", "0"},
                    // {"--refinement-lable-propagation-faraj20-max-iterations", "25"},

                    // Refinement Quotient graph Faraj20
                    // {"--refinement-quotient-graph-faraj20-enable", "0"},
                    //{"--refinement-quotient-graph-faraj20-max-iterations", "1"},
                    //{"--refinement-quotient-graph-faraj20-max-moves-without-max", "500"},

                    // Refinement k-way FM Faraj20
                    // {"--refinement-k-way-fm-faraj20-enable", "0"},
                    //{"--refinement-k-way-fm-faraj20-max-iterations", "1"},

                    // Refinement Multi-Try FM Faraj20
                    // {"--refinement-multi-try-fm-faraj20-enable", "0"},
                    // {"--refinement-multi-try-fm-faraj20-max-iterations", "1"},

                    // Refinement Label Propagation
                    // {"--refinement-lable-propagation-enable", "0"},
                    // {"--refinement-lable-propagation-max-iterations", "25"},

                    // Refinement Quotient graph
                    // {"--refinement-quotient-graph-enable", "0"},
                    //{"--refinement-quotient-graph-max-iterations", "1"},
                    //{"--refinement-quotient-graph-max-moves-without-max", "500"},

                    // Refinement k-way FM
                    // {"--refinement-k-way-fm-enable", "0"},
                    //{"--refinement-k-way-fm-max-iterations", "1"},

                    // Refinement Multi-Try FM
                    // {"--refinement-multi-try-fm-enable", "0"},
                    // {"--refinement-multi-try-fm-max-iterations", "1"},

                    // Refinement Hierarchy aware cycles
                    // {"--refinement-hierarchy-aware-cycles-enable", "0"},
                };

            std::vector<std::string> args = {"DeepHeiProMap"};
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

            DeepAlgorithmConfiguration ac(argc, argv);

            DeepSolver solver(ac);
            solver.solve();

            for (int i = 0; i < argc; ++i) { delete[] argv[i]; }
            delete[] argv;
        }
    } else {
        DeepAlgorithmConfiguration ac(argc, argv);

        DeepSolver solver(ac);
        solver.solve();
    }

    return 0;
}
