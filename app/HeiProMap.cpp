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

#include "src/datastructures/HeiProMap_solver.h"

int main(const int argc, char *argv[]) {
    auto sp = HeiProMap::get_time_point();

    if (argc == 1) {
        {
            HEIPROMAP_PROFILE_SCOPE("io", "main", "read_args");
            std::vector<std::pair<std::string, std::string> > input = {
                // {"--graph", "../../ProMapRepo/data/mapping/2cubes_sphere.mtx.graph"}, // fast 0.753s 8,704,035 comm cost // eco 2.559s 7,469,493 // strong 22.498s 7,223,358
                // {"--mapping", "../data/out/partition/2cubes_sphere.mtx.txt"},
                {"--graph", "../../ProMapRepo/data/mapping/del23.graph"}, // fast 9.72s, 4,285,098.0 comm cost, // eco 17.53s ,3,889,700.0 comm cost, // strong 63.94s 3,838,630.0 comm cost
                {"--mapping", "../data/out/partition/del23.txt"},
                // {"--graph", "../../ProMapRepo/data/mapping/afshell9.graph"}, // fast 1.180s, 11,526,672 comm cost, // eco 3.635s 10,133,094 comm cost, // strong 61.415s 9,595,692 comm cost
                // {"--mapping", "../data/out/partition/afshell9.txt"},
                // {"--graph", "../../ProMapRepo/data/mapping/nlr.graph"}, // fast 6.01s, 3,668,266 comm cost, // eco 24.56s 3,150,316 comm cost, // strong 283.43s 3,047,453 comm cost
                // {"--mapping", "../data/out/partition/nlr.txt"},
                // {"--graph", "../../ProMapRepo/data/mapping/Spielman_k500.graph"}, // fast 6.01s, 3,668,266 comm cost, // eco 24.56s 3,150,316 comm cost, // strong 283.43s 3,047,453 comm cost
                // {"--mapping", "../data/out/partition/Spielman_k500.txt"},
                // {"--graph", "../../ProMapRepo/data/mapping/cont-300.mtx.graph"}, // fast 6.01s, 3,668,266 comm cost, // eco 24.56s 3,150,316 comm cost, // strong 283.43s 3,047,453 comm cost
                // {"--mapping", "../data/out/partition/cont-300.mtx.graph.txt"},
                // {"--graph", "../../ProMapRepo/data/mapping/G2_circuit.mtx.graph"}, // fast 6.01s, 3,668,266 comm cost, // eco 24.56s 3,150,316 comm cost, // strong 283.43s 3,047,453 comm cost
                // {"--mapping", "../data/out/partition/G2_circuit.mtx.graph.txt"},
                {"--hierarchy", "4:8:6"},
                {"--distance", "1:10:100"},
                {"--imbalance", "0.03"},
                {"--config", "strong"},
                {"--seed", "0"},
                {"--threads", "16"},
                {"--hm-level", "0"},
            };

            std::vector<std::string> args = {"HeiProMap"};
            for (const auto &[key, val]: input) {
                args.push_back(key);
                args.push_back(val);
            }

            // Step 3: Prepare argc and argv.
            int temp_argc = (int) args.size();
            if (temp_argc < 0) {
                std::cerr << "Error: Invalid argc size" << std::endl;
                exit(EXIT_FAILURE);
            }

            // Allocate an array of char* for argv.
            char **temp_argv = new char *[temp_argc];

            for (int i = 0; i < temp_argc; ++i) {
                // Allocate enough space for the string plus the null terminator.
                temp_argv[i] = new char[args[i].size() + 1];
                std::strcpy(temp_argv[i], args[i].c_str());
            }

            HeiProMap::AlgorithmConfiguration ac(temp_argc, temp_argv);

            HeiProMap::HeiProMapSolver solver(ac);
            solver.solve();

            for (int i = 0; i < temp_argc; ++i) { delete[] temp_argv[i]; }
            delete[] temp_argv;

            HeiProMap::Profiler::instance().print_table_ascii_colored(std::cout);
        }
    } else {
        HeiProMap::AlgorithmConfiguration ac(argc, argv);

        std::cout << "Seed " << ac.seed << std::endl;

        HeiProMap::HeiProMapSolver solver(ac);
        solver.solve();
    }

    auto ep = HeiProMap::get_time_point();
    std::cout << "Total Time in HeiProMap.cpp: " << HeiProMap::get_seconds(sp, ep) << std::endl;

    return 0;
}
