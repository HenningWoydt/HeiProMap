/*******************************************************************************
 * MIT License
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
 #include <string>
 #include <vector>

 #include "src/datastructures/DynHeiProMap_solver.h"
#include "src/utility/DynHeiProMap_configuration.h"

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cout.tie(nullptr);

    HeiProMap::DynConfiguration config;
    if (argc > 1) {
        config = HeiProMap::DynConfiguration(argc, argv);
    } else {
        HEIPROMAP_PROFILE_SCOPE("io", "main", "parse_args");
        std::vector<std::pair<std::string, std::string> > input = {
            {"--hierarchy", "4:8:6"},
            {"--distance", "1:10:100"},
            {"--imbalance", "0.03"},
            {"--threads", "1"},
            {"--seed", "1"}
        };

        std::vector<std::string> args = {"Dyn-HeiProMap"};
        for (const auto &[key, val]: input) {
            args.push_back(key);
            args.push_back(val);
        }

        // Step 3: Prepare argc and argv.
        int argc_temp = (int) args.size();
        if (argc_temp < 0) {
            std::cerr << "Error: Invalid argc size" << std::endl;
            exit(EXIT_FAILURE);
        }

        // Allocate an array of char* for argv.
        char **argv_temp = new char *[(size_t) argc_temp];

        for (size_t i = 0; i < (size_t) argc_temp; ++i) {
            // Allocate enough space for the string plus the null terminator.
            argv_temp[i] = new char[args[i].size() + 1];
            std::strcpy(argv_temp[i], args[i].c_str());
        }

        config = HeiProMap::DynConfiguration(argc_temp, argv_temp);

        for (int i = 0; i < argc_temp; ++i) { delete[] argv_temp[i]; }
        delete[] argv_temp;

        HeiProMap::DynSolver solver(config);

        std::vector<std::string> commands;
        commands.emplace_back("exec ../../Dyn-HeiProMap-Experiments/data/batches_100/2cubes_sphere.mtx/batch_0.my_seq");
        commands.emplace_back("partition fast");
        commands.emplace_back("save_partition mapping_0.txt");
        for (size_t i = 1; i < 20; ++i) {
            commands.emplace_back("exec ../../Dyn-HeiProMap-Experiments/data/batches_100/2cubes_sphere.mtx/batch_" + std::to_string(i) + ".my_seq");
            commands.emplace_back("autorefine 40");
            if (i % 2 == 1) {
                // commands.emplace_back("refine-fast 5");
            } else {
                // commands.emplace_back("partition fast");
            }
            commands.emplace_back("stats");
            commands.emplace_back("save_partition mapping_" + std::to_string(i) + ".txt");
        }

        // HeiProMap::interactive_mode(solver);
        HeiProMap::execute_commands(solver, commands);

        HeiProMap::Profiler::instance().print_table_ascii_colored(std::cout);

        return 0;
    }

    HeiProMap::DynSolver solver(config);

    HeiProMap::interactive_mode(solver);

    return 0;
}
