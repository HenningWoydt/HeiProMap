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

#include "src/datastructures/HeiPa_solver.h"
#include "src/utility/HeiPa_configuration.h"

int main(const int argc, char *argv[]) {
    auto sp = HeiProMap::get_time_point();

    if (argc == 1) {
        {
            HEIPROMAP_PROFILE_SCOPE("io", "main", "read_args");
            std::vector<std::pair<std::string, std::string> > input = {
                {"--graph", "../../ProMapRepo/data/mapping/del23.graph"},
                {"--mapping", "../data/out/partition/del23.txt"},
                {"--k", "16"},
                {"--imbalance", "0.03"},
                {"--config", "fast"},
                {"--seed", "0"},
                {"--threads", "1"},
                // {"--partitioning-algorithm", "recursive-bisection"},
                {"--partitioning-algorithm-recursive-bisection-kappa", "5"},
                {"--partitioning-algorithm", "kaffpa"},
                {"--partitioning-algorithm-kaffpa-partitioning-mode", "strong"},
                {"--partitioning-algorithm-kaffpa-partitioning-method", "bisection"},
            };

            std::vector<std::string> args = {"HeiPa"};
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

            HeiProMap::HeiPaConfiguration ac(temp_argc, temp_argv);

            HeiProMap::HeiPaSolver solver(ac);
            solver.solve();

            for (int i = 0; i < temp_argc; ++i) { delete[] temp_argv[i]; }
            delete[] temp_argv;

            HeiProMap::Profiler::instance().print_table_ascii_colored(std::cout);
        }
    } else {
        HeiProMap::HeiPaConfiguration ac(argc, argv);

        HeiProMap::HeiPaSolver solver(ac);
        solver.solve();
    }

    auto ep = HeiProMap::get_time_point();
    std::cout << "Total Time in HeiPa.cpp: " << HeiProMap::get_seconds(sp, ep) << std::endl;

    return 0;
}
