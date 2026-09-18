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
#include <vector>

#include "../src/solver/Deep_HeiProMap_solver.h"
#include "../src/configuration/Deep_HeiProMap_configuration.h"
#include "src/datastructures/distance_oracle.h"
#include "src/datastructures/binary_distance_oracle.h"
#include "src/utility/profiler.h"
#include "src/utility/utils.h"

int main(const int argc, char *argv[]) {
    auto sp = HeiProMap::get_time_point();

    if (argc == 1) {
        {
            HEIPROMAP_PROFILE_SCOPE("io", "main", "read_args");
            std::vector<std::pair<std::string, std::string> > input = {
                {"--graph", "../../ProMapRepo/data/mapping/europe_osm.graph"},
                {"--mapping", "../data/out/partition/europe_osm.txt"},
                {"--statistics", "../data/out/statistics/europe_osm.JSON"},
                // {"--hierarchy", "8:2"},
                // {"--distance", "1:10"},
                {"--hierarchy", "32:10:10:32"},
                {"--distance", "1:10:50:100"},
                {"--imbalance", "0.03"},
                {"--config", "fast"},
                {"--threads", "16"},
                {"--seed", "5"},
                {"--distance-oracle", "binary-based"}
            };

            std::vector<std::string> args = {"DeepHeiProMap"};
            for (const auto &[key, val]: input) {
                args.push_back(key);
                args.push_back(val);
            }

            int temp_argc = (int) args.size();
            if (temp_argc < 0) {
                std::cerr << "Error: Invalid argc size" << std::endl;
                exit(EXIT_FAILURE);
            }

            char **temp_argv = new char *[temp_argc];
            for (int i = 0; i < temp_argc; ++i) {
                temp_argv[i] = new char[args[i].size() + 1];
                std::strcpy(temp_argv[i], args[i].c_str());
            }

            HeiProMap::DeepHeiProMapConfiguration ac(temp_argc, temp_argv);

            if (ac.use_binary_oracle() || ac.distance_oracle_algorithm_id == HeiProMap::DISTANCE_ORACLE_ALGS_BINARY) {
                HeiProMap::DeepHeiProMapSolver<HeiProMap::BinaryDistanceOracle> solver(ac);
                solver.solve();
            } else if (ac.distance_oracle_algorithm_id == HeiProMap::DISTANCE_ORACLE_ALGS_DIVISION ||
                       ac.distance_oracle_algorithm_id == HeiProMap::DISTANCE_ORACLE_ALGS_STORE_DIVISION) {
                HeiProMap::DeepHeiProMapSolver<HeiProMap::DistanceOracle> solver(ac);
                solver.solve();
            } else {
                std::cerr << "Distance Oracle Algorithm not recognized : " << ac.distance_oracle_algorithm_string << std::endl;
                std::exit(EXIT_FAILURE);
            }

            for (int i = 0; i < temp_argc; ++i) { delete[] temp_argv[i]; }
            delete[] temp_argv;

            HeiProMap::Profiler::instance().print_table_ascii_colored(std::cout);
        }
    } else {
        HeiProMap::DeepHeiProMapConfiguration ac(argc, argv);

        if (ac.use_binary_oracle() || ac.distance_oracle_algorithm_id == HeiProMap::DISTANCE_ORACLE_ALGS_BINARY) {
            HeiProMap::DeepHeiProMapSolver<HeiProMap::BinaryDistanceOracle> solver(ac);
            solver.solve();
        } else if (ac.distance_oracle_algorithm_id == HeiProMap::DISTANCE_ORACLE_ALGS_DIVISION ||
                   ac.distance_oracle_algorithm_id == HeiProMap::DISTANCE_ORACLE_ALGS_STORE_DIVISION) {
            HeiProMap::DeepHeiProMapSolver<HeiProMap::DistanceOracle> solver(ac);
            solver.solve();
        } else {
            std::cerr << "Distance Oracle Algorithm not recognized : " << ac.coarsening_algorithm_id << " " << ac.distance_oracle_algorithm_string << std::endl;
            std::exit(EXIT_FAILURE);
        }

        HeiProMap::Profiler::instance().print_table_ascii_colored(std::cout);
    }

    auto ep = HeiProMap::get_time_point();
    std::cout << "Total Time in Deep-HeiProMap.cpp (s): " << HeiProMap::get_seconds(sp, ep) << std::endl;

    return 0;
}
