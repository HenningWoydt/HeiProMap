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

#include <csignal>
#include <cstring>
#include <iostream>
#include <vector>

#include "../src/solver/DeepHeiProMap_solver.h"
#include "../src/configuration/DeepHeiProMap_configuration.h"
#include "src/distance_oracles/distance_oracle.h"
#include "../src/distance_oracles/binary_distance_oracle.h"
#include "src/utility/profiler.h"
#include "src/utility/utils.h"

static void signal_handler(int sig) {
    std::cout << "\n[Deep-HeiProMap] Caught signal " << sig << ", printing profiler before exit:\n";
    HeiProMap::Profiler::instance().print_table_ascii_colored(std::cout);
    std::exit(128 + sig);
}

int main(const int argc, char *argv[]) {
    std::signal(SIGTERM, signal_handler);
    std::signal(SIGINT, signal_handler);
    std::signal(SIGHUP, signal_handler);
    std::signal(SIGABRT, signal_handler);

    auto sp = HeiProMap::get_time_point();

    int run_argc = argc;
    char **run_argv = argv;

    std::vector<std::string> default_args_storage;
    std::vector<char *> default_argv_ptrs;

    if (argc == 1) {
        HEIPROMAP_PROFILE_SCOPE("io", "main", "read_args");
        std::vector<std::pair<std::string, std::string> > input = {
            // {"--graph", "../../ProMapRepo/data/mapping/europe_osm.graph"},
            {"--graph", "../../ProMapRepo/data/mapping_backup/rgg27.graph"},
            {"--mapping", "../data/out/partition/europe_osm.txt"},
            {"--statistics", "../data/out/statistics/europe_osm.JSON"},
            {"--hierarchy", "32:10:10:32"},
            {"--distance", "1:10:50:100"},
            {"--imbalance", "0.03"},
            {"--config", "fast"},
            {"--threads", "2"},
            {"--seed", "5"},
            {"--distance-oracle", "binary-based"},
        };

        default_args_storage.push_back("DeepHeiProMap");
        for (const auto &[key, val]: input) {
            default_args_storage.push_back(key);
            default_args_storage.push_back(val);
        }

        for (auto &arg: default_args_storage) {
            default_argv_ptrs.push_back(arg.data());
        }

        run_argc = static_cast<int>(default_argv_ptrs.size());
        run_argv = default_argv_ptrs.data();
    }

    HeiProMap::DeepHeiProMapConfiguration ac(run_argc, run_argv);

    auto t_read_start = HeiProMap::get_time_point();
    HeiProMap::CSRGraph g(ac.graph_in);
    auto t_read_end = HeiProMap::get_time_point();
    HeiProMap::f64 read_io_ms = HeiProMap::get_milli_seconds(t_read_start, t_read_end);

    std::vector<HeiProMap::partition_t> partition;

    if (ac.use_binary_oracle() || ac.distance_oracle_algorithm_id == HeiProMap::DISTANCE_ORACLE_ALGS_BINARY) {
        if (ac.k >= 1024) {
            HeiProMap::DeepHeiProMapSolver<true, HeiProMap::BinaryDistanceOracle, HeiProMap::LargeQuotientGraph> solver(std::move(g), ac);
            partition = solver.solve();
        } else {
            HeiProMap::DeepHeiProMapSolver<false, HeiProMap::BinaryDistanceOracle, HeiProMap::LargeQuotientGraph> solver(std::move(g), ac);
            partition = solver.solve();
        }
    } else if (ac.distance_oracle_algorithm_id == HeiProMap::DISTANCE_ORACLE_ALGS_DIVISION || ac.distance_oracle_algorithm_id == HeiProMap::DISTANCE_ORACLE_ALGS_STORE_DIVISION) {
        if (ac.k >= 1024) {
            HeiProMap::DeepHeiProMapSolver<true, HeiProMap::DistanceOracle, HeiProMap::LargeQuotientGraph> solver(std::move(g), ac);
            partition = solver.solve();
        } else {
            HeiProMap::DeepHeiProMapSolver<false, HeiProMap::DistanceOracle, HeiProMap::LargeQuotientGraph> solver(std::move(g), ac);
            partition = solver.solve();
        }
    } else {
        std::cerr << "Distance Oracle Algorithm not recognized : " << ac.coarsening_algorithm_id << " " << ac.distance_oracle_algorithm_string << std::endl;
        std::exit(EXIT_FAILURE);
    }

    auto t_write_start = HeiProMap::get_time_point();
    HeiProMap::write_partition(partition, ac.mapping_out);
    auto t_write_end = HeiProMap::get_time_point();
    HeiProMap::f64 write_io_ms = HeiProMap::get_milli_seconds(t_write_start, t_write_end);

    std::cout << "IO Read Time (ms)     : " << read_io_ms << std::endl;
    std::cout << "IO Write Time (ms)    : " << write_io_ms << std::endl;

    auto ep = HeiProMap::get_time_point();
    std::cout << "Total Time in Deep-HeiProMap.cpp (s): " << HeiProMap::get_seconds(sp, ep) << std::endl;

    HeiProMap::Profiler::instance().print_table_ascii_colored(std::cout);

    return 0;
}
