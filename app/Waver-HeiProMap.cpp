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

#include "../src/configuration/WaverHeiProMap_configuration.h"
#include "../src/solver/WaverHeiProMap_solver.h"
#include "../src/utility/profiler.h"
#include "../src/utility/utils.h"

static void signal_handler(int sig) {
    std::cout << "\n[Waver-HeiProMap] Caught signal " << sig << ", printing profiler before exit:\n";
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
        std::vector<std::pair<std::string, std::string>> input = {
            {"--topology", "../data/topologies/256_PEs.graph"},
            {"--graph", "../data/heipromap_input.graph"},
            {"--mapping", "../data/out/heipromap_input.txt"},
            {"--imbalance", "0.03"},
            {"--config", "fast"},
            {"--seed", "0"},
            {"--threads", "1"},
        };

        default_args_storage.push_back("Waver-HeiProMap");
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

    HeiProMap::WaverHeiProMapConfiguration ac(run_argc, run_argv);

    if (ac.topology_in.empty()) {
        std::cerr << "Topology required for Waver-HeiProMap!" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    auto t_read_start = HeiProMap::get_time_point();
    HeiProMap::CSRGraph g(ac.graph_in);
    HeiProMap::CSRGraph topology(ac.topology_in);
    auto t_read_end = HeiProMap::get_time_point();
    HeiProMap::f64 read_io_ms = HeiProMap::get_milli_seconds(t_read_start, t_read_end);

    HeiProMap::WaverHeiProMapSolver solver(std::move(g), std::move(topology), ac);
    std::vector<HeiProMap::partition_t> partition = solver.solve();

    auto t_write_start = HeiProMap::get_time_point();
    HeiProMap::write_partition(partition, ac.mapping_out);
    auto t_write_end = HeiProMap::get_time_point();
    HeiProMap::f64 write_io_ms = HeiProMap::get_milli_seconds(t_write_start, t_write_end);

    std::cout << "IO Read Time (ms)     : " << read_io_ms << std::endl;
    std::cout << "IO Write Time (ms)    : " << write_io_ms << std::endl;

    auto ep = HeiProMap::get_time_point();
    std::cout << "Total Time in Waver-HeiProMap.cpp (s): " << HeiProMap::get_seconds(sp, ep) << std::endl;

    HeiProMap::Profiler::instance().print_table_ascii_colored(std::cout);

    return 0;
}
