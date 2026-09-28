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

#include "../src/solver/HeiProMap_solver.h"

static void signal_handler(int sig) {
    std::cout << "\n[HeiProMap] Caught signal " << sig << ", printing profiler before exit:\n";
    HeiProMap::Profiler::instance().print_table_ascii_colored(std::cout);
    std::exit(128 + sig);
}

int main(const int argc, char *argv[]) {
    std::signal(SIGTERM, signal_handler);
    std::signal(SIGINT, signal_handler);
    std::signal(SIGHUP, signal_handler);
    std::signal(SIGABRT, signal_handler);

    auto sp = HeiProMap::get_time_point();

    std::vector<std::string> default_arg_strs;
    std::vector<char *> default_argv_ptrs;
    int run_argc = argc;
    char **run_argv = argv;

    if (argc == 1) {
        HEIPROMAP_PROFILE_SCOPE("io", "main", "read_args");
        std::vector<std::pair<std::string, std::string> > input = {
            {"--graph", "../../ProMapRepo/data/mapping/rgg_n_2_22_s0.graph"},
            {"--mapping", "../data/out/partition/rgg_n_2_22_s0.graph.txt"},
            {"--hierarchy", "4:8:6"},
            {"--distance", "1:10:100"},
            {"--imbalance", "0.03"},
            {"--config", "fast"},
            {"--seed", "0"},
            {"--threads", "2"},
            {"--hm-level", "0"},
        };

        default_arg_strs.push_back("HeiProMap");
        for (const auto &[key, val]: input) {
            default_arg_strs.push_back(key);
            default_arg_strs.push_back(val);
        }

        for (auto &s: default_arg_strs) {
            default_argv_ptrs.push_back(s.data());
        }

        run_argc = static_cast<int>(default_argv_ptrs.size());
        run_argv = default_argv_ptrs.data();
    }

    HeiProMap::AlgorithmConfiguration ac(run_argc, run_argv);

    auto t_read_start = HeiProMap::get_time_point();
    HeiProMap::CSRGraph g(ac.graph_in);
    auto t_read_end = HeiProMap::get_time_point();
    HeiProMap::f64 read_io_ms = HeiProMap::get_milli_seconds(t_read_start, t_read_end);

    std::vector<HeiProMap::partition_t> partition;

    if (ac.use_binary_oracle()) {
        HeiProMap::HeiProMapSolver<HeiProMap::BinaryDistanceOracle> solver(std::move(g), ac);
        partition = solver.solve();
    } else {
        HeiProMap::HeiProMapSolver<HeiProMap::DistanceOracle> solver(std::move(g), ac);
        partition = solver.solve();
    }

    auto t_write_start = HeiProMap::get_time_point();
    HeiProMap::write_partition(partition, ac.mapping_out);
    auto t_write_end = HeiProMap::get_time_point();
    HeiProMap::f64 write_io_ms = HeiProMap::get_milli_seconds(t_write_start, t_write_end);

    if (argc == 1) {
        HeiProMap::Profiler::instance().print_table_ascii_colored(std::cout);
    }

    std::cout << "IO Read Time (ms)      : " << read_io_ms << std::endl;
    std::cout << "IO Write Time (ms)     : " << write_io_ms << std::endl;

    auto ep = HeiProMap::get_time_point();
    std::cout << "Total Time in HeiProMap.cpp (s): " << HeiProMap::get_seconds(sp, ep) << std::endl;

    return 0;
}
