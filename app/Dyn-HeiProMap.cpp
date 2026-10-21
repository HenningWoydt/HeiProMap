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

#include <algorithm>
#include <csignal>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "../src/solver/DynHeiProMap_solver.h"
#include "../src/configuration/DynHeiProMap_configuration.h"

static void signal_handler(int sig) {
    std::cout << "\n[Dyn-HeiProMap] Caught signal " << sig << ", printing stack trace before exit:\n";
    HeiProMap::print_stack_trace();
    HeiProMap::Profiler::instance().print_table_ascii_colored(std::cout);
    std::exit(128 + sig);
}

int main(int argc, char **argv) {
    std::signal(SIGSEGV, signal_handler);
    std::signal(SIGTERM, signal_handler);
    std::signal(SIGINT, signal_handler);
    std::signal(SIGHUP, signal_handler);
    std::signal(SIGABRT, signal_handler);

    std::ios::sync_with_stdio(false);
    std::cout.tie(nullptr);

    if (argc > 1) {
        HeiProMap::DynConfiguration config(argc, argv);
        HeiProMap::DynHeiProMapSolver solver(config);
        HeiProMap::interactive_mode(solver);
    } else {
        HeiProMap::DynConfiguration config;
        config.hierarchy = {4, 8, 6};
        config.distance = {1, 10, 100};
        config.imbalance = 0.03;
        config.n_threads = 1;
        config.seed = 1;

        HeiProMap::DynHeiProMapSolver solver(config);

        // Replay the advection_cellular frames in order: the first frame is partitioned
        // from scratch, every following frame is applied as a batch and refined.
        std::vector<std::string> frames;
        for (const auto &entry: std::filesystem::directory_iterator("../data/advection_cellular")) {
            if (entry.path().extension() == ".my_seq") {
                frames.push_back(entry.path().string());
            }
        }
        std::sort(frames.begin(), frames.end());

        std::vector<std::string> commands;
        for (size_t i = 0; i < frames.size(); ++i) {
            commands.emplace_back("exec " + frames[i]);
            if (i == 0) {
                commands.emplace_back("partition fast");
            } else {
                commands.emplace_back("autorefine 40");
                commands.emplace_back("stats");
            }
        }

        HeiProMap::execute_commands(solver, commands);
    }

    HeiProMap::Profiler::instance().print_table_ascii_colored(std::cout);
    return 0;
}
