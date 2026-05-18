/*******************************************************************************
 * MIT License
 *
 * This file is part of Dyn-HeiProMap.
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

#ifndef DYN_HEIPROMAP_CONFIGURATION_H
#define DYN_HEIPROMAP_CONFIGURATION_H

#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <sstream>

#include "../../definitions.h"

namespace HeiProMap {
    class Configuration {
        std::vector<CommandLineOption> options = {
            {"--help", "", "Produces the help message", "", "", false},
            {"--hierarchy", "-h", "Hierarchy (e.g., 2:2:2).", "", "", false},
            {"--distance", "-d", "Distance (e.g., 1:10:100).", "", "", false},
            {"--imbalance", "-e", "Allowed imbalance (e.g., 0.03).", "0.03", "", false},
            {"--threads", "-t", "Number of threads.", "1", "", false},
            {"--seed", "-s", "Seed for randomness.", "0", "", false},
        };

    public:
        std::vector<partition_t> hierarchy;
        std::vector<vertex_t> distance;
        f64 imbalance = 0.03;
        u64 n_threads = 1;
        u64 seed = 0;

        Configuration() = default;

        Configuration(int argc, char *argv[]) {
            std::vector<std::string> args(argv, argv + argc);

            for (size_t i = 1; i < (size_t) argc; ++i) {
                if (args[i] == "--help") {
                    print_help_message();
                    exit(EXIT_SUCCESS);
                }
            }

            for (size_t i = 1; i < (size_t) argc; ++i) {
                for (auto &[large_key, small_key, description, default_val, input, is_set]: options) {
                    if (large_key == args[i] || small_key == args[i]) {
                        input = args[i + 1];
                        is_set = true;
                        i += 1;
                        break;
                    }
                }
            }

            if (is_set("--hierarchy")) {
                std::string h_str = get("--hierarchy");
                std::istringstream h_iss(h_str);
                std::string token;
                while (std::getline(h_iss, token, ':')) {
                    hierarchy.push_back(std::stoull(token));
                }
            }

            if (is_set("--distance")) {
                std::string d_str = get("--distance");
                std::istringstream d_iss(d_str);
                std::string token;
                while (std::getline(d_iss, token, ':')) {
                    distance.push_back(std::stoull(token));
                }
            }

            imbalance = std::stod(get("--imbalance"));
            n_threads = std::stoull(get("--threads"));

            if (is_set("--seed")) {
                seed = std::stoull(get("--seed"));
            } else {
                seed = std::random_device{}();
            }
        }

        /**
         * Returns whether the option was entered.
         */
        bool is_set(const std::string &var) {
            for (const auto &[large_key, small_key, description, default_val, input, is_set]: options) {
                if (large_key == var || small_key == var) {
                    return is_set;
                }
            }
            return false;
        }

        /**
         * Gets the entered input as a string.
         */
        std::string get(const std::string &var) {
            for (const auto &[large_key, small_key, description, default_val, input, is_set]: options) {
                if (large_key == var || small_key == var) {
                    if (input.empty() && default_val.empty()) {
                        std::cout << "Command Line \"" << var << "\" not set!" << std::endl;
                        exit(EXIT_FAILURE);
                    } else if (input.empty()) {
                        return default_val;
                    }
                    return input;
                }
            }
            std::cout << "Command Line \"" << var << "\" is not an allowed name!" << std::endl;
            exit(EXIT_FAILURE);
        }

        /**
         * Prints the help message.
         */
        void print_help_message() {
            std::cout << "Dyn-HeiProMap - Dynamic Graph Partitioning Tool\n\n";
            std::cout << "Usage: Dyn-HeiProMap [options]\n\n";
            std::cout << "Options:\n";
            for (const auto &[large_key, small_key, description, default_val, input, is_set]: options) {
                if (small_key.empty()) {
                    std::cout << "  " << large_key << " - " << description;
                } else {
                    std::cout << "  " << large_key << ", " << small_key << " - " << description;
                }
                if (!default_val.empty()) {
                    std::cout << " (default: " << default_val << ")";
                }
                std::cout << "\n";
            }
            std::cout << "\nInteractive Commands:\n";
            std::cout << "  +v <v> [weight]                - Add vertex\n";
            std::cout << "  +e <u> <v> [weight]            - Add edge\n";
            std::cout << "  -v <v> [weight]                - Remove vertex\n";
            std::cout << "  -e <u> <v> [weight]            - Remove edge\n";
            std::cout << "  vw <v> <weight>                - Set vertex weight\n";
            std::cout << "  ew <u> <v> <weight>            - Set edge weight\n";
            std::cout << "  exec <file>                    - Execute commands from file\n";
            std::cout << "  set-hierarchy <h> <d>          - Set hierarchy (e.g., 2:2:2) and distance (e.g., 1:10:100)\n";
            std::cout << "  partition [config]             - Partition graph (config: fast|eco|strong)\n";
            std::cout << "  refine-fast [iter]             - Fast refinement (dirty vertices only)\n";
            std::cout << "  autorefine [threshold] [config]- Partition if dirty % > threshold, else refine-fast. config: fast|eco|strong\n";
            std::cout << "  no-refinement                  - Dummy command, does nothing\n";
            std::cout << "  save-partition <file>          - Save partition to file\n";
            std::cout << "  save-graph <file>              - Save graph to METIS format\n";
            std::cout << "  stats                          - Print graph statistics\n";
            std::cout << "  profile                        - Print profiling information\n";
            std::cout << "  quit                           - Exit program\n";
        }
    };
}

#endif //DYN_HEIPROMAP_CONFIGURATION_H
