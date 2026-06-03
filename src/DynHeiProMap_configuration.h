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

#ifndef HEIPROMAP_DYN_CONFIGURATION_H
#define HEIPROMAP_DYN_CONFIGURATION_H

#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>

#include "definitions.h"
#include "HeiProMap_configuration.h"

namespace HeiProMap {
    class DynConfiguration {
    public:
        std::vector<partition_t> hierarchy;
        std::vector<weight_t> distance;
        f64 imbalance = 0.03;
        u64 n_threads = 1;
        u64 seed = 0;

        AlgorithmConfiguration ac;

        DynConfiguration() = default;

        DynConfiguration(int argc, char *argv[]) : ac(argc, argv) {
            std::vector<std::string> args(argv, argv + argc);

            for (size_t i = 1; i < (size_t) argc; ++i) {
                if (args[i] == "--help") {
                    print_help_message();
                    exit(EXIT_SUCCESS);
                }
            }

            set_hierarchy();
            set_distance();
            set_imbalance();
            set_threads();
            set_seed();
        }

        void set_hierarchy() {
            if (is_set("--hierarchy")) {
                hierarchy.clear();
                std::vector<std::string> tokens = split(get("--hierarchy"), ':');
                for (auto &t: tokens) {
                    hierarchy.push_back(std::stoull(t));
                }
            } else {
                hierarchy = {2, 2, 2};
            }
            ac.hierarchy = hierarchy;
            ac.k = 1;
            for (auto h : hierarchy) ac.k *= h;
        }

        void set_distance() {
            if (is_set("--distance")) {
                distance.clear();
                std::vector<std::string> tokens = split(get("--distance"), ':');
                for (auto &t: tokens) {
                    distance.push_back(std::stoull(t));
                }
            } else {
                distance = {1, 10, 100};
            }
            ac.distance = distance;
        }

        void set_imbalance() {
            if (is_set("--imbalance")) {
                imbalance = std::stod(get("--imbalance"));
            }
            ac.imbalance = imbalance;
        }

        void set_threads() {
            if (is_set("--threads")) {
                n_threads = std::stoull(get("--threads"));
            }
            ac.threads = (u32)n_threads;
        }

        void set_seed() {
            if (is_set("--seed")) {
                seed = std::stoull(get("--seed"));
            }
            ac.seed = (u32)seed;
        }

        /**
         * Returns whether the option was entered.
         */
        bool is_set(const std::string &var) {
            for (const auto &o : ac.options) {
                if (o.large_key == var || o.small_key == var) {
                    return o.is_set;
                }
            }
            return false;
        }

        /**
         * Gets the entered input as a string.
         */
        std::string get(const std::string &var) {
            for (const auto &o : ac.options) {
                if (o.large_key == var || o.small_key == var) {
                    if (o.input.empty()) return o.default_val;
                    return o.input;
                }
            }
            return "";
        }

        /**
         * Prints the help message.
         */
        void print_help_message() {
            std::cout << "Dyn-HeiProMap - Dynamic Graph Partitioning Tool\n\n";
            std::cout << "Usage: Dyn-HeiProMap [options]\n\n";
            std::cout << "Options:\n";
            for (const auto &o : ac.options) {
                if (o.small_key.empty()) {
                    std::cout << "  " << std::left << std::setw(30) << o.large_key << " - " << o.description;
                } else {
                    std::cout << "  " << std::left << std::setw(30) << (o.large_key + ", " + o.small_key) << " - " << o.description;
                }
                if (!o.default_val.empty()) {
                    std::cout << " (default: " << o.default_val << ")";
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

#endif //HEIPROMAP_DYN_CONFIGURATION_H
