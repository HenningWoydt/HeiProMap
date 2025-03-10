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

#ifndef HEIPROMAP_PARALLEL_ALGORITHM_CONFIGURATION_H
#define HEIPROMAP_PARALLEL_ALGORITHM_CONFIGURATION_H

#include <string>
#include <vector>
#include <random>

#include "../../commons.h"
#include "../../commons/utils.h"
#include "parallel_utils.h"
#include "../../definitions.h"
#include "../coarsening/parallel_heavy_edge_matcher.h"
#include "../partitioning/parallel_kaffpa_partitioner.h"

namespace HeiProMap {
    enum PARALLEL_COARSENING_ALGS {
        PARALLEL_COARSENING_ALG_UNDEFINED,
        PARALLEL_COARSENING_ALG_HEAVY_MATCHING,
    };

    inline PARALLEL_COARSENING_ALGS string_to_coarsening_algorithm(const std::string &str) {
        if (str == "UNDEFINED") return PARALLEL_COARSENING_ALG_UNDEFINED;
        if (str == "heavy-matching") return PARALLEL_COARSENING_ALG_HEAVY_MATCHING;
        return PARALLEL_COARSENING_ALG_UNDEFINED;
    }

    inline std::string coarsening_algorithm_to_string(PARALLEL_COARSENING_ALGS alg) {
        switch (alg) {
            case PARALLEL_COARSENING_ALG_UNDEFINED:
                return "UNDEFINED";
            case PARALLEL_COARSENING_ALG_HEAVY_MATCHING:
                return "heavy-matching";
            default:
                return "UNDEFINED";
        }
    }

    enum PARALLEL_PARTITIONING_ALGS {
        PARALLEL_PARTITIONING_ALG_UNDEFINED,
        PARALLEL_PARTITIONING_ALG_KAFFPA
    };

    inline PARALLEL_PARTITIONING_ALGS parallel_string_to_partitioning_algorithm(const std::string &str) {
        if (str == "UNDEFINED") return PARALLEL_PARTITIONING_ALG_UNDEFINED;
        if (str == "kaffpa") return PARALLEL_PARTITIONING_ALG_KAFFPA;
        return PARALLEL_PARTITIONING_ALG_UNDEFINED;
    }

    inline std::string parallel_partitioning_algorithm_to_string(PARALLEL_PARTITIONING_ALGS alg) {
        switch (alg) {
            case PARALLEL_PARTITIONING_ALG_UNDEFINED:
                return "UNDEFINED";
            case PARALLEL_PARTITIONING_ALG_KAFFPA:
                return "kaffpa";
            default:
                return "UNDEFINED";
        }
    }

    class ParallelAlgorithmConfiguration {
    private:
        std::vector<CommandLineOption> options = {
                {"--help",                                                       "",   "Produces the help message",                                                                                                                        "",                     "", false},
                {"--graph",                                                      "-g", "Filepath to the graph.",                                                                                                                           "",                     "", false},
                {"--mapping",                                                    "-m", "Output filepath to the generated mapping.",                                                                                                        "",                     "", false},
                {"--statistics",                                                 "",   "Output filepath to the statistics file.",                                                                                                          "HeiProMap_stats.JSON", "", false},
                {"--hierarchy",                                                  "-h", "Hierarchy in the form a1:a2:...:al .",                                                                                                             "",                     "", false},
                {"--distance",                                                   "-d", "Distance in the form d1:d2:...:dl .",                                                                                                              "",                     "", false},
                {"--imbalance",                                                  "-e", "Allowed imbalance (for example 0.03).",                                                                                                            "0.03",                 "", false},
                {"--threads",                                                    "-t", "Number of threads to use.",                                                                                                                        "16",                   "", false},
                {"--config",                                                     "-c", "The configuration.",                                                                                                                               "",                     "", false},
                {"--seed",                                                       "",   "Seed for diversifying results.",                                                                                                                   "",                     "", false},

                /** Coarsening */
                {"--parallel-coarsening-algorithm",                              "",   "Which parallel coarsening algorithm to use. Allowed values are {heavy-matching}.",                                                                 "heavy-matching",       "", false},

                // Coarsening heavy matching
                {"--parallel-coarsening-algorithm-heavy-matching-pendant-first", "",   "Whether the parallel heavy matching algorithm should handle pendant vertices first. 1 enables first matching pendant vertices, while 0 does not.", "1",                    "", false},

                /** Partitioning */
                {"--parallel-partitioning-algorithm",                            "",   "Which parallel partitioning algorithm to use. Allowed values are {kaffpa}.",                                                                       "kaffpa",               "", false},

                // Partitioning kaffpa
                {"--parallel-partitioning-algorithm-kaffpa-partitioning-mode",   "",   "Which mode {strong, eco, fast} to use.",                                                                                                           "strong",               "", false},
                {"--parallel-partitioning-algorithm-kaffpa-partitioning-method", "",   "Which mode {bisection, multisection} to use.",                                                                                                     "multisection",         "", false},
        };

    public:
        // graph information
        std::string graph_in;
        std::string mapping_out;
        std::string statistics_out;

        // hierarchy information
        std::string              hierarchy_string;
        std::vector<partition_t> hierarchy;
        partition_t              k = 0;

        // distance information
        std::string           distance_string;
        std::vector<weight_t> distance;

        // thread information
        u64 threads = 1;

        // balancing information
        f64 imbalance = -1.0;

        // random initialization
        u64 seed = 0;

        // coarsening algorithm
        std::string              parallel_coarsening_algorithm_string;
        PARALLEL_COARSENING_ALGS parallel_coarsening_algorithm_id = PARALLEL_COARSENING_ALG_UNDEFINED;

        ParallelHeavyEdgeMatcherConfiguration parallel_heavy_edge_matcher_config;

        // partitioning algorithm
        std::string                parallel_partitioning_algorithm_string;
        PARALLEL_PARTITIONING_ALGS parallel_partitioning_algorithm_id = PARALLEL_PARTITIONING_ALG_UNDEFINED;

        ParallelKaffpaPartitionerConfiguration parallel_kaffpa_partitioner_config;

        ParallelAlgorithmConfiguration() = default;

        void set_hierarchy() {
            hierarchy_string = get("--hierarchy");
            hierarchy        = convert<partition_t>(split(hierarchy_string, ':'));
            k                = prod<partition_t>(hierarchy);
        }

        void set_distance() {
            distance_string = get("--distance");
            distance        = convert<weight_t>(split(distance_string, ':'));
        }

        void set_imbalance() {
            imbalance = std::stod(get("--imbalance"));
        }

        void set_threads() {
            threads = std::stoi(get("--threads"));
        }

        void set_seed() {
            if (is_set("--seed")) {
                seed = std::stoi(get("--seed"));
            } else {
                seed = std::random_device{}();
            }
        }

        void set_coarsening_algorithm(const bool use_default = false) {
            // initialize heavy matching config
            if (use_default || is_set("--parallel-coarsening-algorithm-heavy-matching-pendant-first")) {
                parallel_heavy_edge_matcher_config.match_pendant_vertices_first = get("--parallel-coarsening-algorithm-heavy-matching-pendant-first") == "1";
            }

            // actually set which algorithm to use
            if (use_default || is_set("--parallel-coarsening-algorithm")) {
                parallel_coarsening_algorithm_string = get("--parallel-coarsening-algorithm");
                parallel_coarsening_algorithm_id     = string_to_coarsening_algorithm(parallel_coarsening_algorithm_string);
            }
        }

        void set_partitioning_algorithm(const bool use_default = false) {
            // initialize kaffpa partitioning config
            if (use_default || is_set("--parallel-partitioning-algorithm-kaffpa-partitioning-mode")) {
                parallel_kaffpa_partitioner_config.mode_string = get("--parallel-partitioning-algorithm-kaffpa-partitioning-mode");
                parallel_kaffpa_partitioner_config.mode        = parallel_string_to_kaffpa_partitioner_mode(parallel_kaffpa_partitioner_config.mode_string);
            }
            if (use_default || is_set("--parallel-partitioning-algorithm-kaffpa-partitioning-method")) {
                parallel_kaffpa_partitioner_config.method_string = get("--parallel-partitioning-algorithm-kaffpa-partitioning-method");
                parallel_kaffpa_partitioner_config.method        = parallel_string_to_kaffpa_partitioner_method(parallel_kaffpa_partitioner_config.method_string);
            }

            // actually set which algorithm to use
            if (use_default || is_set("--parallel-partitioning-algorithm")) {
                parallel_partitioning_algorithm_string = get("--parallel-partitioning-algorithm");
                parallel_partitioning_algorithm_id     = parallel_string_to_partitioning_algorithm(parallel_partitioning_algorithm_string);
            }
        }

        void set_rebalancing_algorithm(const bool use_default = false) {

        }

        void enable_refinement_algorithms(const bool use_default = false) {

        }

        ParallelAlgorithmConfiguration(int argc, char *argv[]) {
            // read command lines into vector
            std::vector<std::string> args(argv, argv + argc);

            // check for a help message
            for (int i = 1; i < argc; ++i) {
                if (args[i] == "--help") {
                    print_help_message();
                    exit(EXIT_SUCCESS);
                }
            }

            // read all command line args
            for (int i = 1; i < argc; ++i) {
                for (auto &[large_key, small_key, description, default_val, input, is_set]: options) {
                    if (large_key == args[i] || small_key == args[i]) {
                        input  = args[i + 1];
                        is_set = true;
                        i += 1;
                        break;
                    }
                }
            }

            graph_in       = get("--graph");
            mapping_out    = get("--mapping");
            statistics_out = get("--statistics");

            set_hierarchy();
            set_distance();
            set_imbalance();
            set_seed();

            set_coarsening_algorithm(true);
            set_partitioning_algorithm(true);
            set_rebalancing_algorithm(true);
            enable_refinement_algorithms(true);

            if (get("--config") == "experimental") {
                set_experimental();
            } else {
                std::cout << "Config " << get("--config") << " not recognized!" << std::endl;
                exit(EXIT_FAILURE);
            }

            set_coarsening_algorithm();
            set_partitioning_algorithm();
            set_rebalancing_algorithm();
            enable_refinement_algorithms();
        }


        void set_experimental() {

        }

        /**
         * Gets the entered input as a string.
         *
         * @param var The option in interest.
         * @return The input.
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
         * Returns whether the option was entered.
         *
         * @param var The option in interest.
         * @return True if the option was entered, false else.
         */
        bool is_set(const std::string &var) {
            for (const auto &[large_key, small_key, description, default_val, input, is_set]: options) {
                if (large_key == var || small_key == var) {
                    return is_set;
                }
            }
            std::cout << "Command Line \"" << var << "\" is not an allowed name!" << std::endl;
            exit(EXIT_FAILURE);
        }

        /**
         * Prints the help message.
         */
        void print_help_message() {
            for (const auto &[large_key, small_key, description, default_val, input, is_set]: options) {
                if (small_key.empty()) {
                    std::cout << "[ " << large_key << "] - " << description << std::endl;
                } else {
                    std::cout << "[ " << large_key << ", " << small_key << "] - " << description << std::endl;
                }
            }
        }
    };
}

#endif //HEIPROMAP_PARALLEL_ALGORITHM_CONFIGURATION_H
