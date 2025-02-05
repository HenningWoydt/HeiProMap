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

#ifndef HEIPROMAP_ALGORITHMCONFIGURATION_H
#define HEIPROMAP_ALGORITHMCONFIGURATION_H

#include <string>
#include <vector>

#include "utils.h"
#include "../../definitions.h"
#include "../coarsening/greedy_edge_matcher.h"
#include "../partitioning/global_multisection.h"

namespace HeiProMap {
    enum COARSENING_ALGS {
        COARSENING_ALG_UNDEFINED,
        COARSENING_ALG_GREEDY_MATCHING,
        COARSENING_ALG_HEAVY_MATCHING,
        COARSENING_ALG_GLOBAL_PATHS,
    };

    inline COARSENING_ALGS string_to_coarsening_algorithm(const std::string& str) {
        if (str == "UNDEFINED") return COARSENING_ALG_UNDEFINED;
        if (str == "greedy-matching") return COARSENING_ALG_GREEDY_MATCHING;
        if (str == "heavy-matching") return COARSENING_ALG_HEAVY_MATCHING;
        if (str == "global-paths") return COARSENING_ALG_GLOBAL_PATHS;
        return COARSENING_ALG_UNDEFINED;
    }

    inline std::string coarsening_algorithm_to_string(COARSENING_ALGS alg) {
        switch (alg) {
        case COARSENING_ALG_UNDEFINED:
            return "UNDEFINED";
        case COARSENING_ALG_GREEDY_MATCHING:
            return "greedy-matching";
        case COARSENING_ALG_HEAVY_MATCHING:
            return "heavy-matching";
        case COARSENING_ALG_GLOBAL_PATHS:
            return "global-paths";
        default:
            return "UNDEFINED";
        }
    }

    enum PARTITIONING_ALGS {
        PARTITIONING_ALG_UNDEFINED,
        PARTITIONING_ALG_KAFFPA_MULTISECTION,
        PARTITIONING_ALG_MULTISECTION,
    };

    inline PARTITIONING_ALGS string_to_partitioning_algorithm(const std::string& str) {
        if (str == "UNDEFINED") return PARTITIONING_ALG_UNDEFINED;
        if (str == "kaffpa-multisection") return PARTITIONING_ALG_KAFFPA_MULTISECTION;
        if (str == "multisection") return PARTITIONING_ALG_MULTISECTION;
        return PARTITIONING_ALG_UNDEFINED;
    }

    inline std::string partitioning_algorithm_to_string(PARTITIONING_ALGS alg) {
        switch (alg) {
        case PARTITIONING_ALG_UNDEFINED:
            return "UNDEFINED";
        case PARTITIONING_ALG_KAFFPA_MULTISECTION:
            return "kaffpa-multisection";
        case PARTITIONING_ALG_MULTISECTION:
            return "multisection";
        default:
            return "UNDEFINED";
        }
    }

    struct CommandLineOption {
        std::string large_key;
        std::string small_key;
        std::string description;
        std::string default_val;
        std::string input;
        bool is_set;
    };

    class AlgorithmConfiguration {
    private:
        std::vector<CommandLineOption> options = {
                {"--help", "", "Produces the help message", "", "", false},
                {"--graph", "-g", "Filepath to the graph.", "", "", false},
                {"--mapping", "-m", "Output filepath to the generated mapping.", "", "", false},
                {"--statistics", "", "Output filepath to the statistics file.", "HeiProMap_stats.JSON", "", false},
                {"--hierarchy", "-h", "Hierarchy in the form a1:a2:...:al .", "", "", false},
                {"--distance", "-d", "Distance in the form d1:d2:...:dl .", "", "", false},
                {"--imbalance", "-e", "Allowed imbalance (for example 0.03).", "0.03", "", false},
                {"--config", "-c", "The configuration.", "", "", false},
                {"--seed", "", "Seed for diversifying results.", "", "", false},

                /** Coarsening */
                {"--coarsening-algorithm", "", "Which coarsening algorithm to use. Allowed values are {greedy-matching, heavy-matching-global-paths}.", "global-paths", "", false},

                // Coarsening greedy matching
                {"--coarsening-algorithm-greedy-matching-pendant-first", "", "Whether the greedy matching algorithm should handle pendant vertices first. 1 enables first matching pendant vertices, while 0 does not.", "1", "", false},
                {"--coarsening-algorithm-greedy-matching-no-overload", "", "Whether the greedy matching algorithm can match vertices, if it would lead to overload. 1 equals no overload possible, while 0 enables overloads.", "1", "", false},

                /** Partitioning */
                {"--partitioning-algorithm", "", "Which partitioning algorithm to use. Allowed values are {kaffpa-multisection, multisection}.", "multisection", "", false},

                // Partitioning multisection
                {"--partitioning-algorithm-multisection-mode", "", "Which mode {strong, eco, fast} to use.", "strong", "", false},

                // Refinement
                {"--enable-refinement-lable-propagation-faraj20", "", "Enables the label propagation by Faraj20.", "1", "", false},
                {"--enable-refinement-quotient-graph-faraj20", "", "Enables the quotient graph by Faraj20.", "1", "", false}
            };

    public:
        // graph information
        std::string graph_in;
        std::string mapping_out;
        std::string statistics_out;

        // hierarchy information
        std::string hierarchy_string;
        std::vector<partition_t> hierarchy;
        partition_t k = 0;

        // distance information
        std::string distance_string;
        std::vector<weight_t> distance;

        // balancing information
        f64 imbalance = -1.0;

        // random initialization
        u64 seed = 0;

        // coarsening algorithm
        std::string coarsening_algorithm_string;
        COARSENING_ALGS coarsening_algorithm_id = COARSENING_ALG_UNDEFINED;

        GreedyEdgeMatcherConfiguration greedy_edge_matcher_config;

        // partitioning algorithm
        std::string partitioning_algorithm_string;
        PARTITIONING_ALGS partitioning_algorithm_id = PARTITIONING_ALG_UNDEFINED;

        GlobalMultisectionConfiguration global_multisection_config;

        // refinement algorithms
        bool do_refinement_label_propagation_faraj20 = false;
        bool do_refinement_quotient_graph_faraj20    = false;

        AlgorithmConfiguration() = default;

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

        void set_seed() {
            if (is_set("--seed")) {
                seed = std::stoi(get("--seed"));
            } else {
                seed = std::random_device{}();
            }
        }

        void set_coarsening_algorithm() {
            // initialize greedy matching config
            greedy_edge_matcher_config.match_pendant_vertices_first = get("--coarsening-algorithm-greedy-matching-pendant-first") == "1";
            greedy_edge_matcher_config.no_overload                  = get("--coarsening-algorithm-greedy-matching-no-overload") == "1";

            // actually set which algorithm to use
            coarsening_algorithm_string = get("--coarsening-algorithm");
            coarsening_algorithm_id     = string_to_coarsening_algorithm(coarsening_algorithm_string);
        }

        void set_partitioning_algorithm() {
            // initialize global multisection config
            global_multisection_config.mode_string = get("--partitioning-algorithm-multisection-mode");
            global_multisection_config.mode = string_to_global_multisection_mode(global_multisection_config.mode_string);

            // actually set which algorithm to use
            partitioning_algorithm_string = get("--partitioning-algorithm");
            partitioning_algorithm_id     = string_to_partitioning_algorithm(partitioning_algorithm_string);
        }

        void enable_refinement_algorithms() {
            // initialize label propagation faraj20 configuration

            // initialize quotient graph faraj20 configuration

            // actually enable the configurations
            do_refinement_label_propagation_faraj20 = get("--enable-refinement-lable-propagation-faraj20") == "1";
            do_refinement_quotient_graph_faraj20    = get("--enable-refinement-quotient-graph-faraj20") == "1";
        }

        AlgorithmConfiguration(int argc, char* argv[]) {
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
                for (auto& [large_key, small_key, description, default_val, input, is_set] : options) {
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

            set_coarsening_algorithm();
            set_partitioning_algorithm();
            enable_refinement_algorithms();
        }

        /**
         * Gets the entered input as a string.
         *
         * @param var The option in interest.
         * @return The input.
         */
        std::string get(const std::string& var) {
            for (const auto& [large_key, small_key, description, default_val, input, is_set] : options) {
                if (large_key == var || small_key == var) {
                    if (input.empty()) {
                        std::cout << "Command Line \"" << var << "\" not set!" << std::endl;
                        exit(EXIT_FAILURE);
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
        bool is_set(const std::string& var) {
            for (const auto& [large_key, small_key, description, default_val, input, is_set] : options) {
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
            for (const auto& [large_key, small_key, description, default_val, input, is_set] : options) {
                if (small_key.empty()) {
                    std::cout << "[ " << large_key << "] - " << description << std::endl;
                } else {
                    std::cout << "[ " << large_key << ", " << small_key << "] - " << description << std::endl;
                }
            }
        }
    };
}

#endif //HEIPROMAP_ALGORITHMCONFIGURATION_H
