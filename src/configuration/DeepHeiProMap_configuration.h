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

#ifndef HEIPROMAP_DEEP_HEIPROMAP_CONFIGURATION_H
#define HEIPROMAP_DEEP_HEIPROMAP_CONFIGURATION_H

#include <algorithm>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "../definitions.h"
#include "../utility/profiler.h"
#include "../utility/utils.h"
#include "../coarsening/heavy_edge_matching.h"
#include "../coarsening/global_path_algorithm.h"
#include "../coarsening/size_constrained_lp.h"
#include "../refinement/label_propagation_refinement.h"
#include "../refinement/simple_label_propagation_refinement.h"
#include "../refinement/quotient_graph_refinement.h"
#include "../refinement/simple_quotient_graph_refinement.h"
#include "../refinement/flow_based_refinement.h"

namespace HeiProMap {
    #ifndef HEIPROMAP_COARSENING_ALGS_DEFINED
    #define HEIPROMAP_COARSENING_ALGS_DEFINED

    enum COARSENING_ALGS {
        COARSENING_ALG_UNDEFINED,
        COARSENING_ALG_GLOBAL_PATHS,
        COARSENING_ALG_SIZE_CONSTRAINED_LP,
        COARSENING_ALG_HEAVY_EDGE,
        COARSENING_ALG_HEAVY_MATCHING = COARSENING_ALG_HEAVY_EDGE
    };

    inline COARSENING_ALGS string_to_coarsening_algorithm(const std::string &str) {
        if (str == "UNDEFINED") return COARSENING_ALG_UNDEFINED;
        if (str == "global-paths") return COARSENING_ALG_GLOBAL_PATHS;
        if (str == "size-constrained-lp") return COARSENING_ALG_SIZE_CONSTRAINED_LP;
        if (str == "heavy-edge" || str == "heavy-matching") return COARSENING_ALG_HEAVY_EDGE;
        return COARSENING_ALG_UNDEFINED;
    }

    inline std::string coarsening_algorithm_to_string(COARSENING_ALGS alg) {
        switch (alg) {
            case COARSENING_ALG_UNDEFINED:
                return "UNDEFINED";
            case COARSENING_ALG_GLOBAL_PATHS:
                return "global-paths";
            case COARSENING_ALG_SIZE_CONSTRAINED_LP:
                return "size-constrained-lp";
            case COARSENING_ALG_HEAVY_EDGE:
                return "heavy-edge";
            default:
                return "UNDEFINED";
        }
    }
    #endif

    #ifndef HEIPROMAP_REBALANCING_ALGS_DEFINED
    #define HEIPROMAP_REBALANCING_ALGS_DEFINED

    enum REBALANCING_ALGS {
        REBALANCING_ALG_UNDEFINED,
        REBALANCING_ALG_SIMPLE
    };

    inline REBALANCING_ALGS string_to_rebalancing_algorithm(const std::string &str) {
        if (str == "UNDEFINED") return REBALANCING_ALG_UNDEFINED;
        if (str == "simple") return REBALANCING_ALG_SIMPLE;
        return REBALANCING_ALG_UNDEFINED;
    }

    inline std::string rebalancing_algorithm_to_string(REBALANCING_ALGS alg) {
        switch (alg) {
            case REBALANCING_ALG_UNDEFINED:
                return "UNDEFINED";
            case REBALANCING_ALG_SIMPLE:
                return "simple";
            default:
                return "UNDEFINED";
        }
    }
    #endif

    enum DISTANCE_ORACLE_ALGS {
        DISTANCE_ORACLE_ALGS_UNDEFINED,
        DISTANCE_ORACLE_ALGS_DIVISION,
        DISTANCE_ORACLE_ALGS_STORE_DIVISION,
        DISTANCE_ORACLE_ALGS_BINARY
    };

    inline DISTANCE_ORACLE_ALGS string_to_distance_oracle_algorithm(const std::string &str) {
        if (str == "UNDEFINED") return DISTANCE_ORACLE_ALGS_UNDEFINED;
        if (str == "division-based") return DISTANCE_ORACLE_ALGS_DIVISION;
        if (str == "store-division-based") return DISTANCE_ORACLE_ALGS_STORE_DIVISION;
        if (str == "binary-based") return DISTANCE_ORACLE_ALGS_BINARY;
        return DISTANCE_ORACLE_ALGS_UNDEFINED;
    }

    inline std::string distance_oracle_algorithm_to_string(DISTANCE_ORACLE_ALGS alg) {
        switch (alg) {
            case DISTANCE_ORACLE_ALGS_UNDEFINED:
                return "UNDEFINED";
            case DISTANCE_ORACLE_ALGS_DIVISION:
                return "division-based";
            case DISTANCE_ORACLE_ALGS_STORE_DIVISION:
                return "store-division-based";
            case DISTANCE_ORACLE_ALGS_BINARY:
                return "binary-based";
            default:
                return "UNDEFINED";
        }
    }

    class DeepHeiProMapConfiguration {
    public:
        std::vector<CommandLineOption> options = {
            {"--help", "", "Produces the help message", "", "", false},
            {"--graph", "-g", "Filepath to the graph.", "", "", false},
            {"--mapping", "-m", "Output filepath to the generated mapping.", "", "", false},
            {"--statistics", "", "Output filepath to the statistics file.", "HeiProMap_stats.JSON", "", false},
            {"--hierarchy", "-h", "Hierarchy in the form a1:a2:...:al .", "", "", false},
            {"--distance", "-d", "Distance in the form d1:d2:...:dl .", "", "", false},
            {"--imbalance", "-e", "Allowed imbalance (for example 0.03).", "0.03", "", false},
            {"--config", "-c", "The configuration.", "", "", false},
            {"--threads", "-t", "The number of threads.", "1", "", false},
            {"--seed", "", "Seed for diversifying results.", "", "", false},
            {"--distance-oracle", "", "Which Distance Oracle to use. {division-based, store-division-based, binary-based}", "binary-based", "", false},
            {"--coarsening-alg", "", "Which coarsening algorithm to use. {global-paths, size-constrained-lp}", "size-constrained-lp", "", false},
            {"--hierarchical-coarsening", "", "Use staged hierarchical coarsening (true) or aggressive flat coarsening (false).", "true", "", false},
        };

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

        // configuration
        std::string config_string;

        // balancing information
        f64 imbalance = -1.0;

        // random initialization
        u64 seed = 0;

        // threads
        u64 threads = 1;

        std::string distance_oracle_algorithm_string;
        DISTANCE_ORACLE_ALGS distance_oracle_algorithm_id = DISTANCE_ORACLE_ALGS_UNDEFINED;

        // coarsening algorithm
        std::string coarsening_algorithm_string;
        COARSENING_ALGS coarsening_algorithm_id = COARSENING_ALG_UNDEFINED;

        GlobalPathAlgorithmConfiguration global_path_algorithm_config;
        HeavyEdgeMatchingConfiguration parallel_heavy_edge_matching_configuration;
        SizeConstrainedLPConfiguration size_constrained_lp_clustering_configuration;

        // hierarchy contraction thresholds
        bool hierarchical_coarsening = true;
        vertex_t initial_C = 32;
        u64 initial_kappa = 10;
        vertex_t intermediate_C = 32;
        u64 intermediate_kappa = 1;

        // datastructure flags
        bool enable_q_graph = true;
        bool enable_block_conn = true;

        // refinement algorithms
        SimpleLabelPropagationConfiguration deep_label_propagation_refinement_config = SimpleLabelPropagationConfiguration("Deep Label Propagation Refinement");
        QuotientGraphRefinementConfiguration deep_quotient_graph_refinement_config = QuotientGraphRefinementConfiguration("Deep Quotient Graph Refinement");
        SimpleQuotientGraphRefinementConfiguration deep_simple_qg_refinement_config = SimpleQuotientGraphRefinementConfiguration("Deep Simple QG Refinement");
        FlowBasedRefinementConfiguration deep_flow_based_refinement_config = FlowBasedRefinementConfiguration("Deep Flow Based Refinement");

        bool use_binary_oracle() const {
            return distance_oracle_algorithm_id == DISTANCE_ORACLE_ALGS_BINARY;
        }

        DeepHeiProMapConfiguration() = default;

        DeepHeiProMapConfiguration(int argc, char *argv[]) {
            HEIPROMAP_PROFILE_SCOPE("io", "DeepHeiProMapConfiguration", "parse_command_line");
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
                        input = args[i + 1];
                        is_set = true;
                        i += 1;
                        break;
                    }
                }
            }

            graph_in = get("--graph");
            mapping_out = get("--mapping");
            statistics_out = get("--statistics");

            hierarchy_string = get("--hierarchy");
            hierarchy = convert<partition_t>(split(hierarchy_string, ':'));
            k = prod<partition_t>(hierarchy);

            distance_string = get("--distance");
            distance = convert<weight_t>(split(distance_string, ':'));

            imbalance = std::stod(get("--imbalance"));

            if (is_set("--threads")) {
                threads = std::stoi(get("--threads"));
            }

            if (is_set("--seed")) {
                seed = std::stoi(get("--seed"));
            } else {
                seed = std::random_device{}();
            }

            config_string = get("--config");
            if (config_string == "fast") {
                set_fast();
            } else if (config_string == "eco") {
                set_eco();
            } else if (config_string == "strong") {
                set_strong();
            } else if (config_string == "experimental") {
                set_experimental();
            } else {
                std::cout << "Config " << config_string << " not recognized!" << std::endl;
                exit(EXIT_FAILURE);
            }

            if (is_set("--distance-oracle")) {
                distance_oracle_algorithm_string = get("--distance-oracle");
                distance_oracle_algorithm_id = string_to_distance_oracle_algorithm(distance_oracle_algorithm_string);
            }

            if (is_set("--coarsening-alg")) {
                coarsening_algorithm_string = get("--coarsening-alg");
                coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);
            }

            if (is_set("--hierarchical-coarsening")) {
                std::string val = get("--hierarchical-coarsening");
                hierarchical_coarsening = (val == "true" || val == "1" || val == "yes");
            }
        }

        void set_fast() {
            coarsening_algorithm_string = "size-constrained-lp";
            // coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            distance_oracle_algorithm_string = "binary-based";
            distance_oracle_algorithm_id = string_to_distance_oracle_algorithm(distance_oracle_algorithm_string);

            // refinement
            deep_label_propagation_refinement_config.enabled = true;
            deep_label_propagation_refinement_config.max_iteration = 5;
            deep_label_propagation_refinement_config.force_parallel_alg = false;

            deep_quotient_graph_refinement_config.enabled = false;
            deep_quotient_graph_refinement_config.max_iteration = 2;
            deep_quotient_graph_refinement_config.alpha = 5.0;
            deep_quotient_graph_refinement_config.min_n_steps = 3;
            deep_quotient_graph_refinement_config.use_preemptive_exit = true;

            deep_flow_based_refinement_config.enabled = false;
            deep_flow_based_refinement_config.max_global_iteration = 1;
            deep_flow_based_refinement_config.max_local_iteration = 1;
            deep_flow_based_refinement_config.alpha = 1.0;
            deep_flow_based_refinement_config.alpha_upper_bound = 64.0;
            deep_flow_based_refinement_config.alpha_modifier = 2.0;
            deep_flow_based_refinement_config.use_closed_vertex_set = true;
            deep_flow_based_refinement_config.always_include_boundary = true;
            deep_flow_based_refinement_config.closed_vertex_sets_repeats = 500;

            size_constrained_lp_clustering_configuration.rating_function = EdgeRatingFunction::EXPANSIONSTAR;
            size_constrained_lp_clustering_configuration.f = 32;
            size_constrained_lp_clustering_configuration.max_rounds = 5;

            enable_q_graph = false;
            enable_block_conn = false;
            deep_label_propagation_refinement_config.enable_q_graph = false;
            deep_label_propagation_refinement_config.enable_block_conn = false;
        }

        void set_eco() {
            coarsening_algorithm_string = "size-constrained-lp";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            distance_oracle_algorithm_string = "binary-based";
            distance_oracle_algorithm_id = string_to_distance_oracle_algorithm(distance_oracle_algorithm_string);

            // refinement
            deep_label_propagation_refinement_config.enabled = true;
            deep_label_propagation_refinement_config.max_iteration = 5;
            deep_label_propagation_refinement_config.force_parallel_alg = false;

            deep_quotient_graph_refinement_config.enabled = false;

            deep_simple_qg_refinement_config.enabled = true;
            deep_simple_qg_refinement_config.max_iteration = 1;
            deep_simple_qg_refinement_config.alpha = 5.0;
            deep_simple_qg_refinement_config.min_n_steps = 3;
            deep_simple_qg_refinement_config.use_preemptive_exit = true;

            deep_flow_based_refinement_config.enabled = false;

            enable_q_graph = false;
            enable_block_conn = false;
            deep_label_propagation_refinement_config.enable_q_graph = false;
            deep_label_propagation_refinement_config.enable_block_conn = false;
            deep_flow_based_refinement_config.use_active_block_scheduling = true;
            deep_flow_based_refinement_config.max_global_iteration = 5;
            deep_flow_based_refinement_config.max_local_iteration = 5;
            deep_flow_based_refinement_config.alpha = 2.0;
            deep_flow_based_refinement_config.alpha_upper_bound = 64.0;
            deep_flow_based_refinement_config.alpha_modifier = 2.0;
            deep_flow_based_refinement_config.use_closed_vertex_set = false;
            deep_flow_based_refinement_config.closed_vertex_sets_repeats = 500;
            deep_flow_based_refinement_config.always_include_boundary = true;
            deep_flow_based_refinement_config.growth_strategy = GrowthStrategy::BFS;

            size_constrained_lp_clustering_configuration.rating_function = EdgeRatingFunction::EXPANSION;
            size_constrained_lp_clustering_configuration.f = 32;
            size_constrained_lp_clustering_configuration.max_rounds = 10;
        }

        void set_strong() {
            coarsening_algorithm_string = "size-constrained-lp";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            distance_oracle_algorithm_string = "binary-based";
            distance_oracle_algorithm_id = string_to_distance_oracle_algorithm(distance_oracle_algorithm_string);

            // refinement
            deep_label_propagation_refinement_config.enabled = true;
            deep_label_propagation_refinement_config.max_iteration = 5;
            deep_label_propagation_refinement_config.force_parallel_alg = false;

            deep_quotient_graph_refinement_config.enabled = true;
            deep_quotient_graph_refinement_config.max_iteration = 2;
            deep_quotient_graph_refinement_config.alpha = 5.0;
            deep_quotient_graph_refinement_config.min_n_steps = 3;
            deep_quotient_graph_refinement_config.use_preemptive_exit = true;

            deep_flow_based_refinement_config.enabled = true;
            deep_flow_based_refinement_config.use_active_block_scheduling = true;
            deep_flow_based_refinement_config.max_global_iteration = 5;
            deep_flow_based_refinement_config.max_local_iteration = 5;
            deep_flow_based_refinement_config.alpha = 2.0;
            deep_flow_based_refinement_config.alpha_upper_bound = 64.0;
            deep_flow_based_refinement_config.alpha_modifier = 2.0;
            deep_flow_based_refinement_config.use_closed_vertex_set = false;
            deep_flow_based_refinement_config.closed_vertex_sets_repeats = 500;
            deep_flow_based_refinement_config.always_include_boundary = true;
            deep_flow_based_refinement_config.growth_strategy = GrowthStrategy::BFS;

            size_constrained_lp_clustering_configuration.rating_function = EdgeRatingFunction::EXPANSION;
            size_constrained_lp_clustering_configuration.f = 32;
            size_constrained_lp_clustering_configuration.max_rounds = 10;
        }

        void set_experimental() {
            coarsening_algorithm_string = "size-constrained-lp";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            distance_oracle_algorithm_string = "binary-based";
            distance_oracle_algorithm_id = string_to_distance_oracle_algorithm(distance_oracle_algorithm_string);

            // refinement
            deep_label_propagation_refinement_config.enabled = true;
            deep_label_propagation_refinement_config.max_iteration = 5;
            deep_label_propagation_refinement_config.force_parallel_alg = false;

            deep_quotient_graph_refinement_config.enabled = true;
            deep_quotient_graph_refinement_config.max_iteration = 2;
            deep_quotient_graph_refinement_config.alpha = 5.0;
            deep_quotient_graph_refinement_config.min_n_steps = 3;
            deep_quotient_graph_refinement_config.use_preemptive_exit = true;

            deep_flow_based_refinement_config.enabled = true;
            deep_flow_based_refinement_config.use_active_block_scheduling = true;
            deep_flow_based_refinement_config.max_global_iteration = 5;
            deep_flow_based_refinement_config.max_local_iteration = 5;
            deep_flow_based_refinement_config.alpha = 2.0;
            deep_flow_based_refinement_config.alpha_upper_bound = 64.0;
            deep_flow_based_refinement_config.alpha_modifier = 2.0;
            deep_flow_based_refinement_config.use_closed_vertex_set = false;
            deep_flow_based_refinement_config.closed_vertex_sets_repeats = 500;
            deep_flow_based_refinement_config.always_include_boundary = true;
            deep_flow_based_refinement_config.growth_strategy = GrowthStrategy::BFS;
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

        static std::string to_JSON() {
            return "{}";
        }
    };

    using AlgorithmConfiguration = DeepHeiProMapConfiguration;
}

#endif //HEIPROMAP_DEEP_HEIPROMAP_CONFIGURATION_H
