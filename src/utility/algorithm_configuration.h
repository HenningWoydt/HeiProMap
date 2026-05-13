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

#include "../definitions.h"
#include "utils.h"
#include "../partitioning/global_multisection.h"
#include "../refinement/label_propagation_refinement.h"
#include "../refinement/quotient_graph_refinement.h"
#include "../coarsening/global_path_algorithm.h"
#include "../coarsening/size_constrained_lp.h"
#include "../partitioning/kaffpa_partitioner.h"


namespace HeiProMap {
    enum COARSENING_ALGS {
        COARSENING_ALG_UNDEFINED,
        COARSENING_ALG_GLOBAL_PATHS,
        COARSENING_ALG_SIZE_CONSTRAINED_LP
    };

    inline COARSENING_ALGS string_to_coarsening_algorithm(const std::string &str) {
        if (str == "UNDEFINED") return COARSENING_ALG_UNDEFINED;
        if (str == "global-paths") return COARSENING_ALG_GLOBAL_PATHS;
        if (str == "size-constrained-lp") return COARSENING_ALG_SIZE_CONSTRAINED_LP;
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
            default:
                return "UNDEFINED";
        }
    }

    enum PARTITIONING_ALGS {
        PARTITIONING_ALG_UNDEFINED,
        PARTITIONING_ALG_KAFFPA,
        PARTITIONING_ALG_MULTISECTION,
    };

    inline PARTITIONING_ALGS string_to_partitioning_algorithm(const std::string &str) {
        if (str == "UNDEFINED") return PARTITIONING_ALG_UNDEFINED;
        if (str == "kaffpa") return PARTITIONING_ALG_KAFFPA;
        if (str == "multisection") return PARTITIONING_ALG_MULTISECTION;
        return PARTITIONING_ALG_UNDEFINED;
    }

    inline std::string partitioning_algorithm_to_string(PARTITIONING_ALGS alg) {
        switch (alg) {
            case PARTITIONING_ALG_UNDEFINED:
                return "UNDEFINED";
            case PARTITIONING_ALG_KAFFPA:
                return "kaffpa";
            case PARTITIONING_ALG_MULTISECTION:
                return "multisection";
            default:
                return "UNDEFINED";
        }
    }

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

    class AlgorithmConfiguration {
    private:
        std::vector<CommandLineOption> options = {
            {"--help", "", "Produces the help message", "", "", false},
            {"--graph", "-g", "Filepath to the graph.", "", "", false},
            {"--mapping", "-m", "Output filepath to the generated mapping.", "", "", false},
            {"--hierarchy", "-h", "Hierarchy in the form a1:a2:...:al .", "", "", false},
            {"--distance", "-d", "Distance in the form d1:d2:...:dl .", "", "", false},
            {"--imbalance", "-e", "Allowed imbalance (for example 0.03).", "0.03", "", false},
            {"--config", "-c", "The configuration.", "", "", false},
            {"--threads", "-t", "Number of threads.", "1", "", false},
            {"--seed", "", "Seed for diversifying results.", "", "", false},
            {"--hm-level", "", "Level of hierarchical multisection.", "0", "", false},

            /** Coarsening */
            {"--coarsening-algorithm", "", "Which coarsening algorithm to use. Allowed values are {global-paths, size-constrained-lp}.", "global-paths", "", false},

            // Coarsening global-path
            {"--coarsening-algorithm-global-paths-random-level", "", "On which levels to run random if global-paths is chosen. Smaller-equal than use random, greater than use GPA.", "4", "", false},
            {"--coarsening-algorithm-global-paths-rating-function", "", "Which rating function to use. Allowed values are {weight, expansion, heavy-edge, greedy}.", "greedy", "", false},

            /** Partitioning */
            {"--partitioning-algorithm", "", "Which partitioning algorithm to use. Allowed values are {kaffpa, multisection}.", "multisection", "", false},

            // Partitioning kaffpa
            {"--partitioning-algorithm-kaffpa-partitioning-mode", "", "Which mode {strong, eco, fast} to use.", "strong", "", false},
            {"--partitioning-algorithm-kaffpa-partitioning-method", "", "Which mode {bisection, multisection} to use.", "multisection", "", false},

            // Partitioning multisection
            {"--partitioning-algorithm-multisection-mode", "", "Which mode {strong, eco, fast} to use.", "strong", "", false},

            /** Rebalancing */
            {"--rebalancing-algorithm", "", "Which rebalancing algorithm to use. Allowed values are {simple}.", "simple", "", false},

            /** Refinement */
            // Refinement label propagation
            {"--refinement-label-propagation-enable", "", "Enables the label propagation refinement.", "0", "", false},
            {"--refinement-label-propagation-max-iterations", "", "For how many iterations to run label propagation refinement.", "25", "", false},

            // Refinement quotient graph
            {"--refinement-quotient-graph-enable", "", "Enables the quotient graph refinement.", "0", "", false},

            // Refinement Flow Based
            {"--refinement-flow-enable", "", "Enables the flow based refinement.", "0", "", false},
        };

    public:
        // graph information
        std::string graph_in;
        std::string mapping_out;

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

        u64 threads = 1;

        u64 hm_level = 0;

        u64 initial_c = 8;

        // coarsening algorithm
        std::string coarsening_algorithm_string;
        COARSENING_ALGS coarsening_algorithm_id = COARSENING_ALG_UNDEFINED;

        GlobalPathAlgorithmConfiguration global_path_algorithm_config;
        SizeConstrainedLPConfiguration size_constrained_lp_config;

        // partitioning algorithm
        std::string partitioning_algorithm_string;
        PARTITIONING_ALGS partitioning_algorithm_id = PARTITIONING_ALG_UNDEFINED;

        GlobalMultisectionConfiguration global_multisection_config;

        // rebalance algorithm
        std::string rebalancing_algorithm_string;
        REBALANCING_ALGS rebalancing_algorithm_id = REBALANCING_ALG_UNDEFINED;

        u64 n_refinement_iterations = 1;

        // refinement algorithms
        LabelPropagationConfiguration label_propagation_config = LabelPropagationConfiguration("Label Propagation");
        QuotientGraphRefinementConfiguration quotient_graph_refinement_config = QuotientGraphRefinementConfiguration("Quotient Graph");
        FlowBasedRefinementConfiguration flow_based_refinement_config = FlowBasedRefinementConfiguration("Flow Based");

        AlgorithmConfiguration() = default;

        void set_hierarchy() {
            hierarchy_string = get("--hierarchy");
            hierarchy = convert<partition_t>(split(hierarchy_string, ':'));
            k = prod<partition_t>(hierarchy);
        }

        void set_distance() {
            distance_string = get("--distance");
            distance = convert<weight_t>(split(distance_string, ':'));
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

        void set_threads() {
            if (is_set("--threads")) {
                threads = std::stoi(get("--threads"));
            } else {
                threads = 1;
            }
        }

        void set_hm_level() {
            if (is_set("--hm-level")) {
                hm_level = std::stoi(get("--hm-level"));
            }
        }

        static EdgeRatingFunction string_to_rating_function(const std::string &rating_function) {
            if (rating_function == "weight") {
                return EdgeRatingFunction::WEIGHT;
            } else if (rating_function == "expansion") {
                return EdgeRatingFunction::EXPANSION;
            } else if (rating_function == "heavy-edge") {
                return EdgeRatingFunction::HEAVY_EDGE;
            } else if (rating_function == "greedy") {
                return EdgeRatingFunction::GREEDY;
            } else {
                std::cerr << "Unknown rating function: " << rating_function << std::endl;
                exit(1);
            }
        }

        void set_coarsening_algorithm(const bool use_default = false) {
            // initialize global-paths config
            if (use_default || is_set("--coarsening-algorithm-global-paths-random-level")) {
                global_path_algorithm_config.random_level = std::stoi(get("--coarsening-algorithm-global-paths-random-level"));
            }

            if (use_default || is_set("--coarsening-algorithm-global-paths-rating-function")) {
                global_path_algorithm_config.rating_function = string_to_rating_function(get("--coarsening-algorithm-global-paths-rating-function"));
            }

            // actually set which algorithm to use
            if (use_default || is_set("--coarsening-algorithm")) {
                coarsening_algorithm_string = get("--coarsening-algorithm");
                coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);
            }
        }

        void set_partitioning_algorithm(const bool use_default = false) {
            // initialize global multisection config
            if (use_default || is_set("--partitioning-algorithm-multisection-mode")) {
                global_multisection_config.mode_string = get("--partitioning-algorithm-multisection-mode");
                global_multisection_config.mode = string_to_global_multisection_mode(global_multisection_config.mode_string);
            }

            // actually set which algorithm to use
            if (use_default || is_set("--partitioning-algorithm")) {
                partitioning_algorithm_string = get("--partitioning-algorithm");
                partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            }
        }

        void set_rebalancing_algorithm(const bool use_default = false) {
            if (use_default || is_set("--rebalancing-algorithm")) {
                rebalancing_algorithm_string = get("--rebalancing-algorithm");
                rebalancing_algorithm_id = string_to_rebalancing_algorithm(rebalancing_algorithm_string);
            }
        }

        void enable_refinement_algorithms(const bool use_default = false) {
            // initialize label propagation configuration
            if (use_default || is_set("--refinement-label-propagation-max-iterations")) {
                label_propagation_config.max_iteration = std::stoi(get("--refinement-label-propagation-max-iterations"));
            }
            if (use_default || is_set("--refinement-label-propagation-enable")) {
                label_propagation_config.enabled = get("--refinement-label-propagation-enable") == "1";
            }

            // initialize quotient graph configuration
            if (use_default || is_set("--refinement-quotient-graph-enable")) {
                quotient_graph_refinement_config.enabled = get("--refinement-quotient-graph-enable") == "1";
            }

            // initialize flow based refinement
            if (use_default || is_set("--refinement-flow-enable")) {
                flow_based_refinement_config.enabled = get("--refinement-flow-enable") == "1";
            }
        }

        AlgorithmConfiguration(int argc, char *argv[]) {
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

            set_hierarchy();
            set_distance();
            set_imbalance();
            set_seed();
            set_threads();
            set_hm_level();

            set_coarsening_algorithm(true);
            set_partitioning_algorithm(true);
            set_rebalancing_algorithm(true);
            enable_refinement_algorithms(true);

            if (get("--config") == "fast") {
                set_fast();
            } else if (get("--config") == "eco") {
                set_eco();
            } else if (get("--config") == "strong") {
                set_strong();
            } else if (get("--config") == "super-strong") {
                set_super_strong();
            } else if (get("--config") == "experimental") {
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

        void set_fast() {
            initial_c = 8;

            // set GPA matching algorithm
            coarsening_algorithm_string = "size-constrained-lp";
            // coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.random_level = 4;

            size_constrained_lp_config.max_rounds = 1;
            size_constrained_lp_config.min_threshold = 0.10;
            size_constrained_lp_config.multiplier = 16;

            // set multisection
            partitioning_algorithm_string = "multisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            global_multisection_config.mode_string = "metis-kway";
            global_multisection_config.mode = string_to_global_multisection_mode(global_multisection_config.mode_string);
            global_multisection_config.kappa = 3;

            // enable label propagation
            label_propagation_config.enabled = true;
            label_propagation_config.max_iteration = 25;

            // enable quotient graph refinement
            quotient_graph_refinement_config.enabled = true;
            quotient_graph_refinement_config.max_iteration = 5;
            quotient_graph_refinement_config.alpha = 1000.0;
            quotient_graph_refinement_config.min_n_steps = 4;
            quotient_graph_refinement_config.use_preemptive_exit = true;
        }

        void set_eco() {
            initial_c = 8;

            n_refinement_iterations = 1;

            // set GPA matching algorithm
            coarsening_algorithm_string = "global-paths";
            // coarsening_algorithm_string = "size-constrained-lp";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.random_level = 0;

            size_constrained_lp_config.max_rounds = 1;
            size_constrained_lp_config.min_threshold = 0.10;
            size_constrained_lp_config.multiplier = 4;

            // set multisection
            partitioning_algorithm_string = "multisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            global_multisection_config.mode_string = "kaffpa-eco";
            // global_multisection_config.mode_string = "metis-kway";
            global_multisection_config.mode = string_to_global_multisection_mode(global_multisection_config.mode_string);
            global_multisection_config.kappa = 3;

            // enable label propagation
            label_propagation_config.enabled = true;
            label_propagation_config.max_iteration = 25;

            // enable quotient graph refinement
            quotient_graph_refinement_config.enabled = true;
            quotient_graph_refinement_config.max_iteration = 10;
            quotient_graph_refinement_config.alpha = 1000.0;
            quotient_graph_refinement_config.min_n_steps = 10;
            quotient_graph_refinement_config.use_active_scheduling = true;

            // enable flow based refinement
            flow_based_refinement_config.enabled = true;
            flow_based_refinement_config.max_global_iteration = 1;
            flow_based_refinement_config.max_local_iteration = 1;
            flow_based_refinement_config.alpha = 1.0;
            flow_based_refinement_config.alpha_upper_bound = 64.0;
            flow_based_refinement_config.alpha_modifier = 2.0;
            flow_based_refinement_config.use_closed_vertex_set = true;
            flow_based_refinement_config.closed_vertex_sets_repeats = 500;
        }

        void set_strong() {
            initial_c = 16;

            n_refinement_iterations = 1;

            // set GPA matching algorithm
            coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.random_level = 0;

            size_constrained_lp_config.max_rounds = 5;
            size_constrained_lp_config.min_threshold = 0.10;
            size_constrained_lp_config.multiplier = 2;

            // set multisection
            partitioning_algorithm_string = "multisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            global_multisection_config.mode_string = "kaffpa-strong";
            global_multisection_config.mode = string_to_global_multisection_mode(global_multisection_config.mode_string);
            global_multisection_config.kappa = 1;

            // enable label propagation
            label_propagation_config.enabled = true;
            label_propagation_config.max_iteration = 25;

            // enable quotient graph refinement
            quotient_graph_refinement_config.enabled = true;
            quotient_graph_refinement_config.max_iteration = 2;
            quotient_graph_refinement_config.alpha = 1000.0;
            quotient_graph_refinement_config.min_n_steps = 10;

            // enable flow based refinement
            flow_based_refinement_config.enabled = true;
            flow_based_refinement_config.use_active_block_scheduling = true;
            flow_based_refinement_config.max_global_iteration = 3;
            flow_based_refinement_config.max_local_iteration = 3;
            flow_based_refinement_config.alpha = 1.0;
            flow_based_refinement_config.alpha_upper_bound = 16.0;
            flow_based_refinement_config.alpha_modifier = 2.0;
            flow_based_refinement_config.use_closed_vertex_set = true;
            flow_based_refinement_config.closed_vertex_sets_repeats = 500;
        }

        void set_super_strong() {
            initial_c = 16;

            n_refinement_iterations = 1;

            // set GPA matching algorithm
            coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.random_level = 0;

            size_constrained_lp_config.max_rounds = 5;
            size_constrained_lp_config.min_threshold = 0.10;
            size_constrained_lp_config.multiplier = 2;

            // set multisection
            partitioning_algorithm_string = "multisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            global_multisection_config.mode_string = "kaffpa-strong";
            global_multisection_config.mode = string_to_global_multisection_mode(global_multisection_config.mode_string);
            global_multisection_config.kappa = 1;

            // enable label propagation
            label_propagation_config.enabled = true;
            label_propagation_config.max_iteration = 25;

            // enable quotient graph refinement
            quotient_graph_refinement_config.enabled = true;
            quotient_graph_refinement_config.max_iteration = 2;
            quotient_graph_refinement_config.alpha = 1000.0;
            quotient_graph_refinement_config.min_n_steps = 10;

            // enable flow based refinement
            flow_based_refinement_config.enabled = true;
            flow_based_refinement_config.use_active_block_scheduling = true;
            flow_based_refinement_config.max_global_iteration = 5;
            flow_based_refinement_config.max_local_iteration = 10;
            flow_based_refinement_config.alpha = 8.0;
            flow_based_refinement_config.alpha_upper_bound = 16.0;
            flow_based_refinement_config.alpha_modifier = 2.0;
            flow_based_refinement_config.use_closed_vertex_set = true;
            flow_based_refinement_config.closed_vertex_sets_repeats = 500;
        }

        void set_experimental() {
            // set GPA matching algorithm
            coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.random_level = 0;

            // enable flow based refinement
            flow_based_refinement_config.enabled = false;
            flow_based_refinement_config.max_global_iteration = 1;
            flow_based_refinement_config.max_local_iteration = 5;
            flow_based_refinement_config.alpha = 8.0;
            flow_based_refinement_config.alpha_upper_bound = 32.0;
            flow_based_refinement_config.alpha_modifier = 2.0;
            flow_based_refinement_config.use_closed_vertex_set = true;
            flow_based_refinement_config.closed_vertex_sets_repeats = 100;
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
                        std::cout << "Command Line "" << var << "" not set!" << std::endl;
                        exit(EXIT_FAILURE);
                    } else if (input.empty()) {
                        return default_val;
                    }
                    return input;
                }
            }
            std::cout << "Command Line "" << var << "" is not an allowed name!" << std::endl;
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
            std::cout << "Command Line "" << var << "" is not an allowed name!" << std::endl;
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

#endif //HEIPROMAP_ALGORITHMCONFIGURATION_H
