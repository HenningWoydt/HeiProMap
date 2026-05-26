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
 *
 ******************************************************************************/

#ifndef HEIPROMAP_HEIPA_CONFIGURATION_H
#define HEIPROMAP_HEIPA_CONFIGURATION_H

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
#include "HeiProMap_configuration.h"

namespace HeiProMap {

    class HeiPaConfiguration {
    private:
        std::vector<CommandLineOption> options = {
            {"--help", "", "Produces the help message", "", "", false},
            {"--graph", "-g", "Filepath to the graph.", "", "", false},
            {"--mapping", "-m", "Output filepath to the generated mapping.", "", "", false},
            {"--imbalance", "-e", "Allowed imbalance (for example 0.03).", "0.03", "", false},
            {"--config", "-c", "The configuration.", "", "", false},
            {"--threads", "-t", "Number of threads.", "1", "", false},
            {"--seed", "", "Seed for diversifying results.", "", "", false},
            {"--k", "-k", "Number of partitions.", "0", "", false},

            /** Coarsening */
            {"--coarsening-algorithm", "", "Which coarsening algorithm to use. Allowed values are {global-paths, size-constrained-lp}.", "global-paths", "", false},

            // Coarsening global-path
            {"--coarsening-algorithm-global-paths-random-level", "", "On which levels to run random if global-paths is chosen. Smaller-equal than use random, greater than use GPA.", "4", "", false},
            {"--coarsening-algorithm-global-paths-rating-function", "", "Which rating function to use. Allowed values are {weight, expansion, heavy-edge, greedy}.", "greedy", "", false},

            /** Partitioning */
            {"--partitioning-algorithm", "", "Which partitioning algorithm to use. Allowed values are {kaffpa, multisection, recursive-bisection}.", "recursive-bisection", "", false},

            // Partitioning kaffpa
            {"--partitioning-algorithm-kaffpa-partitioning-mode", "", "Which mode {strong, eco, fast} to use.", "strong", "", false},
            {"--partitioning-algorithm-kaffpa-partitioning-method", "", "Which mode {bisection, multisection} to use.", "multisection", "", false},

            // Partitioning multisection
            {"--partitioning-algorithm-multisection-mode", "", "Which mode {strong, eco, fast} to use.", "strong", "", false},

            // Partitioning recursive-bisection
            {"--partitioning-algorithm-recursive-bisection-kappa", "", "Number of independent trials for recursive bisection (best edge-cut is kept).", "1", "", false},

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

        partition_t k = 0;

        // balancing information
        f64 imbalance = -1.0;

        // random initialization
        u64 seed = 0;

        u64 threads = 1;

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
        KaffpaPartitionMode kaffpa_partition_mode = KAFFPA_PARTITION_STRONG;
        u64 rb_kappa = 1;

        // rebalance algorithm
        std::string rebalancing_algorithm_string;
        REBALANCING_ALGS rebalancing_algorithm_id = REBALANCING_ALG_UNDEFINED;

        // refinement algorithms
        LabelPropagationConfiguration label_propagation_config = LabelPropagationConfiguration("Label Propagation");
        QuotientGraphRefinementConfiguration quotient_graph_refinement_config = QuotientGraphRefinementConfiguration("Quotient Graph");
        FlowBasedRefinementConfiguration flow_based_refinement_config = FlowBasedRefinementConfiguration("Flow Based");

        HeiPaConfiguration() = default;

        void set_imbalance() { imbalance = std::stod(get("--imbalance")); }

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

        static std::string rating_function_to_string(const EdgeRatingFunction rating_function) {
            switch (rating_function) {
                case EdgeRatingFunction::WEIGHT:
                    return "weight";
                case EdgeRatingFunction::EXPANSION:
                    return "expansion";
                case EdgeRatingFunction::HEAVY_EDGE:
                    return "heavy-edge";
                case EdgeRatingFunction::GREEDY:
                    return "greedy";
                default:
                    return "UNKNOWN";
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

            if (use_default || is_set("--partitioning-algorithm-kaffpa-partitioning-mode")) {
                kaffpa_partition_mode = string_to_kaffpa_partition_mode(get("--partitioning-algorithm-kaffpa-partitioning-mode"));
            }

            if (use_default || is_set("--partitioning-algorithm-recursive-bisection-kappa")) {
                rb_kappa = std::stoull(get("--partitioning-algorithm-recursive-bisection-kappa"));
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

        HeiPaConfiguration(int argc, char *argv[]) {
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

            if (is_set("--k")) {
                k = std::stoull(get("--k"));
            }

            set_imbalance();
            set_seed();
            set_threads();

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

            coarsening_algorithm_string = "global-paths";
            global_path_algorithm_config.rating_function = EdgeRatingFunction::HEAVY_EDGE;
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            global_path_algorithm_config.random_level = 4;

            size_constrained_lp_config.max_rounds = 1;
            size_constrained_lp_config.min_threshold = 0.10;
            size_constrained_lp_config.multiplier = 16;

            partitioning_algorithm_string = "recursive-bisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            rb_kappa = 1;
            kaffpa_partition_mode = KAFFPA_PARTITION_FAST;

            label_propagation_config.enabled = true;
            label_propagation_config.max_iteration = 5;

            quotient_graph_refinement_config.enabled = true;
            quotient_graph_refinement_config.max_iteration = 2;
            quotient_graph_refinement_config.alpha = 5.0;
            quotient_graph_refinement_config.min_n_steps = 8;
            quotient_graph_refinement_config.use_preemptive_exit = true;
        }

        void set_eco() {
            initial_c = 8;

            coarsening_algorithm_string = "global-paths";
            global_path_algorithm_config.rating_function = EdgeRatingFunction::HEAVY_EDGE;
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            global_path_algorithm_config.random_level = 4;

            size_constrained_lp_config.max_rounds = 1;
            size_constrained_lp_config.min_threshold = 0.10;
            size_constrained_lp_config.multiplier = 4;

            partitioning_algorithm_string = "recursive-bisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            rb_kappa = 3;
            kaffpa_partition_mode = KAFFPA_PARTITION_ECO;

            label_propagation_config.enabled = true;
            label_propagation_config.max_iteration = 5;

            quotient_graph_refinement_config.enabled = true;
            quotient_graph_refinement_config.max_iteration = 2;
            quotient_graph_refinement_config.alpha = 5.0;
            quotient_graph_refinement_config.min_n_steps = 8;
            quotient_graph_refinement_config.use_preemptive_exit = true;

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

            coarsening_algorithm_string = "global-paths";
            global_path_algorithm_config.rating_function = EdgeRatingFunction::HEAVY_EDGE;
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            global_path_algorithm_config.random_level = 4;

            size_constrained_lp_config.max_rounds = 5;
            size_constrained_lp_config.min_threshold = 0.10;
            size_constrained_lp_config.multiplier = 2;

            partitioning_algorithm_string = "recursive-bisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            rb_kappa = 5;
            kaffpa_partition_mode = KAFFPA_PARTITION_STRONG;

            label_propagation_config.enabled = true;
            label_propagation_config.max_iteration = 5;

            quotient_graph_refinement_config.enabled = true;
            quotient_graph_refinement_config.max_iteration = 2;
            quotient_graph_refinement_config.alpha = 5.0;
            quotient_graph_refinement_config.min_n_steps = 8;
            quotient_graph_refinement_config.use_preemptive_exit = true;

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

            coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            global_path_algorithm_config.random_level = 0;

            size_constrained_lp_config.max_rounds = 5;
            size_constrained_lp_config.min_threshold = 0.10;
            size_constrained_lp_config.multiplier = 2;

            partitioning_algorithm_string = "recursive-bisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            rb_kappa = 10;
            kaffpa_partition_mode = KAFFPA_PARTITION_STRONG;

            label_propagation_config.enabled = true;
            label_propagation_config.max_iteration = 25;

            quotient_graph_refinement_config.enabled = true;
            quotient_graph_refinement_config.max_iteration = 2;
            quotient_graph_refinement_config.alpha = 1000.0;
            quotient_graph_refinement_config.min_n_steps = 10;

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
            coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            global_path_algorithm_config.random_level = 0;

            flow_based_refinement_config.enabled = false;
            flow_based_refinement_config.max_global_iteration = 1;
            flow_based_refinement_config.max_local_iteration = 5;
            flow_based_refinement_config.alpha = 8.0;
            flow_based_refinement_config.alpha_upper_bound = 32.0;
            flow_based_refinement_config.alpha_modifier = 2.0;
            flow_based_refinement_config.use_closed_vertex_set = true;
            flow_based_refinement_config.closed_vertex_sets_repeats = 100;
        }

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

        bool is_set(const std::string &var) {
            for (const auto &[large_key, small_key, description, default_val, input, is_set]: options) {
                if (large_key == var || small_key == var) {
                    return is_set;
                }
            }
            std::cout << "Command Line \"" << var << "\" is not an allowed name!" << std::endl;
            exit(EXIT_FAILURE);
        }

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

#endif //HEIPROMAP_HEIPA_CONFIGURATION_H
