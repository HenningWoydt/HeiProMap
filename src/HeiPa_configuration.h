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

#include "definitions.h"
#include "utility/utils.h"
#include "partitioning/global_multisection.h"
#include "refinement/label_propagation_refinement.h"
#include "refinement/quotient_graph_refinement.h"
#include "coarsening/global_path_algorithm.h"
#include "coarsening/size_constrained_lp.h"
#include "partitioning/kaffpa_partitioner.h"
#include "partitioning/recursive_bisection.h"
#include "refinement/flow_based_refinement.h"
#include "refinement/negative_cycle_detection.h"
#include "HeiProMap_configuration.h"

namespace HeiProMap {
    class HeiPaConfiguration {
    private:
        std::vector<CommandLineOption> options = {
            // General Options
            {"--help", "", "Produces the help message", "", "", false},
            {"--graph", "-g", "Filepath to the graph.", "", "", false},
            {"--mapping", "-m", "Output filepath to the generated mapping.", "", "", false},
            {"--imbalance", "-e", "Allowed imbalance (for example 0.03).", "0.03", "", false},
            {"--config", "-c", "The configuration.", "", "", false},
            {"--threads", "-t", "Number of threads.", "1", "", false},
            {"--seed", "", "Seed for diversifying results.", "", "", false},
            {"--k", "-k", "Number of partitions.", "0", "", false},
            {"--json-output", "", "If true, prints a JSON summary of the result at the end.", "0", "", false},
            {"--use-parallel-contraction", "", "Use parallel contraction (true/false).", "", "", false},

            // Label Propagation (LP) Refinement
            {"--lp-refinement-enabled", "", "Enable label propagation refinement (true/false).", "", "", false},
            {"--lp-refinement-max-iteration", "", "Max iteration for label propagation refinement.", "", "", false},
            {"--lp-refinement-use-parallel-version", "", "Label propagation use parallel version (true/false).", "", "", false},
            {"--lp-refinement-use-edge-cut", "", "Label propagation use edge-cut on last level (true/false).", "", "", false},

            // Quotient Graph (QG) Refinement
            {"--qg-refinement-enabled", "", "Enable quotient graph refinement (true/false).", "", "", false},
            {"--qg-refinement-max-iteration", "", "Max iteration for quotient graph refinement.", "", "", false},
            {"--qg-refinement-min-n-steps", "", "Min n steps for quotient graph refinement.", "", "", false},
            {"--qg-refinement-alpha", "", "Alpha parameter for quotient graph refinement.", "", "", false},
            {"--qg-refinement-beta", "", "Beta parameter for quotient graph refinement.", "", "", false},
            {"--qg-refinement-use-active-scheduling", "", "Use active scheduling in quotient graph refinement (true/false).", "", "", false},
            {"--qg-refinement-use-preemptive-exit", "", "Use preemptive exit in quotient graph refinement (true/false).", "", "", false},
            {"--qg-refinement-use-edge-cut", "", "Quotient graph refinement use edge-cut on last level (true/false).", "", "", false},
        };

    public:
        // graph information
        std::string graph_in;
        std::string mapping_out;

        partition_t k = 0;

        // balancing information
        f64 imbalance = -1.0;
        f64 per_level_imb_add = 0.005;

        // random initialization
        u64 seed = 0;

        u64 threads = 1;

        bool collect_dataset = false;
        std::string data_dir = std::string(HEIPROMAP_SOURCE_DIR) + "/data";
        bool json_output = false;

        u64 initial_c = 8;
        bool use_parallel_contraction = false;

        // coarsening algorithm
        std::string coarsening_algorithm_string;
        COARSENING_ALGS coarsening_algorithm_id = COARSENING_ALG_UNDEFINED;

        GlobalPathAlgorithmConfiguration global_path_algorithm_config;
        SizeConstrainedLPConfiguration size_constrained_lp_config;
        HeavyEdgeMatchingConfiguration heavy_edge_matching_configuration;

        // partitioning algorithm
        std::string partitioning_algorithm_string;
        PARTITIONING_ALGS partitioning_algorithm_id = PARTITIONING_ALG_UNDEFINED;

        RecursiveBisectionConfiguration recursive_bisection_config;

        GlobalMultisectionConfiguration global_multisection_config;
        KaffpaPartitionMode kaffpa_partition_mode = KAFFPA_PARTITION_STRONG;

        // rebalance algorithm
        std::string rebalancing_algorithm_string;
        REBALANCING_ALGS rebalancing_algorithm_id = REBALANCING_ALG_UNDEFINED;

        // refinement algorithms
        LabelPropagationConfiguration label_propagation_config = LabelPropagationConfiguration("Label Propagation");
        QuotientGraphRefinementConfiguration quotient_graph_refinement_config = QuotientGraphRefinementConfiguration("Quotient Graph");
        FlowBasedRefinementConfiguration flow_based_refinement_config = FlowBasedRefinementConfiguration("Flow Based");
        NegativeCycleConfiguration negative_cycle_refinement_config = NegativeCycleConfiguration("Negative Cycle");

        HeiPaConfiguration() = default;

        void set_imbalance() {
            imbalance = std::stod(get("--imbalance"));
        }

        void set_seed() {
            if (is_set("--seed")) {
                seed = std::stoull(get("--seed"));
            } else {
                seed = std::random_device{}();
            }
        }

        void set_threads() {
            if (is_set("--threads")) {
                threads = std::stoull(get("--threads"));
            } else {
                threads = 1;
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
                bool found = false;
                for (auto &[large_key, small_key, description, default_val, input, is_set]: options) {
                    if (large_key == args[i] || (small_key != "" && small_key == args[i])) {
                        if (i + 1 < argc) {
                            input = args[i + 1];
                            is_set = true;
                            i += 1;
                            found = true;
                            break;
                        } else {
                            std::cout << "Error: Value missing for parameter " << args[i] << std::endl;
                            exit(EXIT_FAILURE);
                        }
                    }
                }
                if (!found) {
                    std::cout << "Error: Unknown parameter " << args[i] << std::endl;
                    exit(EXIT_FAILURE);
                }
            }

            graph_in = get("--graph");
            if (is_set("--mapping")) {
                mapping_out = get("--mapping");
            }

            if (is_set("--k")) {
                k = std::stoull(get("--k"));
            }

            set_imbalance();
            set_seed();
            set_threads();

            if (is_set("--json-output")) {
                json_output = get("--json-output") == "1";
            }

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

            size_constrained_lp_config.use_parallel_version = (threads > 1);
            use_parallel_contraction = (threads > 1);
            if (is_set("--use-parallel-contraction")) {
                std::string val = get("--use-parallel-contraction");
                use_parallel_contraction = (val == "true" || val == "1");
            }

            label_propagation_config.use_parallel_alg = (threads > 1);
            if (is_set("--lp-refinement-enabled")) {
                std::string val = get("--lp-refinement-enabled");
                label_propagation_config.enabled = (val == "true" || val == "1");
            }
            if (is_set("--lp-refinement-max-iteration")) {
                label_propagation_config.max_iteration = std::stoull(get("--lp-refinement-max-iteration"));
            }
            if (is_set("--lp-refinement-use-parallel-version")) {
                std::string val = get("--lp-refinement-use-parallel-version");
                label_propagation_config.use_parallel_alg = (val == "true" || val == "1");
            }
            if (is_set("--lp-refinement-use-edge-cut")) {
                std::string val = get("--lp-refinement-use-edge-cut");
                label_propagation_config.use_edge_cut = (val == "true" || val == "1");
            }

            if (is_set("--qg-refinement-enabled")) {
                std::string val = get("--qg-refinement-enabled");
                quotient_graph_refinement_config.enabled = (val == "true" || val == "1");
            }
            if (is_set("--qg-refinement-max-iteration")) {
                quotient_graph_refinement_config.max_iteration = std::stoull(get("--qg-refinement-max-iteration"));
            }
            if (is_set("--qg-refinement-min-n-steps")) {
                quotient_graph_refinement_config.min_n_steps = std::stoull(get("--qg-refinement-min-n-steps"));
            }
            if (is_set("--qg-refinement-alpha")) {
                quotient_graph_refinement_config.alpha = std::stod(get("--qg-refinement-alpha"));
            }
            if (is_set("--qg-refinement-beta")) {
                quotient_graph_refinement_config.beta = std::stod(get("--qg-refinement-beta"));
            }
            if (is_set("--qg-refinement-use-active-scheduling")) {
                std::string val = get("--qg-refinement-use-active-scheduling");
                quotient_graph_refinement_config.use_active_scheduling = (val == "true" || val == "1");
            }
            if (is_set("--qg-refinement-use-preemptive-exit")) {
                std::string val = get("--qg-refinement-use-preemptive-exit");
                quotient_graph_refinement_config.use_preemptive_exit = (val == "true" || val == "1");
            }
            if (is_set("--qg-refinement-use-edge-cut")) {
                std::string val = get("--qg-refinement-use-edge-cut");
                quotient_graph_refinement_config.use_edge_cut = (val == "true" || val == "1");
            }
        }

        void set_fast() {
            initial_c = 32;

            // set GPA matching algorithm
            // coarsening_algorithm_string = "size-constrained-lp";
            coarsening_algorithm_string = "global-paths";
            // coarsening_algorithm_string = "heavy-edge";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;
            global_path_algorithm_config.random_level = 0;
            global_path_algorithm_config.use_edge_rating_tiebreaking = false;

            // configurate global-paths algorithm
            heavy_edge_matching_configuration.rating_function = EdgeRatingFunction::EXPANSIONSTAR;

            size_constrained_lp_config.max_rounds = 1;
            size_constrained_lp_config.min_threshold = 0.10;
            size_constrained_lp_config.f = 16;

            partitioning_algorithm_string = "recursive-bisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            recursive_bisection_config.kappa = 64;
            recursive_bisection_config.use_full_refine = true;
            recursive_bisection_config.method = BisectionMethod::HYBRID;

            recursive_bisection_config.swap_config.enabled = true;
            recursive_bisection_config.swap_config.max_iteration = 3;

            label_propagation_config.enabled = true;
            label_propagation_config.max_iteration = 5;

            quotient_graph_refinement_config.enabled = true;
            quotient_graph_refinement_config.max_iteration = 5;
            quotient_graph_refinement_config.alpha = 5.0;
            quotient_graph_refinement_config.min_n_steps = 3;
            quotient_graph_refinement_config.use_preemptive_exit = true;

            negative_cycle_refinement_config.enabled = false;
            negative_cycle_refinement_config.max_iterations = 10;
            negative_cycle_refinement_config.random_tries = 10;
            negative_cycle_refinement_config.max_path_length = 8;
            negative_cycle_refinement_config.threshold = 100;
        }

        void set_eco() {
            initial_c = 8;

            coarsening_algorithm_string = "global-paths";
            global_path_algorithm_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            global_path_algorithm_config.random_level = 4;

            size_constrained_lp_config.max_rounds = 1;
            size_constrained_lp_config.min_threshold = 0.10;
            size_constrained_lp_config.f = 16;

            partitioning_algorithm_string = "recursive-bisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            recursive_bisection_config.kappa = 3;
            recursive_bisection_config.use_full_refine = true;
            recursive_bisection_config.method = BisectionMethod::HYBRID;
            kaffpa_partition_mode = KAFFPA_PARTITION_ECO;

            label_propagation_config.enabled = true;
            label_propagation_config.max_iteration = 5;

            quotient_graph_refinement_config.enabled = true;
            quotient_graph_refinement_config.max_iteration = 2;
            quotient_graph_refinement_config.alpha = 5.0;
            quotient_graph_refinement_config.min_n_steps = 25;
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
            global_path_algorithm_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            global_path_algorithm_config.random_level = 4;

            size_constrained_lp_config.max_rounds = 5;
            size_constrained_lp_config.min_threshold = 0.10;
            size_constrained_lp_config.f = 16;

            partitioning_algorithm_string = "recursive-bisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            recursive_bisection_config.kappa = 10;
            recursive_bisection_config.use_full_refine = true;
            recursive_bisection_config.method = BisectionMethod::HYBRID;
            kaffpa_partition_mode = KAFFPA_PARTITION_STRONG;
            recursive_bisection_config.swap_config.enabled = true;
            recursive_bisection_config.swap_config.max_iteration = 5;

            label_propagation_config.enabled = true;
            label_propagation_config.max_iteration = 5;

            quotient_graph_refinement_config.enabled = true;
            quotient_graph_refinement_config.max_iteration = 2;
            quotient_graph_refinement_config.alpha = 5.0;
            quotient_graph_refinement_config.min_n_steps = 50;
            quotient_graph_refinement_config.use_preemptive_exit = true;

            flow_based_refinement_config.enabled = true;
            flow_based_refinement_config.use_active_block_scheduling = true;
            flow_based_refinement_config.max_global_iteration = 10;
            flow_based_refinement_config.max_local_iteration = 5;
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
            size_constrained_lp_config.f = 16;

            partitioning_algorithm_string = "recursive-bisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            recursive_bisection_config.kappa = 10;
            recursive_bisection_config.use_full_refine = true;
            recursive_bisection_config.method = BisectionMethod::HYBRID;
            kaffpa_partition_mode = KAFFPA_PARTITION_STRONG;
            recursive_bisection_config.swap_config.enabled = true;
            recursive_bisection_config.swap_config.max_iteration = 5;

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
