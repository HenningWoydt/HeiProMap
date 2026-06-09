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

#include "definitions.h"

namespace HeiProMap {
    enum COARSENING_ALGS {
        COARSENING_ALG_UNDEFINED,
        COARSENING_ALG_GLOBAL_PATHS,
        COARSENING_ALG_SIZE_CONSTRAINED_LP,
        COARSENING_ALG_HEAVY_EDGE
    };

    inline COARSENING_ALGS string_to_coarsening_algorithm(const std::string &str) {
        if (str == "UNDEFINED") return COARSENING_ALG_UNDEFINED;
        if (str == "global-paths") return COARSENING_ALG_GLOBAL_PATHS;
        if (str == "size-constrained-lp") return COARSENING_ALG_SIZE_CONSTRAINED_LP;
        if (str == "heavy-edge") return COARSENING_ALG_HEAVY_EDGE;
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

    enum PARTITIONING_ALGS {
        PARTITIONING_ALG_UNDEFINED,
        PARTITIONING_ALG_MULTISECTION,
        PARTITIONING_ALG_RECURSIVE_BISECTION,
        PARTITIONING_ALG_HEIPA,
        PARTITIONING_ALG_GREEDY,
        PARTITIONING_ALG_GREEDY_GRAPH_GROWING,
        PARTITIONING_ALG_HYBRID,
    };

    inline PARTITIONING_ALGS string_to_partitioning_algorithm(const std::string &str) {
        if (str == "UNDEFINED") return PARTITIONING_ALG_UNDEFINED;
        if (str == "multisection") return PARTITIONING_ALG_MULTISECTION;
        if (str == "recursive-bisection") return PARTITIONING_ALG_RECURSIVE_BISECTION;
        if (str == "heipa") return PARTITIONING_ALG_HEIPA;
        if (str == "greedy") return PARTITIONING_ALG_GREEDY;
        if (str == "greedy-graph-growing") return PARTITIONING_ALG_GREEDY_GRAPH_GROWING;
        if (str == "hybrid") return PARTITIONING_ALG_HYBRID;
        return PARTITIONING_ALG_UNDEFINED;
    }

    inline std::string partitioning_algorithm_to_string(PARTITIONING_ALGS alg) {
        switch (alg) {
            case PARTITIONING_ALG_UNDEFINED: return "UNDEFINED";
            case PARTITIONING_ALG_MULTISECTION: return "multisection";
            case PARTITIONING_ALG_RECURSIVE_BISECTION: return "recursive-bisection";
            case PARTITIONING_ALG_HEIPA: return "heipa";
            case PARTITIONING_ALG_GREEDY: return "greedy";
            case PARTITIONING_ALG_GREEDY_GRAPH_GROWING: return "greedy-graph-growing";
            case PARTITIONING_ALG_HYBRID: return "hybrid";
            default: return "UNDEFINED";
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
}

#include "utility/utils.h"
#include "partitioning/global_multisection.h"
#include "refinement/label_propagation_refinement.h"
#include "refinement/quotient_graph_refinement.h"
#include "coarsening/global_path_algorithm.h"
#include "coarsening/size_constrained_lp.h"
#include "coarsening/heavy_edge_matching.h"
#include "partitioning/kaffpa_partitioner.h"


namespace HeiProMap {
    class AlgorithmConfiguration {
    public:
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
        f64 per_level_imb_add = 0.005;

        // random initialization
        u64 seed = 0;

        u64 threads = 1;

        bool collect_dataset = false;
        std::string data_dir = std::string(HEIPROMAP_SOURCE_DIR) + "/data";

        u64 hm_level = 0;

        u64 initial_c = 8;

        // coarsening algorithm
        std::string coarsening_algorithm_string;
        COARSENING_ALGS coarsening_algorithm_id = COARSENING_ALG_UNDEFINED;

        GlobalPathAlgorithmConfiguration global_path_algorithm_config;
        SizeConstrainedLPConfiguration size_constrained_lp_config;
        HeavyEdgeMatchingConfiguration heavy_edge_matching_config;

        // partitioning algorithm
        std::string partitioning_algorithm_string;
        PARTITIONING_ALGS partitioning_algorithm_id = PARTITIONING_ALG_UNDEFINED;

        u64 rb_kappa = 1;
        bool rb_use_full_refine = false;

        GlobalMultisectionConfiguration global_multisection_config;
        std::string heipa_config_string = "strong";

        // rebalance algorithm
        std::string rebalancing_algorithm_string;
        REBALANCING_ALGS rebalancing_algorithm_id = REBALANCING_ALG_UNDEFINED;

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

        void set_hm_level() {
            if (is_set("--hm-level")) {
                hm_level = std::stoull(get("--hm-level"));
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
            if (is_set("--mapping")) {
                mapping_out = get("--mapping");
            }

            set_hierarchy();
            set_distance();
            set_imbalance();
            set_seed();
            set_threads();
            set_hm_level();

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
        }

        void set_fast() {
            collect_dataset = false;

            initial_c = 32;

            // set GPA matching algorithm
            // coarsening_algorithm_string = "size-constrained-lp";
            coarsening_algorithm_string = "global-paths";
            // coarsening_algorithm_string = "heavy-edge";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;
            global_path_algorithm_config.random_level = 0;
            global_path_algorithm_config.use_adaptive_max_vertex_weight = false;
            global_path_algorithm_config.use_edge_rating_tiebreaking = false;

            heavy_edge_matching_config.rating_function = EdgeRatingFunction::EXPANSION;

            size_constrained_lp_config.max_rounds = 3;
            size_constrained_lp_config.min_threshold = 0.10;
            size_constrained_lp_config.multiplier = 8;
            size_constrained_lp_config.rating_function = EdgeRatingFunction::EXPANSIONSTARSTAR;

            // set multisection
            partitioning_algorithm_string = "multisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            global_multisection_config.mode_string = "heipa-fast";
            global_multisection_config.mode = string_to_global_multisection_mode(global_multisection_config.mode_string);
            global_multisection_config.kappa = 1;
            global_multisection_config.refine = true;
            global_multisection_config.v_cycles = 2;
            global_multisection_config.v_cycle_depth = 5;

            // enable label propagation
            label_propagation_config.enabled = true;
            label_propagation_config.max_iteration = 5;

            // enable quotient graph refinement
            quotient_graph_refinement_config.enabled = true;
            quotient_graph_refinement_config.max_iteration = 2;
            quotient_graph_refinement_config.alpha = 5.0;
            quotient_graph_refinement_config.min_n_steps = 3;
            quotient_graph_refinement_config.use_preemptive_exit = true;

            // enable flow based refinement
            flow_based_refinement_config.enabled = false;
            flow_based_refinement_config.max_global_iteration = 1;
            flow_based_refinement_config.max_local_iteration = 1;
            flow_based_refinement_config.alpha = 1.0;
            flow_based_refinement_config.alpha_upper_bound = 64.0;
            flow_based_refinement_config.alpha_modifier = 2.0;
            flow_based_refinement_config.use_closed_vertex_set = true;
            flow_based_refinement_config.closed_vertex_sets_repeats = 500;
        }

        void set_eco() {
            initial_c = 32;

            // set GPA matching algorithm
            // coarsening_algorithm_string = "size-constrained-lp";
            coarsening_algorithm_string = "global-paths";
            // coarsening_algorithm_string = "heavy-edge";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;
            global_path_algorithm_config.random_level = 0;
            global_path_algorithm_config.use_adaptive_max_vertex_weight = false;
            global_path_algorithm_config.use_edge_rating_tiebreaking = false;
            global_path_algorithm_config.two_hop_threshold = 0.75;

            heavy_edge_matching_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;

            size_constrained_lp_config.max_rounds = 1;
            size_constrained_lp_config.min_threshold = 0.10;
            size_constrained_lp_config.multiplier = 4;

            // set multisection
            partitioning_algorithm_string = "multisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            global_multisection_config.mode_string = "heipa-fast";
            global_multisection_config.mode = string_to_global_multisection_mode(global_multisection_config.mode_string);
            global_multisection_config.kappa = 1;
            global_multisection_config.refine = true;
            global_multisection_config.v_cycles = 5;
            global_multisection_config.v_cycle_depth = 10;

            // enable label propagation
            label_propagation_config.enabled = true;
            label_propagation_config.max_iteration = 5;

            // enable quotient graph refinement
            quotient_graph_refinement_config.enabled = true;
            quotient_graph_refinement_config.max_iteration = 5;
            quotient_graph_refinement_config.alpha = 5.0;
            quotient_graph_refinement_config.min_n_steps = 3;
            quotient_graph_refinement_config.use_preemptive_exit = true;

            // enable flow based refinement
            flow_based_refinement_config.enabled = true;
            flow_based_refinement_config.max_global_iteration = 1;
            flow_based_refinement_config.max_local_iteration = 2;
            flow_based_refinement_config.alpha = 1.0;
            flow_based_refinement_config.alpha_upper_bound = 64.0;
            flow_based_refinement_config.alpha_modifier = 2.0;
            flow_based_refinement_config.use_closed_vertex_set = true;
            flow_based_refinement_config.closed_vertex_sets_repeats = 500;
        }

        void set_strong() {
            initial_c = 32;

            // set GPA matching algorithm
            // coarsening_algorithm_string = "size-constrained-lp";
            coarsening_algorithm_string = "global-paths";
            // coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;
            global_path_algorithm_config.random_level = 0;
            global_path_algorithm_config.use_adaptive_max_vertex_weight = false;
            global_path_algorithm_config.use_edge_rating_tiebreaking = false;
            global_path_algorithm_config.two_hop_threshold = 0.75;

            heavy_edge_matching_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;

            size_constrained_lp_config.max_rounds = 1;
            size_constrained_lp_config.min_threshold = 0.10;
            size_constrained_lp_config.multiplier = 4;

            // set multisection
            partitioning_algorithm_string = "multisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            global_multisection_config.mode_string = "heipa-fast";
            global_multisection_config.mode = string_to_global_multisection_mode(global_multisection_config.mode_string);
            global_multisection_config.kappa = 1;
            global_multisection_config.refine = true;
            global_multisection_config.v_cycles = 5;
            global_multisection_config.v_cycle_depth = 10;

            // enable label propagation
            label_propagation_config.enabled = true;
            label_propagation_config.max_iteration = 5;

            // enable quotient graph refinement
            quotient_graph_refinement_config.enabled = true;
            quotient_graph_refinement_config.max_iteration = 5;
            quotient_graph_refinement_config.alpha = 5.0;
            quotient_graph_refinement_config.min_n_steps = 3;
            quotient_graph_refinement_config.use_preemptive_exit = true;

            // enable flow based refinement
            flow_based_refinement_config.enabled = true;
            flow_based_refinement_config.use_active_block_scheduling = true;
            flow_based_refinement_config.max_global_iteration = 2;
            flow_based_refinement_config.max_local_iteration = 5;
            flow_based_refinement_config.alpha = 1.0;
            flow_based_refinement_config.alpha_upper_bound = 16.0;
            flow_based_refinement_config.alpha_modifier = 2.0;
            flow_based_refinement_config.use_closed_vertex_set = false;
            flow_based_refinement_config.closed_vertex_sets_repeats = 500;
            flow_based_refinement_config.always_include_boundary = true;
            flow_based_refinement_config.growth_strategy = GrowthStrategy::BFS;
        }

        void set_super_strong() {
            initial_c = 16;

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
            global_multisection_config.mode_string = "heipa-super-strong";
            global_multisection_config.mode = string_to_global_multisection_mode(global_multisection_config.mode_string);
            global_multisection_config.kappa = 1;
            global_multisection_config.v_cycles = 3;
            global_multisection_config.v_cycle_depth = 1;

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
            std::cout << "Command Line " << var << " is not an allowed name!" << std::endl;
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
