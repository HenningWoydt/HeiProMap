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

#include "../../commons/utils.h"
#include "../../definitions.h"
#include "../coarsening/greedy_edge_matcher.h"
#include "../partitioning/global_multisection.h"
#include "../refinement/hierarchy_aware_cycle_refinement.h"
#include "../refinement/k_way_fm_refinement.h"
#include "../refinement/k_way_fm_refinement_Faraj20.h"
#include "../refinement/label_propagation_refinement.h"
#include "../refinement/multi_try_fm_refinement.h"
#include "../refinement/multi_try_fm_refinement_Faraj20.h"
#include "../refinement/quotient_graph_refinement.h"
#include "../refinement/two_vertex_label_propagation_refinement.h"

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
        PARTITIONING_ALG_KAFFPA,
        PARTITIONING_ALG_MULTISECTION,
    };

    inline PARTITIONING_ALGS string_to_partitioning_algorithm(const std::string& str) {
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

    inline REBALANCING_ALGS string_to_rebalancing_algorithm(const std::string& str) {
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
                {"--coarsening-algorithm", "", "Which coarsening algorithm to use. Allowed values are {greedy-matching, heavy-matching, global-paths}.", "global-paths", "", false},

                // Coarsening global-path
                {"--coarsening-algorithm-global-paths-random-level", "", "On which levels to run random if global-paths is chosen. Smaller-equal than use random, greater than use GPA.", "4", "", false},

                // Coarsening greedy matching
                {"--coarsening-algorithm-greedy-matching-pendant-first", "", "Whether the greedy matching algorithm should handle pendant vertices first. 1 enables first matching pendant vertices, while 0 does not.", "1", "", false},

                // Coarsening heavy matching
                {"--coarsening-algorithm-heavy-matching-pendant-first", "", "Whether the heavy matching algorithm should handle pendant vertices first. 1 enables first matching pendant vertices, while 0 does not.", "1", "", false},

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
                // Refinement Faraj20 label propagation
                {"--refinement-lable-propagation-faraj20-enable", "", "Enables the label propagation refinement by Faraj20.", "0", "", false},
                {"--refinement-lable-propagation-faraj20-max-iterations", "", "For how many iterations to run label propagation refinement by Faraj20.", "25", "", false},

                // Refinement Faraj20 quotient graph
                {"--refinement-quotient-graph-faraj20-enable", "", "Enables the quotient graph refinement by Faraj20.", "0", "", false},
                {"--refinement-quotient-graph-faraj20-max-iterations", "", "How many iterations to run quotient graph refinement by Faraj20 at most.", "1", "", false},
                {"--refinement-quotient-graph-faraj20-max-moves-without-max", "", "After x moves without discovering a new maximum quotient graph refinement by Faraj20 will stop the fm call.", "5000", "", false},

                // Refinement Faraj20 k-Way FM
                {"--refinement-k-way-fm-faraj20-enable", "", "Enables the K-Way FM refinement by Faraj20.", "0", "", false},
                {"--refinement-k-way-fm-faraj20-max-iterations", "", "How many iterations to run K-Way FM by Faraj20 refinement at most.", "1", "", false},

                // Refinement Faraj20 Multi-Try FM
                {"--refinement-multi-try-fm-faraj20-enable", "", "Enables the Multi-Try FM refinement by Faraj20.", "0", "", false},
                {"--refinement-multi-try-fm-faraj20-max-iterations", "", "How many iterations to run Multi-Try FM refinement by Faraj20 at most.", "1", "", false},

                // Refinement label propagation
                {"--refinement-lable-propagation-enable", "", "Enables the label propagation refinement.", "0", "", false},
                {"--refinement-lable-propagation-max-iterations", "", "For how many iterations to run label propagation refinement.", "25", "", false},

                // Refinement quotient graph
                {"--refinement-quotient-graph-enable", "", "Enables the quotient graph refinement.", "0", "", false},
                {"--refinement-quotient-graph-max-iterations", "", "How many iterations to run quotient graph refinement at most.", "1", "", false},
                {"--refinement-quotient-graph-max-moves-without-max", "", "After x moves without discovering a new maximum quotient graph refinement will stop the fm call.", "5000", "", false},

                // Refinement k-Way FM
                {"--refinement-k-way-fm-enable", "", "Enables the K-Way FM refinement.", "0", "", false},
                {"--refinement-k-way-fm-max-iterations", "", "How many iterations to run K-Way FM refinement at most.", "1", "", false},

                // Refinement Multi-Try FM
                {"--refinement-multi-try-fm-enable", "", "Enables the Multi-Try FM refinement.", "0", "", false},
                {"--refinement-multi-try-fm-max-iterations", "", "How many iterations to run Multi-Try FM refinement at most.", "1", "", false},

                // Refinement Hierarchy Aware Cycles
                {"--refinement-hierarchy-aware-cycles-enable", "", "Enables the hierarchy aware cycles refinement.", "0", "", false},

                // Refinement Two Vertex Label Propagation
                {"--refinement-two-vertex-label-propagation-enable", "", "Enables Label Propagation with two vertices.", "0", "", false},
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
        HeavyEdgeMatcherConfiguration heavy_edge_matcher_config;
        GlobalPathAlgorithmConfiguration global_path_algorithm_config;

        // partitioning algorithm
        std::string partitioning_algorithm_string;
        PARTITIONING_ALGS partitioning_algorithm_id = PARTITIONING_ALG_UNDEFINED;

        GlobalMultisectionConfiguration global_multisection_config;
        KaffpaPartitionerConfiguration kaffpa_partitioner_config;

        // rebalance algorithm
        std::string rebalancing_algorithm_string;
        REBALANCING_ALGS rebalancing_algorithm_id = REBALANCING_ALG_UNDEFINED;

        SimpleRebalancerConfiguration simple_rebalancer_configuration;

        // refinement algorithms
        bool do_refinement_label_propagation_faraj20 = false;
        LabelPropagationFaraj20Configuration label_propagation_faraj20_config;

        bool do_refinement_quotient_graph_faraj20 = false;
        QuotientGraphRefinementFaraj20Configuration quotient_graph_refinement_faraj20_config;

        bool do_refinement_k_way_fm_faraj20 = false;
        KWayFMRefinementFaraj20Configuration k_way_fm_refinement_faraj20_config;

        bool do_refinement_multi_try_fm_faraj20 = false;
        MultiTryFmRefinementFaraj20Configuration multi_try_fm_refinement_faraj20_config;

        bool do_refinement_label_propagation = false;
        LabelPropagationConfiguration label_propagation_config;

        bool do_refinement_quotient_graph = false;
        QuotientGraphRefinementConfiguration quotient_graph_refinement_config;

        bool do_refinement_k_way_fm = false;
        KWayFMRefinementConfiguration k_way_fm_refinement_config;

        bool do_refinement_multi_try_fm = false;
        MultiTryFmRefinementConfiguration multi_try_fm_refinement_config;

        bool do_refinement_hierarchy_aware_cycles_enable = false;
        HierarchyAwareCyclesConfiguration hierarchy_aware_cycles_config;

        bool do_refinement_two_vertex_label_propagation_enable = false;
        TwoVertexLabelPropagationConfiguration two_vertex_label_propagation_config;

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

        void set_coarsening_algorithm(const bool use_default = false) {
            // initialize greedy matching config
            if (use_default || is_set("--coarsening-algorithm-greedy-matching-pendant-first")) {
                greedy_edge_matcher_config.match_pendant_vertices_first = get("--coarsening-algorithm-greedy-matching-pendant-first") == "1";
            }

            // initialize heavy matching config
            if (use_default || is_set("--coarsening-algorithm-heavy-matching-pendant-first")) {
                heavy_edge_matcher_config.match_pendant_vertices_first = get("--coarsening-algorithm-heavy-matching-pendant-first") == "1";
            }

            // initialize global-paths config
            if (use_default || is_set("--coarsening-algorithm-global-paths-random-level")) {
                global_path_algorithm_config.random_level = std::stoi(get("--coarsening-algorithm-global-paths-random-level"));
            }

            // actually set which algorithm to use
            if (use_default || is_set("--coarsening-algorithm")) {
                coarsening_algorithm_string = get("--coarsening-algorithm");
                coarsening_algorithm_id     = string_to_coarsening_algorithm(coarsening_algorithm_string);
            }
        }

        void set_partitioning_algorithm(const bool use_default = false) {
            // initialize global multisection config
            if (use_default || is_set("--partitioning-algorithm-multisection-mode")) {
                global_multisection_config.mode_string = get("--partitioning-algorithm-multisection-mode");
                global_multisection_config.mode        = string_to_global_multisection_mode(global_multisection_config.mode_string);
            }

            // initialize kaffpa partitioning config
            if (use_default || is_set("--partitioning-algorithm-kaffpa-partitioning-mode")) {
                kaffpa_partitioner_config.mode_string = get("--partitioning-algorithm-kaffpa-partitioning-mode");
                kaffpa_partitioner_config.mode        = string_to_kaffpa_partitioner_mode(kaffpa_partitioner_config.mode_string);
            }
            if (use_default || is_set("--partitioning-algorithm-kaffpa-partitioning-method")) {
                kaffpa_partitioner_config.method_string = get("--partitioning-algorithm-kaffpa-partitioning-method");
                kaffpa_partitioner_config.method        = string_to_kaffpa_partitioner_method(kaffpa_partitioner_config.method_string);
            }

            // actually set which algorithm to use
            if (use_default || is_set("--partitioning-algorithm")) {
                partitioning_algorithm_string = get("--partitioning-algorithm");
                partitioning_algorithm_id     = string_to_partitioning_algorithm(partitioning_algorithm_string);
            }
        }

        void set_rebalancing_algorithm(const bool use_default = false) {
            if (use_default || is_set("--rebalancing-algorithm")) {
                rebalancing_algorithm_string = get("--rebalancing-algorithm");
                rebalancing_algorithm_id     = string_to_rebalancing_algorithm(rebalancing_algorithm_string);
            }
        }

        void enable_refinement_algorithms(const bool use_default = false) {
            // initialize label propagation faraj20 configuration
            if (use_default || is_set("--refinement-lable-propagation-faraj20-max-iterations")) {
                label_propagation_faraj20_config.max_iteration = std::stoi(get("--refinement-lable-propagation-faraj20-max-iterations"));
            }
            if (use_default || is_set("--refinement-lable-propagation-faraj20-enable")) {
                do_refinement_label_propagation_faraj20 = get("--refinement-lable-propagation-faraj20-enable") == "1";
            }

            // initialize quotient graph faraj20 configuration
            if (use_default || is_set("--refinement-quotient-graph-faraj20-max-iterations")) {
                quotient_graph_refinement_faraj20_config.max_iteration = std::stoi(get("--refinement-quotient-graph-faraj20-max-iterations"));
            }
            if (use_default || is_set("--refinement-quotient-graph-faraj20-max-moves-without-max")) {
                quotient_graph_refinement_faraj20_config.max_moves_without_max = std::stoi(get("--refinement-quotient-graph-faraj20-max-moves-without-max"));
            }
            if (use_default || is_set("--refinement-quotient-graph-faraj20-enable")) {
                do_refinement_quotient_graph_faraj20 = get("--refinement-quotient-graph-faraj20-enable") == "1";
            }

            // initialize K-Way FM refinement faraj20 configuration
            if (use_default || is_set("--refinement-k-way-fm-faraj20-max-iterations")) {
                k_way_fm_refinement_faraj20_config.max_iteration = std::stoi(get("--refinement-k-way-fm-faraj20-max-iterations"));
            }
            if (use_default || is_set("--refinement-k-way-fm-faraj20-enable")) {
                do_refinement_k_way_fm_faraj20 = get("--refinement-k-way-fm-faraj20-enable") == "1";
            }

            // initialize Multi-Try FM refinement faraj20 configuration
            if (use_default || is_set("--refinement-multi-try-fm-faraj20-max-iterations")) {
                multi_try_fm_refinement_faraj20_config.max_iteration = std::stoi(get("--refinement-multi-try-fm-faraj20-max-iterations"));
            }
            if (use_default || is_set("--refinement-multi-try-fm-faraj20-enable")) {
                do_refinement_multi_try_fm_faraj20 = get("--refinement-multi-try-fm-faraj20-enable") == "1";
            }

            // initialize label propagation configuration
            if (use_default || is_set("--refinement-lable-propagation-max-iterations")) {
                label_propagation_config.max_iteration = std::stoi(get("--refinement-lable-propagation-max-iterations"));
            }
            if (use_default || is_set("--refinement-lable-propagation-enable")) {
                do_refinement_label_propagation = get("--refinement-lable-propagation-enable") == "1";
            }

            // initialize quotient graph configuration
            if (use_default || is_set("--refinement-quotient-graph-max-iterations")) {
                quotient_graph_refinement_config.max_iteration = std::stoi(get("--refinement-quotient-graph-max-iterations"));
            }
            if (use_default || is_set("--refinement-quotient-graph-max-moves-without-max")) {
                quotient_graph_refinement_config.max_moves_without_max = std::stoi(get("--refinement-quotient-graph-max-moves-without-max"));
            }
            if (use_default || is_set("--refinement-quotient-graph-enable")) {
                do_refinement_quotient_graph = get("--refinement-quotient-graph-enable") == "1";
            }

            // initialize K-Way FM refinement configuration
            if (use_default || is_set("--refinement-k-way-fm-max-iterations")) {
                k_way_fm_refinement_config.max_iteration = std::stoi(get("--refinement-k-way-fm-max-iterations"));
            }
            if (use_default || is_set("--refinement-k-way-fm-enable")) {
                do_refinement_k_way_fm = get("--refinement-k-way-fm-enable") == "1";
            }

            // initialize Multi-Try FM refinement configuration
            if (use_default || is_set("--refinement-multi-try-fm-max-iterations")) {
                multi_try_fm_refinement_config.max_iteration = std::stoi(get("--refinement-multi-try-fm-max-iterations"));
            }
            if (use_default || is_set("--refinement-multi-try-fm-enable")) {
                do_refinement_multi_try_fm = get("--refinement-multi-try-fm-enable") == "1";
            }

            // initialize hierarchy aware cycle refinement
            if (use_default || is_set("--refinement-hierarchy-aware-cycles-enable")) {
                do_refinement_hierarchy_aware_cycles_enable = get("--refinement-hierarchy-aware-cycles-enable") == "1";
            }

            // initialize two vertex label propagation
            if (use_default || is_set("--refinement-two-vertex-label-propagation-enable")) {
                do_refinement_hierarchy_aware_cycles_enable = get("--refinement-two-vertex-label-propagation-enable") == "1";
            }
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

            set_coarsening_algorithm(true);
            set_partitioning_algorithm(true);
            set_rebalancing_algorithm(true);
            enable_refinement_algorithms(true);

            // set FARAJ20 config
            if (get("--config") == "Faraj20-fastest") {
                set_Faraj20_fastest();
            } else if (get("--config") == "Faraj20-fast") {
                set_Faraj20_fast();
            } else if (get("--config") == "Faraj20-eco") {
                set_Faraj20_eco();
            } else if (get("--config") == "Faraj20-strong") {
                set_Faraj20_strong();
            } else if (get("--config") == "fastest") {
                set_fastest();
            } else if (get("--config") == "fast") {
                set_fast();
            } else if (get("--config") == "eco") {
                set_eco();
            } else if (get("--config") == "strong") {
                set_strong();
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

        void set_Faraj20_fastest() {
            // set GPA matching algorithm
            coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id     = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.random_level = 4;

            // set kaffpa multisection partitioning
            partitioning_algorithm_string           = "kaffpa";
            partitioning_algorithm_id               = string_to_partitioning_algorithm(partitioning_algorithm_string);
            kaffpa_partitioner_config.mode_string   = "fast";
            kaffpa_partitioner_config.mode          = string_to_kaffpa_partitioner_mode(kaffpa_partitioner_config.mode_string);
            kaffpa_partitioner_config.method_string = "bisection"; // TODO: this should be multisection, but there is a bug
            kaffpa_partitioner_config.method        = string_to_kaffpa_partitioner_method(kaffpa_partitioner_config.method_string);

            // disable all own refinements
            do_refinement_label_propagation = false;
            do_refinement_quotient_graph    = false;
            do_refinement_k_way_fm          = false;
            do_refinement_multi_try_fm      = false;

            // disable all refinements
            do_refinement_label_propagation_faraj20 = false;
            do_refinement_quotient_graph_faraj20    = false;
            do_refinement_k_way_fm_faraj20          = false;
            do_refinement_multi_try_fm_faraj20      = false;
        }

        void set_Faraj20_fast() {
            // set GPA matching algorithm
            coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id     = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.random_level = 4;

            // set kaffpa multisection partitioning
            partitioning_algorithm_string           = "kaffpa";
            partitioning_algorithm_id               = string_to_partitioning_algorithm(partitioning_algorithm_string);
            kaffpa_partitioner_config.mode_string   = "fast";
            kaffpa_partitioner_config.mode          = string_to_kaffpa_partitioner_mode(kaffpa_partitioner_config.mode_string);
            kaffpa_partitioner_config.method_string = "bisection"; // TODO: this should be multisection, but there is a bug
            kaffpa_partitioner_config.method        = string_to_kaffpa_partitioner_method(kaffpa_partitioner_config.method_string);

            // disable all own refinements
            do_refinement_label_propagation = false;
            do_refinement_quotient_graph    = false;
            do_refinement_k_way_fm          = false;
            do_refinement_multi_try_fm      = false;

            // only enable label propagation
            do_refinement_quotient_graph_faraj20 = false;
            do_refinement_k_way_fm_faraj20       = false;
            do_refinement_multi_try_fm_faraj20   = false;

            label_propagation_faraj20_config.max_iteration = 25;
            do_refinement_label_propagation_faraj20        = true;
        }

        void set_Faraj20_eco() {
            // set GPA matching algorithm
            coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id     = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.random_level = 4;

            // set kaffpa multisection partitioning
            partitioning_algorithm_string           = "kaffpa";
            partitioning_algorithm_id               = string_to_partitioning_algorithm(partitioning_algorithm_string);
            kaffpa_partitioner_config.mode_string   = "eco";
            kaffpa_partitioner_config.mode          = string_to_kaffpa_partitioner_mode(kaffpa_partitioner_config.mode_string);
            kaffpa_partitioner_config.method_string = "bisection"; // TODO: this should be multisection, but there is a bug
            kaffpa_partitioner_config.method        = string_to_kaffpa_partitioner_method(kaffpa_partitioner_config.method_string);

            // disable all own refinements
            do_refinement_label_propagation = false;
            do_refinement_quotient_graph    = false;
            do_refinement_k_way_fm          = false;
            do_refinement_multi_try_fm      = false;

            // disable multi try fm
            do_refinement_multi_try_fm_faraj20 = false;

            // enable label propagation
            label_propagation_faraj20_config.max_iteration = 25;
            do_refinement_label_propagation_faraj20        = true;

            // enable quotient graph refinement
            quotient_graph_refinement_faraj20_config.max_iteration         = std::numeric_limits<u64>::max();
            quotient_graph_refinement_faraj20_config.max_moves_without_max = std::numeric_limits<u64>::max();
            do_refinement_quotient_graph_faraj20                           = true;

            // enable k-way fm refinement
            do_refinement_k_way_fm_faraj20 = true;
        }

        void set_Faraj20_strong() {
            // set GPA matching algorithm
            coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id     = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.random_level = 4;

            // set kaffpa multisection partitioning
            partitioning_algorithm_string           = "kaffpa";
            partitioning_algorithm_id               = string_to_partitioning_algorithm(partitioning_algorithm_string);
            kaffpa_partitioner_config.mode_string   = "strong";
            kaffpa_partitioner_config.mode          = string_to_kaffpa_partitioner_mode(kaffpa_partitioner_config.mode_string);
            kaffpa_partitioner_config.method_string = "bisection"; // TODO: this should be multisection, but there is a bug
            kaffpa_partitioner_config.method        = string_to_kaffpa_partitioner_method(kaffpa_partitioner_config.method_string);

            // disable all own refinements
            do_refinement_label_propagation = false;
            do_refinement_quotient_graph    = false;
            do_refinement_k_way_fm          = false;
            do_refinement_multi_try_fm      = false;

            // enable label propagation
            label_propagation_faraj20_config.max_iteration = 25;
            do_refinement_label_propagation_faraj20        = true;

            // enable quotient graph refinement
            quotient_graph_refinement_faraj20_config.max_iteration         = std::numeric_limits<u64>::max();
            quotient_graph_refinement_faraj20_config.max_moves_without_max = std::numeric_limits<u64>::max();
            do_refinement_quotient_graph_faraj20                           = true;

            // enable k-way fm refinement
            do_refinement_k_way_fm_faraj20 = true;

            // enable multi-try fm
            do_refinement_multi_try_fm_faraj20 = true;
        }

        void set_fastest() {
            // set GPA matching algorithm
            coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id     = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.random_level = 4;

            // set multisection
            partitioning_algorithm_string          = "multisection";
            partitioning_algorithm_id              = string_to_partitioning_algorithm(partitioning_algorithm_string);
            global_multisection_config.mode_string = "fast";
            global_multisection_config.mode        = string_to_global_multisection_mode(global_multisection_config.mode_string);

            // set kaffpa multisection partitioning
            partitioning_algorithm_string           = "kaffpa";
            partitioning_algorithm_id               = string_to_partitioning_algorithm(partitioning_algorithm_string);
            kaffpa_partitioner_config.mode_string   = "fast";
            kaffpa_partitioner_config.mode          = string_to_kaffpa_partitioner_mode(kaffpa_partitioner_config.mode_string);
            kaffpa_partitioner_config.method_string = "bisection"; // TODO: this should be multisection, but there is a bug
            kaffpa_partitioner_config.method        = string_to_kaffpa_partitioner_method(kaffpa_partitioner_config.method_string);

            // disable all Faraj20 refinements
            do_refinement_label_propagation_faraj20 = false;
            do_refinement_quotient_graph_faraj20    = false;
            do_refinement_k_way_fm_faraj20          = false;
            do_refinement_multi_try_fm_faraj20      = false;

            // disable all refinements
            do_refinement_label_propagation = false;
            do_refinement_quotient_graph    = false;
            do_refinement_k_way_fm          = false;
            do_refinement_multi_try_fm      = false;
        }

        void set_fast() {
            // set GPA matching algorithm
            coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id     = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.random_level = 4;

            // set multisection
            partitioning_algorithm_string          = "multisection";
            partitioning_algorithm_id              = string_to_partitioning_algorithm(partitioning_algorithm_string);
            global_multisection_config.mode_string = "fast";
            global_multisection_config.mode        = string_to_global_multisection_mode(global_multisection_config.mode_string);

            // set kaffpa multisection partitioning
            partitioning_algorithm_string           = "kaffpa";
            partitioning_algorithm_id               = string_to_partitioning_algorithm(partitioning_algorithm_string);
            kaffpa_partitioner_config.mode_string   = "fast";
            kaffpa_partitioner_config.mode          = string_to_kaffpa_partitioner_mode(kaffpa_partitioner_config.mode_string);
            kaffpa_partitioner_config.method_string = "bisection"; // TODO: this should be multisection, but there is a bug
            kaffpa_partitioner_config.method        = string_to_kaffpa_partitioner_method(kaffpa_partitioner_config.method_string);

            // enable label propagation
            do_refinement_label_propagation        = true;
            label_propagation_config.max_iteration = 25;

            do_refinement_two_vertex_label_propagation_enable = false;
            two_vertex_label_propagation_config.max_iteration = 1;
        }

        void set_eco() {
            // set GPA matching algorithm
            coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id     = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.random_level = 4;

            // set multisection
            partitioning_algorithm_string          = "multisection";
            partitioning_algorithm_id              = string_to_partitioning_algorithm(partitioning_algorithm_string);
            global_multisection_config.mode_string = "eco";
            global_multisection_config.mode        = string_to_global_multisection_mode(global_multisection_config.mode_string);

            // set kaffpa multisection partitioning
            partitioning_algorithm_string           = "kaffpa";
            partitioning_algorithm_id               = string_to_partitioning_algorithm(partitioning_algorithm_string);
            kaffpa_partitioner_config.mode_string   = "fast";
            kaffpa_partitioner_config.mode          = string_to_kaffpa_partitioner_mode(kaffpa_partitioner_config.mode_string);
            kaffpa_partitioner_config.method_string = "bisection"; // TODO: this should be multisection, but there is a bug
            kaffpa_partitioner_config.method        = string_to_kaffpa_partitioner_method(kaffpa_partitioner_config.method_string);

            // enable label propagation
            do_refinement_label_propagation        = true;
            label_propagation_config.max_iteration = 25;

            // enable two vertex label propagation
            do_refinement_two_vertex_label_propagation_enable = false;
            two_vertex_label_propagation_config.max_iteration = 25;

            // enable quotient graph refinement
            do_refinement_quotient_graph                           = true;
            quotient_graph_refinement_config.max_iteration         = std::numeric_limits<u64>::max();
            quotient_graph_refinement_config.max_moves_without_max = std::numeric_limits<u64>::max();

            // enable k-way fm
            do_refinement_k_way_fm = true;
        }

        void set_strong() {
            // set GPA matching algorithm
            coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id     = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.random_level = 4;

            // set multisection
            partitioning_algorithm_string          = "multisection";
            partitioning_algorithm_id              = string_to_partitioning_algorithm(partitioning_algorithm_string);
            global_multisection_config.mode_string = "strong";
            global_multisection_config.mode        = string_to_global_multisection_mode(global_multisection_config.mode_string);

            // set kaffpa multisection partitioning
            partitioning_algorithm_string           = "kaffpa";
            partitioning_algorithm_id               = string_to_partitioning_algorithm(partitioning_algorithm_string);
            kaffpa_partitioner_config.mode_string   = "fast";
            kaffpa_partitioner_config.mode          = string_to_kaffpa_partitioner_mode(kaffpa_partitioner_config.mode_string);
            kaffpa_partitioner_config.method_string = "bisection"; // TODO: this should be multisection, but there is a bug
            kaffpa_partitioner_config.method        = string_to_kaffpa_partitioner_method(kaffpa_partitioner_config.method_string);

            // enable label propagation
            do_refinement_label_propagation        = false;
            label_propagation_config.max_iteration = 25;

            do_refinement_two_vertex_label_propagation_enable = false;
            two_vertex_label_propagation_config.max_iteration = 1;

            // enable quotient graph refinement
            do_refinement_quotient_graph                           = false;
            quotient_graph_refinement_config.max_iteration         = std::numeric_limits<u64>::max();
            quotient_graph_refinement_config.max_moves_without_max = 10;

            // enable k-way fm
            do_refinement_k_way_fm = false;

            // enable multi-try fm
            do_refinement_multi_try_fm = true;
        }

        void set_experimental() {
            // set GPA matching algorithm
            coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id     = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.random_level = 4;

            // set multisection
            partitioning_algorithm_string          = "multisection";
            partitioning_algorithm_id              = string_to_partitioning_algorithm(partitioning_algorithm_string);
            global_multisection_config.mode_string = "eco";
            global_multisection_config.mode        = string_to_global_multisection_mode(global_multisection_config.mode_string);

            // disable multi-try fm
            do_refinement_multi_try_fm = false;

            // enable label propagation
            do_refinement_label_propagation        = true;
            label_propagation_config.max_iteration = 25;

            // enable two vertex label propagation
            do_refinement_two_vertex_label_propagation_enable = false;
            two_vertex_label_propagation_config.max_iteration = 25;

            // enable quotient graph refinement
            do_refinement_quotient_graph                           = true;
            quotient_graph_refinement_config.max_iteration         = std::numeric_limits<u64>::max();
            quotient_graph_refinement_config.max_moves_without_max = std::numeric_limits<u64>::max();

            // enable k-way fm
            do_refinement_k_way_fm = true;

            // enable experimental refinement
            do_refinement_hierarchy_aware_cycles_enable       = false;
            do_refinement_two_vertex_label_propagation_enable = false;
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
