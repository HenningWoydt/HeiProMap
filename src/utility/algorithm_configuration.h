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
#include "../coarsening/greedy_edge_matcher.h"
#include "../partitioning/global_multisection.h"
#include "../refinement/ILP_refinement.h"
#include "../refinement/k_way_fm_refinement.h"
#include "../refinement/label_propagation_refinement.h"
#include "../refinement/multi_try_fm_refinement.h"
#include "../refinement/quotient_graph_refinement.h"
#include "../refinement/two_vertex_label_propagation_refinement.h"
#include "../coarsening/heavy_edge_matcher.h"
#include "../coarsening/global_path_algorithm.h"
#include "../coarsening/random_edge_matcher.h"
#include "../partitioning/kaffpa_partitioner.h"
#include "../refinement/lightning_refinement.h"
#include "../refinement/three_vertex_label_propagation_refinement.h"
#include "../refinement/wave_refinement.h"


namespace HeiProMap {
    enum COARSENING_ALGS {
        COARSENING_ALG_UNDEFINED,
        COARSENING_ALG_GREEDY_MATCHING,
        COARSENING_ALG_HEAVY_MATCHING,
        COARSENING_ALG_RANDOM_MATCHING,
        COARSENING_ALG_GLOBAL_PATHS,
        COARSENING_ALG_SIZE_CONSTRAINED_LP
    };

    inline COARSENING_ALGS string_to_coarsening_algorithm(const std::string &str) {
        if (str == "UNDEFINED") return COARSENING_ALG_UNDEFINED;
        if (str == "greedy-matching") return COARSENING_ALG_GREEDY_MATCHING;
        if (str == "heavy-matching") return COARSENING_ALG_HEAVY_MATCHING;
        if (str == "random-matching") return COARSENING_ALG_RANDOM_MATCHING;
        if (str == "global-paths") return COARSENING_ALG_GLOBAL_PATHS;
        if (str == "size-constrained-lp") return COARSENING_ALG_SIZE_CONSTRAINED_LP;
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
            case COARSENING_ALG_RANDOM_MATCHING:
                return "random-matching";
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

            /** Coarsening */
            {"--coarsening-algorithm", "", "Which coarsening algorithm to use. Allowed values are {greedy-matching, heavy-matching, random-matching, global-paths}.", "global-paths", "", false},

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
            // Refinement label propagation
            {"--refinement-label-propagation-enable", "", "Enables the label propagation refinement.", "0", "", false},
            {"--refinement-label-propagation-max-iterations", "", "For how many iterations to run label propagation refinement.", "25", "", false},

            // Refinement quotient graph
            {"--refinement-quotient-graph-enable", "", "Enables the quotient graph refinement.", "0", "", false},

            // Refinement k-Way FM
            {"--refinement-k-way-fm-enable", "", "Enables the K-Way FM refinement.", "0", "", false},
            {"--refinement-k-way-fm-max-iterations", "", "How many iterations to run K-Way FM refinement at most.", "1", "", false},

            // Refinement Multi-Try FM
            {"--refinement-multi-try-fm-enable", "", "Enables the Multi-Try FM refinement.", "0", "", false},
            {"--refinement-multi-try-fm-max-iterations", "", "How many iterations to run Multi-Try FM refinement at most.", "1", "", false},

            // Refinement Hierarchy Aware Cycles
            {"--refinement-hierarchy-aware-multi-way-fm-enable", "", "Enables the Hierarchy Aware Multi-Way FM refinement.", "0", "", false},

            // Refinement Two Vertex Label Propagation
            {"--refinement-two-vertex-label-propagation-enable", "", "Enables Label Propagation with two vertices.", "0", "", false},
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

        u64 n_max_partitions = 1;
        u64 n_v_cycle = 0;
        bool first_fast_v_cycle = false;
        u64 v_cycle_max_depth = 100;
        u64 initial_c = 8;

        // coarsening algorithm
        std::string coarsening_algorithm_string;
        COARSENING_ALGS coarsening_algorithm_id = COARSENING_ALG_UNDEFINED;

        GreedyEdgeMatcherConfiguration greedy_edge_matcher_config;
        HeavyEdgeMatcherConfiguration heavy_edge_matcher_config;
        GlobalPathAlgorithmConfiguration global_path_algorithm_config;
        RandomEdgeMatcherConfiguration random_edge_matcher_config;
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
        KWayFMRefinementConfiguration k_way_fm_refinement_config = KWayFMRefinementConfiguration("K-Way-FM");
        MultiTryFmRefinementConfiguration multi_try_fm_refinement_config = MultiTryFmRefinementConfiguration("Multi Try-FM");
        // TwoVertexLabelPropagationConfiguration two_vertex_label_propagation_config = TwoVertexLabelPropagationConfiguration("Two Vertex Label");
        // ThreeVertexLabelPropagationConfiguration three_vertex_label_propagation_config = ThreeVertexLabelPropagationConfiguration("Three Vertex Label");
        FlowBasedRefinementConfiguration flow_based_refinement_config = FlowBasedRefinementConfiguration("Flow Based");
        // ILPRefinementConfiguration ilp_refinement_configuration                        = ILPRefinementConfiguration("ILP Refinement");

        // WaveRefinementConfiguration wave_refinement_configuration = WaveRefinementConfiguration("Wave Refinement");
        // LightningRefinementConfiguration lightning_refinement_configuration = LightningRefinementConfiguration("Lightning Refinement");

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

            // initialize K-Way FM refinement configuration
            if (use_default || is_set("--refinement-k-way-fm-max-iterations")) {
                k_way_fm_refinement_config.max_iteration = std::stoi(get("--refinement-k-way-fm-max-iterations"));
            }
            if (use_default || is_set("--refinement-k-way-fm-enable")) {
                k_way_fm_refinement_config.enabled = get("--refinement-k-way-fm-enable") == "1";
            }

            // initialize Multi-Try FM refinement configuration
            if (use_default || is_set("--refinement-multi-try-fm-max-iterations")) {
                multi_try_fm_refinement_config.max_iteration = std::stoi(get("--refinement-multi-try-fm-max-iterations"));
            }
            if (use_default || is_set("--refinement-multi-try-fm-enable")) {
                multi_try_fm_refinement_config.enabled = get("--refinement-multi-try-fm-enable") == "1";
            }

            // initialize two vertex label propagation
            if (use_default || is_set("--refinement-two-vertex-label-propagation-enable")) {
                // two_vertex_label_propagation_config.enabled = get("--refinement-two-vertex-label-propagation-enable") == "1";
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
            coarsening_algorithm_string = "global-paths";
            // coarsening_algorithm_string = "size-constrained-lp";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.random_level = 4;

            size_constrained_lp_config.max_rounds = 3;
            size_constrained_lp_config.min_threshold = 0.10;

            // set multisection
            partitioning_algorithm_string = "multisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            global_multisection_config.mode_string = "kaffpa-fast";
            // global_multisection_config.mode_string = "metis-recursive";
            // global_multisection_config.mode_string = "metis-kway";
            // global_multisection_config.mode_string = "mtkahypar-default";
            // global_multisection_config.mode_string = "mtkahypar-quality";
            // global_multisection_config.mode_string = "mtkahypar-highestquality";
            global_multisection_config.mode = string_to_global_multisection_mode(global_multisection_config.mode_string);
            global_multisection_config.kappa = 3;

            // enable label propagation
            label_propagation_config.enabled = true;
            label_propagation_config.max_iteration = 5;

            // enable multi-try fm
            multi_try_fm_refinement_config.enabled = false;
            multi_try_fm_refinement_config.max_iteration = 2;
            multi_try_fm_refinement_config.alpha = 100.0;

            // enable quotient graph refinement
            quotient_graph_refinement_config.enabled = false;
            quotient_graph_refinement_config.max_iteration = 5;
            quotient_graph_refinement_config.alpha = 1000.0;
            quotient_graph_refinement_config.min_n_steps = 8;

            // enable k-way fm
            k_way_fm_refinement_config.enabled = false;
            k_way_fm_refinement_config.max_iteration = 1;
            k_way_fm_refinement_config.alpha = 100.0;
        }

        void set_eco() {
            initial_c = 8;

            // set GPA matching algorithm
            coarsening_algorithm_string = "size-constrained-lp";
            // coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.random_level = 0;

            size_constrained_lp_config.max_rounds = 2;
            size_constrained_lp_config.min_threshold = 0.10;

            // set multisection
            partitioning_algorithm_string = "multisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            // global_multisection_config.mode_string = "kaffpa-strong";
            global_multisection_config.mode_string = "kaffpa-eco";
            // global_multisection_config.mode_string = "kaffpa-fast";
            // global_multisection_config.mode_string = "metis-recursive";
            // global_multisection_config.mode_string = "metis-kway";
            // global_multisection_config.mode_string = "mtkahypar-default";
            // global_multisection_config.mode_string = "mtkahypar-quality";
            // global_multisection_config.mode_string = "mtkahypar-highestquality";
            global_multisection_config.mode = string_to_global_multisection_mode(global_multisection_config.mode_string);
            global_multisection_config.kappa = 3;

            // enable quotient graph refinement
            quotient_graph_refinement_config.enabled = true;
            quotient_graph_refinement_config.max_iteration = 2;
            quotient_graph_refinement_config.alpha = 10000.0;
            quotient_graph_refinement_config.min_n_steps = 32;

            // enable multi-try fm
            multi_try_fm_refinement_config.enabled = true;
            multi_try_fm_refinement_config.max_iteration = 3;
            multi_try_fm_refinement_config.alpha = 10000.0;
            multi_try_fm_refinement_config.min_n_steps = 8;

            // enable flow based refinement
            flow_based_refinement_config.enabled = true;
            flow_based_refinement_config.max_global_iteration = 2;
            flow_based_refinement_config.max_local_iteration = 2;
            flow_based_refinement_config.alpha = 1.0;
            flow_based_refinement_config.alpha_upper_bound = 64.0;
            flow_based_refinement_config.alpha_modifier = 2.0;
            flow_based_refinement_config.use_closed_vertex_set = true;
            flow_based_refinement_config.closed_vertex_sets_repeats = 500;
        }

        void set_strong() {
            initial_c = 16;

            n_max_partitions = 1;
            n_refinement_iterations = 1;
            n_v_cycle = 2;
            first_fast_v_cycle = false;
            v_cycle_max_depth = 5;

            // set GPA matching algorithm
            coarsening_algorithm_string = "size-constrained-lp";
            // coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.random_level = 0;

            size_constrained_lp_config.max_rounds = 5;
            size_constrained_lp_config.min_threshold = 0.10;

            // set multisection
            partitioning_algorithm_string = "multisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            // global_multisection_config.mode_string = "kaffpa-fast";
            global_multisection_config.mode_string = "kaffpa-strong";
            // global_multisection_config.mode_string = "kaffpa-eco";
            // global_multisection_config.mode_string = "metis-recursive";
            // global_multisection_config.mode_string = "metis-kway";
            // global_multisection_config.mode_string = "mtkahypar-default";
            // global_multisection_config.mode_string = "mtkahypar-quality";
            // global_multisection_config.mode_string = "mtkahypar-highestquality";
            global_multisection_config.mode = string_to_global_multisection_mode(global_multisection_config.mode_string);
            global_multisection_config.kappa = 1;

            // enable label propagation
            label_propagation_config.enabled = false;
            label_propagation_config.max_iteration = 5;

            // enable quotient graph refinement
            quotient_graph_refinement_config.enabled = false;
            quotient_graph_refinement_config.max_iteration = 2;
            quotient_graph_refinement_config.alpha = 1000.0;
            quotient_graph_refinement_config.min_n_steps = 10;

            // enable k-way fm
            k_way_fm_refinement_config.enabled = false;

            // enable multi-try fm
            multi_try_fm_refinement_config.enabled = true;
            multi_try_fm_refinement_config.max_iteration = 2;
            multi_try_fm_refinement_config.alpha = 1000.0;
            multi_try_fm_refinement_config.min_n_steps = 10;

            // enable flow based refinement
            flow_based_refinement_config.enabled = true;
            flow_based_refinement_config.max_global_iteration = 1;
            flow_based_refinement_config.max_local_iteration = 3;
            flow_based_refinement_config.alpha = 2.0;
            flow_based_refinement_config.alpha_upper_bound = 16.0;
            flow_based_refinement_config.alpha_modifier = 2.0;
            flow_based_refinement_config.use_closed_vertex_set = true;
            flow_based_refinement_config.closed_vertex_sets_repeats = 500;

            // wave_refinement_configuration.enabled = false;

            // lightning_refinement_configuration.enabled = false;
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

#endif //HEIPROMAP_ALGORITHMCONFIGURATION_H
