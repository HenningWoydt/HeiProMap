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

    #ifndef HEIPROMAP_PARTITIONING_ALGS_DEFINED
    #define HEIPROMAP_PARTITIONING_ALGS_DEFINED

    enum PARTITIONING_ALGS {
        PARTITIONING_ALG_UNDEFINED,
        PARTITIONING_ALG_MULTISECTION,
        PARTITIONING_ALG_RECURSIVE_BISECTION,
        PARTITIONING_ALG_HEIPA,
        PARTITIONING_ALG_GREEDY,
        PARTITIONING_ALG_GREEDY_GRAPH_GROWING,
        PARTITIONING_ALG_HYBRID,
        PARTITIONING_ALG_FLOW_BISECTION,
    };

    inline PARTITIONING_ALGS string_to_partitioning_algorithm(const std::string &str) {
        if (str == "UNDEFINED") return PARTITIONING_ALG_UNDEFINED;
        if (str == "multisection") return PARTITIONING_ALG_MULTISECTION;
        if (str == "recursive-bisection") return PARTITIONING_ALG_RECURSIVE_BISECTION;
        if (str == "heipa") return PARTITIONING_ALG_HEIPA;
        if (str == "greedy") return PARTITIONING_ALG_GREEDY;
        if (str == "greedy-graph-growing") return PARTITIONING_ALG_GREEDY_GRAPH_GROWING;
        if (str == "hybrid") return PARTITIONING_ALG_HYBRID;
        if (str == "flow-bisection") return PARTITIONING_ALG_FLOW_BISECTION;
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
            case PARTITIONING_ALG_FLOW_BISECTION: return "flow-bisection";
            default: return "UNDEFINED";
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

    enum ORACLE_TYPE {
        ORACLE_AUTOMATIC,
        ORACLE_MATRIX,
        ORACLE_BINARY
    };

    inline ORACLE_TYPE string_to_oracle_type(const std::string &str) {
        if (str == "automatic") return ORACLE_AUTOMATIC;
        if (str == "matrix") return ORACLE_MATRIX;
        if (str == "binary") return ORACLE_BINARY;
        std::cout << "Error: Invalid oracle type '" << str << "'. Allowed values: matrix, binary, automatic" << std::endl;
        exit(EXIT_FAILURE);
    }

    inline std::string oracle_type_to_string(ORACLE_TYPE type) {
        switch (type) {
            case ORACLE_AUTOMATIC: return "automatic";
            case ORACLE_MATRIX: return "matrix";
            case ORACLE_BINARY: return "binary";
            default: return "automatic";
        }
    }

    inline EdgeRatingFunction string_to_edge_rating_function(const std::string &str) {
        if (str == "weight") return EdgeRatingFunction::WEIGHT;
        if (str == "expansion") return EdgeRatingFunction::EXPANSION;
        if (str == "expansionstar") return EdgeRatingFunction::EXPANSIONSTAR;
        if (str == "expansionstarstar") return EdgeRatingFunction::EXPANSIONSTARSTAR;
        if (str == "innerouter") return EdgeRatingFunction::INNEROUTER;
        std::cout << "Error: Invalid rating function '" << str << "'. Allowed values: weight, expansion, expansionstar, expansionstarstar, innerouter" << std::endl;
        exit(EXIT_FAILURE);
    }
}

#include "../utility/utils.h"
#include "../init_partitioning/global_multisection.h"
#include "../init_partitioning/flow_bisection.h"
#include "../init_partitioning/ggg_multisection.h"
#include "../refinement/label_propagation_refinement.h"
#include "../refinement/quotient_graph_refinement.h"
#include "../refinement/flow_based_refinement.h"
#include "../refinement/simple_label_propagation_refinement.h"
#include "../refinement/simple_quotient_graph_refinement.h"
#include "../refinement/simple_flow_based_refinement.h"
#include "../coarsening/global_path_algorithm.h"
#include "../coarsening/size_constrained_lp.h"
#include "../coarsening/heavy_edge_matching.h"
#include "../init_partitioning/kaffpa_partitioner.h"


namespace HeiProMap {
    inline GrowthStrategy string_to_growth_strategy(const std::string &str) {
        if (str == "bfs" || str == "BFS") return GrowthStrategy::BFS;
        if (str == "heavy-first" || str == "heavy_first" || str == "HEAVY_FIRST") return GrowthStrategy::HEAVY_FIRST;
        if (str == "light-first" || str == "light_first" || str == "LIGHT_FIRST") return GrowthStrategy::LIGHT_FIRST;
        std::cout << "Error: Invalid growth strategy '" << str << "'. Allowed values: bfs, heavy-first, light-first" << std::endl;
        exit(EXIT_FAILURE);
    }

    class AlgorithmConfiguration {
    public:
        std::vector<CommandLineOption> options = {
            // General Options
            {"--help", "", "Produces the help message", "", "", false},
            {"--graph", "-g", "Filepath to the graph.", "", "", false},
            {"--mapping", "-m", "Output filepath to the generated mapping.", "", "", false},
            {"--hierarchy", "-h", "Hierarchy in the form a1:a2:...:al .", "", "", false},
            {"--distance", "-d", "Distance in the form d1:d2:...:dl .", "", "", false},
            {"--oracle", "-o", "Distance oracle type (matrix, binary, automatic).", "automatic", "", false},
            {"--imbalance", "-e", "Allowed imbalance (for example 0.03).", "", "", false},
            {"--config", "-c", "The configuration.", "", "", false},
            {"--threads", "-t", "Number of threads.", "1", "", false},
            {"--seed", "", "Seed for diversifying results.", "", "", false},

            // optional
            {"--hm-level", "", "Level of hierarchical multisection.", "0", "", false},
            {"--initial-c", "", "Initial contracting limit.", "64", "", false},
            {"--dump-multisection", "", "After every multisection item (one kaffpa call) write its graph, partition and stats to this directory.", "", "", false},
            {"--save-coarsest-graph", "", "Write the coarsest graph (METIS format, with weights) to this file before initial partitioning.", "", "", false},
            {"--partitioning", "", "Initial partitioning algorithm (multisection, greedy, greedy-graph-growing, flow-bisection). Default: flow-bisection in the fast, eco and strong profiles.", "", "", false},
            {"--flow-bisection-seed-pool", "", "Flow bisection: every seed is drawn at random from the N farthest vertices (1 = the farthest, deterministic).", "1", "", false},
            {"--flow-bisection-kway-rounds", "", "Flow bisection: passes of pairwise flow refinement between adjacent blocks of every level, 0 = off.", "3", "", false},
            {"--flow-bisection-low-kway-rounds", "", "Flow bisection: kway rounds for the levels below the top, -1 = same as the top level.", "-1", "", false},
            {"--flow-bisection-low-small-seeds", "", "Flow bisection: small seeds for the levels below the top, -1 = same as the top level.", "2", "", false},
            {"--flow-bisection-low-refine-rounds", "", "Flow bisection: refine rounds for the levels below the top, -1 = same as the top level.", "-1", "", false},
            {"--flow-bisection-ml-n", "", "Flow bisection: multilevel top level, contracts to this many vertices, 0 = off.", "1000", "", false},
            {"--flow-bisection-qap-rounds", "", "Flow bisection: passes of QAP aware flow refinement between adjacent PEs after the multisection, 0 = off.", "1", "", false},
            {"--flow-bisection-qap-top-levels", "", "Flow bisection: QAP refinement only for PE pairs whose lowest common level is one of this many highest hierarchy levels, 0 = all pairs.", "1", "", false},
            {"--flow-bisection-qap-bands", "", "Flow bisection: number of corridor widths (1-5) of the QAP refinement.", "5", "", false},
            {"--flow-bisection-ml-cap", "", "Flow bisection: cluster weight cap divisor of the multilevel top level.", "20", "", false},
            {"--flow-bisection-small-n", "", "Flow bisection: bisections with fewer vertices use many seeds and no pruning heuristics, 0 = off.", "300", "", false},
            {"--flow-bisection-small-seeds", "", "Flow bisection: number of seeds of the small bisections.", "8", "", false},
            {"--flow-bisection-small-exact", "", "Flow bisection: small bisections continue every job and prune exactly (0/1).", "1", "", false},
            {"--flow-bisection-top-budget", "", "Flow bisection: seed budget of the top hierarchy level, 0 = no limit.", "0", "", false},
            {"--flow-bisection-pair-fraction", "", "Flow bisection: fraction of the seed pairs (not containing the first seed) that are tried.", "1.0", "", false},
            {"--flow-bisection-pair-proxy", "", "Flow bisection: ranking of the seed pairs for pair-fraction: 0 = hop distance, 1 = sum of the seed ball cuts.", "0", "", false},
            {"--flow-bisection-repair", "", "Flow bisection: move vertices out of overloaded blocks after every k'-partition (0/1).", "1", "", false},
            {"--flow-bisection-ibfs", "", "Flow bisection: use IBFS instead of BK for the first max flow of every seed pair (0/1).", "0", "", false},
            {"--flow-bisection-refine-rounds", "", "Flow bisection: max number of improving flow refinements (corridor around the cut) of every bisection, 0 = off.", "100", "", false},
            {"--flow-bisection-refine-band", "", "Flow bisection: weight fraction of the corridor (the refinement tries 0.2, 0.5, 1, 2, 4 times this).", "0.1", "", false},
            {"--flow-bisection-adaptive", "", "Flow bisection: every hierarchy level gets its own (adaptive) imbalance like the multisection, instead of the global one (true/false).", "true", "", false},
            {"--flow-bisection-prune", "", "Flow bisection: stop a seed pair once its cut cannot beat the best feasible cut (true/false).", "true", "", false},
            {"--flow-bisection-init-fraction", "", "Flow bisection: both terminals are first grown to this fraction of the sub graph weight (0 = single vertex).", "0.1", "", false},
            {"--flow-bisection-prune-factor", "", "Flow bisection: drop a pair when its flow reaches this factor times the best feasible cut (1 = exact).", "0.95", "", false},
            {"--flow-bisection-keep-fraction", "", "Flow bisection: fraction of the (pair, orientation) jobs, best lower bound first, that is continued after the first solve.", "0.5", "", false},
            {"--flow-bisection-bound-jobs", "", "Flow bisection: number of jobs whose first solve is complete; the best one is run first and bounds the first solves of the others (0 = off).", "4", "", false},
            {"--flow-bisection-seeds", "", "Number of seed vertices (all pairs are tried) for flow-bisection.", "16", "", false},
            {"--flow-bisection-seed-budget", "", "Flow bisection: a bisection with A arcs uses at most 1 + budget / A seeds (large sub problems get fewer seeds). 0 = no limit.", "300000", "", false},
            {"--flow-bisection-seeds-mid", "", "Flow bisection: seeds on the hierarchy levels between the top and the bottom one (the cuts there weigh less in the QAP).", "2", "", false},
            {"--flow-bisection-seeds-bottom", "", "Flow bisection: seeds on the bottom hierarchy level (cheapest distance).", "2", "", false},
            {"--use-parallel-contraction", "", "Use parallel contraction (true/false).", "", "", false},
            {"--simple-accept-zero-gain", "", "With the simple refiners: the simple flow refinement applies zero-gain cuts like the regular one (true/false).", "true", "", false},
            {"--force-parallel-suite", "", "Always use the parallel suite (simple refiners without quotient graph / block connectivity, parallel LP, parallel contraction), even with 1 thread (true/false). It is enabled automatically when threads > 1.", "false", "", false},
            {"--simple-refinement", "", "Explicitly enable/disable only the simple refiners (true/false). Default: follows the parallel suite.", "", "", false},

            // Coarsening options
            {"--coarsening", "", "Coarsening algorithm (global-paths, size-constrained-lp, heavy-edge).", "", "", false},

            // Global Path Algorithm (GPA) Coarsening
            {"--gpa-coarsening-rating-function", "", "GPA rating function (weight, expansion, expansionstar, expansionstarstar, innerouter).", "", "", false},
            {"--gpa-coarsening-random-level", "", "GPA random level.", "", "", false},
            {"--gpa-coarsening-edge-rating-tiebreaking", "", "GPA use edge rating tiebreaking (true/false).", "", "", false},
            {"--gpa-coarsening-two-hop-threshold", "", "GPA two-hop threshold.", "", "", false},

            // Size-Constrained Label Propagation (SCLP) Coarsening
            {"--sclp-coarsening-max-rounds", "", "SCLP max rounds.", "", "", false},
            {"--sclp-coarsening-min-threshold", "", "SCLP min threshold.", "", "", false},
            {"--sclp-coarsening-f", "", "SCLP tuning parameter f (max cluster weight lmax / f).", "18.0", "", false},
            {"--sclp-coarsening-weight-limit", "", "SCLP cluster weight limit: f (lmax / f) or epsilon-block-weight (KaMinPar: eps * c(V) / clamp(n / C, 2, k)).", "f", "", false},
            {"--sclp-coarsening-rating-function", "", "SCLP rating function (weight, expansion, expansionstar, expansionstarstar, innerouter).", "", "", false},
            {"--sclp-coarsening-use-degree-ordering", "", "SCLP use degree ordering (true/false).", "true", "", false},
            {"--sclp-coarsening-use-parallel-version", "", "SCLP use parallel version (true/false).", "", "", false},

            // Heavy Edge Matching (HEM) Coarsening
            {"--hem-coarsening-rating-function", "", "HEM rating function (weight, expansion, expansionstar, expansionstarstar, innerouter).", "", "", false},

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
            {"--measure-qg-edge-cut", "", "Measure QG edge cut optimization time (true/false).", "false", "", false},

            // Flow-Based Refinement (FB)
            {"--fb-refinement-enabled", "", "Enable flow based refinement (true/false).", "", "", false},
            {"--fb-refinement-use-active-block-scheduling", "", "Use active block scheduling in flow based refinement (true/false).", "", "", false},
            {"--fb-refinement-max-global-iteration", "", "Max global iteration for flow based refinement.", "", "", false},
            {"--fb-refinement-max-local-iteration", "", "Max local iteration for flow based refinement.", "", "", false},
            {"--fb-refinement-alpha", "", "Alpha parameter for flow based refinement.", "", "", false},
            {"--fb-refinement-alpha-upper-bound", "", "Alpha upper bound for flow based refinement.", "", "", false},
            {"--fb-refinement-alpha-modifier", "", "Alpha modifier for flow based refinement.", "", "", false},
            {"--fb-refinement-use-closed-vertex-set", "", "Use closed vertex set in flow based refinement (true/false).", "", "", false},
            {"--fb-refinement-closed-vertex-set-repeats", "", "Number of closed vertex set repeats for flow based refinement.", "", "", false},
            {"--fb-refinement-use-edge-cut", "", "Flow refinement use edge-cut on lowest level (true/false).", "true", "", false},
            {"--fb-refinement-always-include-boundary", "", "Flow refinement always include boundary (true/false).", "", "", false},
            {"--fb-refinement-growth-strategy", "", "Flow refinement growth strategy (bfs, heavy-first, light-first).", "", "", false},
            {"--measure-flow-edge-cut", "", "Measure time for flow edge cut optimization (true/false).", "false", "", false},
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

        u64 initial_c = 64;
        std::string coarsest_graph_out;
        bool use_parallel_contraction = false;

        // if true: no quotient graph / block connectivity datastructures, simple refiners are used
        // parallel suite: automatically enabled when threads > 1 (or forced with --force-parallel-suite)
        bool parallel_suite = false;
        bool simple_refinement = false;
        bool simple_accept_zero_gain = true;

        std::string config_name = "undefined";

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
        FlowBisectionConfiguration flow_bisection_config;
        std::string heipa_config_string = "strong";

        // rebalance algorithm
        std::string rebalancing_algorithm_string;
        REBALANCING_ALGS rebalancing_algorithm_id = REBALANCING_ALG_UNDEFINED;

        // refinement algorithms
        LabelPropagationConfiguration label_propagation_config = LabelPropagationConfiguration("Label Propagation");
        QuotientGraphRefinementConfiguration quotient_graph_refinement_config = QuotientGraphRefinementConfiguration("Quotient Graph");
        FlowBasedRefinementConfiguration flow_based_refinement_config = FlowBasedRefinementConfiguration("Flow Based");
        // distance oracle
        std::string oracle_string = "automatic";
        ORACLE_TYPE oracle_type = ORACLE_AUTOMATIC;

        bool use_binary_oracle() const {
            if (oracle_type == ORACLE_BINARY) return true;
            if (oracle_type == ORACLE_MATRIX) return false;
            return k > 1024;
        }

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

        void set_oracle() {
            if (is_set("--oracle")) {
                oracle_string = get("--oracle");
            } else {
                oracle_string = "automatic";
            }
            oracle_type = string_to_oracle_type(oracle_string);
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

            set_hierarchy();
            set_distance();
            set_oracle();
            set_imbalance();
            set_seed();
            set_threads();
            set_hm_level();

            parallel_suite = (threads > 1);
            if (is_set("--force-parallel-suite")) {
                const std::string val = get("--force-parallel-suite");
                parallel_suite = parallel_suite || (val == "true" || val == "1" || val == "yes");
            }
            simple_refinement = parallel_suite;
            if (is_set("--simple-refinement")) {
                const std::string val = get("--simple-refinement");
                simple_refinement = (val == "true" || val == "1" || val == "yes");
            }
            if (is_set("--simple-accept-zero-gain")) {
                const std::string val = get("--simple-accept-zero-gain");
                simple_accept_zero_gain = (val == "true" || val == "1" || val == "yes");
            }

            if (!is_set("--config")) {
                std::cout << "Error: --config (or -c) must be specified." << std::endl;
                exit(EXIT_FAILURE);
            }
            config_name = get("--config");
            if (config_name == "fast") {
                if (threads == 1) {
                    set_fast();
                } else {
                    set_fast_parallel();
                }
            } else if (config_name == "eco") {
                if (threads == 1) {
                    set_eco();
                } else {
                    set_eco_parallel();
                }
            } else if (config_name == "strong") {
                if (threads == 1) {
                    set_strong();
                } else {
                    set_strong_parallel();
                }
            } else if (config_name == "super-strong") {
                if (threads == 1) {
                    set_super_strong();
                } else {
                    set_super_strong_parallel();
                }
            } else {
                std::cout << "Config " << config_name << " not recognized!" << std::endl;
                exit(EXIT_FAILURE);
            }

            if (is_set("--initial-c")) {
                initial_c = std::stoull(get("--initial-c"));
            }
            if (is_set("--dump-multisection")) {
                global_multisection_config.dump_dir = get("--dump-multisection");
            }
            if (is_set("--save-coarsest-graph")) {
                coarsest_graph_out = get("--save-coarsest-graph");
            }
            if (is_set("--partitioning")) {
                partitioning_algorithm_string = get("--partitioning");
                partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            }
            if (is_set("--flow-bisection-seed-pool")) {
                flow_bisection_config.seed_pool = std::stoull(get("--flow-bisection-seed-pool"));
            }
            if (is_set("--flow-bisection-kway-rounds")) {
                flow_bisection_config.kway_rounds = std::stoull(get("--flow-bisection-kway-rounds"));
            }
            if (is_set("--flow-bisection-low-kway-rounds")) {
                flow_bisection_config.low_kway_rounds = std::stoll(get("--flow-bisection-low-kway-rounds"));
            }
                        if (is_set("--flow-bisection-low-small-seeds")) {
                flow_bisection_config.low_small_seeds = std::stoll(get("--flow-bisection-low-small-seeds"));
            }
                        if (is_set("--flow-bisection-low-refine-rounds")) {
                flow_bisection_config.low_refine_rounds = std::stoll(get("--flow-bisection-low-refine-rounds"));
            }
            if (is_set("--flow-bisection-ml-n")) {
                flow_bisection_config.ml_n = std::stoull(get("--flow-bisection-ml-n"));
            }
            if (is_set("--flow-bisection-qap-rounds")) {
                flow_bisection_config.qap_rounds = std::stoull(get("--flow-bisection-qap-rounds"));
            }
            if (is_set("--flow-bisection-qap-top-levels")) {
                flow_bisection_config.qap_top_levels = std::stoull(get("--flow-bisection-qap-top-levels"));
            }
            if (is_set("--flow-bisection-qap-bands")) {
                flow_bisection_config.qap_bands = std::stoull(get("--flow-bisection-qap-bands"));
            }
            if (is_set("--flow-bisection-ml-cap")) {
                flow_bisection_config.ml_cap = std::stod(get("--flow-bisection-ml-cap"));
            }
            if (is_set("--flow-bisection-small-n")) {
                flow_bisection_config.small_n = std::stoull(get("--flow-bisection-small-n"));
            }
            if (is_set("--flow-bisection-small-seeds")) {
                flow_bisection_config.small_seeds = std::stoull(get("--flow-bisection-small-seeds"));
            }
            if (is_set("--flow-bisection-small-exact")) {
                flow_bisection_config.small_exact = std::stoi(get("--flow-bisection-small-exact")) != 0;
            }
            if (is_set("--flow-bisection-top-budget")) {
                flow_bisection_config.top_seed_budget = std::stoull(get("--flow-bisection-top-budget"));
            }
            if (is_set("--flow-bisection-pair-proxy")) {
                flow_bisection_config.pair_proxy = std::stoi(get("--flow-bisection-pair-proxy"));
            }
            if (is_set("--flow-bisection-pair-fraction")) {
                flow_bisection_config.pair_fraction = std::stod(get("--flow-bisection-pair-fraction"));
            }
            if (is_set("--flow-bisection-repair")) {
                flow_bisection_config.repair = std::stoi(get("--flow-bisection-repair")) != 0;
            }
            if (is_set("--flow-bisection-ibfs")) {
                flow_bisection_config.ibfs = std::stoi(get("--flow-bisection-ibfs")) != 0;
            }
            if (is_set("--flow-bisection-refine-rounds")) {
                flow_bisection_config.refine_rounds = std::stoull(get("--flow-bisection-refine-rounds"));
            }
            if (is_set("--flow-bisection-refine-band")) {
                flow_bisection_config.refine_band = std::stod(get("--flow-bisection-refine-band"));
            }
            if (is_set("--flow-bisection-adaptive")) {
                const std::string val = get("--flow-bisection-adaptive");
                flow_bisection_config.adaptive_imbalance = (val == "true" || val == "1");
            }
            if (is_set("--flow-bisection-prune")) {
                const std::string val = get("--flow-bisection-prune");
                flow_bisection_config.prune = (val == "true" || val == "1");
            }
            if (is_set("--flow-bisection-init-fraction")) {
                flow_bisection_config.init_fraction = std::stod(get("--flow-bisection-init-fraction"));
            }
            if (is_set("--flow-bisection-prune-factor")) {
                flow_bisection_config.prune_factor = std::stod(get("--flow-bisection-prune-factor"));
            }
            if (is_set("--flow-bisection-keep-fraction")) {
                flow_bisection_config.keep_fraction = std::stod(get("--flow-bisection-keep-fraction"));
            }
            if (is_set("--flow-bisection-seeds")) {
                flow_bisection_config.num_seeds = std::stoull(get("--flow-bisection-seeds"));
            }
            if (is_set("--flow-bisection-bound-jobs")) {
                flow_bisection_config.bound_jobs = std::stoull(get("--flow-bisection-bound-jobs"));
            }
            if (is_set("--flow-bisection-seed-budget")) {
                flow_bisection_config.seed_budget = std::stoull(get("--flow-bisection-seed-budget"));
            }
            if (is_set("--flow-bisection-seeds-mid")) {
                flow_bisection_config.num_seeds_mid = std::stoull(get("--flow-bisection-seeds-mid"));
            }
            if (is_set("--flow-bisection-seeds-bottom")) {
                flow_bisection_config.num_seeds_bottom = std::stoull(get("--flow-bisection-seeds-bottom"));
            }

            use_parallel_contraction = parallel_suite;
            if (is_set("--use-parallel-contraction")) {
                std::string val = get("--use-parallel-contraction");
                use_parallel_contraction = (val == "true" || val == "1");
            }

            label_propagation_config.force_parallel_alg = parallel_suite;
            if (is_set("--lp-refinement-enabled")) {
                std::string val = get("--lp-refinement-enabled");
                label_propagation_config.enabled = (val == "true" || val == "1");
            }
            if (is_set("--lp-refinement-max-iteration")) {
                label_propagation_config.max_iteration = std::stoull(get("--lp-refinement-max-iteration"));
            }
            if (is_set("--lp-refinement-use-parallel-version")) {
                std::string val = get("--lp-refinement-use-parallel-version");
                label_propagation_config.force_parallel_alg = (val == "true" || val == "1");
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
            if (is_set("--measure-qg-edge-cut")) {
                std::string val = get("--measure-qg-edge-cut");
                quotient_graph_refinement_config.measure_qg_edge_cut = (val == "true" || val == "1");
            }

            // Override coarsening settings from CLI if they are set
            if (is_set("--coarsening")) {
                coarsening_algorithm_string = get("--coarsening");
                coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);
            }

            // Override GPA configs
            if (is_set("--gpa-coarsening-rating-function")) {
                global_path_algorithm_config.rating_function = string_to_edge_rating_function(get("--gpa-coarsening-rating-function"));
            }
            if (is_set("--gpa-coarsening-random-level")) {
                global_path_algorithm_config.random_level = std::stoull(get("--gpa-coarsening-random-level"));
            }
            if (is_set("--gpa-coarsening-edge-rating-tiebreaking")) {
                std::string val = get("--gpa-coarsening-edge-rating-tiebreaking");
                global_path_algorithm_config.use_edge_rating_tiebreaking = (val == "true" || val == "1");
            }
            if (is_set("--gpa-coarsening-two-hop-threshold")) {
                global_path_algorithm_config.two_hop_threshold = std::stod(get("--gpa-coarsening-two-hop-threshold"));
            }

            size_constrained_lp_config.force_parallel_alg = parallel_suite;
            size_constrained_lp_config.contraction_limit = initial_c;

            // Override SCLP configs
            if (is_set("--sclp-coarsening-max-rounds")) {
                size_constrained_lp_config.max_rounds = std::stoull(get("--sclp-coarsening-max-rounds"));
            }
            if (is_set("--sclp-coarsening-min-threshold")) {
                size_constrained_lp_config.min_threshold = std::stod(get("--sclp-coarsening-min-threshold"));
            }
            if (is_set("--sclp-coarsening-f")) {
                size_constrained_lp_config.f = std::stod(get("--sclp-coarsening-f"));
            }
            if (is_set("--sclp-coarsening-weight-limit")) {
                const std::string val = get("--sclp-coarsening-weight-limit");
                if (val == "epsilon-block-weight") { size_constrained_lp_config.use_epsilon_block_weight = true; }
                else if (val == "f") { size_constrained_lp_config.use_epsilon_block_weight = false; }
                else { std::cerr << "Unknown --sclp-coarsening-weight-limit " << val << " (f, epsilon-block-weight)" << std::endl; exit(EXIT_FAILURE); }
            }
            if (is_set("--sclp-coarsening-rating-function")) {
                size_constrained_lp_config.rating_function = string_to_edge_rating_function(get("--sclp-coarsening-rating-function"));
            }
            if (is_set("--sclp-coarsening-use-degree-ordering")) {
                std::string val = get("--sclp-coarsening-use-degree-ordering");
                size_constrained_lp_config.use_degree_ordering = (val == "true" || val == "1");
            }
            if (is_set("--sclp-coarsening-use-parallel-version")) {
                std::string val = get("--sclp-coarsening-use-parallel-version");
                size_constrained_lp_config.force_parallel_alg = (val == "true" || val == "1");
            }

            // Override HEM configs
            if (is_set("--hem-coarsening-rating-function")) {
                heavy_edge_matching_config.rating_function = string_to_edge_rating_function(get("--hem-coarsening-rating-function"));
            }

            // Override flow refinement configs
            if (is_set("--fb-refinement-enabled")) {
                std::string val = get("--fb-refinement-enabled");
                flow_based_refinement_config.enabled = (val == "true" || val == "1");
            }
            if (is_set("--fb-refinement-use-active-block-scheduling")) {
                std::string val = get("--fb-refinement-use-active-block-scheduling");
                flow_based_refinement_config.use_active_block_scheduling = (val == "true" || val == "1");
            }
            if (is_set("--fb-refinement-max-global-iteration")) {
                flow_based_refinement_config.max_global_iteration = std::stoull(get("--fb-refinement-max-global-iteration"));
            }
            if (is_set("--fb-refinement-max-local-iteration")) {
                flow_based_refinement_config.max_local_iteration = std::stoull(get("--fb-refinement-max-local-iteration"));
            }
            if (is_set("--fb-refinement-alpha")) {
                flow_based_refinement_config.alpha = std::stod(get("--fb-refinement-alpha"));
            }
            if (is_set("--fb-refinement-alpha-upper-bound")) {
                flow_based_refinement_config.alpha_upper_bound = std::stod(get("--fb-refinement-alpha-upper-bound"));
            }
            if (is_set("--fb-refinement-alpha-modifier")) {
                flow_based_refinement_config.alpha_modifier = std::stod(get("--fb-refinement-alpha-modifier"));
            }
            if (is_set("--fb-refinement-use-closed-vertex-set")) {
                std::string val = get("--fb-refinement-use-closed-vertex-set");
                flow_based_refinement_config.use_closed_vertex_set = (val == "true" || val == "1");
            }
            if (is_set("--fb-refinement-closed-vertex-set-repeats")) {
                flow_based_refinement_config.closed_vertex_sets_repeats = std::stoull(get("--fb-refinement-closed-vertex-set-repeats"));
            }
            if (is_set("--fb-refinement-use-edge-cut")) {
                std::string val = get("--fb-refinement-use-edge-cut");
                flow_based_refinement_config.use_edge_cut = (val == "true" || val == "1");
            }
            if (is_set("--fb-refinement-always-include-boundary")) {
                std::string val = get("--fb-refinement-always-include-boundary");
                flow_based_refinement_config.always_include_boundary = (val == "true" || val == "1");
            }
            if (is_set("--fb-refinement-growth-strategy")) {
                flow_based_refinement_config.growth_strategy = string_to_growth_strategy(get("--fb-refinement-growth-strategy"));
            }
            if (is_set("--measure-flow-edge-cut")) {
                std::string val = get("--measure-flow-edge-cut");
                flow_based_refinement_config.measure_flow_edge_cut = (val == "true" || val == "1");
            }
        }

        void set_fast() {
            config_name = "fast";
            collect_dataset = false;

            initial_c = 8; // flow-bisection: tiny coarse graphs are enough (c=8 is 45% faster than 32 for ~1% QAP)

            // set GPA matching algorithm
            // coarsening_algorithm_string = "size-constrained-lp";
            coarsening_algorithm_string = "global-paths";
            // coarsening_algorithm_string = "heavy-edge";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;
            global_path_algorithm_config.random_level = 0;
            global_path_algorithm_config.use_edge_rating_tiebreaking = false;

            heavy_edge_matching_config.rating_function = EdgeRatingFunction::EXPANSION;

            size_constrained_lp_config.max_rounds = 5;
            size_constrained_lp_config.min_threshold = 0.05;
            size_constrained_lp_config.f = 18.0;
            size_constrained_lp_config.rating_function = EdgeRatingFunction::WEIGHT;
            size_constrained_lp_config.use_degree_ordering = true;

            // set multisection
            partitioning_algorithm_string = "flow-bisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            flow_bisection_config.num_seeds = 6;
            flow_bisection_config.ml_n = 500;
            flow_bisection_config.qap_rounds = 0;
            flow_bisection_config.kway_rounds = 1;
            flow_bisection_config.low_kway_rounds = 1;
            global_multisection_config.mode_string = "heipa-fast";
            global_multisection_config.mode = GLOBAL_MULTISECTION_HEIPA_FAST;
            // global_multisection_config.mode_string = "kaffpa-fast";
            // global_multisection_config.mode = GLOBAL_MULTISECTION_KAFFPA_FAST;
            global_multisection_config.kappa = 64;
            global_multisection_config.refine = true;
            global_multisection_config.v_cycles = 0;
            global_multisection_config.v_cycle_depth = 100;

            // enable label propagation
            label_propagation_config.enabled = true;
            label_propagation_config.max_iteration = 5;
            label_propagation_config.force_parallel_alg = false;

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
            flow_based_refinement_config.always_include_boundary = true;
            flow_based_refinement_config.closed_vertex_sets_repeats = 500;
        }

        void set_eco() {
            config_name = "eco";
            initial_c = 32; // flow-bisection: small coarse graphs are enough

            // set GPA matching algorithm
            // coarsening_algorithm_string = "size-constrained-lp";
            coarsening_algorithm_string = "global-paths";
            // coarsening_algorithm_string = "heavy-edge";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;
            global_path_algorithm_config.random_level = 0;
            global_path_algorithm_config.use_edge_rating_tiebreaking = false;
            global_path_algorithm_config.two_hop_threshold = 0.75;

            heavy_edge_matching_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;

            size_constrained_lp_config.max_rounds = 5;
            size_constrained_lp_config.min_threshold = 0.05;
            size_constrained_lp_config.f = 18.0;
            size_constrained_lp_config.rating_function = EdgeRatingFunction::WEIGHT;
            size_constrained_lp_config.use_degree_ordering = true;

            // set multisection
            partitioning_algorithm_string = "flow-bisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            flow_bisection_config.seed_budget = 0; // eco can afford all seeds on large and dense graphs
            flow_bisection_config.num_seeds = 8;
            flow_bisection_config.ml_n = 500;
            flow_bisection_config.kway_rounds = 1;
            flow_bisection_config.low_kway_rounds = 3;
            global_multisection_config.mode_string = "heipa-fast";
            global_multisection_config.mode = GLOBAL_MULTISECTION_HEIPA_FAST;
            // global_multisection_config.mode_string = "kaffpa-fast";
            // global_multisection_config.mode = GLOBAL_MULTISECTION_KAFFPA_FAST;
            global_multisection_config.kappa = 64;
            global_multisection_config.refine = true;
            global_multisection_config.v_cycles = 0;
            global_multisection_config.v_cycle_depth = 100;

            // enable label propagation
            label_propagation_config.enabled = true;
            label_propagation_config.max_iteration = 5;
            label_propagation_config.force_parallel_alg = false;

            // enable quotient graph refinement
            quotient_graph_refinement_config.enabled = true;
            quotient_graph_refinement_config.max_iteration = 2;
            quotient_graph_refinement_config.alpha = 5.0;
            quotient_graph_refinement_config.min_n_steps = 3;
            quotient_graph_refinement_config.use_preemptive_exit = true;

            // enable flow based refinement
            flow_based_refinement_config.enabled = true;
            flow_based_refinement_config.use_active_block_scheduling = true;
            flow_based_refinement_config.max_global_iteration = 5;
            flow_based_refinement_config.max_local_iteration = 5;
            flow_based_refinement_config.alpha = 2.0;
            flow_based_refinement_config.alpha_upper_bound = 64.0;
            flow_based_refinement_config.alpha_modifier = 2.0;
            flow_based_refinement_config.use_closed_vertex_set = false;
            flow_based_refinement_config.closed_vertex_sets_repeats = 500;
            flow_based_refinement_config.always_include_boundary = true;
            flow_based_refinement_config.growth_strategy = GrowthStrategy::BFS;
        }

        void set_strong() {
            config_name = "strong";
            initial_c = 64;

            // set GPA matching algorithm
            // coarsening_algorithm_string = "size-constrained-lp";
            coarsening_algorithm_string = "global-paths";
            // coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;
            global_path_algorithm_config.random_level = 0;
            global_path_algorithm_config.use_edge_rating_tiebreaking = false;
            global_path_algorithm_config.two_hop_threshold = 0.75;

            heavy_edge_matching_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;

            size_constrained_lp_config.max_rounds = 5;
            size_constrained_lp_config.min_threshold = 0.05;
            size_constrained_lp_config.f = 18.0;
            size_constrained_lp_config.rating_function = EdgeRatingFunction::WEIGHT;
            size_constrained_lp_config.use_degree_ordering = true;

            // set multisection
            partitioning_algorithm_string = "flow-bisection"; // the multisection (kaffpa-strong below) is still available with --partitioning multisection
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            // global_multisection_config.mode_string = "heipa-fast";
            // global_multisection_config.mode = GLOBAL_MULTISECTION_HEIPA_FAST;
            global_multisection_config.mode_string = "kaffpa-strong";
            global_multisection_config.mode = GLOBAL_MULTISECTION_KAFFPA_STRONG;
            global_multisection_config.kappa = 3;
            global_multisection_config.refine = true;
            global_multisection_config.v_cycles = 0;
            global_multisection_config.v_cycle_depth = 100;

            // enable label propagation
            label_propagation_config.enabled = true;
            label_propagation_config.max_iteration = 5;
            label_propagation_config.force_parallel_alg = false;

            // enable quotient graph refinement
            quotient_graph_refinement_config.enabled = true;
            quotient_graph_refinement_config.max_iteration = 2;
            quotient_graph_refinement_config.alpha = 5.0;
            quotient_graph_refinement_config.min_n_steps = 3;
            quotient_graph_refinement_config.use_preemptive_exit = true;

            // enable flow based refinement
            flow_based_refinement_config.enabled = true;
            flow_based_refinement_config.use_active_block_scheduling = true;
            flow_based_refinement_config.max_global_iteration = 5;
            flow_based_refinement_config.max_local_iteration = 5;
            flow_based_refinement_config.alpha = 2.0;
            flow_based_refinement_config.alpha_upper_bound = 64.0;
            flow_based_refinement_config.alpha_modifier = 2.0;
            flow_based_refinement_config.use_closed_vertex_set = false;
            flow_based_refinement_config.closed_vertex_sets_repeats = 500;
            flow_based_refinement_config.always_include_boundary = true;
            flow_based_refinement_config.growth_strategy = GrowthStrategy::BFS;
        }

        void set_super_strong() {
            config_name = "super-strong";
            initial_c = 64;

            // set GPA matching algorithm
            coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.random_level = 0;

            size_constrained_lp_config.max_rounds = 5;
            size_constrained_lp_config.min_threshold = 0.10;
            size_constrained_lp_config.f = 8.0;

            // set multisection
            partitioning_algorithm_string = "multisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);
            global_multisection_config.mode_string = "heipa-fast";
            global_multisection_config.mode = GLOBAL_MULTISECTION_HEIPA_FAST;
            global_multisection_config.kappa = 1;
            global_multisection_config.refine = true;
            global_multisection_config.v_cycles = 2;
            global_multisection_config.v_cycle_depth = 100;

            // enable label propagation
            label_propagation_config.enabled = true;
            label_propagation_config.max_iteration = 5;
            label_propagation_config.force_parallel_alg = false;

            // enable quotient graph refinement
            quotient_graph_refinement_config.enabled = true;
            quotient_graph_refinement_config.max_iteration = 2;
            quotient_graph_refinement_config.alpha = 5.0;
            quotient_graph_refinement_config.min_n_steps = 3;
            quotient_graph_refinement_config.use_preemptive_exit = true;

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

        void set_fast_parallel() {
            set_fast();
        }

        void set_eco_parallel() {
            set_eco();
        }

        void set_strong_parallel() {
            set_strong();
        }

        void set_super_strong_parallel() {
            set_super_strong();
        }

        /**
         * Simple refiner configurations, derived from the regular refinement
         * configurations. Only meaningful if simple_refinement is set.
         */
        SimpleLabelPropagationConfiguration get_simple_label_propagation_config() const {
            SimpleLabelPropagationConfiguration c("Simple Label Propagation");
            c.enabled = label_propagation_config.enabled;
            c.max_iteration = label_propagation_config.max_iteration;
            c.force_parallel_alg = label_propagation_config.force_parallel_alg;
            c.use_edge_cut = label_propagation_config.use_edge_cut;
            // parallel suite does not use q_graph / block_conn, skip rebuilding them
            c.enable_q_graph = false;
            c.enable_block_conn = false;
            return c;
        }

        SimpleQuotientGraphRefinementConfiguration get_simple_qg_refinement_config() const {
            SimpleQuotientGraphRefinementConfiguration c("Simple Quotient Graph");
            c.enabled = quotient_graph_refinement_config.enabled;
            c.max_iteration = quotient_graph_refinement_config.max_iteration;
            c.min_n_steps = quotient_graph_refinement_config.min_n_steps;
            c.alpha = quotient_graph_refinement_config.alpha;
            c.use_preemptive_exit = quotient_graph_refinement_config.use_preemptive_exit;
            c.use_edge_cut = quotient_graph_refinement_config.use_edge_cut;
            return c;
        }

        SimpleFlowBasedRefinementConfiguration get_simple_flow_based_refinement_config() const {
            SimpleFlowBasedRefinementConfiguration c("Simple Flow Based");
            c.enabled = flow_based_refinement_config.enabled;
            c.max_global_iteration = flow_based_refinement_config.max_global_iteration;
            c.max_local_iteration = flow_based_refinement_config.max_local_iteration;
            c.alpha = flow_based_refinement_config.alpha;
            c.alpha_upper_bound = flow_based_refinement_config.alpha_upper_bound;
            c.alpha_modifier = flow_based_refinement_config.alpha_modifier;
            c.use_closed_vertex_set = flow_based_refinement_config.use_closed_vertex_set;
            c.closed_vertex_sets_repeats = flow_based_refinement_config.closed_vertex_sets_repeats;
            c.use_edge_cut = flow_based_refinement_config.use_edge_cut;
            c.always_include_boundary = flow_based_refinement_config.always_include_boundary;
            c.accept_zero_gain = simple_accept_zero_gain;
            return c;
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
