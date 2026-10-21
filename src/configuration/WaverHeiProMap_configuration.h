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

#ifndef HEIPROMAP_WAVERHEIPROMAP_CONFIGURATION_H
#define HEIPROMAP_WAVERHEIPROMAP_CONFIGURATION_H

#include <string>
#include <vector>

#include "../definitions.h"
#include "../refinement/simple_label_propagation_refinement.h"
#include "../refinement/waver_label_propagation_refinement.h"
#include "../refinement/waver_message_refinement.h"
namespace HeiProMap {

    class NegativeCycleConfiguration final {
    public:
        NegativeCycleConfiguration() = default;
        explicit NegativeCycleConfiguration(const std::string &t_name) { name = t_name; }
        std::string name = "Negative Cycle Refinement";
        bool enabled = false;
        u64 random_tries = 10;
        u64 max_iterations = 10;
        u64 max_path_length = 8;
        u64 threshold = 10;
    };

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

    inline EdgeRatingFunction string_to_edge_rating_function(const std::string &str) {
        if (str == "weight") return EdgeRatingFunction::WEIGHT;
        if (str == "expansion") return EdgeRatingFunction::EXPANSION;
        if (str == "expansion*") return EdgeRatingFunction::EXPANSIONSTAR;
        if (str == "expansion**") return EdgeRatingFunction::EXPANSIONSTARSTAR;
        if (str == "innerouter") return EdgeRatingFunction::INNEROUTER;
        std::cout << "Error: Invalid rating function '" << str << "'. Allowed values: weight, expansion, expansion*, expansion**, innerouter" << std::endl;
        exit(EXIT_FAILURE);
    }
}

#include "../utility/utils.h"
// #include "../init_partitioning/global_multisection.h"
#include "../coarsening/global_path_algorithm.h"
#include "../coarsening/size_constrained_lp.h"
#include "../coarsening/heavy_edge_matching.h"
#include "../init_partitioning/kaffpa_partitioner.h"


namespace HeiProMap {
    class WaverHeiProMapConfiguration {
    public:
        std::vector<CommandLineOption> options = {
            {"--help", "", "Produces the help message", "", "", false},
            {"--graph", "-g", "Filepath to the graph.", "", "", false},
            {"--mapping", "-m", "Output filepath to the generated mapping.", "", "", false},
            {"--topology", "-tp", "Filepath to the topology graph for the distance oracle.", "", "", false},
            {"--imbalance", "-e", "Allowed imbalance of the partition of the input graph into PEs, used by KaMinPar (for example 0.03). The mapping itself is always a bijection.", "0.03", "", false},
            {"--swap-radius", "", "Neighborhood radius for swap refinement (1 = direct neighbors, 0 = disabled).", "0", "", false},
            {"--swap-extended", "", "Also consider neighborhoods of communication partners for swaps.", "0", "", false},
            {"--long-msg-padding", "", "Rectangle padding for long-message swap refinement (0 = disabled).", "1", "", false},
            {"--sa-iterations", "", "Number of simulated annealing swap iterations (0 = disabled).", "1000", "", false},
            {"--sa-radius", "", "Max swap radius for simulated annealing.", "10", "", false},
            {"--kaminpar-threads", "", "Threads for the KaMinPar partition (0 = all hardware threads). Independent of --threads.", "0", "", false},
            {"--kaminpar-preset", "", "KaMinPar preset used for that partition.", "largek", "", false},
            {"--prepartitioner", "", "Partitioner of the input graph into k blocks if it has more vertices than PEs: kaminpar or spectral3d (multilevel three-eigenvector spectral embedding cut by recursive coordinate bisection).", "kaminpar", "", false},
            {"--mode", "", "How the mapping is found: spectral (return the spectral Laplacian-RCB placement, then refine only the finished mapping), multilevel (coarsen, partition, uncoarsen with refinement), or both (multilevel, adopting the spectral placement at a level whenever it is cheaper).", "spectral", "", false},
            {"--vcycle", "", "After the bijection is found, run a V-cycle: coarsen topology and mapped vertices together, then uncoarsen with swap refinement at every level (0 = off).", "0", "", false},
            {"--vcycle-swap-radius", "", "Neighborhood radius of the swap refinement used at every V-cycle level.", "2", "", false},
            {"--vertex-refine", "", "After the mapping is propagated back to the input vertices, refine the individual vertex placement with message-group refinement within the allowed --imbalance (0 = off).", "1", "", false},
            {"--vertex-refine-iters", "", "Maximum iterations of the input-vertex refinement.", "10", "", false},
            {"--topo-split-factor", "", "Split a coarse topology cell only if its block weight is at least this factor times its number of children.", "4", "", false},
            {"--vmove-radius", "", "Neighborhood radius for vertex move refinement.", "5", "", false},
            {"--vmove-extended", "", "Also consider neighborhoods of communication partners for vertex moves.", "1", "", false},
            {"--config", "-c", "The configuration.", "", "", false},
            {"--threads", "-t", "Number of threads.", "1", "", false},
            {"--seed", "", "Seed for diversifying results.", "", "", false},

            {"--alpha-send", "", "Fixed send-side overhead per message (alpha_send).", "100.0", "", false},
            {"--alpha-recv", "", "Fixed receive-side overhead per message (alpha_recv).", "100.0", "", false},
            {"--alpha-link", "", "Fixed per-hop router overhead per message (alpha_link).", "2.0", "", false},
            {"--beta-send", "", "Volume-dependent send overhead per byte (beta_send).", "0.0", "", false},
            {"--beta-recv", "", "Volume-dependent receive overhead per byte (beta_recv).", "8.0", "", false},
            {"--beta-link", "", "Volume-dependent per-hop link cost per byte (beta_link).", "1.0", "", false},

            {"--coarsening", "", "Coarsening algorithm (global-paths, size-constrained-lp, heavy-edge).", "", "", false},
            {"--gpa-rating-function", "", "GPA rating function (weight, expansion, expansion*, expansion**, innerouter).", "", "", false},
            {"--gpa-random-level", "", "GPA random level.", "", "", false},
            {"--gpa-edge-rating-tiebreaking", "", "GPA use edge rating tiebreaking (true/false).", "", "", false},
            {"--gpa-two-hop-threshold", "", "GPA two-hop threshold.", "", "", false},
            {"--sclp-max-rounds", "", "SCLP max rounds.", "", "", false},
            {"--sclp-min-threshold", "", "SCLP min threshold.", "", "", false},
            {"--sclp-f", "", "SCLP tuning parameter f.", "8.0", "", false},
            {"--sclp-rating-function", "", "SCLP rating function (weight, expansion, expansion*, expansion**, innerouter).", "", "", false},
            {"--sclp-use-degree-ordering", "", "SCLP use degree ordering (true/false).", "true", "", false},
            {"--hem-rating-function", "", "HEM rating function (weight, expansion, expansion*, expansion**, innerouter).", "", "", false},

            {"--flow-refinement-enabled", "", "Enable flow based refinement (true/false).", "", "", false},
            {"--flow-refinement-use-active-block-scheduling", "", "Use active block scheduling in flow based refinement (true/false).", "", "", false},
            {"--flow-refinement-max-global-iteration", "", "Max global iteration for flow based refinement.", "", "", false},
            {"--flow-refinement-max-local-iteration", "", "Max local iteration for flow based refinement.", "", "", false},
            {"--flow-refinement-alpha", "", "Alpha parameter for flow based refinement.", "", "", false},
            {"--flow-refinement-alpha-upper-bound", "", "Alpha upper bound for flow based refinement.", "", "", false},
            {"--flow-refinement-alpha-modifier", "", "Alpha modifier for flow based refinement.", "", "", false},
            {"--flow-refinement-use-closed-vertex-set", "", "Use closed vertex set in flow based refinement (true/false).", "", "", false},
            {"--flow-refinement-closed-vertex-set-repeats", "", "Number of closed vertex set repeats for flow based refinement.", "", "", false},
            {"--flow-refinement-enabled-edge-cut-optimization", "", "Enable lowest level edge cut optimization (true/false).", "true", "", false},

            {"--initial-kappa", "", "Initial partitioning: number of attempts per thread, best is kept.", "10", "", false},
        };

    public:
        // graph information
        std::string graph_in;
        std::string topology_in;
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
        u64 swap_radius = 0;
        bool swap_extended = false;
        u64 vmove_radius = 2;
        u64 long_msg_padding = 3;
        u64 sa_iterations = 0;
        u64 sa_radius = 10;
        std::string mode = "spectral"; // spectral | multilevel | both
        bool do_spectral = true;     // derived from mode: solve the spectral placement (spectral, both)
        bool do_multilevel = false;  // derived from mode: partition and refine multilevel (multilevel, both)
        bool vcycle = true;
        bool vertex_refine = true;
        u64 vertex_refine_iters = 10;
        u64 vcycle_swap_radius = 2;
        std::string kaminpar_preset = "default";
        u64 kaminpar_threads = 0;
        std::string prepartitioner = "kaminpar";
        u64 topo_split_factor = 4;
        bool vmove_extended = false;
        f64 per_level_imb_add = 0.005;

        // random initialization
        u64 seed = 0;

        u64 threads = 1;

        bool collect_dataset = false;
        std::string data_dir = std::string(HEIPROMAP_SOURCE_DIR) + "/data";

        u64 hm_level = 0;

        u64 initial_c = 8;

        // communication cost model: Cost(m) = (alpha_send + alpha_recv) + (beta_send + beta_recv)*S(m)
        //                                    + h(m) * alpha_link + h(m) * beta_link * S(m)
        f64 alpha_send = 100.0;
        f64 alpha_recv = 100.0;
        f64 alpha_link =   2.0;
        f64 beta_send  =   0.0;
        f64 beta_recv  =   8.0;
        f64 beta_link  =   1.0;

        // derived convenience accessors
        f64 alpha_end() const { return alpha_send + alpha_recv; }
        f64 beta_end()  const { return beta_send  + beta_recv;  }

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

        u64 initial_kappa = 100;
        u64 rb_kappa = 1;
        bool rb_use_full_refine = false;

        std::string heipa_config_string = "strong";

        // rebalance algorithm
        std::string rebalancing_algorithm_string;
        REBALANCING_ALGS rebalancing_algorithm_id = REBALANCING_ALG_UNDEFINED;

        // refinement algorithms
        SimpleLabelPropagationConfiguration simple_label_propagation_config = SimpleLabelPropagationConfiguration("Label Propagation");
        WaverLabelPropagationConfiguration   waver_label_propagation_config   = WaverLabelPropagationConfiguration("Waver Label Propagation");
        WaverMessageRefinementConfiguration  waver_message_refinement_config  = WaverMessageRefinementConfiguration("Waver Message Refinement");
        NegativeCycleConfiguration negative_cycle_config = NegativeCycleConfiguration("Negative Cycle");

        WaverHeiProMapConfiguration() = default;

        void set_hierarchy() {
            hierarchy = convert<partition_t>(split(hierarchy_string, ':'));
            k = prod<partition_t>(hierarchy);
        }

        void set_distance() {
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


        WaverHeiProMapConfiguration(int argc, char *argv[]) {
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
            if (is_set("--topology")) { topology_in = get("--topology"); }
            if (is_set("--mapping")) {
                mapping_out = get("--mapping");
            }

            set_hierarchy();
            set_distance();
            set_imbalance();
            swap_radius = std::stoull(get("--swap-radius"));
            swap_extended = std::stoull(get("--swap-extended")) != 0;
            long_msg_padding = std::stoull(get("--long-msg-padding"));
            sa_iterations = std::stoull(get("--sa-iterations"));
            sa_radius = std::stoull(get("--sa-radius"));
            mode = get("--mode");
            if (mode != "spectral" && mode != "multilevel" && mode != "both") {
                std::cout << "Error: Invalid --mode '" << mode
                          << "'. Allowed values: spectral, multilevel, both" << std::endl;
                exit(EXIT_FAILURE);
            }
            do_spectral = mode != "multilevel";
            do_multilevel = mode != "spectral";
            vcycle = std::stoull(get("--vcycle")) != 0;
            vertex_refine = std::stoull(get("--vertex-refine")) != 0;
            vertex_refine_iters = std::stoull(get("--vertex-refine-iters"));
            vcycle_swap_radius = std::stoull(get("--vcycle-swap-radius"));
            kaminpar_preset = get("--kaminpar-preset");
            kaminpar_threads = std::stoull(get("--kaminpar-threads"));
            prepartitioner = get("--prepartitioner");
            if (prepartitioner != "kaminpar" && prepartitioner != "spectral3d") {
                std::cout << "Error: Invalid --prepartitioner '" << prepartitioner
                          << "'. Allowed values: kaminpar, spectral3d" << std::endl;
                exit(EXIT_FAILURE);
            }
            topo_split_factor = std::stoull(get("--topo-split-factor"));
            vmove_radius = std::stoull(get("--vmove-radius"));
            vmove_extended = std::stoull(get("--vmove-extended")) != 0;
            set_seed();
            set_threads();

            if (is_set("--initial-kappa")) initial_kappa = std::stoull(get("--initial-kappa"));

            if (is_set("--alpha-send")) alpha_send = std::stod(get("--alpha-send"));
            if (is_set("--alpha-recv")) alpha_recv = std::stod(get("--alpha-recv"));
            if (is_set("--alpha-link")) alpha_link = std::stod(get("--alpha-link"));
            if (is_set("--beta-send"))  beta_send  = std::stod(get("--beta-send"));
            if (is_set("--beta-recv"))  beta_recv  = std::stod(get("--beta-recv"));
            if (is_set("--beta-link"))  beta_link  = std::stod(get("--beta-link"));

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

            // Propagate cost-model parameters into the refinement configs
            waver_label_propagation_config.alpha_end  = alpha_end();
            waver_label_propagation_config.beta_end   = beta_end();
            waver_label_propagation_config.alpha_link = alpha_link;
            waver_label_propagation_config.beta_link  = beta_link;

            // Override coarsening settings from CLI if they are set
            if (is_set("--coarsening")) {
                coarsening_algorithm_string = get("--coarsening");
                coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);
            }

            // Override GPA configs
            if (is_set("--gpa-rating-function")) {
                global_path_algorithm_config.rating_function = string_to_edge_rating_function(get("--gpa-rating-function"));
            }
            if (is_set("--gpa-random-level")) {
                global_path_algorithm_config.random_level = std::stoull(get("--gpa-random-level"));
            }
            if (is_set("--gpa-edge-rating-tiebreaking")) {
                std::string val = get("--gpa-edge-rating-tiebreaking");
                global_path_algorithm_config.use_edge_rating_tiebreaking = (val == "true" || val == "1");
            }
            if (is_set("--gpa-two-hop-threshold")) {
                global_path_algorithm_config.two_hop_threshold = std::stod(get("--gpa-two-hop-threshold"));
            }

            // Override SCLP configs
            if (is_set("--sclp-max-rounds")) {
                size_constrained_lp_config.max_rounds = std::stoull(get("--sclp-max-rounds"));
            }
            if (is_set("--sclp-min-threshold")) {
                size_constrained_lp_config.min_threshold = std::stod(get("--sclp-min-threshold"));
            }
            if (is_set("--sclp-f")) {
                size_constrained_lp_config.f = std::stod(get("--sclp-f"));
            }
            if (is_set("--sclp-rating-function")) {
                size_constrained_lp_config.rating_function = string_to_edge_rating_function(get("--sclp-rating-function"));
            }
            if (is_set("--sclp-use-degree-ordering")) {
                std::string val = get("--sclp-use-degree-ordering");
                size_constrained_lp_config.use_degree_ordering = (val == "true" || val == "1");
            }

            // Override HEM configs
            if (is_set("--hem-rating-function")) {
                heavy_edge_matching_config.rating_function = string_to_edge_rating_function(get("--hem-rating-function"));
            }

        }

        void set_fast() {
            config_name = "fast";
            collect_dataset = false;

            initial_c = 64;

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
            size_constrained_lp_config.f = 32.0;
            size_constrained_lp_config.rating_function = EdgeRatingFunction::WEIGHT;
            size_constrained_lp_config.use_degree_ordering = true;

            // set multisection
            partitioning_algorithm_string = "multisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);


            simple_label_propagation_config.enabled = false;

            waver_label_propagation_config.enabled = true;
            waver_label_propagation_config.max_iteration = 25;
            waver_label_propagation_config.force_parallel_alg = false;

            waver_message_refinement_config.enabled   = false;
            waver_message_refinement_config.max_iters = 5;

            negative_cycle_config.enabled = false;
        }

        void set_eco() {
            config_name = "eco";
            initial_c = 64;

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
            size_constrained_lp_config.f = 16.0;

            // set multisection
            partitioning_algorithm_string = "multisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);


            simple_label_propagation_config.enabled = false;

            waver_label_propagation_config.enabled = true;
            waver_label_propagation_config.max_iteration = 5;

            negative_cycle_config.enabled = false;
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
            size_constrained_lp_config.f = 16.0;

            // set multisection
            partitioning_algorithm_string = "multisection";
            partitioning_algorithm_id = string_to_partitioning_algorithm(partitioning_algorithm_string);


            simple_label_propagation_config.enabled = false;

            waver_label_propagation_config.enabled = true;
            waver_label_propagation_config.max_iteration = 5;

            negative_cycle_config.enabled = false;
        }

        void set_super_strong() {
            config_name = "super-strong";
            initial_c = 16;

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


            simple_label_propagation_config.enabled = false;

            waver_label_propagation_config.enabled = true;
            waver_label_propagation_config.max_iteration = 5;

            negative_cycle_config.enabled = false;
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

#endif //HEIPROMAP_WAVERHEIPROMAP_CONFIGURATION_H
