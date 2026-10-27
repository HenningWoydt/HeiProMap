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

#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "../definitions.h"
#include "../coarsening/global_path_algorithm.h"
#include "../coarsening/heavy_edge_matching.h"
#include "../coarsening/size_constrained_lp.h"
#include "../refinement/waver_label_propagation_refinement.h"
#include "../refinement/waver_long_message_refinement.h"
#include "../refinement/waver_message_refinement.h"
#include "../refinement/waver_sa_refinement.h"
#include "../refinement/waver_swap_refinement.h"
#include "../utility/spectral_tracker.h"
#include "../utility/utils.h"

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

    inline EdgeRatingFunction string_to_edge_rating_function(const std::string &str) {
        if (str == "weight") return EdgeRatingFunction::WEIGHT;
        if (str == "expansion") return EdgeRatingFunction::EXPANSION;
        if (str == "expansion*") return EdgeRatingFunction::EXPANSIONSTAR;
        if (str == "expansion**") return EdgeRatingFunction::EXPANSIONSTARSTAR;
        if (str == "innerouter") return EdgeRatingFunction::INNEROUTER;
        std::cout << "Error: Invalid rating function '" << str << "'. Allowed values: weight, expansion, expansion*, expansion**, innerouter" << std::endl;
        exit(EXIT_FAILURE);
    }

    // All settings of Waver-HeiProMap. --config selects a preset (fast, eco, strong, super-strong) that sets the effort
    // related settings; the options of the same names given on the command line override the preset.
    class WaverHeiProMapConfiguration {
    public:
        std::vector<CommandLineOption> options = {
            {"--help", "", "Produces the help message", "", "", false},
            {"--graph", "-g", "Filepath to the graph.", "", "", false},
            {"--mapping", "-m", "Output filepath to the generated mapping.", "", "", false},
            {"--topology", "-tp", "Filepath to the topology graph for the distance oracle.", "", "", false},
            {"--config", "-c", "The preset: fast, eco, strong or super-strong (increasing effort). The options marked (preset) take their default from it.", "", "", false},
            {"--mode", "", "How the mapping is found; the uncoarsening with refinement is the same for all: spectral (spectral Laplacian-RCB placement of the coarsest level), multilevel (greedy partition of the coarsest level), or both (greedy partition of the coarsest level). In spectral and both mode, the spectral placement of the finest level is computed independently and replaces the result if it is cheaper.", "spectral", "", false},
            {"--imbalance", "-e", "Allowed imbalance of the PE weights (for example 0.03). It applies to the partition of the input graph into PEs (KaMinPar) and to the final level of the individual input vertices; all coarser levels map one vertex to each PE.", "0.03", "", false},
            {"--prepartitioner", "", "Partitioner of the input graph into one block per PE if it has more vertices than PEs: kaminpar or spectral3d (multilevel three-eigenvector spectral embedding cut by recursive coordinate bisection).", "kaminpar", "", false},
            {"--kaminpar-threads", "", "Threads for the KaMinPar partition (0 = all hardware threads). Independent of --threads.", "0", "", false},
            {"--kaminpar-preset", "", "KaMinPar preset used for that partition.", "largek", "", false},
            {"--threads", "-t", "Number of threads of the coarsening and the prepartitioner (the rest of the pipeline is single-threaded).", "1", "", false},
            {"--seed", "", "Seed for diversifying results (random if not given).", "", "", false},

            {"--initial-kappa", "", "(preset) Greedy partitions tried on the coarsest level; the best is kept (multilevel and both mode).", "", "", false},
            {"--topo-split-factor", "", "(preset) While graph levels are left, the topology is only uncoarsened if the graph has at least this many vertices per finer PE on average.", "", "", false},
            {"--swap-radius", "", "(preset) Neighborhood radius of the swap refinement (1 = direct neighbors, 0 = off).", "", "", false},
            {"--swap-extended", "", "(preset) Also consider the neighborhoods of communication partners for swaps (0/1).", "", "", false},
            {"--long-msg-padding", "", "(preset) Rectangle padding of the long-message swap refinement (0 = off).", "", "", false},
            {"--sa-iterations", "", "(preset) Number of simulated annealing swap iterations on the finest level (0 = off).", "", "", false},
            {"--sa-radius", "", "(preset) Max swap radius of the simulated annealing.", "", "", false},
            {"--vertex-refine", "", "(preset) Run the message-group refinement on every level, including the individual input vertices (0/1).", "", "", false},
            {"--vertex-refine-iters", "", "(preset) Maximum iterations of the message-group refinement on one level.", "", "", false},
            {"--refine-verbose", "", "Print the statistics of every message-group refinement iteration (0/1).", "", "", false},
            {"--vcycle", "", "(preset) After the mapping is found, run a V-cycle: coarsen topology and graph together (one vertex per PE on every level), then uncoarsen with the same legalization and refinements as the main uncoarsening on every level (0/1).", "", "", false},

            {"--spectral-rounds", "", "Re-smoothing rounds of the spectral embedding after each prolongation (levels up to --spectral-big-n vertices).", "3", "", false},
            {"--spectral-depth", "", "CG iterations per re-smoothing round (levels up to --spectral-big-n vertices).", "24", "", false},
            {"--spectral-big-rounds", "", "Re-smoothing rounds on levels with more than --spectral-big-n vertices (the prepartitioner of a large input graph).", "1", "", false},
            {"--spectral-big-depth", "", "CG iterations per round on levels with more than --spectral-big-n vertices.", "8", "", false},
            {"--spectral-big-n", "", "Vertex count above which the --spectral-big-* settings apply.", "500000", "", false},

            {"--alpha-send", "", "Fixed send-side overhead per message (alpha_send).", "100.0", "", false},
            {"--alpha-recv", "", "Fixed receive-side overhead per message (alpha_recv).", "100.0", "", false},
            {"--alpha-link", "", "Fixed per-hop router overhead per message (alpha_link).", "2.0", "", false},
            {"--beta-send", "", "Volume-dependent send overhead per byte (beta_send).", "0.0", "", false},
            {"--beta-recv", "", "Volume-dependent receive overhead per byte (beta_recv).", "8.0", "", false},
            {"--beta-link", "", "Volume-dependent per-hop link cost per byte (beta_link).", "1.0", "", false},

            {"--coarsening", "", "(preset) Coarsening algorithm (global-paths, size-constrained-lp, heavy-edge).", "", "", false},
            {"--gpa-rating-function", "", "(preset) GPA rating function (weight, expansion, expansion*, expansion**, innerouter).", "", "", false},
            {"--gpa-random-level", "", "(preset) GPA random level.", "", "", false},
            {"--gpa-edge-rating-tiebreaking", "", "(preset) GPA use edge rating tiebreaking (true/false).", "", "", false},
            {"--gpa-two-hop-threshold", "", "(preset) GPA two-hop threshold.", "", "", false},
            {"--sclp-max-rounds", "", "(preset) SCLP max rounds.", "", "", false},
            {"--sclp-min-threshold", "", "(preset) SCLP min threshold.", "", "", false},
            {"--sclp-f", "", "(preset) SCLP tuning parameter f.", "", "", false},
            {"--sclp-rating-function", "", "(preset) SCLP rating function (weight, expansion, expansion*, expansion**, innerouter).", "", "", false},
            {"--sclp-use-degree-ordering", "", "(preset) SCLP use degree ordering (true/false).", "", "", false},
            {"--hem-rating-function", "", "(preset) HEM rating function (weight, expansion, expansion*, expansion**, innerouter).", "", "", false},
        };

    public:
        // graph information
        std::string graph_in;
        std::string topology_in;
        std::string mapping_out;
        partition_t k = 0; // run state: the number of PEs of the current topology level (set by the solver)

        f64 imbalance = 0.03;
        std::string mode = "spectral"; // spectral | multilevel | both
        bool do_spectral = true;       // derived from mode: spectral placement (spectral, both)
        bool do_multilevel = false;    // derived from mode: greedy partition of the coarsest level (multilevel, both)
        std::string prepartitioner = "kaminpar";
        std::string kaminpar_preset = "largek";
        u64 kaminpar_threads = 0;
        u64 seed = 0;
        u64 threads = 1;
        SpectralSmoothing spectral_smoothing;

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

        // ---- set by the preset ----
        std::string config_name = "undefined";
        u64 initial_c = 64;          // the graph is coarsened to about initial_c vertices per (coarse) PE
        u64 initial_kappa = 10;      // greedy partitions of the coarsest level
        u64 topo_split_factor = 4;
        bool vcycle = false;

        // coarsening algorithm
        std::string coarsening_algorithm_string;
        COARSENING_ALGS coarsening_algorithm_id = COARSENING_ALG_UNDEFINED;
        GlobalPathAlgorithmConfiguration global_path_algorithm_config;
        SizeConstrainedLPConfiguration size_constrained_lp_config;
        HeavyEdgeMatchingConfiguration heavy_edge_matching_config;

        // refinement algorithms
        WaverLabelPropagationConfiguration      waver_label_propagation_config      = WaverLabelPropagationConfiguration("Waver Label Propagation");
        WaverMessageRefinementConfiguration     waver_message_refinement_config     = WaverMessageRefinementConfiguration("Waver Message Refinement");
        WaverSwapRefinementConfiguration        waver_swap_refinement_config        = WaverSwapRefinementConfiguration("Waver Swap Refinement");
        WaverLongMessageRefinementConfiguration waver_long_message_refinement_config = WaverLongMessageRefinementConfiguration("Waver Long Message Refinement");
        WaverSARefinementConfiguration          waver_sa_refinement_config          = WaverSARefinementConfiguration("Waver Simulated Annealing");

        WaverHeiProMapConfiguration() = default;

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
                        if (i + 1 >= argc) {
                            std::cout << "Error: Command Line \"" << args[i] << "\" needs a value!" << std::endl;
                            exit(EXIT_FAILURE);
                        }
                        input = args[i + 1];
                        is_set = true;
                        i += 1;
                        break;
                    }
                }
            }

            graph_in = get("--graph");
            if (is_set("--topology")) { topology_in = get("--topology"); }
            if (is_set("--mapping")) { mapping_out = get("--mapping"); }

            imbalance = std::stod(get("--imbalance"));
            mode = get("--mode");
            if (mode != "spectral" && mode != "multilevel" && mode != "both") {
                std::cout << "Error: Invalid --mode '" << mode
                          << "'. Allowed values: spectral, multilevel, both" << std::endl;
                exit(EXIT_FAILURE);
            }
            do_spectral = mode != "multilevel";
            do_multilevel = mode != "spectral";
            prepartitioner = get("--prepartitioner");
            if (prepartitioner != "kaminpar" && prepartitioner != "spectral3d") {
                std::cout << "Error: Invalid --prepartitioner '" << prepartitioner
                          << "'. Allowed values: kaminpar, spectral3d" << std::endl;
                exit(EXIT_FAILURE);
            }
            kaminpar_preset = get("--kaminpar-preset");
            kaminpar_threads = std::stoull(get("--kaminpar-threads"));
            seed = is_set("--seed") ? std::stoull(get("--seed")) : std::random_device{}();
            threads = std::stoull(get("--threads"));

            spectral_smoothing.rounds = std::stoi(get("--spectral-rounds"));
            spectral_smoothing.depth = std::stoi(get("--spectral-depth"));
            spectral_smoothing.big_rounds = std::stoi(get("--spectral-big-rounds"));
            spectral_smoothing.big_depth = std::stoi(get("--spectral-big-depth"));
            spectral_smoothing.big_n = std::stoull(get("--spectral-big-n"));
            if (spectral_smoothing.rounds < 0 || spectral_smoothing.depth < 1
                || spectral_smoothing.big_rounds < 0 || spectral_smoothing.big_depth < 1) {
                std::cout << "Error: --spectral-rounds and --spectral-big-rounds must be >= 0, "
                          << "--spectral-depth and --spectral-big-depth >= 1" << std::endl;
                exit(EXIT_FAILURE);
            }

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
                set_fast();
            } else if (config_name == "eco") {
                set_eco();
            } else if (config_name == "strong") {
                set_strong();
            } else if (config_name == "super-strong") {
                set_super_strong();
            } else {
                std::cout << "Config " << config_name << " not recognized!" << std::endl;
                exit(EXIT_FAILURE);
            }

            apply_command_line_overrides();

            // Propagate cost-model parameters into the refinement configs
            waver_label_propagation_config.alpha_end  = alpha_end();
            waver_label_propagation_config.beta_end   = beta_end();
            waver_label_propagation_config.alpha_link = alpha_link;
            waver_label_propagation_config.beta_link  = beta_link;
        }

        // The presets differ in how much effort they put in; the settings of the effort scale are explained in
        // set_effort. Everything else (the coarsening algorithms) is chosen per preset below.

        void set_fast() {
            config_name = "fast";
            initial_c = 64;
            set_effort(/*kappa*/ 10, /*swap radius*/ 0, /*swap extended*/ false, /*SA iterations*/ 1000,
                       /*message refinement iterations*/ 10, /*V-cycle*/ false);
            set_global_paths_coarsening();

            global_path_algorithm_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;
            global_path_algorithm_config.random_level = 0;
            global_path_algorithm_config.use_edge_rating_tiebreaking = false;
            heavy_edge_matching_config.rating_function = EdgeRatingFunction::EXPANSION;
            size_constrained_lp_config.max_rounds = 5;
            size_constrained_lp_config.min_threshold = 0.05;
            size_constrained_lp_config.f = 32.0;
            size_constrained_lp_config.rating_function = EdgeRatingFunction::WEIGHT;
            size_constrained_lp_config.use_degree_ordering = true;

            waver_label_propagation_config.enabled = true;
            waver_label_propagation_config.max_iteration = 25;
        }

        void set_eco() {
            config_name = "eco";
            initial_c = 64;
            set_effort(25, 1, false, 2000, 10, false);
            set_global_paths_coarsening();

            global_path_algorithm_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;
            global_path_algorithm_config.random_level = 0;
            global_path_algorithm_config.use_edge_rating_tiebreaking = false;
            global_path_algorithm_config.two_hop_threshold = 0.75;
            heavy_edge_matching_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;
            size_constrained_lp_config.max_rounds = 5;
            size_constrained_lp_config.min_threshold = 0.05;
            size_constrained_lp_config.f = 16.0;

            waver_label_propagation_config.enabled = true;
            waver_label_propagation_config.max_iteration = 5;
        }

        void set_strong() {
            config_name = "strong";
            initial_c = 64;
            set_effort(50, 2, false, 5000, 20, true);
            set_global_paths_coarsening();

            global_path_algorithm_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;
            global_path_algorithm_config.random_level = 0;
            global_path_algorithm_config.use_edge_rating_tiebreaking = false;
            global_path_algorithm_config.two_hop_threshold = 0.75;
            heavy_edge_matching_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;
            size_constrained_lp_config.max_rounds = 5;
            size_constrained_lp_config.min_threshold = 0.05;
            size_constrained_lp_config.f = 16.0;

            waver_label_propagation_config.enabled = true;
            waver_label_propagation_config.max_iteration = 5;
        }

        void set_super_strong() {
            config_name = "super-strong";
            initial_c = 16;
            set_effort(100, 2, true, 10000, 50, true);
            set_global_paths_coarsening();

            global_path_algorithm_config.random_level = 0;
            size_constrained_lp_config.max_rounds = 5;
            size_constrained_lp_config.min_threshold = 0.10;
            size_constrained_lp_config.f = 8.0;

            waver_label_propagation_config.enabled = true;
            waver_label_propagation_config.max_iteration = 5;
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

    private:
        void set_global_paths_coarsening() {
            coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);
        }

        // The effort related settings of a preset: the greedy attempts on the coarsest level, the swap refinement, the
        // simulated annealing and the message-group refinement of every level, and the V-cycle. All presets keep the
        // long-message refinement (padding 1) and message refinement on, because they matter most for the quality.
        void set_effort(const u64 kappa, const u64 swap_radius, const bool swap_extended, const u64 sa_iterations,
                        const u64 message_iterations, const bool use_vcycle) {
            initial_kappa = kappa;
            topo_split_factor = 4;

            waver_swap_refinement_config.radius = swap_radius;
            waver_swap_refinement_config.enabled = swap_radius > 0;
            waver_swap_refinement_config.extended = swap_extended;

            waver_long_message_refinement_config.padding = 1;
            waver_long_message_refinement_config.enabled = true;

            waver_sa_refinement_config.iterations = sa_iterations;
            waver_sa_refinement_config.enabled = sa_iterations > 0;
            waver_sa_refinement_config.radius = 10;

            waver_message_refinement_config.enabled = true;
            waver_message_refinement_config.max_iters = message_iterations;

            vcycle = use_vcycle;
        }

        // The options of the same names as the preset's settings override the preset.
        void apply_command_line_overrides() {
            if (is_set("--initial-kappa")) initial_kappa = std::stoull(get("--initial-kappa"));
            if (is_set("--topo-split-factor")) topo_split_factor = std::stoull(get("--topo-split-factor"));
            if (is_set("--swap-radius")) {
                waver_swap_refinement_config.radius = std::stoull(get("--swap-radius"));
                waver_swap_refinement_config.enabled = waver_swap_refinement_config.radius > 0;
            }
            if (is_set("--swap-extended")) waver_swap_refinement_config.extended = std::stoull(get("--swap-extended")) != 0;
            if (is_set("--long-msg-padding")) {
                waver_long_message_refinement_config.padding = std::stoull(get("--long-msg-padding"));
                waver_long_message_refinement_config.enabled = waver_long_message_refinement_config.padding > 0;
            }
            if (is_set("--sa-iterations")) {
                waver_sa_refinement_config.iterations = std::stoull(get("--sa-iterations"));
                waver_sa_refinement_config.enabled = waver_sa_refinement_config.iterations > 0;
            }
            if (is_set("--sa-radius")) waver_sa_refinement_config.radius = std::stoull(get("--sa-radius"));
            if (is_set("--vertex-refine")) waver_message_refinement_config.enabled = std::stoull(get("--vertex-refine")) != 0;
            if (is_set("--vertex-refine-iters")) waver_message_refinement_config.max_iters = std::stoull(get("--vertex-refine-iters"));
            if (is_set("--refine-verbose")) waver_message_refinement_config.verbose = std::stoull(get("--refine-verbose")) != 0;
            if (is_set("--vcycle")) vcycle = std::stoull(get("--vcycle")) != 0;

            if (is_set("--coarsening")) {
                coarsening_algorithm_string = get("--coarsening");
                coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);
            }

            // GPA
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

            // SCLP
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

            // HEM
            if (is_set("--hem-rating-function")) {
                heavy_edge_matching_config.rating_function = string_to_edge_rating_function(get("--hem-rating-function"));
            }
        }
    };
}

#endif //HEIPROMAP_WAVERHEIPROMAP_CONFIGURATION_H
