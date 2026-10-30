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
#include "../refinement/waver_long_message_refinement.h"
#include "../refinement/waver_merge_refinement.h"
#include "../refinement/waver_busy_refinement.h"
#include "../refinement/waver_cycle_refinement.h"
#include "../refinement/waver_crowd_refinement.h"
#include "../refinement/waver_shift_refinement.h"
#include "../refinement/waver_smooth_refinement.h"
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
    // related settings; the options of the sections marked "overrides the preset" replace its settings.
    class WaverHeiProMapConfiguration {
    public:
        std::vector<CommandLineOption> options = {
            {"# Input and output", "", "", "", "", false},
            {"--help", "", "Produces the help message", "", "", false},
            {"--graph", "-g", "Filepath to the graph.", "", "", false},
            {"--topology", "-tp", "Filepath to the topology graph for the distance oracle.", "", "", false},
            {"--mapping", "-m", "Output filepath to the generated mapping.", "", "", false},
            {"--config", "-c", "The preset: fast, eco, strong or super-strong (increasing effort). The options of the sections marked \"overrides the preset\" replace its settings.", "", "", false},
            {"--mode", "", "How the mapping is found; the uncoarsening with refinement is the same for both: spectral (spectral Laplacian-RCB placement of the coarsest level; the spectral placement of the finest level is computed independently and replaces the result if it is cheaper) or multilevel (greedy partition of the coarsest level).", "spectral", "", false},

            {"# General", "", "", "", "", false},
            {"--imbalance", "-e", "Allowed imbalance of the PE weights (for example 0.03). It applies to the partition of the input graph into PEs (KaMinPar) and to the final level of the individual input vertices; all coarser levels map one vertex to each PE.", "0.03", "", false},
            {"--min-fill", "", "Multilevel mode: the lowest weight a PE may be refined down to, as a fraction of the average PE weight (for example 0.5). A coarse cell stands for several PEs and keeps that fraction of their total. Without a lower limit the refinement drains PEs to cut messages, and the topology can no longer be uncoarsened in time.", "0.5", "", false},
            {"--seed", "", "Seed for diversifying results (random if not given).", "", "", false},
            {"--threads", "-t", "Number of threads of the coarsening and the prepartitioner (the rest of the pipeline is single-threaded).", "1", "", false},

            {"# Prepartition of the input graph", "", "", "", "", false},
            {"--prepartitioner", "", "Partitioner of the input graph into one block per PE if it has more vertices than PEs: kaminpar or spectral3d (multilevel three-eigenvector spectral embedding cut by recursive coordinate bisection).", "spectral3d", "", false},
            {"--kaminpar-preset", "", "KaMinPar preset used for that partition.", "largek", "", false},
            {"--kaminpar-threads", "", "Threads for the KaMinPar partition (0 = all hardware threads). Independent of --threads.", "0", "", false},

            {"# Spectral embedding", "", "", "", "", false},
            {"--spectral-rounds", "", "Re-smoothing rounds of the spectral embedding after each prolongation (levels up to --spectral-big-n vertices).", "3", "", false},
            {"--spectral-depth", "", "CG iterations per re-smoothing round (levels up to --spectral-big-n vertices).", "24", "", false},
            {"--spectral-big-rounds", "", "Re-smoothing rounds on levels with more than --spectral-big-n vertices (the prepartitioner of a large input graph).", "1", "", false},
            {"--spectral-big-depth", "", "CG iterations per round on levels with more than --spectral-big-n vertices.", "8", "", false},
            {"--spectral-big-n", "", "Vertex count above which the --spectral-big-* settings apply.", "500000", "", false},

            {"# Communication cost model", "", "", "", "", false},
            {"--alpha-send", "", "Fixed send-side overhead per message (alpha_send).", "100.0", "", false},
            {"--alpha-recv", "", "Fixed receive-side overhead per message (alpha_recv).", "100.0", "", false},
            {"--alpha-link", "", "Fixed per-hop router overhead per message (alpha_link).", "2.0", "", false},
            {"--beta-send", "", "Volume-dependent send overhead per byte (beta_send).", "0.0", "", false},
            {"--beta-recv", "", "Volume-dependent receive overhead per byte (beta_recv).", "8.0", "", false},
            {"--beta-link", "", "Volume-dependent per-hop link cost per byte (beta_link).", "1.0", "", false},

            {"# Coarsening (overrides the preset)", "", "", "", "", false},
            {"--coarsening-vertices-per-pe", "", "The graph is coarsened to about this many vertices per (coarse) PE.", "", "", false},
            {"--coarsening", "", "Coarsening algorithm (global-paths, size-constrained-lp, heavy-edge).", "", "", false},
            {"--gpa-rating-function", "", "Global-paths coarsening: rating function (weight, expansion, expansion*, expansion**, innerouter).", "", "", false},
            {"--gpa-random-level", "", "Global-paths coarsening: the coarsening levels below this one use a random matching instead of the rating-based one.", "", "", false},
            {"--gpa-edge-rating-tiebreaking", "", "Global-paths coarsening: also sort the edges of a thread when all its ratings are equal (true/false).", "", "", false},
            {"--gpa-two-hop-threshold", "", "Global-paths coarsening: minimum matched fraction below which two-hop matching is tried.", "", "", false},
            {"--sclp-max-rounds", "", "Size-constrained label propagation coarsening: maximum rounds.", "", "", false},
            {"--sclp-min-threshold", "", "Size-constrained label propagation coarsening: stop when fewer than this fraction of the vertices move.", "", "", false},
            {"--sclp-f", "", "Size-constrained label propagation coarsening: the cluster weight limit is the PE weight limit divided by f.", "", "", false},
            {"--sclp-rating-function", "", "Size-constrained label propagation coarsening: rating function (weight, expansion, expansion*, expansion**, innerouter).", "", "", false},
            {"--sclp-use-degree-ordering", "", "Size-constrained label propagation coarsening: visit vertices by increasing degree (true/false).", "", "", false},
            {"--hem-rating-function", "", "Heavy-edge matching coarsening: rating function (weight, expansion, expansion*, expansion**, innerouter).", "", "", false},

            {"# Initial mapping (overrides the preset)", "", "", "", "", false},
            {"--greedy-attempts", "", "Greedy partitions tried on the coarsest level; the best is kept (multilevel mode).", "", "", false},

            {"# Uncoarsening and refinement (overrides the preset)", "", "", "", "", false},
            {"--topology-uncoarsening-factor", "", "While graph levels are left, the topology is only uncoarsened if the graph has at least this many vertices per finer PE on average.", "", "", false},
            {"--message-refinement", "", "Run the message-group refinement on every level, including the individual input vertices (0/1).", "", "", false},
            {"--message-refinement-iterations", "", "Maximum iterations of the message-group refinement on one level.", "", "", false},
            {"--long-message-padding", "", "Rectangle padding of the long-message swap refinement (0 = off).", "", "", false},
            {"--merge-refinement", "", "Multilevel mode: run the merge refinement, which moves a sender group into a block without room for it and passes the overload on along a chain of moves (0/1).", "", "", false},
            {"--merge-depth", "", "Maximum number of moves in one chain of the merge refinement.", "", "", false},
            {"--merge-rounds", "", "Maximum rounds of the merge refinement on one level.", "", "", false},
            {"--merge-partial", "", "Merge refinement: an overloaded block may also give up only part of a sender group, the vertices tied most closely to the receiving block (0/1).", "", "", false},
            {"--merge-split", "", "Merge refinement: a sender group may also be split among up to this many blocks, each vertex going to the block it has the most edge weight to (below 2: never).", "", "", false},
            {"--merge-scatter", "", "Merge refinement without chains: the sender group is spread over the receiver and all its neighbors, within their room, and no block is overloaded (0/1).", "", "", false},
            {"--shift-refinement", "", "Multilevel mode: run the shift refinement, which places the blocks of a long message next to each other and shifts the blocks on the way by one cell (0/1).", "", "", false},
            {"--shift-rounds", "", "Maximum rounds of the shift refinement on one level.", "", "", false},
            {"--shift-messages", "", "Shift refinement: the number of longest messages worked on per round.", "", "", false},
            {"--smooth-refinement", "", "Multilevel mode: run the smooth refinement, which moves single boundary vertices to lower the neighbor cells, regions and hops (0/1).", "", "", false},
            {"--smooth-rounds", "", "Maximum sweeps of the smooth refinement on one level.", "", "", false},
            {"--smooth-hop-weight", "", "Smooth refinement: cost of one hop of one cut edge.", "", "", false},
            {"--smooth-verbose", "", "Print the statistics of every sweep of the smooth refinement (0/1).", "", "", false},
            {"--shift-verbose", "", "Print the statistics of every round of the shift refinement (0/1).", "", "", false},
            {"--crowd-refinement", "", "Multilevel mode: run the crowd refinement, which shortens and separates the bounding boxes of the messages (0/1).", "", "", false},
            {"--crowd-rounds", "", "Maximum rounds of the crowd refinement on one level.", "", "", false},
            {"--crowd-top", "", "Crowd refinement: the share of messages, most crowded first, worked on per round.", "", "", false},
            {"--crowd-verbose", "", "Print the statistics of every round of the crowd refinement (0/1).", "", "", false},
            {"--cycle-refinement", "", "Multilevel mode: run the cycle refinement, which moves sender groups around a cycle of blocks so that every block gives one and receives one (0/1).", "", "", false},
            {"--cycle-length", "", "Maximum number of blocks on a cycle of the cycle refinement.", "", "", false},
            {"--cycle-rounds", "", "Maximum rounds of the cycle refinement on one level.", "", "", false},
            {"--cycle-tries", "", "Cycle refinement: first moves tried per block and round.", "", "", false},
            {"--cycle-verbose", "", "Print the statistics of every round of the cycle refinement (0/1).", "", "", false},
            {"--busy-refinement", "", "Multilevel mode: run the busy-PE refinement, which moves sender groups away from the blocks with the highest engine load (0/1).", "", "", false},
            {"--busy-rounds", "", "Maximum rounds of the busy-PE refinement on one level.", "", "", false},
            {"--busy-top", "", "Busy-PE refinement: the share of blocks, busiest first, worked on per round.", "", "", false},
            {"--busy-power", "", "Busy-PE refinement: the loads touched by a move are compared as sum(load^power).", "", "", false},
            {"--busy-alpha", "", "Busy-PE refinement: load per message (0: the alpha of the cost model).", "", "", false},
            {"--busy-beta", "", "Busy-PE refinement: load per unit of volume (0: the beta of the cost model).", "", "", false},
            {"--busy-exact", "", "Busy-PE refinement: the load of a block counts its send and receive regions and its values as well as its messages (0/1).", "", "", false},
            {"--busy-tries", "", "Busy-PE refinement: sender groups tried per busy block and round.", "", "", false},
            {"--busy-verbose", "", "Print the statistics of every round of the busy-PE refinement (0/1).", "", "", false},
            {"--merge-verbose", "", "Print the statistics of every round of the merge refinement (0/1).", "", "", false},
            {"--swap-radius", "", "Neighborhood radius of the swap refinement (1 = direct neighbors, 0 = off).", "", "", false},
            {"--swap-extended-neighborhood", "", "Also consider the neighborhoods of communication partners for swaps (0/1).", "", "", false},
            {"--annealing-iterations", "", "Number of simulated annealing swap iterations on every level (0 = off).", "0", "", false},
            {"--annealing-radius", "", "Max swap radius of the simulated annealing.", "", "", false},
            {"--smash-rounds", "", "Multilevel mode: rounds of smash and repair on the finished mapping: the windows around the busiest PEs are placed again from scratch and the mapping is repaired; kept if the busiest PE gets less load (0: off).", "", "", false},
            {"--smash-count", "", "Smash and repair: number of busiest PEs smashed per round.", "", "", false},
            {"--smash-radius", "", "Smash and repair: the window reaches this many PEs on every side of the busiest PE.", "", "", false},
            {"--level-stats", "", "Multilevel mode: print metrics of the mapping on every level, after the level change and after the refinement (0/1).", "", "", false},
            {"--diffusion-equalize", "", "Multilevel mode: after a level change also pull vertices from above-average into below-average PEs, without regard to the boundary (0/1, default 1).", "", "", false},
            {"--intense-max-n", "", "Multilevel mode: on the finest topology, graphs up to this many vertices use the --intense-* merge settings (0: off).", "", "", false},
            {"--intense-depth", "", "Merge depth on the intense levels.", "", "", false},
            {"--intense-rounds", "", "Merge rounds on the intense levels.", "", "", false},
            {"--intense-partial", "", "Partial moves on the intense levels (0/1).", "", "", false},
            {"--intense-split", "", "Merge split parts on the intense levels.", "", "", false},
            {"--window-rounds", "", "Multilevel mode: rounds of the window re-embedding on every level of the finest topology (0: off).", "", "", false},
            {"--window-size", "", "Window re-embedding: cells per side of a window.", "", "", false},
            {"--window-min-n", "", "Window re-embedding: only on graphs with at least this many vertices.", "", "", false},
            {"--window-max-n", "", "Window re-embedding: only on graphs with at most this many vertices (0: no limit).", "", "", false},
            {"--relax-rounds", "", "Multilevel mode: rounds of the warm-started relaxation on every level of the finest topology (0: off).", "", "", false},
            {"--relax-sweeps", "", "Relaxation: smoothing sweeps of the vertex positions per round.", "", "", false},
            {"--relax-anchor", "", "Relaxation: weight of the pull back to the current cell, relative to the pull of the neighbors.", "", "", false},
            {"--relax-min-n", "", "Relaxation: only on graphs with at least this many vertices.", "", "", false},
            {"--diffusion-rebalance", "", "Multilevel mode: before the rebalancer, move the excess weight of overloaded PEs along the grid by a diffusion flow, one hop at a time, taking the boundary vertices that fit the receiving PE best (0/1).", "", "", false},
            {"--dump-hierarchy", "", "Multilevel mode: write the cluster of every input vertex on every graph level to <prefix>_<level>.txt (level 1 is the first coarse graph).", "", "", false},
            {"--spectral-rotations", "", "Multilevel mode: turns of the two smoothest vectors tried in the spectral placement of the coarsest level, evenly over a quarter turn (1: no turn).", "", "", false},
            {"--global-split", "", "Multilevel mode: at every topology split also place the whole current graph spectrally onto the finer grid, and keep it if its message cost is lower than that of the local splits (0/1).", "", "", false},
            {"--terminal-propagation", "", "Multilevel mode: when a coarse PE is split, pull the placement of its vertices towards their neighbors in the neighboring cells (0/1).", "", "", false},
            {"--vcycle-vertices-per-pe", "", "Multilevel mode: the V-cycle coarsens the graph to about this many vertices per PE.", "", "", false},
            {"--vcycle", "", "Multilevel mode: after the mapping is found, run a V-cycle: coarsen the graph again within the PEs of the current mapping, then uncoarsen with legalization and refinement on every level; the old mapping is kept if the cost does not drop (0/1).", "", "", false},
        };

    public:
        // graph information
        std::string graph_in;
        std::string topology_in;
        std::string mapping_out;
        partition_t k = 0; // run state: the number of PEs of the current topology level (set by the solver)

        f64 imbalance = 0.03;
        f64 min_fill = 0.5;
        std::string mode = "spectral"; // spectral | multilevel
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
        u64 smash_rounds = 0;   // multilevel: rounds of smash and repair on the finished mapping (0: off)
        u64 smash_count = 10;   // busiest PEs smashed per round
        u64 smash_radius = 1;   // a smashed window reaches this many PEs on every side of the busiest one
        bool level_stats = false;          // multilevel: print metrics of the mapping on every level
        bool diffusion_equalize = true;    // multilevel: the diffusion also equalizes the weights, not only fills empty PEs
        u64 intense_max_n = 0;             // multilevel: graphs up to this size, on the finest topology, merge more intensely
        u64 intense_depth = 3, intense_rounds = 6, intense_split = 3;
        bool intense_partial = true;
        u64 window_rounds = 0, window_size = 3, window_min_n = 0, window_max_n = 0; // multilevel: window re-embedding
        u64 relax_rounds = 0, relax_sweeps = 8, relax_min_n = 0; // multilevel: warm-started relaxation
        f64 relax_anchor = 0.5;
        bool diffusion_rebalance = false;  // multilevel: a diffusion flow on the grid moves the excess weight before the rebalancer
        std::string dump_hierarchy;        // multilevel: prefix of the files the clusters of every level are written to
        u64 spectral_rotations = 1;        // multilevel: turns of the two smoothest vectors tried in the coarsest spectral placement
        bool global_split = false;         // multilevel: a global spectral placement competes with the local splits
        bool terminal_propagation = false; // multilevel: the split of a coarse cell is pulled towards the neighbors outside it
        u64 vcycle_c = 8;                 // the V-cycle coarsens the graph to about this many vertices per PE

        // coarsening algorithm
        std::string coarsening_algorithm_string;
        COARSENING_ALGS coarsening_algorithm_id = COARSENING_ALG_UNDEFINED;
        GlobalPathAlgorithmConfiguration global_path_algorithm_config;
        SizeConstrainedLPConfiguration size_constrained_lp_config;
        HeavyEdgeMatchingConfiguration heavy_edge_matching_config;

        // refinement algorithms
        WaverMessageRefinementConfiguration     waver_message_refinement_config     = WaverMessageRefinementConfiguration("Waver Message Refinement");
        WaverSwapRefinementConfiguration        waver_swap_refinement_config        = WaverSwapRefinementConfiguration("Waver Swap Refinement");
        WaverLongMessageRefinementConfiguration waver_long_message_refinement_config = WaverLongMessageRefinementConfiguration("Waver Long Message Refinement");
        WaverShiftRefinementConfiguration       waver_shift_refinement_config       = WaverShiftRefinementConfiguration("Waver Shift Refinement");
        WaverSmoothRefinementConfiguration      waver_smooth_refinement_config      = WaverSmoothRefinementConfiguration("Waver Smooth Refinement");
        WaverCrowdRefinementConfiguration       waver_crowd_refinement_config       = WaverCrowdRefinementConfiguration("Waver Crowd Refinement");
        WaverCycleRefinementConfiguration       waver_cycle_refinement_config       = WaverCycleRefinementConfiguration("Waver Cycle Refinement");
        WaverBusyRefinementConfiguration        waver_busy_refinement_config        = WaverBusyRefinementConfiguration("Waver Busy Refinement");
        WaverMergeRefinementConfiguration       waver_merge_refinement_config       = WaverMergeRefinementConfiguration("Waver Merge Refinement");
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
            min_fill = std::stod(get("--min-fill"));
            mode = get("--mode");
            if (mode != "spectral" && mode != "multilevel") {
                std::cout << "Error: Invalid --mode '" << mode
                          << "'. Allowed values: spectral, multilevel" << std::endl;
                exit(EXIT_FAILURE);
            }
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

        }

        // The presets differ in how much effort they put in (greedy attempts, refinements, V-cycle) and in the coarsening
        // algorithm settings.

        void set_fast() {
            config_name = "fast";
            initial_c = 64;
            initial_kappa = 10;
            topo_split_factor = 4;
            vcycle = false;

            coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id = COARSENING_ALG_GLOBAL_PATHS;

            global_path_algorithm_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;
            global_path_algorithm_config.random_level = 0;
            global_path_algorithm_config.use_edge_rating_tiebreaking = false;
            heavy_edge_matching_config.rating_function = EdgeRatingFunction::EXPANSION;
            size_constrained_lp_config.max_rounds = 5;
            size_constrained_lp_config.min_threshold = 0.05;
            size_constrained_lp_config.f = 32.0;
            size_constrained_lp_config.rating_function = EdgeRatingFunction::WEIGHT;
            size_constrained_lp_config.use_degree_ordering = true;

            waver_message_refinement_config.enabled = true;
            waver_message_refinement_config.max_iters = 1;

            waver_merge_refinement_config.enabled = true;
            waver_merge_refinement_config.max_depth = 2;

            waver_long_message_refinement_config.enabled = false;
            waver_long_message_refinement_config.padding = 2;

            waver_swap_refinement_config.enabled = false;
            waver_swap_refinement_config.radius = 3;
            waver_swap_refinement_config.max_rounds = 3;
            waver_swap_refinement_config.extended = false;

            waver_sa_refinement_config.enabled = false;
            waver_sa_refinement_config.iterations = 1;
            waver_sa_refinement_config.radius = 10;
        }

        void set_eco() {
            config_name = "eco";
            initial_c = 64;
            initial_kappa = 25;
            topo_split_factor = 4;
            vcycle = false;

            waver_message_refinement_config.enabled = true;
            waver_message_refinement_config.max_iters = 10;
            waver_long_message_refinement_config.padding = 1;
            waver_long_message_refinement_config.enabled = false;
            waver_swap_refinement_config.radius = 1;
            waver_swap_refinement_config.enabled = true;
            waver_swap_refinement_config.extended = false;
            waver_sa_refinement_config.iterations = 2000;
            waver_sa_refinement_config.enabled = true;
            waver_sa_refinement_config.radius = 10;

            coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id = COARSENING_ALG_GLOBAL_PATHS;

            global_path_algorithm_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;
            global_path_algorithm_config.random_level = 0;
            global_path_algorithm_config.use_edge_rating_tiebreaking = false;
            global_path_algorithm_config.two_hop_threshold = 0.75;
            heavy_edge_matching_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;
            size_constrained_lp_config.max_rounds = 5;
            size_constrained_lp_config.min_threshold = 0.05;
            size_constrained_lp_config.f = 16.0;
        }

        void set_strong() {
            config_name = "strong";
            initial_c = 64;
            initial_kappa = 50;
            topo_split_factor = 4;
            vcycle = false;

            waver_message_refinement_config.enabled = true;
            waver_message_refinement_config.max_iters = 20;

            waver_merge_refinement_config.enabled = true;
            waver_merge_refinement_config.max_depth = 2;
            waver_busy_refinement_config.enabled = true;
            waver_busy_refinement_config.max_rounds = 20;
            waver_busy_refinement_config.top_fraction = 0.5;
            waver_busy_refinement_config.power = 8.0;
            waver_busy_refinement_config.max_tries = 30;
            waver_long_message_refinement_config.padding = 1;
            waver_long_message_refinement_config.enabled = false;
            waver_swap_refinement_config.radius = 2;
            waver_swap_refinement_config.enabled = true;
            waver_swap_refinement_config.extended = false;
            waver_sa_refinement_config.iterations = 5000;
            waver_sa_refinement_config.enabled = true;
            waver_sa_refinement_config.radius = 10;

            coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id = COARSENING_ALG_GLOBAL_PATHS;

            global_path_algorithm_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;
            global_path_algorithm_config.random_level = 0;
            global_path_algorithm_config.use_edge_rating_tiebreaking = false;
            global_path_algorithm_config.two_hop_threshold = 0.75;
            heavy_edge_matching_config.rating_function = EdgeRatingFunction::EXPANSIONSTAR;
            size_constrained_lp_config.max_rounds = 5;
            size_constrained_lp_config.min_threshold = 0.05;
            size_constrained_lp_config.f = 16.0;
        }

        void set_super_strong() {
            config_name = "super-strong";
            initial_c = 16;
            initial_kappa = 100;
            topo_split_factor = 4;
            vcycle = true;

            waver_message_refinement_config.enabled = true;
            waver_message_refinement_config.max_iters = 50;
            waver_long_message_refinement_config.padding = 1;
            waver_long_message_refinement_config.enabled = false;
            waver_swap_refinement_config.radius = 2;
            waver_swap_refinement_config.enabled = true;
            waver_swap_refinement_config.extended = true;
            waver_sa_refinement_config.iterations = 10000;
            waver_sa_refinement_config.enabled = true;
            waver_sa_refinement_config.radius = 10;

            coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id = COARSENING_ALG_GLOBAL_PATHS;

            global_path_algorithm_config.random_level = 0;
            size_constrained_lp_config.max_rounds = 5;
            size_constrained_lp_config.min_threshold = 0.10;
            size_constrained_lp_config.f = 8.0;
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
                if (large_key[0] == '#') {
                    std::cout << "\n" << large_key.substr(2) << ":" << std::endl;
                } else if (small_key.empty()) {
                    std::cout << "  " << large_key << " - " << description << std::endl;
                } else {
                    std::cout << "  " << large_key << ", " << small_key << " - " << description << std::endl;
                }
            }
        }

    private:
        // The options of the same names as the preset's settings override the preset.
        void apply_command_line_overrides() {
            if (is_set("--coarsening-vertices-per-pe")) initial_c = std::stoull(get("--coarsening-vertices-per-pe"));
            if (is_set("--greedy-attempts")) initial_kappa = std::stoull(get("--greedy-attempts"));
            if (is_set("--topology-uncoarsening-factor")) topo_split_factor = std::stoull(get("--topology-uncoarsening-factor"));
            if (is_set("--swap-radius")) {
                waver_swap_refinement_config.radius = std::stoull(get("--swap-radius"));
                waver_swap_refinement_config.enabled = waver_swap_refinement_config.radius > 0;
            }
            if (is_set("--swap-extended-neighborhood")) waver_swap_refinement_config.extended = std::stoull(get("--swap-extended-neighborhood")) != 0;
            if (is_set("--merge-refinement")) waver_merge_refinement_config.enabled = std::stoull(get("--merge-refinement")) != 0;
            if (is_set("--merge-depth")) waver_merge_refinement_config.max_depth = std::stoull(get("--merge-depth"));
            if (is_set("--merge-rounds")) waver_merge_refinement_config.max_rounds = std::stoull(get("--merge-rounds"));
            if (is_set("--merge-partial")) waver_merge_refinement_config.partial = std::stoull(get("--merge-partial")) != 0;
            if (is_set("--merge-split")) waver_merge_refinement_config.split_parts = std::stoull(get("--merge-split"));
            if (is_set("--merge-scatter")) waver_merge_refinement_config.scatter = std::stoull(get("--merge-scatter")) != 0;
            if (is_set("--smooth-refinement")) waver_smooth_refinement_config.enabled = std::stoull(get("--smooth-refinement")) != 0;
            if (is_set("--smooth-rounds")) waver_smooth_refinement_config.max_rounds = std::stoull(get("--smooth-rounds"));
            if (is_set("--smooth-hop-weight")) waver_smooth_refinement_config.hop_weight = std::stod(get("--smooth-hop-weight"));
            if (is_set("--smooth-verbose")) waver_smooth_refinement_config.verbose = std::stoull(get("--smooth-verbose")) != 0;
            if (is_set("--shift-refinement")) waver_shift_refinement_config.enabled = std::stoull(get("--shift-refinement")) != 0;
            if (is_set("--shift-rounds")) waver_shift_refinement_config.max_rounds = std::stoull(get("--shift-rounds"));
            if (is_set("--shift-messages")) waver_shift_refinement_config.max_messages = std::stoull(get("--shift-messages"));
            if (is_set("--shift-verbose")) waver_shift_refinement_config.verbose = std::stoull(get("--shift-verbose")) != 0;
            if (is_set("--crowd-refinement")) waver_crowd_refinement_config.enabled = std::stoull(get("--crowd-refinement")) != 0;
            if (is_set("--crowd-rounds")) waver_crowd_refinement_config.max_rounds = std::stoull(get("--crowd-rounds"));
            if (is_set("--crowd-top")) waver_crowd_refinement_config.top_fraction = std::stod(get("--crowd-top"));
            if (is_set("--crowd-verbose")) waver_crowd_refinement_config.verbose = std::stoull(get("--crowd-verbose")) != 0;
            if (is_set("--cycle-refinement")) waver_cycle_refinement_config.enabled = std::stoull(get("--cycle-refinement")) != 0;
            if (is_set("--cycle-length")) waver_cycle_refinement_config.max_length = std::stoull(get("--cycle-length"));
            if (is_set("--cycle-rounds")) waver_cycle_refinement_config.max_rounds = std::stoull(get("--cycle-rounds"));
            if (is_set("--cycle-tries")) waver_cycle_refinement_config.max_tries = std::stoull(get("--cycle-tries"));
            if (is_set("--cycle-verbose")) waver_cycle_refinement_config.verbose = std::stoull(get("--cycle-verbose")) != 0;
            if (is_set("--busy-refinement")) waver_busy_refinement_config.enabled = std::stoull(get("--busy-refinement")) != 0;
            if (is_set("--busy-rounds")) waver_busy_refinement_config.max_rounds = std::stoull(get("--busy-rounds"));
            if (is_set("--busy-top")) waver_busy_refinement_config.top_fraction = std::stod(get("--busy-top"));
            if (is_set("--busy-power")) waver_busy_refinement_config.power = std::stod(get("--busy-power"));
            if (is_set("--busy-alpha")) waver_busy_refinement_config.alpha = std::stod(get("--busy-alpha"));
            if (is_set("--busy-beta")) waver_busy_refinement_config.beta = std::stod(get("--busy-beta"));
            if (is_set("--busy-exact")) waver_busy_refinement_config.exact = std::stoull(get("--busy-exact")) != 0;
            if (is_set("--busy-tries")) waver_busy_refinement_config.max_tries = std::stoull(get("--busy-tries"));
            if (is_set("--busy-verbose")) waver_busy_refinement_config.verbose = std::stoull(get("--busy-verbose")) != 0;
            if (is_set("--merge-verbose")) waver_merge_refinement_config.verbose = std::stoull(get("--merge-verbose")) != 0;
            if (is_set("--long-message-padding")) {
                waver_long_message_refinement_config.padding = std::stoull(get("--long-message-padding"));
                waver_long_message_refinement_config.enabled = waver_long_message_refinement_config.padding > 0;
            }
            if (is_set("--annealing-iterations")) {
                waver_sa_refinement_config.iterations = std::stoull(get("--annealing-iterations"));
                waver_sa_refinement_config.enabled = waver_sa_refinement_config.iterations > 0;
            }
            if (is_set("--annealing-radius")) waver_sa_refinement_config.radius = std::stoull(get("--annealing-radius"));
            if (is_set("--message-refinement")) waver_message_refinement_config.enabled = std::stoull(get("--message-refinement")) != 0;
            if (is_set("--message-refinement-iterations")) waver_message_refinement_config.max_iters = std::stoull(get("--message-refinement-iterations"));
            if (is_set("--smash-rounds")) smash_rounds = std::stoull(get("--smash-rounds"));
            if (is_set("--smash-count")) smash_count = std::stoull(get("--smash-count"));
            if (is_set("--smash-radius")) smash_radius = std::stoull(get("--smash-radius"));
            if (is_set("--level-stats")) level_stats = std::stoull(get("--level-stats")) != 0;
            if (is_set("--diffusion-equalize")) diffusion_equalize = std::stoull(get("--diffusion-equalize")) != 0;
            if (is_set("--intense-max-n")) intense_max_n = std::stoull(get("--intense-max-n"));
            if (is_set("--intense-depth")) intense_depth = std::stoull(get("--intense-depth"));
            if (is_set("--intense-rounds")) intense_rounds = std::stoull(get("--intense-rounds"));
            if (is_set("--intense-partial")) intense_partial = std::stoull(get("--intense-partial")) != 0;
            if (is_set("--intense-split")) intense_split = std::stoull(get("--intense-split"));
            if (is_set("--window-rounds")) window_rounds = std::stoull(get("--window-rounds"));
            if (is_set("--window-size")) window_size = std::stoull(get("--window-size"));
            if (is_set("--window-min-n")) window_min_n = std::stoull(get("--window-min-n"));
            if (is_set("--window-max-n")) window_max_n = std::stoull(get("--window-max-n"));
            if (is_set("--relax-rounds")) relax_rounds = std::stoull(get("--relax-rounds"));
            if (is_set("--relax-sweeps")) relax_sweeps = std::stoull(get("--relax-sweeps"));
            if (is_set("--relax-anchor")) relax_anchor = std::stod(get("--relax-anchor"));
            if (is_set("--relax-min-n")) relax_min_n = std::stoull(get("--relax-min-n"));
            if (is_set("--diffusion-rebalance")) diffusion_rebalance = std::stoull(get("--diffusion-rebalance")) != 0;
            if (is_set("--dump-hierarchy")) dump_hierarchy = get("--dump-hierarchy");
            if (is_set("--spectral-rotations")) spectral_rotations = std::stoull(get("--spectral-rotations"));
            if (is_set("--global-split")) global_split = std::stoull(get("--global-split")) != 0;
            if (is_set("--terminal-propagation")) terminal_propagation = std::stoull(get("--terminal-propagation")) != 0;
            if (is_set("--vcycle-vertices-per-pe")) vcycle_c = std::stoull(get("--vcycle-vertices-per-pe"));
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
