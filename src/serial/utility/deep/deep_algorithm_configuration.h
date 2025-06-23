#ifndef HEIPROMAP_DEEP_ALGORITHM_CONFIGURATION_H
#define HEIPROMAP_DEEP_ALGORITHM_CONFIGURATION_H

#include <vector>

#include "../../../commons/definitions.h"
#include "../algorithm_configuration.h"
#include "../../refinement/deep/deep_quotient_graph_refinement.h"
#include "../../refinement/deep/deep_flow_based_refinement.h"

namespace HeiProMap {

    class DeepAlgorithmConfiguration {

        std::vector<CommandLineOption> options = {
                {"--help",       "",   "Produces the help message",                 "",                     "", false},
                {"--graph",      "-g", "Filepath to the graph.",                    "",                     "", false},
                {"--mapping",    "-m", "Output filepath to the generated mapping.", "",                     "", false},
                {"--statistics", "",   "Output filepath to the statistics file.",   "HeiProMap_stats.JSON", "", false},
                {"--hierarchy",  "-h", "Hierarchy in the form a1:a2:...:al .",      "",                     "", false},
                {"--distance",   "-d", "Distance in the form d1:d2:...:dl .",       "",                     "", false},
                {"--imbalance",  "-e", "Allowed imbalance (for example 0.03).",     "0.03",                 "", false},
                {"--config",     "-c", "The configuration.",                        "",                     "", false},
                {"--threads",    "-t", "The number of threads.",                    "1",                    "", false},
                {"--seed",       "",   "Seed for diversifying results.",            "",                     "", false}
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

        // threads
        u64 threads = 1;

        // coarsening algorithm
        std::string coarsening_algorithm_string;
        COARSENING_ALGS coarsening_algorithm_id = COARSENING_ALG_UNDEFINED;

        GlobalPathAlgorithmConfiguration global_path_algorithm_config;

        // partitioning algorithm
        std::string partitioning_algorithm_string;
        PARTITIONING_ALGS partitioning_algorithm_id = PARTITIONING_ALG_UNDEFINED;

        KaffpaKWayPartitionerConfiguration kaffpa_kway_partitioner_config;

        // refinement algorithms
        DeepQuotientGraphRefinementConfiguration deep_quotient_graph_refinement_config = DeepQuotientGraphRefinementConfiguration("Deep Quotient Graph Refinement");
        DeepFlowBasedRefinementConfiguration deep_flow_based_refinement_config = DeepFlowBasedRefinementConfiguration("Deep Flow Based Refinement");

        DeepAlgorithmConfiguration() = default;

        DeepAlgorithmConfiguration(int argc, char *argv[]) {
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
            statistics_out = get("--statistics");

            hierarchy_string = get("--hierarchy");
            hierarchy = convert<partition_t>(split(hierarchy_string, ':'));
            k = prod<partition_t>(hierarchy);

            distance_string = get("--distance");
            distance = convert<weight_t>(split(distance_string, ':'));

            imbalance = std::stod(get("--imbalance"));

            if (is_set("--threads")) {
                threads = std::stoi(get("--threads"));
            }

            if (is_set("--seed")) {
                seed = std::stoi(get("--seed"));
            } else {
                seed = std::random_device{}();
            }

            // set FARAJ20 config
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
        }

        void set_fast() {

        }

        void set_eco() {

        }

        void set_strong() {

        }

        void set_experimental() {
            // set GPA matching algorithm
            coarsening_algorithm_string = "global-paths";
            coarsening_algorithm_id = string_to_coarsening_algorithm(coarsening_algorithm_string);

            // configurate global-paths algorithm
            global_path_algorithm_config.random_level = 0;

            // refinement
            deep_quotient_graph_refinement_config.enabled = false;
            deep_quotient_graph_refinement_config.max_iteration = 3;
            deep_quotient_graph_refinement_config.alpha = 1000.0;
            deep_quotient_graph_refinement_config.beta_factor = 1.0;
            deep_quotient_graph_refinement_config.alpha_edge_cut = 5.0;
            deep_quotient_graph_refinement_config.beta_factor_edge_cut = 1.0;

            deep_flow_based_refinement_config.enabled = false;
            deep_flow_based_refinement_config.min_level = 0;
            deep_flow_based_refinement_config.max_level = 100;
            deep_flow_based_refinement_config.max_global_iteration = 2;
            deep_flow_based_refinement_config.max_local_iteration = 5;
            deep_flow_based_refinement_config.alpha = 2.0;
            deep_flow_based_refinement_config.alpha_upper_bound = 16.0;
            deep_flow_based_refinement_config.alpha_modifier = 2.0;
            deep_flow_based_refinement_config.use_closed_vertex_set = true;
            deep_flow_based_refinement_config.closed_vertex_sets_repeats = 100;
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

#endif //HEIPROMAP_DEEP_ALGORITHM_CONFIGURATION_H
