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

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "test_utils.h"
#include "../../src/serial/coarsening/global_path_algorithm_arrays.h"
#include "../../src/serial/coarsening/heavy_edge_matcher.h"
#include "../../src/serial/datastructures/active_vertex_manager.h"
#include "../../src/serial/datastructures/graph_csr.h"
#include "../../src/serial/datastructures/graph_csr_arrays.h"
#include "../../src/serial/datastructures/simple_graph.h"
#include "../../src/serial/datastructures/sorted_graph_csr.h"

namespace HeiProMap {
    std::vector<EdgeUV> read_solution(std::string &solution) {
        std::vector<EdgeUV> edges;
        std::ifstream ifs(solution);

        int u, v;
        // Read two integers per line until we reach the end of the file.
        while (ifs >> u >> v) {
            edges.emplace_back(u, v);
        }

        return edges;
    }

    void compare(std::string &graph_in, std::string &solution) {
        std::cout << graph_in << " " << solution << std::endl;

        // calculate solution
        GraphCSRArrays g(graph_in);

        ActiveVertexManager av_manager;
        av_manager.initialize(g.get_n());

        GlobalPathAlgorithmArraysMatcher gpa_matcher;
        GlobalPathAlgorithmConfiguration gpa_config;
        gpa_matcher.initialize(g.get_n(), g.get_m(), 2, std::numeric_limits<weight_t>::max(), 0);

        size_t n_64 = round_up_64(g.get_n() / 2);
        EdgeUV* matches = (EdgeUV*) aligned_alloc(64, n_64 * sizeof(EdgeUV));
        size_t matches_size = 0;
        gpa_matcher.match(gpa_config, g, av_manager, matches, matches_size);

        // read in solution
        std::vector<EdgeUV> solution_edges = read_solution(solution);

        // matchings must have the same weight, but not necessarily the same edges
        weight_t own_weight = 0;
        for (size_t i = 0; i < matches_size; i++) {
            std::cout << matches[i].u << " A " << matches[i].v << std::endl;
            for (auto [v, w] : g[matches[i].u]) {
                if (v == matches[i].v) {
                    own_weight += w;
                    break;
                }
            }
        }

        weight_t solution_weight = 0;
        for (size_t i = 0; i < solution_edges.size(); i++) {
            std::cout << solution_edges[i].u << " B " << solution_edges[i].v << std::endl;
            for (auto [v, w] : g[solution_edges[i].u]) {
                if (v == solution_edges[i].v) {
                    solution_weight += w;
                    break;
                }
            }
        }

        EXPECT_EQ(own_weight, solution_weight);
    }

    TEST(GlobalPathAlgorithm, chain_10_graph) {
        const std::string graph_in = "../data/test/chain_10.graph";

        GraphCSRArrays g(graph_in);

        ActiveVertexManager av_manager;
        av_manager.initialize(g.get_n());

        GlobalPathAlgorithmArraysMatcher gpa_matcher;
        GlobalPathAlgorithmConfiguration gpa_config;
        gpa_matcher.initialize(g.get_n(), g.get_m(), 2, std::numeric_limits<weight_t>::max(), 0);

        size_t n_64 = round_up_64(g.get_n() / 2);
        EdgeUV* matches = (EdgeUV*) aligned_alloc(64, n_64 * sizeof(EdgeUV));
        size_t matches_size = 0;
        gpa_matcher.match(gpa_config, g, av_manager, matches, matches_size);

        EXPECT_EQ(matches_size, 4);
        EdgeUV e1(0, 0);
        EdgeUV e2(0, 0);
        bool found_e1;
        bool found_e2;

        e1 = EdgeUV(8, 7);
        e2 = EdgeUV(7, 8);
        found_e1 = std::find(matches, matches + matches_size, e1) != matches + matches_size;
        found_e2 = std::find(matches, matches + matches_size, e2) != matches + matches_size;
        EXPECT_TRUE(found_e1 || found_e2);

        e1 = EdgeUV(6, 5);
        e2 = EdgeUV(5, 6);
        found_e1 = std::find(matches, matches + matches_size, e1) != matches + matches_size;
        found_e2 = std::find(matches, matches + matches_size, e2) != matches + matches_size;
        EXPECT_TRUE(found_e1 || found_e2);

        e1 = EdgeUV(4, 3);
        e2 = EdgeUV(3, 4);
        found_e1 = std::find(matches, matches + matches_size, e1) != matches + matches_size;
        found_e2 = std::find(matches, matches + matches_size, e2) != matches + matches_size;
        EXPECT_TRUE(found_e1 || found_e2);

        e1 = EdgeUV(2, 1);
        e2 = EdgeUV(1, 2);
        found_e1 = std::find(matches, matches + matches_size, e1) != matches + matches_size;
        found_e2 = std::find(matches, matches + matches_size, e2) != matches + matches_size;
        EXPECT_TRUE(found_e1 || found_e2);

        free(matches);
    }

    TEST(GlobalPathAlgorithm, gpa_graphs) {
        int n_min = 5;
        int n_max = 5;
        int reps = 12;
        for (int n = n_min; n <= n_max; n++) {
            for (int rep = 11; rep < reps; rep++) {
                std::string graph_in = "../test_data/gpa_graphs/" + std::to_string(rep) + "_" + std::to_string(n) + ".graph";
                std::string solution = "../test_data/gpa_solutions/matching_" + std::to_string(rep) + "_" + std::to_string(n) + ".edges";
                compare(graph_in, solution);
            }
        }
    }
}