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
#include "../../src/serial/coarsening/heavy_edge_matcher.h"
#include "../../src/serial/datastructures/active_vertex_manager.h"
#include "../../src/serial/datastructures/graph_csr.h"
#include "../../src/serial/datastructures/simple_graph.h"

namespace HeiProMap {

    void compare_coarsening(const std::string &graph_in){
        SimpleGraph g(graph_in);
        GraphCSR csr_g(graph_in);

        ActiveVertexManager av_manager;
        av_manager.initialize(g.get_n());

        HeavyEdgeMatcher he_matcher;
        he_matcher.initialize(g.get_n(), g.get_m(), 0, std::numeric_limits<weight_t>::max());

        std::vector<EdgeUV> matches;
        matches.reserve(av_manager.get_n_active() / 2);
        he_matcher.match(g, av_manager, matches);

        auto g_sp = std::chrono::high_resolution_clock::now();
        SimpleGraph g1(g, matches); // coarse the graph
        auto g_ep = std::chrono::high_resolution_clock::now();

        auto csr_g_sp = std::chrono::high_resolution_clock::now();
        GraphCSR csr_g1(csr_g, matches); // coarse the graph
        auto csr_g_ep = std::chrono::high_resolution_clock::now();

        double t_g = get_seconds(g_sp, g_ep);
        double t_csr_g = get_seconds(csr_g_sp, csr_g_ep);
        std::cout << "coarsening " << graph_in << " : g = " << t_g << " csr_g = " << t_csr_g << " speedup = " << t_g / t_csr_g << std::endl;

        graphs_are_equal(g, csr_g);
        graphs_are_equal(g1, csr_g1);
    }

    TEST(GraphCoarsening, graph_csrgraph_PGPgiantcompo_graph) {
        const std::string graph_in = "../data/mapping/PGPgiantcompo.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_cop20k_A_mtx_graph) {
        const std::string graph_in = "../data/mapping/cop20k_A.mtx.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_2cubes_sphere_mtx_graph) {
        const std::string graph_in = "../data/mapping/2cubes_sphere.mtx.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_thermomech_TC_mtx_graph) {
        const std::string graph_in = "../data/mapping/thermomech_TC.mtx.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_598a_graph) {
        const std::string graph_in = "../data/mapping/598a.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_cfd2_mtx_graph) {
        const std::string graph_in = "../data/mapping/cfd2.mtx.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_boneS01_mtx_graph) {
        const std::string graph_in = "../data/mapping/boneS01.mtx.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_fe_ocean_graph) {
        const std::string graph_in = "../data/mapping/fe_ocean.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_144_graph) {
        const std::string graph_in = "../data/mapping/144.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_Dubcova3_mtx_graph) {
        const std::string graph_in = "../data/mapping/Dubcova3.mtx.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_bmwcra_1_mtx_graph) {
        const std::string graph_in = "../data/mapping/bmwcra_1.mtx.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_G2_circuit_mtx_graph) {
        const std::string graph_in = "../data/mapping/G2_circuit.mtx.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_wave_graph) {
        const std::string graph_in = "../data/mapping/wave.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_shipsec5_mtx_graph) {
        const std::string graph_in = "../data/mapping/shipsec5.mtx.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_cont_300_mtx_graph) {
        const std::string graph_in = "../data/mapping/cont-300.mtx.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_m14b_graph) {
        const std::string graph_in = "../data/mapping/m14b.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_auto_graph) {
        const std::string graph_in = "../data/mapping/auto.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_afshell9_graph) {
        const std::string graph_in = "../data/mapping/afshell9.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_thermal2_graph) {
        const std::string graph_in = "../data/mapping/thermal2.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_nlr_graph) {
        const std::string graph_in = "../data/mapping/nlr.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_deu_graph) {
        const std::string graph_in = "../data/mapping/deu.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_del23_graph) {
        const std::string graph_in = "../data/mapping/del23.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_rgg23_graph) {
        const std::string graph_in = "../data/mapping/rgg23.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_del24_graph) {
        const std::string graph_in = "../data/mapping/del24.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_rgg24_graph) {
        const std::string graph_in = "../data/mapping/rgg24.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_eur_graph) {
        const std::string graph_in = "../data/mapping/eur.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_del26_graph) {
        const std::string graph_in = "../data/mapping/del26.graph";
        // compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_rgg_n26_graph) {
        const std::string graph_in = "../data/mapping/rgg_n26.graph";
        // compare_coarsening(graph_in);
    }
}