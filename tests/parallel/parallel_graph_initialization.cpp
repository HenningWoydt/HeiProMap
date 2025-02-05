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

#include <vector>
#include <string>

#include <gtest/gtest.h>
#include "../../src/serial/datastructures/graph.h"
#include "../serial/test_utils.h"
#include "../../src/parallel/datastructures/parallel_graph.h"
#include "../../src/parallel/datastructures/parallel_static_csr_graph.h"

namespace HeiProMap {

    void compare_parallel_initialization(const std::string &graph_in){
        auto sp_g = std::chrono::high_resolution_clock::now();
        Graph g;g.initialize(graph_in);
        auto ep_g = std::chrono::high_resolution_clock::now();

        auto sp_p_g = std::chrono::high_resolution_clock::now();
        ParallelGraph p_g;p_g.initialize(graph_in, 4);
        auto ep_p_g = std::chrono::high_resolution_clock::now();

        auto sp_static_csr_g = std::chrono::high_resolution_clock::now();
        ParallelStaticCSRGraph static_csr_g;static_csr_g.initialize(graph_in, 4);
        auto ep_static_csr_g = std::chrono::high_resolution_clock::now();

        double t_g = get_seconds(sp_g, ep_g);
        double t_p_g = get_seconds(sp_p_g, ep_p_g);
        double t_static_csr_g = get_seconds(sp_static_csr_g, ep_static_csr_g);
        std::cout << "initialization " << graph_in << " : g = " << t_g << " p_g = " << t_p_g << " speedup = " << t_g / t_p_g << std::endl;
        std::cout << "initialization " << graph_in << " : g = " << t_g << " static_csr_g = " << t_static_csr_g << " speedup = " << t_g / t_static_csr_g << std::endl;

        graphs_are_equal(g, p_g);
        graphs_are_equal(g, static_csr_g);
    }

    TEST(GraphParallelInitialization, graph_pgraph_PGPgiantcompo_graph) {
        std::string graph_in = "../data/mapping/PGPgiantcompo.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_cop20k_A_mtx_graph) {
        std::string graph_in = "../data/mapping/cop20k_A.mtx.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_2cubes_sphere_mtx_graph) {
        std::string graph_in = "../data/mapping/2cubes_sphere.mtx.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_thermomech_TC_mtx_graph) {
        std::string graph_in = "../data/mapping/thermomech_TC.mtx.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_598a_graph) {
        std::string graph_in = "../data/mapping/598a.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_cfd2_mtx_graph) {
        std::string graph_in = "../data/mapping/cfd2.mtx.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_boneS01_mtx_graph) {
        std::string graph_in = "../data/mapping/boneS01.mtx.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_fe_ocean_graph) {
        std::string graph_in = "../data/mapping/fe_ocean.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_144_graph) {
        std::string graph_in = "../data/mapping/144.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_Dubcova3_mtx_graph) {
        std::string graph_in = "../data/mapping/Dubcova3.mtx.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_bmwcra_1_mtx_graph) {
        std::string graph_in = "../data/mapping/bmwcra_1.mtx.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_G2_circuit_mtx_graph) {
        std::string graph_in = "../data/mapping/G2_circuit.mtx.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_wave_graph) {
        std::string graph_in = "../data/mapping/wave.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_shipsec5_mtx_graph) {
        std::string graph_in = "../data/mapping/shipsec5.mtx.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_cont_300_mtx_graph) {
        std::string graph_in = "../data/mapping/cont-300.mtx.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_m14b_graph) {
        std::string graph_in = "../data/mapping/m14b.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_auto_graph) {
        std::string graph_in = "../data/mapping/auto.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_afshell9_graph) {
        std::string graph_in = "../data/mapping/afshell9.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_thermal2_graph) {
        std::string graph_in = "../data/mapping/thermal2.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_nlr_graph) {
        std::string graph_in = "../data/mapping/nlr.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_deu_graph) {
        std::string graph_in = "../data/mapping/deu.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_del23_graph) {
        std::string graph_in = "../data/mapping/del23.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_rgg23_graph) {
        std::string graph_in = "../data/mapping/rgg23.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_del24_graph) {
        std::string graph_in = "../data/mapping/del24.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_rgg24_graph) {
        std::string graph_in = "../data/mapping/rgg24.graph";
        compare_parallel_initialization(graph_in);
    }

    TEST(GraphParallelInitialization, graph_pgraph_eur_graph) {
        std::string graph_in = "../data/mapping/eur.graph";
        compare_parallel_initialization(graph_in);
    }
}

