#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "test_utils.h"
#include "../../src/serial/datastructures/graph_csr.h"
#include "../../src/serial/datastructures/simple_graph.h"

namespace HeiProMap {

    void compare_initialization(const std::string &graph_in){
        auto sp_g = std::chrono::high_resolution_clock::now();
        SimpleGraph simple_graph(graph_in);
        auto ep_g = std::chrono::high_resolution_clock::now();

        auto sp_csr_g = std::chrono::high_resolution_clock::now();
        GraphCSR csr_g(graph_in);
        auto ep_csr_g = std::chrono::high_resolution_clock::now();

        double t_g = get_seconds(sp_g, ep_g);
        double t_csr_g = get_seconds(sp_csr_g, ep_csr_g);
        std::cout << "initialization " << graph_in << " : g = " << t_g << " csr_g = " << t_csr_g << " speedup = " << t_g / t_csr_g << std::endl;

        graphs_are_equal(simple_graph, csr_g);
    }

    TEST(GraphInitialization, graph_csrgraph_PGPgiantcompo_graph) {
        const std::string graph_in = "../data/mapping/PGPgiantcompo.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_cop20k_A_mtx_graph) {
        const std::string graph_in = "../data/mapping/cop20k_A.mtx.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_2cubes_sphere_mtx_graph) {
        const std::string graph_in = "../data/mapping/2cubes_sphere.mtx.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_thermomech_TC_mtx_graph) {
        const std::string graph_in = "../data/mapping/thermomech_TC.mtx.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_598a_graph) {
        const std::string graph_in = "../data/mapping/598a.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_cfd2_mtx_graph) {
        const std::string graph_in = "../data/mapping/cfd2.mtx.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_boneS01_mtx_graph) {
        const std::string graph_in = "../data/mapping/boneS01.mtx.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_fe_ocean_graph) {
        const std::string graph_in = "../data/mapping/fe_ocean.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_144_graph) {
        const std::string graph_in = "../data/mapping/144.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_Dubcova3_mtx_graph) {
        const std::string graph_in = "../data/mapping/Dubcova3.mtx.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_bmwcra_1_mtx_graph) {
        const std::string graph_in = "../data/mapping/bmwcra_1.mtx.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_G2_circuit_mtx_graph) {
        const std::string graph_in = "../data/mapping/G2_circuit.mtx.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_wave_graph) {
        const std::string graph_in = "../data/mapping/wave.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_shipsec5_mtx_graph) {
        const std::string graph_in = "../data/mapping/shipsec5.mtx.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_cont_300_mtx_graph) {
        const std::string graph_in = "../data/mapping/cont-300.mtx.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_m14b_graph) {
        const std::string graph_in = "../data/mapping/m14b.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_auto_graph) {
        const std::string graph_in = "../data/mapping/auto.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_afshell9_graph) {
        const std::string graph_in = "../data/mapping/afshell9.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_thermal2_graph) {
        const std::string graph_in = "../data/mapping/thermal2.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_nlr_graph) {
        const std::string graph_in = "../data/mapping/nlr.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_deu_graph) {
        const std::string graph_in = "../data/mapping/deu.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_del23_graph) {
        const std::string graph_in = "../data/mapping/del23.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_rgg23_graph) {
        const std::string graph_in = "../data/mapping/rgg23.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_del24_graph) {
        const std::string graph_in = "../data/mapping/del24.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_rgg24_graph) {
        const std::string graph_in = "../data/mapping/rgg24.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_eur_graph) {
        const std::string graph_in = "../data/mapping/eur.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_del26_graph) {
        const std::string graph_in = "../data/mapping/del26.graph";
        compare_initialization(graph_in);
    }

    TEST(GraphInitialization, graph_csrgraph_rgg_n26_graph) {
        const std::string graph_in = "../data/mapping/rgg_n26.graph";
        compare_initialization(graph_in);
    }
}

