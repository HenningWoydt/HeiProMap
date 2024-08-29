#include <vector>
#include <string>

#include <gtest/gtest.h>
#include "../../src/serial/datastructures/graph.h"
#include "../../src/serial/datastructures/csr_graph.h"
#include "test_utils.h"
#include "../../src/serial/datastructures/active_vertex_manager.h"
#include "../../src/serial/coarsening/heavy_edge_matcher.h"

namespace HeiProMap {

    void compare_coarsening(const std::string &graph_in){
        Graph g;g.initialize(graph_in);
        ActiveVertexManager g_av_manager;g_av_manager.initialize(&g);
        HeavyEdgeMatcher g_he_matcher;g_he_matcher.initialize(&g, &g_av_manager);
        std::vector<EdgeUV> g_matches;
        g_matches.reserve(g_av_manager.get_n_active() / 2);
        g_he_matcher.match(g_matches);

        auto g_sp = std::chrono::high_resolution_clock::now();
        for (auto &e: g_matches) {
            g.contract(e.u, e.v);
            g_av_manager.contract(e.u, e.v);
        }
        auto g_ep = std::chrono::high_resolution_clock::now();

        CSRGraph csr_g;csr_g.initialize(graph_in);
        ActiveVertexManager csr_g_av_manager;csr_g_av_manager.initialize(&csr_g);
        HeavyEdgeMatcher csr_g_he_matcher;csr_g_he_matcher.initialize(&csr_g, &csr_g_av_manager);
        std::vector<EdgeUV> csr_g_matches;
        csr_g_matches.reserve(csr_g_av_manager.get_n_active() / 2);
        csr_g_he_matcher.match(csr_g_matches);
        auto csr_g_sp = std::chrono::high_resolution_clock::now();
        for (auto &e: csr_g_matches) {
            csr_g.contract(e.u, e.v);
            csr_g_av_manager.contract(e.u, e.v);
        }
        auto csr_g_ep = std::chrono::high_resolution_clock::now();

        double t_g = get_seconds(g_sp, g_ep);
        double t_csr_g = get_seconds(csr_g_sp, csr_g_ep);
        std::cout << "coarsening " << graph_in << " : g = " << t_g << " csr_g = " << t_csr_g << " speedup = " << t_g / t_csr_g << std::endl;

        graphs_are_equal(g, g_av_manager, csr_g, csr_g_av_manager);
    }

    TEST(GraphCoarsening, graph_csrgraph_PGPgiantcompo_graph) {
        std::string graph_in = "../data/mapping/PGPgiantcompo.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_cop20k_A_mtx_graph) {
        std::string graph_in = "../data/mapping/cop20k_A.mtx.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_2cubes_sphere_mtx_graph) {
        std::string graph_in = "../data/mapping/2cubes_sphere.mtx.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_thermomech_TC_mtx_graph) {
        std::string graph_in = "../data/mapping/thermomech_TC.mtx.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_598a_graph) {
        std::string graph_in = "../data/mapping/598a.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_cfd2_mtx_graph) {
        std::string graph_in = "../data/mapping/cfd2.mtx.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_boneS01_mtx_graph) {
        std::string graph_in = "../data/mapping/boneS01.mtx.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_fe_ocean_graph) {
        std::string graph_in = "../data/mapping/fe_ocean.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_144_graph) {
        std::string graph_in = "../data/mapping/144.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_Dubcova3_mtx_graph) {
        std::string graph_in = "../data/mapping/Dubcova3.mtx.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_bmwcra_1_mtx_graph) {
        std::string graph_in = "../data/mapping/bmwcra_1.mtx.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_G2_circuit_mtx_graph) {
        std::string graph_in = "../data/mapping/G2_circuit.mtx.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_wave_graph) {
        std::string graph_in = "../data/mapping/wave.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_shipsec5_mtx_graph) {
        std::string graph_in = "../data/mapping/shipsec5.mtx.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_cont_300_mtx_graph) {
        std::string graph_in = "../data/mapping/cont-300.mtx.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_m14b_graph) {
        std::string graph_in = "../data/mapping/m14b.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_auto_graph) {
        std::string graph_in = "../data/mapping/auto.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_afshell9_graph) {
        std::string graph_in = "../data/mapping/afshell9.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_thermal2_graph) {
        std::string graph_in = "../data/mapping/thermal2.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_nlr_graph) {
        std::string graph_in = "../data/mapping/nlr.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_deu_graph) {
        std::string graph_in = "../data/mapping/deu.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_del23_graph) {
        std::string graph_in = "../data/mapping/del23.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_rgg23_graph) {
        std::string graph_in = "../data/mapping/rgg23.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_del24_graph) {
        std::string graph_in = "../data/mapping/del24.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_rgg24_graph) {
        std::string graph_in = "../data/mapping/rgg24.graph";
        compare_coarsening(graph_in);
    }

    TEST(GraphCoarsening, graph_csrgraph_eur_graph) {
        std::string graph_in = "../data/mapping/eur.graph";
        compare_coarsening(graph_in);
    }
}