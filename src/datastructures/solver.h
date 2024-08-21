#ifndef SERIALPROCESSMAPPING_SOLVER_H
#define SERIALPROCESSMAPPING_SOLVER_H

#include "../utility/definitions.h"
#include "../utility/macros.h"
#include "../utility/utils.h"
#include "graph.h"
#include "statistic_collector.h"
#include "../partitioning/kaffpa_partitioner.h"
#include "../utility/qap.h"
#include "../coarsening/simple_edge_matcher.h"
#include "../coarsening/greedy_edge_matcher.h"
#include "../coarsening/heavy_edge_matcher.h"
#include "../refinement/identity_refinement.h"
#include "../refinement/label_propagation_refinement.h"
#include "partition_manager.h"
#include "../refinement/quotient_graph_refinement.h"
#include "../coarsening/simple_clustering.h"

namespace SPM {
    class Solver {
    private:
        // graph
        Graph g;
        std::string graph_in;

        // distance
        std::vector<u64> hierarchy;
        std::vector<u64> distance;
        u64 k;
        f64 imbalance;
        DistanceOracle dist_o;

        // balance
        u64 lmax = 0;

        // partition
        PartitionManager pm;

        // multilevel
        u64 threshold;

        // matching
        GreedyEdgeMatcher gem;
        HeavyEdgeMatcher hem;
        SimpleClustering sc;

        std::vector<std::vector<Edge>> matches_per_level;

        // partitioning
        KaffpaPartitioner kp;

        // refinement
        LabelPropagationRefinement lpr;
        IdentityRefinement ir;
        QuotientGraphRefinement qgr;

        // statistics
        StatisticCollector stat_collect;
        u64 qap = 0;
        vertex_t n_active = 0;
        size_t matches_size = 0;

    public:
        Solver(const std::string &graph_in,
               std::vector<u64> &hierarchy,
               std::vector<u64> &distance,
               u64 k,
               f64 imbalance)
                :
                graph_in(graph_in),
                hierarchy(hierarchy),
                distance(distance),
                k(k),
                imbalance(imbalance) {
            auto sp_graph_io = std::chrono::high_resolution_clock::now();
            g = Graph(graph_in);
            auto ep_graph_io = std::chrono::high_resolution_clock::now();

            auto sp_io = std::chrono::high_resolution_clock::now();

            // distance
            dist_o = DistanceOracle(k, hierarchy, distance);

            // balance
            lmax = get_lmax(imbalance, k, g.get_sum_vertex_weights());

            // partition
            pm.initialize(&g, k);

            // multilevel
            threshold = g.get_n() / 500;

            // matching
            gem.initialize(&g);
            hem.initialize(&g);
            sc.initialize(&g);

            // partitioning
            kp.initialize(&g, hierarchy, distance, k, imbalance);

            // refinement
            lpr.initialize(&g, hierarchy, distance, k, imbalance, lmax, &dist_o);
            qgr.initialize(&g, hierarchy, distance, k, imbalance, lmax, &dist_o);

            auto ep_io = std::chrono::high_resolution_clock::now();
            stat_collect.set_io(get_seconds(sp_graph_io, ep_graph_io), get_seconds(sp_io, ep_io));
        }

        std::vector<vertex_t> solve() {
            recursive_solve(0);

#if STATISTICCOLLECTOR
            qap = get_qap(g, pm, hierarchy, distance);
#endif
            stat_collect.set_final(qap, pm.get_pweights(), lmax);
            stat_collect.finalize();

            std::cout << stat_collect.to_JSON() << std::endl;

            return pm.get_partition();
        }

    private:
        void recursive_solve(s32 level) {

            /*************
             * PARTITION *
             *************/
            if (g.get_n_active() <= threshold || g.get_n_active() <= k * 64) {
                // small enough, now partition
                auto sp_partition = std::chrono::high_resolution_clock::now();
                kp.partition(pm, FAST);
                pm.init_after_partition();
                auto ep_partition = std::chrono::high_resolution_clock::now();

#if STATISTICCOLLECTOR
                qap = get_qap(g, pm, hierarchy, distance);
#endif

                stat_collect.set_partition(get_seconds(sp_partition, ep_partition), qap, pm.get_pweights(), lmax);
                return;
            }

            /***********************
             * MATCHING/CLUSTERING *
             ***********************/
            auto sp_match = std::chrono::high_resolution_clock::now();
            matches_per_level.emplace_back();
            matches_per_level.back().reserve(g.get_n_active() / 2);

            // SimpleEdgeMatcher sem;
            // sem.match(g, matches, marker, level);

            hem.match(matches_per_level.back(), level);
            // gem.match(matches, level);
            // sc.match(matches, level);

            auto ep_match = std::chrono::high_resolution_clock::now();

#if STATISTICCOLLECTOR
            matches_size = matches_per_level.back().size();
#endif
            stat_collect.set_matching(get_seconds(sp_match, ep_match), level, matches_size);

            /************
             * COARSING *
             ************/
            auto sp_coarse = std::chrono::high_resolution_clock::now();
            for (auto &e: matches_per_level.back()) {
                g.contract_edge(e.u, e.v);
            }
            auto ep_coarse = std::chrono::high_resolution_clock::now();

#if STATISTICCOLLECTOR
            n_active = g.get_n_active();
#endif
            stat_collect.set_coarsening(get_seconds(sp_coarse, ep_coarse), level, n_active);

            /**************
             * NEXT LEVEL *
             **************/
            recursive_solve(level + 1);

            /**************
             * UNCOARSING *
             **************/
            auto sp_uncoarse = std::chrono::high_resolution_clock::now();
            for (size_t i = 0; i < matches_per_level.back().size(); ++i) {
                u64 idx = matches_per_level.back().size() - 1 - i;
                vertex_t u = matches_per_level.back()[idx].u;
                vertex_t v = matches_per_level.back()[idx].v;
                pm.uncontract_edge(u, v);
            }
            matches_per_level.pop_back();
            auto ep_uncoarse = std::chrono::high_resolution_clock::now();

#if STATISTICCOLLECTOR
            n_active = g.get_n_active();
#endif
            stat_collect.set_uncoarsening(get_seconds(sp_uncoarse, ep_uncoarse), level, n_active);

            /**************
             * REFINEMENT *
             **************/
            auto sp_refinement = std::chrono::high_resolution_clock::now();

            // ir.refine(partition, pweights);
            lpr.refine(pm);
            // qgr.refine(pm);

            auto ep_refinement = std::chrono::high_resolution_clock::now();

#if STATISTICCOLLECTOR
            qap = get_qap(g, pm, hierarchy, distance);
#endif
            stat_collect.set_refinement(get_seconds(sp_refinement, ep_refinement), level, qap);
        }
    };
}

#endif //SERIALPROCESSMAPPING_SOLVER_H
