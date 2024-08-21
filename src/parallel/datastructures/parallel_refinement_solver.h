#ifndef SERIALPROCESSMAPPING_PARALLEL_REFINEMENT_SOLVER_H
#define SERIALPROCESSMAPPING_PARALLEL_REFINEMENT_SOLVER_H

#include "../../utility/definitions.h"
#include "../../utility/macros.h"
#include "../../utility/utils.h"
#include "../../utility/qap.h"
#include "../../datastructures/graph.h"
#include "../../datastructures/distance_oracle.h"
#include "../../datastructures/partition_manager.h"
#include "../../coarsening/greedy_edge_matcher.h"
#include "../../coarsening/heavy_edge_matcher.h"
#include "../../coarsening/simple_clustering.h"
#include "../../partitioning/kaffpa_partitioner.h"
#include "../../refinement/label_propagation_refinement.h"
#include "../../refinement/identity_refinement.h"
#include "../../refinement/quotient_graph_refinement.h"
#include "../../datastructures/statistic_collector.h"
#include "../refinement/parallel_label_propagation_refinement.h"

namespace SPM {

    class ParallelRefinementSolver {
    private:
        // graph
        Graph g;
        std::string graph_in;

        // partition
        PartitionManager pm;
        std::string mapping_in;

        // distance
        std::vector<u64> hierarchy;
        std::vector<u64> distance;
        u64 k;
        f64 imbalance;
        DistanceOracle dist_o;

        // threads
        u64 n_threads = 1;

        // balance
        u64 lmax = 0;

        // refinement
        ParallelLabelPropagationRefinement mt_label_prop;

        // statistics
        StatisticCollector stat_collect;
        u64 qap = 0;
        vertex_t n_active = 0;
        size_t matches_size = 0;

    public:
        ParallelRefinementSolver(const std::string &graph_in,
                                 const std::string &mapping_in,
                                 std::vector<u64> &hierarchy,
                                 std::vector<u64> &distance,
                                 u64 k,
                                 f64 imbalance,
                                 u64 n_threads)
                :
                graph_in(graph_in),
                mapping_in(mapping_in),
                hierarchy(hierarchy),
                distance(distance),
                k(k),
                imbalance(imbalance),
                n_threads(n_threads) {
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
            pm.set_partition(read_partition(mapping_in));

            // refinement
            mt_label_prop.initialize(&g, hierarchy, distance, k, imbalance, lmax, &dist_o, n_threads);

            auto ep_io = std::chrono::high_resolution_clock::now();
            stat_collect.set_io(get_seconds(sp_graph_io, ep_graph_io), get_seconds(sp_io, ep_io));
        }

        std::vector<vertex_t> solve() {
            partition();
            refine();

            stat_collect.set_final(qap, pm.get_pweights(), lmax);
            stat_collect.finalize();
            std::cout << stat_collect.to_JSON() << std::endl;

            return pm.get_partition();
        }

    private:
        void partition() {
            auto sp_partition = std::chrono::high_resolution_clock::now();

            pm.init_after_partition();

            auto ep_partition = std::chrono::high_resolution_clock::now();
#if STATISTICCOLLECTOR
            qap = get_qap(g, pm, hierarchy, distance);
#endif

            stat_collect.set_partition(get_seconds(sp_partition, ep_partition), qap, pm.get_pweights(), lmax);
        }

        void refine() {
            auto sp_refinement = std::chrono::high_resolution_clock::now();

            mt_label_prop.refine(pm);

            auto ep_refinement = std::chrono::high_resolution_clock::now();

#if STATISTICCOLLECTOR
            qap = get_qap(g, pm, hierarchy, distance);
#endif
            stat_collect.set_refinement(get_seconds(sp_refinement, ep_refinement), 0, qap);

#if STATISTICCOLLECTOR
            qap = get_qap(g, pm, hierarchy, distance);
#endif
        }
    };

}


#endif //SERIALPROCESSMAPPING_PARALLEL_REFINEMENT_SOLVER_H
