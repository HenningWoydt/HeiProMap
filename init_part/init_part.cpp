// Standalone initial partitioning: reads a (coarsest) graph, computes an initial
// partition / mapping onto the hierarchy and reports QAP, balance and time.
//
// Usage: InitPart --graph g.graph --hierarchy 4:8:6 --distance 1:10:100 --imbalance 0.03
//                 --config strong [--partitioning multisection|greedy|flow-bisection]
//                 [--seed 0] [--threads 1] [--mapping out.txt]
// Defaults to the algorithm of the chosen --config (strong: multisection, kaffpa-strong).
#include <iostream>
#include <iomanip>
#include <cmath>
#include <queue>

#include "../src/configuration/HeiProMap_configuration.h"
#include "../src/datastructures/csr_graph.h"
#include "../src/datastructures/partition_manager.h"
#include "../src/distance_oracles/matrix_distance_oracle.h"
#include "../src/init_partitioning/global_multisection.h"
#include "../src/init_partitioning/flow_bisection.h"
#include "../src/init_partitioning/greedy_partitioner.h"
#include "../src/solver/HeiPa_solver.h"
#include "../src/utility/qap.h"
#include "../src/utility/utils.h"

using namespace HeiProMap;

int main(const int argc, char *argv[]) {
    AlgorithmConfiguration ac(argc, argv);

    CSRGraph g(ac.graph_in);
    MemoryPool pool;
    g.set_pool(&pool);

    PartitionManager p_manager;
    p_manager.initialize(g.n, ac.k, g.g_weight, &pool);
    DistanceOracle d_oracle;
    d_oracle.initialize(ac.hierarchy, ac.distance, &pool);

    const f64 imbalance = ac.imbalance;

    const auto sp = get_time_point();
    if (ac.partitioning_algorithm_id == PARTITIONING_ALG_MULTISECTION) {
        GlobalMultisectionPartitioner::partition(g, p_manager, ac.hierarchy, ac.distance, imbalance, ac.global_multisection_config, ac.seed);
    } else if (ac.partitioning_algorithm_id == PARTITIONING_ALG_FLOW_BISECTION) {
        FlowBisectionPartitioner::partition(g, p_manager, ac.hierarchy, imbalance, ac.flow_bisection_config, ac.seed);
    } else if (ac.partitioning_algorithm_id == PARTITIONING_ALG_GREEDY) {
        AlignedArray<vertex_t> vertices;
        vertices.initialize(g.n);
        greedy_partition(g, d_oracle, imbalance, ac.seed, p_manager, vertices);
    } else {
        std::cerr << "Unknown partitioning algorithm" << std::endl;
        return EXIT_FAILURE;
    }
    const auto ep = get_time_point();

    p_manager.reset_weights();
    for (vertex_t u = 0; u < g.n; ++u) {
        p_manager.set(u, g.uniform_v_weights ? 1 : g.v_weight(u), p_manager[u]);
    }

    const weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) ac.k));
    size_t n_empty = 0, n_oload = 0;
    weight_t sum_oload = 0;
    for (partition_t id = 0; id < ac.k; ++id) {
        n_empty += p_manager.get_bweight(id) == 0;
        n_oload += p_manager.get_bweight(id) > lmax;
        sum_oload += std::max((weight_t) 0, p_manager.get_bweight(id) - lmax);
    }

    std::cout << "Graph             : " << ac.graph_in << std::endl;
    std::cout << "#Nodes            : " << g.n << std::endl;
    std::cout << "#Edges            : " << g.m << std::endl;
    std::cout << "k                 : " << ac.k << std::endl;
    std::cout << "Partitioner       : " << ac.partitioning_algorithm_string << std::endl;
    std::cout << "Seed              : " << ac.seed << std::endl;
    std::cout << "Lmax              : " << lmax << std::endl;
    std::cout << "QAP               : " << get_qap(g, p_manager, d_oracle) << std::endl;
    std::cout << "max block w       : " << max(p_manager.get_bweights()) << std::endl;
    std::cout << "#empty partitions : " << n_empty << std::endl;
    std::cout << "#oload partitions : " << n_oload << std::endl;
    std::cout << "Sum oload weights : " << sum_oload << std::endl;
    std::cout << "Time (ms)         : " << get_milli_seconds(sp, ep) << std::endl;

    if (!ac.mapping_out.empty()) {
        std::vector<partition_t> p(g.n);
        for (vertex_t u = 0; u < g.n; ++u) { p[u] = p_manager[u]; }
        write_partition(p, ac.mapping_out);
    }
    if (std::getenv("HEIPROMAP_PROFILE") != nullptr) { HeiProMap::Profiler::instance().print_table_ascii_colored(std::cout); }
    return 0;
}
