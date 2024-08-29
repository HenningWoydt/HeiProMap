#include <iostream>
#include <chrono>

#include "../src/definitions.h"
#include "../src/macros.h"
#include "../src/serial/utility/utils.h"
#include "../src/parallel/datastructures/parallel_solver.h"

using namespace HeiProMap;

int main(int argc, char *argv[]) {
    auto sp = std::chrono::high_resolution_clock::now();
    {
        // std::string graph_in = "../data/mapping/afshell9.graph"; std::string mapping_out = "../data/out/partition/afshell9.txt";
        // std::string graph_in = "../data/mapping/2cubes_sphere.mtx.graph"; std::string mapping_out = "../data/out/partition/2cubes_sphere.txt";
        // std::string graph_in = "../data/mapping/eur.graph"; std::string mapping_out = "../data/out/partition/eur.txt";
        std::string graph_in = "../data/mapping/rgg24.graph"; std::string mapping_out = "../data/out/partition/rgg24.txt";
        // std::string graph_in = "../data/mapping/deu.graph"; std::string mapping_out = "../data/out/partition/deu.txt";
        // std::string graph_in = "../data/mapping/PGPgiantcompo.graph"; std::string mapping_out = "../data/out/partition/PGPgiantcompo.txt";
        // std::string graph_in = "../data/test/manual_graphs/0.graph";
        std::string statistics_out = "statistics.JSON";
        std::string hierarchy_string = "4:8:6";
        std::string distance_string = "1:10:100";
        f64 imbalance = 0.03;
        u64 n_thread = 6;

        std::vector<partition_t> hierarchy = convert<partition_t>(split(hierarchy_string, ':'));
        std::vector<weight_t> distance = convert<weight_t>(split(distance_string, ':'));

        ParallelSolver parallel_solver(graph_in,
                                       hierarchy,
                                       distance,
                                       imbalance,
                                       n_thread);
        std::vector<partition_t> partition = parallel_solver.solve();
        write_partition(partition, mapping_out);
    }
    auto ep = std::chrono::high_resolution_clock::now();
    std::cout << "Total time: " << get_seconds(sp, ep) << std::endl;

    return 0;
}