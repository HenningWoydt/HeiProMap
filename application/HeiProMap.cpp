#include <iostream>

#include "../src/utility/definitions.h"
#include "../src/utility/macros.h"
#include "../src/utility/utils.h"
#include "../src/datastructures/solver.h"
#include "../src/parallel/datastructures/parallel_refinement_solver.h"

using namespace SPM;

int main(int argc, char *argv[]) {
    auto sp = std::chrono::high_resolution_clock::now();
    {
        std::string graph_in = "../data/mapping/afshell9.graph"; std::string mapping_out = "afshell9.txt";
        // std::string graph_in = "../data/mapping/2cubes_sphere.mtx.graph"; std::string mapping_out = "2cubes_sphere.txt";
        // std::string graph_in = "../data/mapping/eur.graph"; std::string mapping_out = "eur.txt";
        // std::string graph_in = "../data/mapping/rgg24.graph"; std::string mapping_out = "rgg24.txt";
        // std::string graph_in = "../data/mapping/deu.graph"; std::string mapping_out = "deu.txt";
        // std::string graph_in = "../data/mapping/PGPgiantcompo.graph"; std::string mapping_out = "PGPgiantcompo.txt";
        // std::string graph_in = "../data/test/manual_graphs/0.graph";
        std::string statistics_out = "statistics.JSON";
        std::string hierarchy_string = "4:8:6";
        std::string distance_string = "1:10:100";
        f64 imbalance = 0.03;

        std::vector<u64> hierarchy = convert<u64>(split(hierarchy_string, ':'));
        std::vector<u64> distance = convert<u64>(split(distance_string, ':'));
        u64 k = prod<u64>(hierarchy);

        Solver solver(graph_in,
                      hierarchy,
                      distance,
                      k,
                      imbalance);
        std::vector<partition_t> partition = solver.solve();
        write_partition(partition, mapping_out);
    }
    auto ep = std::chrono::high_resolution_clock::now();
    std::cout << "Total time: " << get_seconds(sp, ep) << std::endl;

    return 0;
}
