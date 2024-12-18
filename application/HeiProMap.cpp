#include <iostream>

#include "../src/definitions.h"
#include "../src/macros.h"
#include "../src/serial/datastructures/solver.h"
#include "../src/serial/utility/utils.h"

using namespace HeiProMap;

int main() {
    const auto sp = std::chrono::high_resolution_clock::now();
    {
        // std::string graph_in = "../data/mapping/afshell9.graph"; std::string mapping_out = "../data/out/partition/afshell9.txt";
        std::string graph_in = "../data/mapping/2cubes_sphere.mtx.graph"; std::string mapping_out = "../data/out/partition/2cubes_sphere.txt";
        // const std::string graph_in = "../data/mapping/eur.graph"; const std::string mapping_out = "../data/out/partition/eur.txt";
        // const std::string graph_in = "../data/mapping/eur.graph"; const std::string mapping_out = "../data/out/partition/eur.txt";
        // std::string graph_in = "../data/mapping/rgg24.graph"; std::string mapping_out = "../data/out/partition/rgg24.txt";
        // std::string graph_in = "../data/mapping/rgg_n26.graph"; std::string mapping_out = "../data/out/partition/rgg_n26.txt";
        // std::string graph_in = "../data/mapping/deu.graph"; std::string mapping_out = "../data/out/partition/deu.txt";
        // std::string graph_in = "../data/mapping/PGPgiantcompo.graph"; std::string mapping_out = "../data/out/partition/PGPgiantcompo.txt";
        // std::string graph_in = "../data/test/manual_graphs/0.graph";
        const std::string statistics_out   = "statistics.JSON";
        const std::string hierarchy_string = "4:8:6";
        const std::string distance_string  = "1:10:100";
        constexpr f64 imbalance            = 0.03;

        const std::vector<partition_t> hierarchy = convert<partition_t>(split(hierarchy_string, ':'));
        const std::vector<weight_t> distance = convert<weight_t>(split(distance_string, ':'));

        Solver solver(graph_in,
                      hierarchy,
                      distance,
                      imbalance);
        const std::vector<partition_t> partition = solver.solve();
        write_partition(partition, mapping_out);
    }
    const auto ep = std::chrono::high_resolution_clock::now();
    std::cout << "Total time: " << get_seconds(sp, ep) << std::endl;

    return 0;
}
