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

#include <iostream>

#include "../src/definitions.h"
#include "../src/macros.h"
#include "../src/serial/datastructures/solver.h"
#include "../src/serial/utility/utils.h"

using namespace HeiProMap;

int main(const int argc, const char* argv[]) {
    const auto sp = std::chrono::high_resolution_clock::now();
    {
        // std::string graph_in = "../data/mapping/afshell9.graph"; std::string mapping_out = "../data/out/partition/afshell9.txt";
        std::string graph_in    = "../data/mapping/2cubes_sphere.mtx.graph";
        std::string mapping_out = "../data/out/partition/2cubes_sphere.txt"; // 7.135.366
        // const std::string graph_in = "../data/mapping/eur.graph"; const std::string mapping_out = "../data/out/partition/eur.txt";
        // const std::string graph_in = "../data/mapping/as-skitter.graph"; const std::string mapping_out = "../data/out/partition/as-skitter.txt"; // 13.801.154
        // const std::string graph_in = "../data/mapping/wiki-Talk.graph"; const std::string mapping_out = "../data/out/partition/wiki-Talk.txt"; // 80.794.192
        // std::string graph_in = "../data/mapping/rgg24.graph"; std::string mapping_out = "../data/out/partition/rgg24.txt";
        // std::string graph_in = "../data/mapping/rgg_n26.graph"; std::string mapping_out = "../data/out/partition/rgg_n26.txt";
        // std::string graph_in = "../data/mapping/deu.graph"; std::string mapping_out = "../data/out/partition/deu.txt"; // 177.128
        // std::string graph_in = "../data/mapping/PGPgiantcompo.graph"; std::string mapping_out = "../data/out/partition/PGPgiantcompo.txt";
        // std::string graph_in = "../data/test/manual_graphs/0.graph";
        const std::string statistics_out   = "statistics.JSON";
        const std::string hierarchy_string = "4:8:6";
        const std::string distance_string  = "1:10:100";
        constexpr f64 imbalance            = 0.03;

        const std::vector<partition_t> hierarchy = convert<partition_t>(split(hierarchy_string, ':'));
        const std::vector<weight_t> distance     = convert<weight_t>(split(distance_string, ':'));

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
