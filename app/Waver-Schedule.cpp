// Schedules the halo exchange of a mapping with the timed SpMV model of wse-spmv and reports the makespan and the chain
// of messages that determines it.
//
//   Waver-Schedule --graph g.graph --topology t.graph --mapping m.txt [--costs costs.json] [--schedule-seed 0]
//                  [--chain 20]
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../src/datastructures/csr_graph.h"
#include "../src/distance_oracles/grid_distance_oracle.h"
#include "../src/scheduler/waver_scheduler.h"
#include "../src/utility/utils.h"

int main(int argc, char *argv[]) {
    using namespace HeiProMap;
    std::string graph, topology, mapping, costs;
    u64 seed = 0, chain_print = 20, top_pes = 0;
    for (int i = 1; i + 1 < argc; i += 2) {
        const std::string key = argv[i], value = argv[i + 1];
        if (key == "--graph") graph = value;
        else if (key == "--topology") topology = value;
        else if (key == "--mapping") mapping = value;
        else if (key == "--costs") costs = value;
        else if (key == "--schedule-seed") seed = std::stoull(value);
        else if (key == "--chain") chain_print = std::stoull(value);
        else if (key == "--top-pes") top_pes = std::stoull(value);
        else { std::cerr << "unknown option " << key << std::endl; return 1; }
    }
    if (graph.empty() || topology.empty() || mapping.empty()) {
        std::cerr << "usage: Waver-Schedule --graph g --topology t --mapping m [--costs costs.json] [--schedule-seed s] [--chain n]" << std::endl;
        return 1;
    }
    CSRGraph g(graph);
    GridDistanceOracle oracle;
    oracle.initialize(CSRGraph(topology), 1);
    const u32 cols = oracle.current_W(), K = oracle.current_k(), rows = K / cols;
    std::vector<partition_t> part(g.n);
    read_partition(mapping, part);
    SchedulerCosts C;
    if (!costs.empty() && !C.load(costs)) { std::cerr << "cannot read " << costs << std::endl; return 1; }

    const auto r = WaverScheduler::run(g, part, cols, rows, C, seed);
    std::cout << "makespan " << r.makespan << "\nmessages " << r.messages.size() << "\nmax engine load " << r.max_engine_load
              << "\nmax link load " << r.max_link_load << "\ncritical chain " << r.chain.size() << " messages ("
              << r.chain_engine_edges << " engine, " << r.chain_link_edges << " link dependencies)\n";
    const auto pe = [&](u32 p) { return "(" + std::to_string(p % cols) + "," + std::to_string(p / cols) + ")"; };
    for (size_t i = 0; i < r.chain.size() && i < chain_print; ++i) {
        const auto &m = r.messages[r.chain[i]];
        std::cout << "  " << std::setw(3) << i << ": " << pe(m.src) << " -> " << pe(m.dst) << "  values " << m.count << "  hops "
                  << m.hops << "  inject " << m.inject << "  finish " << m.finish << "  waited for "
                  << (m.binding >= 0 ? (m.binding_is_link ? "link" : "engine") : (m.binding < -1 ? "startup" : "nothing")) << "\n";
    }
    if (top_pes > 0) {
        // the PEs with the most engine work: what their messages are
        std::vector<u32> order(K);
        for (u32 p = 0; p < K; ++p) order[p] = p;
        std::sort(order.begin(), order.end(), [&](u32 a, u32 b) { return r.engine_load[a] > r.engine_load[b]; });
        std::vector<u32> sends(K, 0), recvs(K, 0), svals(K, 0), rvals(K, 0), far(K, 0);
        for (const auto &m: r.messages) {
            ++sends[m.src]; ++recvs[m.dst]; svals[m.src] += m.count; rvals[m.dst] += m.count;
            if (m.hops > 3) { ++far[m.src]; ++far[m.dst]; }
        }
        std::cout << "busiest engines (startup included):\n";
        for (u64 i = 0; i < top_pes && i < K; ++i) {
            const u32 p = order[i];
            std::cout << "  " << pe(p) << "  load " << r.engine_load[p] << "  sends " << sends[p] << " (" << svals[p] << " values)  receives "
                      << recvs[p] << " (" << rvals[p] << " values)  messages over >3 hops: " << far[p] << "\n";
        }
    }
    return 0;
}
