// Schedules the halo exchange of a mapping with the timed SpMV model of wse-spmv and reports the makespan and the chain
// of messages that determines it.
//
//   Waver-Schedule --graph g.graph --topology t.graph --mapping m.txt [--costs costs.json] [--schedule-seed 0]
//                  [--chain 20]
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../src/datastructures/csr_graph.h"
#include "../src/distance_oracles/grid_distance_oracle.h"
#include "../src/refinement/waver_corner_refinement.h"
#include "../src/scheduler/waver_scheduler.h"
#include "../src/utility/utils.h"

int main(int argc, char *argv[]) {
    using namespace HeiProMap;
    std::string graph, topology, mapping, costs;
    u64 seed = 0, chain_print = 20, top_pes = 0, drop_hops = 0, drop_values = 0, seeds = 1, remap_rounds = 0, crit_orders = 0, remap_iters = 1;
    f64 imbalance = 0.03, power = 6.0, tau = 2.0, crit_lambda = 4.0;
    bool pair_objective = false;
    std::string remap_out, order_name = "random", analysis_out;
    for (int i = 1; i + 1 < argc; i += 2) {
        const std::string key = argv[i], value = argv[i + 1];
        if (key == "--graph") graph = value;
        else if (key == "--topology") topology = value;
        else if (key == "--mapping") mapping = value;
        else if (key == "--costs") costs = value;
        else if (key == "--schedule-seed") seed = std::stoull(value);
        else if (key == "--chain") chain_print = std::stoull(value);
        else if (key == "--remap-pair") pair_objective = value != "0";
        else if (key == "--analysis-out") analysis_out = value;
        else if (key == "--crit-orders") crit_orders = std::stoull(value);
        else if (key == "--crit-lambda") crit_lambda = std::stod(value);
        else if (key == "--remap-iters") remap_iters = std::stoull(value);
        else if (key == "--order") order_name = value;
        else if (key == "--seeds") seeds = std::stoull(value);
        else if (key == "--remap-rounds") remap_rounds = std::stoull(value);
        else if (key == "--remap-out") remap_out = value;
        else if (key == "--remap-power") power = std::stod(value);
        else if (key == "--remap-tau") tau = std::stod(value);
        else if (key == "--imbalance") imbalance = std::stod(value);
        else if (key == "--drop-hops") drop_hops = std::stoull(value);
        else if (key == "--drop-values") drop_values = std::stoull(value);
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

    const WaverScheduler::Order order_mode = order_name == "rounds" ? WaverScheduler::Order::Rounds
                                           : order_name == "longest" ? WaverScheduler::Order::Longest : WaverScheduler::Order::Random;
    const auto median_makespan = [&](const std::vector<partition_t> &p) {
        std::vector<s64> v;
        for (u64 s = 0; s < seeds; ++s) { v.push_back(WaverScheduler::run(g, p, cols, rows, C, seed + s, 0, 0, order_mode).makespan); }
        std::sort(v.begin(), v.end());
        return v[v.size() / 2];
    };
    if (remap_rounds > 0) {
        WaverCornerRefinementConfiguration rc;
        rc.enabled = true;
        rc.max_rounds = remap_rounds;
        rc.power = power;
        rc.tau = tau;
        rc.pair_objective = pair_objective;
        rc.verbose = false;
        std::cout << "makespan before (median of " << seeds << "): " << median_makespan(part) << std::endl;
        const weight_t lmax = static_cast<weight_t>(std::ceil((1.0 + imbalance) * static_cast<f64>(g.n) / K));
        u64 moves = 0;
        std::vector<f64> crit_weights;
        for (u64 it = 0; it < remap_iters; ++it) {
            if (crit_orders > 0) {
                // how often each PE is on the critical chain, over orders other than the ones the result is measured with
                std::vector<f64> hits(K, 0.0);
                for (u64 o = 0; o < crit_orders; ++o) {
                    const auto sr = WaverScheduler::run(g, part, cols, rows, C, 1000 + o, 0, 0, order_mode);
                    std::vector<char> seen(K, 0);
                    for (const u32 id: sr.chain) { seen[sr.messages[id].src] = 1; seen[sr.messages[id].dst] = 1; }
                    for (u32 p = 0; p < K; ++p) { hits[p] += seen[p]; }
                }
                crit_weights.assign(K, 1.0);
                for (u32 p = 0; p < K; ++p) { crit_weights[p] = 1.0 + crit_lambda * hits[p] / static_cast<f64>(crit_orders); }
                rc.weights = &crit_weights;
                u32 critical = 0;
                for (u32 p = 0; p < K; ++p) { critical += hits[p] > 0; }
                std::cout << "iteration " << it << ": " << critical << " PEs were on a critical chain at least once" << std::endl;
            }
            moves += WaverCornerRefinement::refine(g, part, K, lmax, rc);
        }
        std::cout << "corner refinement: " << moves << " moves\nmakespan after (median of " << seeds << "): "
                  << median_makespan(part) << std::endl;
        if (!remap_out.empty()) { write_partition(part, remap_out); }
    }
    if (!analysis_out.empty()) {
        // per PE: engine load and how often it is on the critical chain over crit_orders random orders; and the messages
        const auto base = WaverScheduler::run(g, part, cols, rows, C, seed, 0, 0, order_mode);
        std::vector<f64> hits(K, 0.0), chain_msgs(K, 0.0);
        for (u64 o = 0; o < crit_orders; ++o) {
            const auto sr = WaverScheduler::run(g, part, cols, rows, C, 1000 + o, 0, 0, order_mode);
            std::vector<char> seen(K, 0);
            for (const u32 id: sr.chain) {
                seen[sr.messages[id].src] = 1; seen[sr.messages[id].dst] = 1;
                chain_msgs[sr.messages[id].src] += 1; chain_msgs[sr.messages[id].dst] += 1;
            }
            for (u32 p = 0; p < K; ++p) { hits[p] += seen[p]; }
        }
        std::ofstream pe_file(analysis_out + ".pe"), msg_file(analysis_out + ".msg");
        for (u32 p = 0; p < K; ++p) { pe_file << p << " " << base.engine_load[p] << " " << hits[p] << " " << chain_msgs[p] << "\n"; }
        for (const auto &m: base.messages) { msg_file << m.src << " " << m.dst << " " << m.count << "\n"; }
    }
    const auto r = WaverScheduler::run(g, part, cols, rows, C, seed, static_cast<u32>(drop_hops), static_cast<u32>(drop_values), order_mode);
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
