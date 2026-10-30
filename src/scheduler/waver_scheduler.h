#ifndef HEIPROMAP_WAVER_SCHEDULER_H
#define HEIPROMAP_WAVER_SCHEDULER_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "../definitions.h"
#include "../datastructures/csr_graph.h"

namespace HeiProMap {

    // Cycle constants of the timed SpMV model (the top-level numbers of wse-spmv's costs.json). The defaults are the
    // values of that file; load() reads another one.
    struct SchedulerCosts {
        s64 injection_early = 0, injection_late = 0, hop = 2;
        f64 hop_bias = 1.0;
        s64 arrival_const = 0, arrival_jitter = 0;
        f64 hop_uncertainty = 0.0;
        s64 launch = 102, matvec_setup = 77, row_empty = 20, row_even = 19, row_odd = 21, nnz_pair = 9, row_single = 1;
        s64 send_first = 196, send_lead = 196, send_setup = 9, send_assemble = 110, send_pair = 8, send_tail = 5;
        s64 send_term = 14;
        f64 wire_value = 1.18;
        s64 recv_prepare = 70, recv_fast_setup = 109, recv_fast_pair = 16, recv_fast_tail = 4;
        s64 recv_general_setup = 127, recv_general_value = 11, recv_general_pair = 14, recv_general_tail = 7;
        s64 dispatch = 64;
        s64 guard = 0; // user guard; every window also carries kEventGuard

        static constexpr s64 kEventGuard = 4;

        // Reads the top-level numeric entries ("name": number on a line with two spaces of indentation) of a costs.json.
        bool load(const std::string &path) {
            std::ifstream in(path);
            if (!in) { return false; }
            std::map<std::string, f64> v;
            std::string line;
            while (std::getline(in, line)) {
                if (line.size() < 5 || line[0] != ' ' || line[1] != ' ' || line[2] != '"') { continue; }
                const size_t q = line.find('"', 3);
                if (q == std::string::npos) { continue; }
                const size_t colon = line.find(':', q);
                if (colon == std::string::npos) { continue; }
                const std::string name = line.substr(3, q - 3);
                const std::string rest = line.substr(colon + 1);
                char *end = nullptr;
                const f64 value = std::strtod(rest.c_str(), &end);
                if (end != rest.c_str()) { v[name] = value; }
            }
            const auto I = [&](const char *k, s64 &field) { if (v.count(k)) { field = static_cast<s64>(std::llround(v[k])); } };
            const auto F = [&](const char *k, f64 &field) { if (v.count(k)) { field = v[k]; } };
            I("injection_early", injection_early); I("injection_late", injection_late); I("hop", hop);
            F("hop_bias", hop_bias); I("arrival_const", arrival_const); I("arrival_jitter", arrival_jitter);
            F("hop_uncertainty", hop_uncertainty);
            I("launch", launch); I("matvec_setup", matvec_setup); I("row_empty", row_empty); I("row_even", row_even);
            I("row_odd", row_odd); I("nnz_pair", nnz_pair); I("row_single", row_single);
            I("send_first", send_first); I("send_lead", send_lead); I("send_setup", send_setup);
            I("send_assemble", send_assemble); I("send_pair", send_pair); I("send_tail", send_tail);
            I("send_term", send_term); F("wire_value", wire_value);
            I("recv_prepare", recv_prepare); I("recv_fast_setup", recv_fast_setup); I("recv_fast_pair", recv_fast_pair);
            I("recv_fast_tail", recv_fast_tail); I("recv_general_setup", recv_general_setup);
            I("recv_general_value", recv_general_value); I("recv_general_pair", recv_general_pair);
            I("recv_general_tail", recv_general_tail); I("dispatch", dispatch);
            return true;
        }

        s64 event_guard() const { return std::max(guard, kEventGuard); }
    };

    // A greedy list scheduler for the halo exchange of one SpMV on a cols x rows PE grid, following the timed model of
    // wse-spmv. Every ordered pair of PEs that exchange values is one coalesced message. A message holds, relative to its
    // injection time r:
    //   - the engine of the sender, from before r (arming and filling the send buffer) until the data is on the wire,
    //   - every directed link of its XY route (along the source row first, then down the destination column) for the
    //     time the data takes to pass,
    //   - the engine of the receiver, from before the header arrives until the values are scattered.
    // Each PE has one engine, which first runs the startup (local SpMV). The messages are visited in a random order and
    // each is injected at the earliest time at which all its windows are free; the makespan is the time the last message
    // is received. For every message the reservation that pushed it to its injection time is recorded, so the chain of
    // messages that determines the makespan can be followed back.
    class WaverScheduler {
    public:
        struct Message {
            u32 src = 0, dst = 0;
            u32 count = 0;     // values in the message
            u32 work = 0;      // matrix entries the receiver updates with them
            u32 odd = 0;       // values with an odd number of such entries
            u32 hops = 0;
            s64 inject = 0, finish = 0;
            s32 binding = -1;  // the message whose reservation pushed this one to `inject`, -1 none, -2 - pe: startup of pe
            bool binding_is_link = false;
        };

        struct Result {
            s64 makespan = 0;
            std::vector<Message> messages;
            std::vector<s64> startup_end;   // per PE: end of the startup (local SpMV) reservation
            std::vector<s64> engine_load;   // per PE: startup plus all message windows on its engine
            std::vector<s64> link_load;     // per directed link (PE * 4 + direction)
            std::vector<u32> chain;         // messages of the critical chain, last one first
            u32 last_message = 0;
            u32 chain_link_edges = 0, chain_engine_edges = 0;
            s64 max_engine_load = 0, max_link_load = 0;
        };

        static constexpr u32 kEast = 0, kWest = 1, kSouth = 2, kNorth = 3;

        static Result run(const CSRGraph &g, const std::vector<partition_t> &part, const u32 cols, const u32 rows,
                          const SchedulerCosts &C, const u64 seed) {
            const u32 K = cols * rows;
            Result res;
            std::vector<Message> &msgs = res.messages;

            // ---- messages and startup features from the graph and the mapping ----
            std::vector<s64> row_even(K, 0), row_odd(K, 0), pairs(K, 0), single(K, 0), nrows(K, 0), empty(K, 0);
            struct Entry {
                u64 key;
                u32 nb;
            };
            std::vector<Entry> entries;
            std::vector<partition_t> others;
            for (vertex_t u = 0; u < g.n; ++u) {
                const partition_t a = part[u];
                ++nrows[a];
                others.clear();
                u32 same = 0;
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const vertex_t v = g.edges_v[j];
                    if (part[v] == a) { ++same; } else { others.push_back(part[v]); }
                }
                const u32 len = same + 1; // the diagonal entry
                row_even[a] += len % 2 == 0;
                row_odd[a] += len % 2;
                pairs[a] += len / 2;
                single[a] += len == 1;
                std::sort(others.begin(), others.end());
                for (size_t i = 0; i < others.size();) {
                    size_t j = i;
                    while (j < others.size() && others[j] == others[i]) { ++j; }
                    entries.push_back({static_cast<u64>(a) * K + others[i], static_cast<u32>(j - i)});
                    i = j;
                }
            }
            std::sort(entries.begin(), entries.end(), [](const Entry &x, const Entry &y) { return x.key < y.key; });
            for (size_t i = 0; i < entries.size();) {
                size_t j = i;
                Message m;
                m.src = static_cast<u32>(entries[i].key / K);
                m.dst = static_cast<u32>(entries[i].key % K);
                while (j < entries.size() && entries[j].key == entries[i].key) {
                    ++m.count;
                    m.work += entries[j].nb;
                    m.odd += entries[j].nb % 2;
                    ++j;
                }
                m.hops = static_cast<u32>(std::abs(static_cast<int>(m.src % cols) - static_cast<int>(m.dst % cols))
                                          + std::abs(static_cast<int>(m.src / cols) - static_cast<int>(m.dst / cols)));
                msgs.push_back(m);
                i = j;
            }

            // ---- reservations ----
            struct Reservation {
                s64 s, e;
                s32 owner; // message index, or -2 - pe for the startup
            };
            const u32 num_links = 4 * K;
            std::vector<std::vector<Reservation> > table(num_links + K);
            res.startup_end.assign(K, 0);
            res.engine_load.assign(K, 0);
            res.link_load.assign(num_links, 0);
            for (u32 p = 0; p < K; ++p) {
                const s64 ready = std::max<s64>(0, C.launch + C.matvec_setup + C.row_even * row_even[p] + C.row_odd * row_odd[p]
                                                       + C.nnz_pair * pairs[p] + C.row_single * single[p]
                                                       + C.row_empty * empty[p]);
                const s64 end = ready + C.send_first + C.event_guard();
                res.startup_end[p] = end;
                res.engine_load[p] = end;
                table[num_links + p].push_back({0, end, -2 - static_cast<s32>(p)});
            }

            struct Window {
                u32 resource;
                s64 s, e; // relative to the injection time
            };
            const auto link_of = [&](const u32 from_pe, const u32 dir) { return from_pe * 4 + dir; };
            std::vector<Window> windows;
            const auto build = [&](const Message &m) {
                windows.clear();
                const s64 wire = std::llround(C.wire_value * static_cast<f64>(m.count));
                const s64 wire_span = wire + C.send_term;
                const s64 fill = C.send_assemble + C.send_pair * (m.count / 2) + C.send_tail * (m.count % 2);
                const s64 blast = C.send_setup + wire;
                s64 scatter;
                if (m.work == m.count) {
                    scatter = C.recv_fast_setup + C.recv_fast_pair * (m.count / 2) + C.recv_fast_tail * (m.count % 2);
                } else {
                    scatter = C.recv_general_setup + C.recv_general_value * m.count
                              + C.recv_general_pair * ((m.work - m.odd) / 2) + C.recv_general_tail * m.odd;
                }
                const s64 recv = wire_span + scatter;
                const auto arrival = [&](const u32 k) -> s64 {
                    return k == 0 ? 0 : C.hop * k + std::llround(C.hop_bias * static_cast<f64>(k));
                };
                // sender engine: arming and filling before r, then the blast and the dispatch
                windows.push_back({num_links + m.src, -std::max<s64>(C.send_lead + fill, C.injection_early),
                                   blast + C.injection_late + C.dispatch + C.event_guard()});
                // links of the XY route
                const s64 duration = wire_span + 2 + C.guard;
                u32 pe = m.src, k = 0;
                const u32 dc = m.dst % cols, dr = m.dst / cols;
                while (pe % cols != dc) {
                    const bool east = dc > pe % cols;
                    windows.push_back({link_of(pe, east ? kEast : kWest), arrival(k), arrival(k) + duration});
                    pe = east ? pe + 1 : pe - 1;
                    ++k;
                }
                while (pe / cols != dr) {
                    const bool south = dr > pe / cols;
                    windows.push_back({link_of(pe, south ? kSouth : kNorth), arrival(k), arrival(k) + duration});
                    pe = south ? pe + cols : pe - cols;
                    ++k;
                }
                // receiver engine: from before the earliest header until the values are scattered
                windows.push_back({num_links + m.dst, arrival(k) - C.recv_prepare,
                                   arrival(k) + C.arrival_const + recv + C.dispatch + C.event_guard()});
            };

            // the first reservation of the resource that overlaps [s, e), or -1
            const auto conflict = [&](const std::vector<Reservation> &v, const s64 s, const s64 e) -> s64 {
                size_t lo = 0, hi = v.size();
                while (lo < hi) {
                    const size_t mid = (lo + hi) / 2;
                    if (v[mid].e > s) { hi = mid; } else { lo = mid + 1; }
                }
                return lo < v.size() && v[lo].s < e ? static_cast<s64>(lo) : -1;
            };

            std::vector<u32> order(msgs.size());
            for (u32 i = 0; i < order.size(); ++i) { order[i] = i; }
            std::mt19937_64 rng(seed * 0x9E3779B97F4A7C15ULL + 12345);
            std::shuffle(order.begin(), order.end(), rng);

            res.makespan = 0;
            for (const u32 id: order) {
                Message &m = msgs[id];
                build(m);
                s64 r = 0;
                s32 binding = -1;
                bool binding_link = false;
                for (bool moved = true; moved;) {
                    moved = false;
                    for (const Window &w: windows) {
                        const std::vector<Reservation> &v = table[w.resource];
                        const s64 c = conflict(v, r + w.s, r + w.e);
                        if (c < 0) { continue; }
                        r = v[static_cast<size_t>(c)].e - w.s;
                        binding = v[static_cast<size_t>(c)].owner;
                        binding_link = w.resource < num_links;
                        moved = true;
                        break;
                    }
                }
                m.inject = r;
                m.binding = binding;
                m.binding_is_link = binding_link;
                m.finish = r + windows.back().e;
                for (const Window &w: windows) {
                    std::vector<Reservation> &v = table[w.resource];
                    const s64 s = r + w.s, e = r + w.e;
                    size_t lo = 0, hi = v.size();
                    while (lo < hi) {
                        const size_t mid = (lo + hi) / 2;
                        if (v[mid].s >= s) { hi = mid; } else { lo = mid + 1; }
                    }
                    v.insert(v.begin() + static_cast<std::ptrdiff_t>(lo), {s, e, static_cast<s32>(id)});
                    if (w.resource < num_links) {
                        res.link_load[w.resource] += e - s;
                    } else {
                        res.engine_load[w.resource - num_links] += e - s;
                    }
                }
                res.makespan = std::max(res.makespan, m.finish);
            }
            for (u32 p = 0; p < K; ++p) { res.makespan = std::max(res.makespan, res.startup_end[p]); }
            for (const s64 l: res.engine_load) { res.max_engine_load = std::max(res.max_engine_load, l); }
            for (const s64 l: res.link_load) { res.max_link_load = std::max(res.max_link_load, l); }

            // the critical chain: from the message that finishes last, follow what pushed each one
            if (!msgs.empty()) {
                u32 last = 0;
                for (u32 i = 1; i < msgs.size(); ++i) { if (msgs[i].finish > msgs[last].finish) { last = i; } }
                res.last_message = last;
                for (s32 cur = static_cast<s32>(last); cur >= 0; cur = msgs[static_cast<size_t>(cur)].binding) {
                    res.chain.push_back(static_cast<u32>(cur));
                    const Message &m = msgs[static_cast<size_t>(cur)];
                    if (m.binding >= 0) { (m.binding_is_link ? res.chain_link_edges : res.chain_engine_edges)++; }
                }
            }
            return res;
        }
    };

}

#endif //HEIPROMAP_WAVER_SCHEDULER_H
