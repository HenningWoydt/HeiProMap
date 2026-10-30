#ifndef HEIPROMAP_WAVER_REGION_LOAD_H
#define HEIPROMAP_WAVER_REGION_LOAD_H

#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../definitions.h"
#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"

namespace HeiProMap {
    // Per-block counts that decide the engine load of a PE, kept up to date while vertices move.
    //
    // A separator vertex of a block x (a vertex with an edge to another block) has an involved set: the foreign blocks
    // its neighbors are in. A send region of x is a maximal group of separator vertices of x with the same involved set,
    // and it is shipped to every block of that set as one message region. So the load of a block depends on
    //   - its neighbor blocks (one coalesced message to and from each),
    //   - its send regions (SR),
    //   - the regions of other blocks that it receives (RR: regions whose involved set contains it),
    //   - the values it sends (SV: its separator vertices) and receives (RV: separator vertices of other blocks with an
    //     edge to it).
    // The cost per unit of each is set by the Weights (the defaults are read off the cost model of the scheduler).
    //
    // For every vertex the number of its edges into each block is kept, from which the involved set follows. Moving a
    // vertex changes its own set and the sets of its neighbors, so only these are re-registered.
    class WaverRegionLoad {
        using graph_t = CSRGraph;
        using p_manager_t = PartitionManager;

    public:
        // The engine cost model of the scheduler's timed model, fitted to costs.json: a message costs the sender
        // send_message and the receiver receive_message, every value in it send_value / receive_value on the engine
        // that sends / receives it (a value sent to two PEs is in two messages), and every local nonzero local_nnz.
        // The terms that do not depend on the mapping (startup, rows) are left out. The region terms are 0: they
        // are not part of the fitted model (they can be set to count the regions of the messages).
        struct Weights {
            f64 message = 655.0;        // per neighbor block: send (383) and receive (272) of one coalesced message
            f64 send_region = 0.0;      // per send region
            f64 receive_region = 0.0;   // per region received
            f64 send_value = 5.1;       // per value sent (once per receiving PE)
            f64 receive_value = 15.3;   // per value received
            f64 local_nnz = 4.5;        // per matrix entry whose row and column are on the PE
        };

        void initialize(const graph_t &t_g, const p_manager_t &t_pm, const partition_t t_k, const Weights &t_w) {
            g = &t_g;
            k = t_k;
            w = t_w;
            part.assign(g->n, 0);
            for (vertex_t u = 0; u < g->n; ++u) { part[u] = t_pm[u]; }
            internal.assign(k, 0);
            sr.assign(k, 0);
            rr.assign(k, 0);
            sv.assign(k, 0);
            rv.assign(k, 0);
            regions.assign(k, {});
            stamp.assign(k, 0);
            epoch = 1;
            dirty.clear();
            zobrist.resize(k);
            u64 s = 0x9E3779B97F4A7C15ULL;
            for (partition_t b = 0; b < k; ++b) {
                s += 0x9E3779B97F4A7C15ULL;
                u64 z = s;
                z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
                z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
                zobrist[b] = z ^ (z >> 31);
            }
            conn.assign(g->n, {});
            for (vertex_t u = 0; u < g->n; ++u) {
                for (size_t j = g->neighborhoods[u]; j < g->neighborhoods[u + 1]; ++j) { add_edge(u, part[g->edges_v[j]], +1); }
            }
            for (vertex_t u = 0; u < g->n; ++u) {
                for (size_t j = g->neighborhoods[u]; j < g->neighborhoods[u + 1]; ++j) {
                    if (part[g->edges_v[j]] == part[u]) { ++internal[part[u]]; }
                }
            }
            for (vertex_t u = 0; u < g->n; ++u) { reg(u, +1); }
            dirty.clear();
            ++epoch;
        }

        // The load of block x, given the number of its neighbor blocks.
        f64 load(const partition_t x, const size_t degree) const {
            return w.local_nnz * static_cast<f64>(internal[x]) + w.message * static_cast<f64>(degree) + w.send_region * static_cast<f64>(sr[x]) +
                   w.receive_region * static_cast<f64>(rr[x]) + w.send_value * static_cast<f64>(sv[x]) +
                   w.receive_value * static_cast<f64>(rv[x]);
        }

        // Moves vertex u from block `from` to block `to`.
        void move(const vertex_t u, const partition_t from, const partition_t to) {
            const size_t begin = g->neighborhoods[u], end = g->neighborhoods[u + 1];
            // the entries between u and the vertices of `from` stop being local, those with the vertices of `to` start
            // to be (an edge is two entries)
            for (const auto &[block, count]: conn[u]) {
                if (block == from) { internal[from] -= 2 * static_cast<int64_t>(count); }
                if (block == to) { internal[to] += 2 * static_cast<int64_t>(count); }
            }
            touch(from);
            touch(to);
            reg(u, -1);
            for (size_t j = begin; j < end; ++j) { reg(g->edges_v[j], -1); }
            for (size_t j = begin; j < end; ++j) {
                const vertex_t v = g->edges_v[j];
                add_edge(v, from, -1);
                add_edge(v, to, +1);
            }
            part[u] = to;
            reg(u, +1);
            for (size_t j = begin; j < end; ++j) { reg(g->edges_v[j], +1); }
        }

        int64_t send_regions(const partition_t x) const { return sr[x]; }
        int64_t receive_regions(const partition_t x) const { return rr[x]; }
        int64_t send_values(const partition_t x) const { return sv[x]; }
        int64_t receive_values(const partition_t x) const { return rv[x]; }

        // The blocks whose counters changed since the last clear_dirty().
        const std::vector<partition_t> &dirty_blocks() const { return dirty; }

        void clear_dirty() {
            dirty.clear();
            ++epoch;
        }

    private:
        const graph_t *g = nullptr;
        partition_t k = 0;
        Weights w;
        std::vector<partition_t> part;
        std::vector<std::vector<std::pair<partition_t, uint32_t> > > conn; // vertex -> (block, number of edges into it)
        struct Region {
            uint32_t count = 0;
        };
        std::vector<std::unordered_map<u64, Region> > regions; // block -> involved set (hash) -> separator vertices
        std::vector<int64_t> internal, sr, rr, sv, rv; // internal: local matrix entries of the block
        std::vector<u64> zobrist;
        std::vector<u32> stamp;
        u32 epoch = 1;
        std::vector<partition_t> dirty;

        void touch(const partition_t x) {
            if (stamp[x] != epoch) {
                stamp[x] = epoch;
                dirty.push_back(x);
            }
        }

        void add_edge(const vertex_t v, const partition_t block, const int delta) {
            auto &list = conn[v];
            for (size_t i = 0; i < list.size(); ++i) {
                if (list[i].first == block) {
                    list[i].second = static_cast<uint32_t>(static_cast<int64_t>(list[i].second) + delta);
                    if (list[i].second == 0) {
                        list[i] = list.back();
                        list.pop_back();
                    }
                    return;
                }
            }
            list.emplace_back(block, static_cast<uint32_t>(delta));
        }

        // Registers (+1) or unregisters (-1) vertex v in the counters of its block and the blocks of its involved set.
        void reg(const vertex_t v, const int sign) {
            const partition_t b = part[v];
            u64 key = 0;
            size_t foreign = 0;
            for (const auto &[block, count]: conn[v]) {
                if (block != b) {
                    key ^= zobrist[block];
                    ++foreign;
                }
            }
            if (foreign == 0) { return; }
            auto &map = regions[b];
            Region &r = map[key];
            const bool region_changes = sign > 0 ? r.count == 0 : r.count == 1;
            if (sign > 0) {
                ++r.count;
            } else {
                --r.count;
            }
            if (r.count == 0) { map.erase(key); }
            touch(b);
            sv[b] += sign * static_cast<int64_t>(foreign); // the value is in one message per receiving PE
            if (region_changes) { sr[b] += sign; }
            for (const auto &[block, count]: conn[v]) {
                if (block == b) { continue; }
                touch(block);
                rv[block] += sign;
                if (region_changes) { rr[block] += sign; }
            }
        }
    };
} // namespace HeiProMap

#endif //HEIPROMAP_WAVER_REGION_LOAD_H
