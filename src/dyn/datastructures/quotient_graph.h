#ifndef DYN_HEIPROMAP_QUOTIENT_GRAPH_H
#define DYN_HEIPROMAP_QUOTIENT_GRAPH_H

#include <vector>
#include <map>
#include "dyn_graph.h"

namespace HeiProMap {
    class QuotientGraph {
    public:
        u64 n = 0;
        std::vector<weight_t> v_weights;
        std::vector<std::map<partition_t, weight_t> > adj;

        std::vector<partition_t> dirty_list;
        std::vector<s64> dirty_indices;
        std::vector<partition_t> dirty_list_partition;
        std::vector<s64> dirty_indices_partition;

        QuotientGraph() = default;

        QuotientGraph(const DynGraph &g, const std::vector<partition_t> &partition, u64 num_blocks)
            : n(num_blocks), v_weights(num_blocks, 0), adj(num_blocks) {
            rebuild(g, partition);
        }

        void mark_dirty(partition_t b) {
            if (b < dirty_indices.size() && dirty_indices[b] == -1) {
                dirty_indices[b] = (s64) dirty_list.size();
                dirty_list.push_back(b);
            }
            if (b < dirty_indices_partition.size() && dirty_indices_partition[b] == -1) {
                dirty_indices_partition[b] = (s64) dirty_list_partition.size();
                dirty_list_partition.push_back(b);
            }
        }

        void clear_dirty_status() {
            for (partition_t b: dirty_list) {
                if (b < dirty_indices.size()) dirty_indices[b] = -1;
            }
            dirty_list.clear();
        }

        void clear_dirty_status_partition() {
            for (partition_t b: dirty_list_partition) {
                if (b < dirty_indices_partition.size()) dirty_indices_partition[b] = -1;
            }
            dirty_list_partition.clear();
        }

        void rebuild(const DynGraph &g, const std::vector<partition_t> &partition) {
            v_weights.assign(n, 0);
            if (adj.size() != n) adj.resize(n);
            for (auto &a: adj) a.clear();

            dirty_indices.assign(n, -1);
            dirty_list.clear();
            dirty_indices_partition.assign(n, -1);
            dirty_list_partition.clear();

            for (vertex_t v = 0; v < g.n; ++v) {
                if (!g.vertex_exists(v)) continue;
                partition_t p_v = partition[v];
                v_weights[p_v] += g.v_weights[v];

                for (const auto &edge: g.neighbors[v]) {
                    if (v < edge.u) {
                        // count each edge only once
                        partition_t p_u = partition[edge.u];
                        if (p_v != p_u) {
                            adj[p_v][p_u] += edge.w;
                            adj[p_u][p_v] += edge.w;
                        }
                    }
                }
            }
        }

        void move_vertex(vertex_t v, partition_t from, partition_t to, const DynGraph &g,
                         const std::vector<partition_t> &partition) {
            if (from == to) return;
            weight_t v_weight = g.v_weights[v];
            v_weights[from] -= v_weight;
            v_weights[to] += v_weight;

            for (const auto &edge: g.neighbors[v]) {
                partition_t nb_block = partition[edge.u];
                if (nb_block == from) {
                    // Edge was internal to 'from', now it's (to, from)
                    adj[to][from] += edge.w;
                    adj[from][to] += edge.w;
                } else if (nb_block == to) {
                    // Edge was (from, to), now it's internal to 'to'
                    adj[from][to] -= edge.w;
                    if (adj[from][to] == 0) adj[from].erase(to);
                    adj[to][from] -= edge.w;
                    if (adj[to][from] == 0) adj[to].erase(from);
                } else {
                    // Edge was (from, nb_block), now it's (to, nb_block)
                    adj[from][nb_block] -= edge.w;
                    if (adj[from][nb_block] == 0) adj[from].erase(nb_block);
                    adj[nb_block][from] -= edge.w;
                    if (adj[nb_block][from] == 0) adj[nb_block].erase(from);

                    adj[to][nb_block] += edge.w;
                    adj[nb_block][to] += edge.w;
                }
            }
        }

        // Returns the parent block ID at a coarser hierarchical level
        partition_t get_parent_id(partition_t block_id, int level, const std::vector<partition_t> &hierarchy) const {
            // level 0: finest, 1: mid, 2: coarsest
            if (level == 0) return block_id / hierarchy[0];                  // Coarser is L1
            if (level == 1) return (block_id / hierarchy[0]) / hierarchy[1]; // Coarser is L2
            return 0;                                                        // Coarsest level
        }

        void contract_blocks(const std::vector<int> &block_mapping) {
            u64 new_n = 0;
            for (int b: block_mapping) if (b >= (int) new_n) new_n = b + 1;

            std::vector<weight_t> new_v_weights(new_n, 0);
            std::vector<std::map<partition_t, weight_t> > new_adj(new_n);

            for (u64 i = 0; i < n; ++i) {
                int target = block_mapping[i];
                if (target == -1) continue;
                new_v_weights[target] += v_weights[i];
                for (auto const &[neighbor, weight]: adj[i]) {
                    int neighbor_target = block_mapping[neighbor];
                    if (neighbor_target != -1) {
                        new_adj[target][neighbor_target] += weight;
                    }
                }
            }
            n = new_n;
            v_weights = new_v_weights;
            adj = new_adj;
        }
    };
}

#endif
