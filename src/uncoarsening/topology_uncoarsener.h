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

#ifndef HEIPROMAP_TOPOLOGY_UNCOARSENER_H
#define HEIPROMAP_TOPOLOGY_UNCOARSENER_H

#include <algorithm>
#include <limits>
#include <numeric>
#include <vector>

#include "../datastructures/csr_graph.h"
#include "../datastructures/partition_manager.h"
#include "../definitions.h"
#include "../distance_oracles/grid_distance_oracle.h"
#include "../rebalance/empty_block_diffusion.h"
#include "../utility/grid_utils.h"
#include "../utility/mapping.h"

namespace HeiProMap {
    // Expands a partition from a coarse topology level to the next finer one. Every graph vertex of a coarse
    // block is first assigned to the cheapest child cell, then the result is repaired in several passes:
    //   1. steal one vertex for each empty cell from its siblings,
    //   2. avalanche: move whole sender groups along the density gradient into the remaining empty cells,
    //   3. evacuate small islands of non-empty cells into the largest connected region,
    //   4. try all assignments of the children's vertex groups to the child cells (only for 2-4 children),
    //   5. diffuse from dense borders into empty regions.
    // The caller has already uncoarsened the oracle and re-initialized the partition manager for the finer k.
    class TopologyUncoarsener {
    public:
        struct CostParams {
            f64 alpha_end;
            f64 beta_end;
            f64 alpha_link;
            f64 beta_link;
        };

        // coarse_to_fine[c]: the child cells of coarse cell c. Must be built before the oracle is uncoarsened.
        static std::vector<std::vector<partition_t> > children(const GridDistanceOracle &oracle) {
            const partition_t coarse_k = oracle.current_k();
            const partition_t fine_k = oracle.fine_k();
            const Mapping &topo_map = oracle.top_mapping();
            std::vector<std::vector<partition_t> > coarse_to_fine(coarse_k);
            for (partition_t f = 0; f < fine_k; ++f) {
                coarse_to_fine[topo_map.get(f)].push_back(f);
            }
            return coarse_to_fine;
        }

        TopologyUncoarsener(const CSRGraph &t_g,
                            PartitionManager &t_p_manager,
                            const GridDistanceOracle &t_oracle,
                            const CostParams &t_costs,
                            const std::vector<std::vector<partition_t> > &t_coarse_to_fine)
            : g(t_g), p_manager(t_p_manager), oracle(t_oracle), costs(t_costs),
              coarse_to_fine(t_coarse_to_fine),
              coarse_k(static_cast<partition_t>(t_coarse_to_fine.size())),
              fine_k(t_oracle.current_k()) {}

        // old_partition[u]: coarse cell of vertex u.
        void run(const std::vector<partition_t> &old_partition) {
            expand(old_partition);
            steal_for_empty_cells();
            avalanche_into_empty_cells();
            evacuate_islands();
            optimize_sibling_permutations();
            EmptyBlockDiffusion::diffuse(g, p_manager, static_cast<size_t>(oracle.current_W()));
        }

    private:
        const CSRGraph &g;
        PartitionManager &p_manager;
        const GridDistanceOracle &oracle;
        CostParams costs;
        const std::vector<std::vector<partition_t> > &coarse_to_fine;
        partition_t coarse_k;
        partition_t fine_k;

        // blk_verts[b]: vertices assigned to fine cell b (may contain stale entries; check p_manager.partition).
        std::vector<std::vector<vertex_t> > blk_verts;

        // Assign each vertex of a coarse block to its cheapest child cell (capacity permitting).
        void expand(const std::vector<partition_t> &old_partition) {
            // mark all unassigned so cost computation can skip them
            for (vertex_t u = 0; u < g.n; ++u) {
                p_manager.partition[u] = NO_ID;
            }

            blk_verts.assign(fine_k, {});

            for (vertex_t u = 0; u < g.n; ++u) {
                const auto &fine_blocks = coarse_to_fine[old_partition[u]];
                const weight_t u_weight = g.v_weight(u);

                // collect distinct assigned neighbor fine blocks (for alpha_end)
                partition_t neighbor_blocks[256];
                size_t n_neighbor_blocks = 0;
                for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                    const partition_t v_block = p_manager.partition[g.edges_v[j]];
                    if (v_block == NO_ID) { continue; }
                    bool found = false;
                    for (size_t x = 0; x < n_neighbor_blocks; ++x) {
                        if (neighbor_blocks[x] == v_block) {
                            found = true;
                            break;
                        }
                    }
                    if (!found && n_neighbor_blocks < 256) { neighbor_blocks[n_neighbor_blocks++] = v_block; }
                }

                partition_t best = NO_ID;
                f64 min_cost = std::numeric_limits<f64>::max();

                for (const partition_t block: fine_blocks) {
                    if (p_manager.bweights[block] + u_weight > p_manager.get_lmax(block)) { continue; }

                    f64 cost = 0.0;
                    for (size_t x = 0; x < n_neighbor_blocks; ++x) {
                        if (neighbor_blocks[x] != block) { cost += costs.alpha_end; }
                    }
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const vertex_t v = g.edges_v[j];
                        const partition_t v_block = p_manager.partition[v];
                        if (v_block == NO_ID) { continue; }
                        const f64 w = static_cast<f64>(g.edge_w(j));
                        const f64 dist = static_cast<f64>(oracle.get(block, v_block));
                        cost += costs.beta_end * w * static_cast<f64>(block != v_block)
                                + costs.alpha_link * dist
                                + costs.beta_link * w * dist;
                    }

                    // Penalize blocks above lmin to spread load across siblings
                    const weight_t block_lmin = p_manager.get_lmin(block);
                    if (block_lmin > 0 && p_manager.bweights[block] >= block_lmin) {
                        f64 over_ratio = static_cast<f64>(p_manager.bweights[block]) / static_cast<f64>(block_lmin);
                        cost += costs.alpha_end * over_ratio;
                    }

                    if (best == NO_ID || cost < min_cost) {
                        min_cost = cost;
                        best = block;
                    }
                }

                // fallback: all fine blocks overloaded — pick least loaded
                if (best == NO_ID) {
                    best = fine_blocks[0];
                    weight_t min_w = p_manager.bweights[fine_blocks[0]];
                    for (size_t i = 1; i < fine_blocks.size(); ++i) {
                        if (p_manager.bweights[fine_blocks[i]] < min_w) {
                            min_w = p_manager.bweights[fine_blocks[i]];
                            best = fine_blocks[i];
                        }
                    }
                }

                p_manager.set(u, u_weight, best);
                blk_verts[best].push_back(u);
            }

        }

        // Guarantee every child cell at least one vertex by stealing the cheapest-to-move vertex from a sibling.
        void steal_for_empty_cells() {
            // Guarantee every fine block gets at least one vertex.
            // For each empty fine block e, steal the cheapest-to-move vertex from
            // its most-loaded sibling (same coarse parent).  "Cheapest" = minimum
            // increase in hop cost when switching from donor block s to e.
            for (partition_t c = 0; c < coarse_k; ++c) {
                const auto &siblings = coarse_to_fine[c];
                if (siblings.size() <= 1) { continue; }

                for (const partition_t e: siblings) {
                    if (p_manager.bweights[e] > 0) { continue; }

                    vertex_t best_u = NO_ID;
                    partition_t best_donor = NO_ID;
                    f64 best_delta = std::numeric_limits<f64>::max();

                    for (const partition_t s: siblings) {
                        if (s == e) { continue; }
                        for (const vertex_t u: blk_verts[s]) {
                            if (p_manager.partition[u] != s) { continue; } // already moved
                            if (g.v_weight(u) > p_manager.get_lmax(e)) { continue; }

                            // hop-cost delta for moving u from s to e
                            f64 delta = 0.0;
                            for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                                const partition_t vb = p_manager.partition[g.edges_v[j]];
                                if (vb == NO_ID || vb == s || vb == e) { continue; }
                                const f64 w = static_cast<f64>(g.edge_w(j));
                                delta += (costs.alpha_link + costs.beta_link * w)
                                        * (static_cast<f64>(oracle.get(e, vb))
                                           - static_cast<f64>(oracle.get(s, vb)));
                            }

                            if (delta < best_delta) {
                                best_delta = delta;
                                best_u = u;
                                best_donor = s;
                            }
                        }
                    }

                    if (best_u == NO_ID) { continue; } // coarse block has zero vertices — nothing to steal

                    // Perform the move (bv_manager/q_graph not built yet, use raw arrays)
                    const weight_t uw = g.v_weight(best_u);
                    p_manager.bweights[best_donor] -= uw;
                    p_manager.partition[best_u] = e;
                    p_manager.bweights[e] += uw;
                    blk_verts[e].push_back(best_u);
                }
            }

        }

        // Fill remaining empty cells by moving sender groups along the density gradient.
        void avalanche_into_empty_cells() {
            // Avalanche steal along the density gradient.
            // For each empty block, find the adjacent PE with highest weight,
            // pick a sender group from it (all vertices sending to the same
            // target), and move the whole group.  If that empties the donor,
            // the "hole" cascades along the gradient until it reaches a block
            // that can absorb the loss.  Each move is one grid hop.
            const size_t W = static_cast<size_t>(oracle.current_W());
            // scratch for building sender groups
            std::vector<size_t> target_idx(fine_k, SIZE_MAX);
            std::vector<partition_t> target_dirty;

            std::vector<partition_t> still_empty;
            for (partition_t e = 0; e < fine_k; ++e) {
                if (p_manager.bweights[e] == 0) {
                    still_empty.push_back(e);
                }
            }

            for (const partition_t e: still_empty) {
                if (p_manager.bweights[e] > 0) { continue; }

                partition_t hole = e;
                constexpr size_t MAX_CHAIN = 64;
                for (size_t step = 0; step < MAX_CHAIN; ++step) {
                    partition_t nbrs[4];
                    const size_t n_nbrs = grid_neighbors(hole, W, fine_k, nbrs);

                    // Pick neighbor with highest block weight
                    partition_t donor = NO_ID;
                    weight_t donor_bw = 0;
                    for (size_t i = 0; i < n_nbrs; ++i) {
                        if (p_manager.bweights[nbrs[i]] > donor_bw) {
                            donor_bw = p_manager.bweights[nbrs[i]];
                            donor = nbrs[i];
                        }
                    }
                    if (donor == NO_ID) { break; }

                    // Build sender groups in donor
                    std::vector<std::vector<vertex_t> > sg_verts;
                    std::vector<partition_t> sg_target;
                    std::vector<weight_t> sg_weight;

                    for (const vertex_t u: blk_verts[donor]) {
                        if (p_manager.partition[u] != donor) { continue; }
                        const weight_t uw = g.v_weight(u);
                        partition_t ext_tgts[256];
                        size_t n_ext = 0;
                        for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                            const partition_t vb = p_manager.partition[g.edges_v[j]];
                            if (vb == NO_ID || vb == donor) { continue; }
                            bool dup = false;
                            for (size_t x = 0; x < n_ext; ++x) {
                                if (ext_tgts[x] == vb) {
                                    dup = true;
                                    break;
                                }
                            }
                            if (!dup && n_ext < 256) { ext_tgts[n_ext++] = vb; }
                        }
                        for (size_t x = 0; x < n_ext; ++x) {
                            const partition_t D = ext_tgts[x];
                            if (target_idx[D] == SIZE_MAX) {
                                target_idx[D] = sg_verts.size();
                                sg_verts.emplace_back();
                                sg_target.push_back(D);
                                sg_weight.push_back(0);
                                target_dirty.push_back(D);
                            }
                            const size_t idx = target_idx[D];
                            sg_verts[idx].push_back(u);
                            sg_weight[idx] += uw;
                        }
                    }
                    for (const partition_t D: target_dirty) { target_idx[D] = SIZE_MAX; }
                    target_dirty.clear();

                    // Pick lightest group that fits in hole
                    const weight_t hole_cap = p_manager.get_lmax(hole);
                    size_t pick = SIZE_MAX;
                    weight_t pick_w = std::numeric_limits<weight_t>::max();
                    for (size_t i = 0; i < sg_verts.size(); ++i) {
                        if (sg_weight[i] > hole_cap) { continue; }
                        if (sg_weight[i] < pick_w) {
                            pick_w = sg_weight[i];
                            pick = i;
                        }
                    }
                    if (pick == SIZE_MAX) { break; }

                    // Move the group
                    for (const vertex_t u: sg_verts[pick]) {
                        if (p_manager.partition[u] != donor) { continue; }
                        const weight_t uw = g.v_weight(u);
                        p_manager.bweights[donor] -= uw;
                        p_manager.partition[u] = hole;
                        p_manager.bweights[hole] += uw;
                        blk_verts[hole].push_back(u);
                    }

                    if (p_manager.bweights[donor] > 0) { break; }
                    hole = donor;
                }
            }
        }

        // Move vertices of small connected regions of non-empty cells into the largest region.
        void evacuate_islands() {
            // Evacuate small islands of non-empty PEs.
            // BFS on the PE grid to find connected components of non-empty
            // blocks.  Small components get their vertices moved to the
            // nearest block in the largest component.
            const size_t W = static_cast<size_t>(oracle.current_W());
            // BFS to label connected components of non-empty blocks
            std::vector<partition_t> comp_id(fine_k, NO_ID);
            std::vector<std::vector<partition_t> > components;
            std::vector<partition_t> bfs_queue;

            for (partition_t b = 0; b < fine_k; ++b) {
                if (p_manager.bweights[b] == 0 || comp_id[b] != NO_ID) { continue; }
                const partition_t cid = static_cast<partition_t>(components.size());
                components.emplace_back();
                bfs_queue.clear();
                bfs_queue.push_back(b);
                comp_id[b] = cid;
                for (size_t head = 0; head < bfs_queue.size(); ++head) {
                    const partition_t cur = bfs_queue[head];
                    components[cid].push_back(cur);
                    partition_t nbrs[4];
                    const size_t n_nbrs = grid_neighbors(cur, W, fine_k, nbrs);
                    for (size_t i = 0; i < n_nbrs; ++i) {
                        if (p_manager.bweights[nbrs[i]] > 0 && comp_id[nbrs[i]] == NO_ID) {
                            comp_id[nbrs[i]] = cid;
                            bfs_queue.push_back(nbrs[i]);
                        }
                    }
                }
            }

            // Find the largest component
            size_t largest_idx = 0;
            for (size_t i = 1; i < components.size(); ++i) {
                if (components[i].size() > components[largest_idx].size()) {
                    largest_idx = i;
                }
            }

            // Evacuate all small components into the largest
            for (size_t ci = 0; ci < components.size(); ++ci) {
                if (ci == largest_idx) { continue; }
                for (const partition_t b: components[ci]) {
                    if (p_manager.bweights[b] == 0) { continue; }
                    // Find nearest block in largest component with capacity
                    const weight_t b_weight = p_manager.bweights[b];
                    partition_t target = NO_ID;
                    weight_t best_dist = std::numeric_limits<weight_t>::max();
                    for (const partition_t s: components[largest_idx]) {
                        if (p_manager.bweights[s] + b_weight > p_manager.get_lmax(s)) { continue; }
                        weight_t d = oracle.get(b, s);
                        if (d < best_dist) {
                            best_dist = d;
                            target = s;
                        }
                    }
                    if (target == NO_ID) { continue; }
                    for (const vertex_t u: blk_verts[b]) {
                        if (p_manager.partition[u] != b) { continue; }
                        const weight_t uw = g.v_weight(u);
                        p_manager.bweights[b] -= uw;
                        p_manager.partition[u] = target;
                        p_manager.bweights[target] += uw;
                        blk_verts[target].push_back(u);
                    }
                }
            }
        }

        // Brute-force the assignment of vertex groups to child cells for parents with 2-4 children.
        void optimize_sibling_permutations() {
            // Brute-force sibling permutation optimization.
            // For each coarse parent with 2-4 children, try all permutations
            // of which vertex group maps to which child PE; keep the best.
            auto perm_cost = [&](const std::vector<vertex_t> &verts, partition_t block) -> f64 {
                f64 cost = 0.0;
                for (const vertex_t u : verts) {
                    for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                        const partition_t nb = p_manager.partition[g.edges_v[j]];
                        if (nb == NO_ID) continue;
                        const f64 w = static_cast<f64>(g.edge_w(j));
                        const f64 dist = static_cast<f64>(oracle.get(block, nb));
                        if (block != nb) {
                            cost += costs.alpha_end + costs.beta_end * w;
                        }
                        cost += costs.alpha_link * dist + costs.beta_link * w * dist;
                    }
                }
                return cost;
            };

            for (partition_t c = 0; c < coarse_k; ++c) {
                const auto &siblings = coarse_to_fine[c];
                const size_t ns = siblings.size();
                if (ns <= 1 || ns > 4) continue;

                // Collect vertex sets for each sibling
                std::vector<std::vector<vertex_t>> sib_verts(ns);
                std::vector<weight_t> sib_weights(ns);
                for (size_t i = 0; i < ns; ++i) {
                    for (const vertex_t u : blk_verts[siblings[i]]) {
                        if (p_manager.partition[u] == siblings[i]) {
                            sib_verts[i].push_back(u);
                        }
                    }
                    sib_weights[i] = p_manager.bweights[siblings[i]];
                }

                // Generate all permutations of indices [0..ns-1]
                std::vector<size_t> perm(ns);
                std::iota(perm.begin(), perm.end(), 0);

                std::vector<size_t> best_perm = perm;
                f64 best_cost = std::numeric_limits<f64>::max();

                // Temporarily assign identity permutation to compute baseline
                // (partition already has identity mapping)
                do {
                    // Check capacity: vertex group i goes to PE siblings[perm[i]]
                    bool feasible = true;
                    for (size_t i = 0; i < ns; ++i) {
                        if (sib_weights[i] > p_manager.get_lmax(siblings[perm[i]])) {
                            feasible = false;
                            break;
                        }
                    }
                    if (!feasible) continue;

                    // Apply permutation temporarily
                    for (size_t i = 0; i < ns; ++i) {
                        for (const vertex_t u : sib_verts[i]) {
                            p_manager.partition[u] = siblings[perm[i]];
                        }
                    }

                    f64 cost = 0.0;
                    for (size_t i = 0; i < ns; ++i) {
                        cost += perm_cost(sib_verts[i], siblings[perm[i]]);
                    }

                    if (cost < best_cost) {
                        best_cost = cost;
                        best_perm = perm;
                    }
                } while (std::next_permutation(perm.begin(), perm.end()));

                // Apply best permutation
                for (size_t i = 0; i < ns; ++i) {
                    for (const vertex_t u : sib_verts[i]) {
                        p_manager.partition[u] = siblings[best_perm[i]];
                    }
                    p_manager.bweights[siblings[best_perm[i]]] = sib_weights[i];
                }

                // Rebuild blk_verts for affected siblings
                for (size_t i = 0; i < ns; ++i) {
                    blk_verts[siblings[i]].clear();
                }
                for (size_t i = 0; i < ns; ++i) {
                    for (const vertex_t u : sib_verts[i]) {
                        blk_verts[siblings[best_perm[i]]].push_back(u);
                    }
                }
            }
        }
    };
}

#endif //HEIPROMAP_TOPOLOGY_UNCOARSENER_H
