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

#ifndef HEIPROMAP_FLOW_BASED_REFINEMENT_H
#define HEIPROMAP_FLOW_BASED_REFINEMENT_H

#include <algorithm>

#include "quotient_graph_refinement.h"
#include "../../commons/indexed_update_heap.h"
#include "../../commons/utils.h"
#include "../datastructures/functions.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/qap.h"

namespace HeiProMap {
    class FlowNetwork {
        vertex_t s;
        partition_t s_id;
        vertex_t t;
        partition_t t_id;
        std::vector<std::vector<EdgeVW>> adj;
        std::vector<partition_t> partition;

    public:
        explicit FlowNetwork(size_t n) {
            s    = 0;
            s_id = 0;
            t    = 0;
            t_id = 0;
            adj.resize(n);
            partition.resize(n);
        }

        void add(vertex_t u, vertex_t v, weight_t w) {
            adj[u].emplace_back(v, w);
        }

        void set_s_t_vertex(vertex_t t_s, partition_t t_s_id, vertex_t t_t, partition_t t_t_id) {
            s    = t_s;
            s_id = t_s_id;
            t    = t_t;
            t_id = t_t_id;
        }

        void solve() {
            // solve the flow network
        }

        partition_t get(vertex_t u) {
            return partition[u];
        }
    };

    class FlowBasedRefinementConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit FlowBasedRefinementConfiguration(const std::string& t_name) : ISerialRefinerConfiguration(t_name) {}
        u64 max_iteration = 1;
    };

    class FlowBasedRefinement final : public ISerialRefiner {
        vertex_t m_n    = 0;
        vertex_t m_m    = 0;
        partition_t m_k = 0;
        weight_t m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;

        u32* vertex_used = nullptr;
        u32 vertex_mark  = 0;

        u32* block_used = nullptr;
        u32 block_mark  = 0;

        vertex_t* curr_boundary   = nullptr;
        size_t curr_boundary_size = 0;

        Move* moves       = nullptr;
        size_t moves_size = 0;

        // active block scheduling
        u8* active_this_round = nullptr;
        u8* active_next_round = nullptr;
        PairWeight* pairs     = nullptr;
        size_t pairs_size     = 0;

        RandomEngine* random_engine                    = nullptr;
        const FlowBasedRefinementConfiguration* config = nullptr;
        StatisticCollector* m_stat_collector           = nullptr;

    public:
        FlowBasedRefinement() = default;

        ~FlowBasedRefinement() override {
            free(vertex_used);
            free(block_used);
            free(curr_boundary);
            free(moves);
        }

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const weight_t t_lmax,
                        const std::vector<partition_t>& t_hierarchy,
                        const std::vector<weight_t>& t_distance,
                        RandomEngine& t_random_engine,
                        const ISerialRefinerConfiguration& i_config,
                        StatisticCollector& t_stat_collect) override {
            m_n         = t_n;
            m_m         = t_m;
            m_k         = t_k;
            m_lmax      = t_lmax;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;

            random_engine    = &t_random_engine;
            config           = dynamic_cast<const FlowBasedRefinementConfiguration*>(&i_config);
            m_stat_collector = &t_stat_collect;

            vertex_t m_n_64        = round_up_64(m_n);
            partition_t m_k_64     = round_up_64(m_k);
            partition_t m_k_m_k_64 = round_up_64(m_k * m_k);

            vertex_used = (u32*)aligned_alloc(64, sizeof(u32) * m_n_64);
            std::fill_n(vertex_used, m_n_64, 0);
            vertex_mark = 0;

            block_used = (u32*)aligned_alloc(64, sizeof(u32) * m_k_64);
            std::fill_n(block_used, m_k_64, 0);
            block_mark = 0;

            curr_boundary      = (vertex_t*)aligned_alloc(64, sizeof(vertex_t) * m_n_64);
            curr_boundary_size = 0;

            moves      = (Move*)aligned_alloc(64, sizeof(Move) * m_n_64);
            moves_size = 0;

            // active block scheduling
            active_this_round = (u8*)aligned_alloc(64, m_k_64 * sizeof(u8));
            active_next_round = (u8*)aligned_alloc(64, m_k_64 * sizeof(u8));
            pairs             = (PairWeight*)aligned_alloc(64, m_k_m_k_64 * sizeof(PairWeight));
            pairs_size        = 0;
        }

        void refine(const u64 level,
                    const u64 max_level,
                    const graph_t& g,
                    const d_oracle_t& d_oracle,
                    bv_manager_t& bv_manager,
                    p_manager_t& p_manager,
                    q_graph_t& q_graph) override {
            for (u64 iteration = 0; iteration < config->max_iteration; ++iteration) {
                // determine all pairs in the quotient graph
                pairs_size = 0;
                for (partition_t u_id = 0; u_id < m_k; ++u_id) {
                    for (partition_t v_id = u_id + 1; v_id < m_k; ++v_id) {
                        if (q_graph.has_edge(u_id, v_id) && (active_this_round[u_id] || active_this_round[v_id])) {
                            pairs[pairs_size++] = {u_id, v_id, d_oracle.get(u_id, v_id)};
                        }
                    }
                }
                std::sort(pairs, pairs + pairs_size, std::greater<>());

                for (size_t i = 0; i < pairs_size; ++i) {
                    partition_t left_id  = pairs[i].id1;
                    partition_t right_id = pairs[i].id2;
                    refine_blocks(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, left_id, right_id);
                }
            }
        }

        void refine_blocks(const u64 level,
                           const u64 max_level,
                           const graph_t& g,
                           const d_oracle_t& d_oracle,
                           bv_manager_t& bv_manager,
                           p_manager_t& p_manager,
                           q_graph_t& q_graph,
                           partition_t left_id,
                           partition_t right_id) {
            // get boundary vertices
            std::vector<vertex_t> left_boundary_vertices  = get_boundary_vertices(g, bv_manager, p_manager, left_id, right_id);
            std::vector<vertex_t> right_boundary_vertices = get_boundary_vertices(g, bv_manager, p_manager, right_id, left_id);

            // calc max weight for each bfs
            weight_t left_bfs_max_weight  = m_lmax - p_manager.get_bweight(right_id);
            weight_t right_bfs_max_weight = m_lmax - p_manager.get_bweight(left_id);

            // get both regions
            std::vector<vertex_t> left_region  = grow_bfs(g, p_manager, left_boundary_vertices, left_id, left_bfs_max_weight);
            std::vector<vertex_t> right_region = grow_bfs(g, p_manager, right_boundary_vertices, right_id, right_bfs_max_weight);

            // determine penalties for all vertices
            std::vector<weight_t> left_penalties  = determine_penalties(g, p_manager, d_oracle, left_boundary_vertices, left_id, right_boundary_vertices);
            std::vector<weight_t> right_penalties = determine_penalties(g, p_manager, d_oracle, right_boundary_vertices, right_id, left_boundary_vertices);

            // build a translation table from graph to flow network
            TranslationTable<vertex_t> translation_table;
            translation_table.reserve(g.get_n(), g.get_n());

            vertex_t new_u = 0;
            for (const vertex_t u : left_region) { translation_table.add(u, new_u++); }
            for (const vertex_t u : right_region) { translation_table.add(u, new_u++); }
            vertex_t s = new_u++;
            vertex_t t = new_u;

            // build flownetwork
            size_t n = left_region.size() + right_region.size() + 2;
            FlowNetwork flow_network(n);
            flow_network.set_s_t_vertex(s, left_id, t, right_id);
            for (const vertex_t u : left_region) {
                forall_guivw(g, u, i, v, w)
                    {
                        partition_t v_id = p_manager[v];
                        if (v_id != left_id && v_id != right_id) { continue; }
                        vertex_t new_u = translation_table.get_n(u);
                        vertex_t new_v = translation_table.get_n(v);

                        flow_network.add(new_u, new_v, w * d_oracle.get(left_id, v_id));
                    }
                endfor
            }
            for (const vertex_t u : right_region) {
                forall_guivw(g, u, i, v, w)
                    {
                        partition_t v_id = p_manager[v];
                        if (v_id != left_id && v_id != right_id) { continue; }
                        vertex_t new_u = translation_table.get_n(u);
                        vertex_t new_v = translation_table.get_n(v);

                        flow_network.add(new_u, new_v, w * d_oracle.get(right_id, v_id));
                    }
                endfor
            }

            // add the penalties
            for (const vertex_t u : left_region) {
                vertex_t new_u         = translation_table.get_n(u);
                weight_t left_penalty  = left_penalties[u];
                weight_t right_penalty = right_penalties[u];
                flow_network.add(s, new_u, left_penalty);
                flow_network.add(t, new_u, right_penalty);
            }
            for (const vertex_t u : right_region) {
                vertex_t new_u         = translation_table.get_n(u);
                weight_t left_penalty  = left_penalties[u];
                weight_t right_penalty = right_penalties[u];
                flow_network.add(s, new_u, left_penalty);
                flow_network.add(t, new_u, right_penalty);
            }

            // solve the flow network
            flow_network.solve();

            // make the changes
            std::vector<vertex_t> vertices_to_move;
            for (const vertex_t u : left_region) {
                vertex_t new_u = translation_table.get_n(u);
                if (left_id != flow_network.get(new_u)) {
                    vertices_to_move.push_back(u);
                }
            }
            for (const vertex_t u : right_region) {
                vertex_t new_u = translation_table.get_n(u);
                if (right_id != flow_network.get(new_u)) {
                    vertices_to_move.push_back(u);
                }
            }

            while (!vertices_to_move.empty()) {
                for (size_t i = 0; i < vertices_to_move.size(); ++i) {
                    vertex_t u          = vertices_to_move[i];
                    weight_t u_weight   = g.weight(u);
                    partition_t u_id    = p_manager[u];
                    partition_t move_id = u_id == left_id ? right_id : left_id;

                    if (bv_manager.is_boundary(u)) {
                        // the vertex is boundary so move it
                        bv_manager.move(g, p_manager, u, u_id, move_id);
                        q_graph.move(g, p_manager, u, u_id, move_id);
                        p_manager.move(u, u_weight, u_id, move_id);

                        std::swap(vertices_to_move[i], vertices_to_move.back());
                        vertices_to_move.pop_back();
                        i -= 1;
                    }
                }
            }
        }

        std::vector<vertex_t> get_boundary_vertices(const graph_t& g,
                                                    bv_manager_t& bv_manager,
                                                    p_manager_t& p_manager,
                                                    partition_t left_id,
                                                    partition_t right_id) {
            std::vector<vertex_t> boundary_vertices;
            forall_bv_id_iu(bv_manager, left_id, i, u)
                {
                    forall_guiv(g, u, j, v)
                        {
                            partition_t v_id = p_manager[v];
                            if (v_id == right_id) {
                                boundary_vertices.push_back(v);
                                break;
                            }
                        }
                    endfor
                }
            endfor
            return boundary_vertices;
        }

        std::vector<vertex_t> grow_bfs(const graph_t& g,
                                       p_manager_t& p_manager,
                                       std::vector<vertex_t>& boundary_vertices,
                                       partition_t id,
                                       weight_t max_weight) {
            std::vector<u8> seen(g.get_n(), 0);
            weight_t curr_weight = 0;

            std::deque<vertex_t> queue;
            for (vertex_t u : boundary_vertices) {
                queue.push_back(u);
            }
            std::vector<vertex_t> region;

            while (!queue.empty()) {
                vertex_t u = queue.front();
                queue.pop_front();
                if (seen[u] == 1) { continue; }
                seen[u] = 1;
                if (curr_weight + g.weight(u) <= max_weight) {
                    region.push_back(u);
                    curr_weight += g.weight(u);
                    forall_guiv(g, u, i, v)
                        {
                            partition_t v_id = p_manager[v];
                            if (v_id == id && seen[v] == 0) {
                                queue.push_back(v);
                            }
                        }
                    endfor
                }
            }

            return region;
        }

        std::vector<weight_t> determine_penalties(const graph_t& g,
                                                  p_manager_t& p_manager,
                                                  const d_oracle_t& d_oracle,
                                                  std::vector<vertex_t>& boundary_vertices_1,
                                                  partition_t id_1,
                                                  std::vector<vertex_t>& boundary_vertices_2) {
            std::vector<weight_t> penalties(g.get_n());
            for (vertex_t u : boundary_vertices_1) {
                penalties[u] = 0;
                forall_guivw(g, u, i, v, w)
                    {
                        if (exists(boundary_vertices_1, v) || exists(boundary_vertices_2, v)) { continue; } // ignore neighbors that are in the region
                        partition_t v_id = p_manager[v];
                        penalties[u] += w * d_oracle.get(id_1, v_id);
                    }
                endfor
            }
            for (vertex_t u : boundary_vertices_2) {
                penalties[u] = 0;
                forall_guivw(g, u, i, v, w)
                    {
                        if (exists(boundary_vertices_1, v) || exists(boundary_vertices_2, v)) { continue; } // ignore neighbors that are in the region
                        partition_t v_id = p_manager[v];
                        penalties[u] += w * d_oracle.get(id_1, v_id);
                    }
                endfor
            }
            return penalties;
        }

        JSONString get_stats() override {
            std::string stats = "{ \n";
#if COLLECT_METRICS

#endif
            stats.pop_back();
            stats.pop_back();
            stats += "\n}";

            JSONString json_stats;
            json_stats.s = stats;
            return json_stats;
        }
    };
}

#endif //HEIPROMAP_FLOW_BASED_REFINEMENT_H
