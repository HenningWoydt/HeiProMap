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

#include "../../extern/maxflow-v3.04.src/graph.h"

#include "quotient_graph_refinement.h"
#include "../../commons/indexed_update_heap.h"
#include "../../commons/random_engine.h"
#include "../../commons/statistic_collector.h"
#include "../../commons/utils.h"
#include "../datastructures/functions.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/qap.h"

namespace HeiProMap {
    class FlowNetwork {
        Graph<weight_t, weight_t, weight_t> g;
        partition_t                         s_id;
        partition_t                         t_id;

    public:
        FlowNetwork() : g(Graph<weight_t, weight_t, weight_t>(0, 0)) {
            s_id = 0;
            t_id = 0;
        }

        void initialize(size_t t_n) {
            g.reset();
            g.add_node(t_n);
            s_id = 0;
            t_id = 0;
        }

        void add(vertex_t u, vertex_t v, weight_t w) {
            g.add_edge(u, v, w, w);
        }

        void add_s_edge(vertex_t v, weight_t w) {
            g.add_tweights(v, w, 0);
        }

        void add_t_edge(vertex_t v, weight_t w) {
            g.add_tweights(v, 0, w);
        }

        void set_s_t_vertex(partition_t t_s_id, partition_t t_t_id) {
            s_id = t_s_id;
            t_id = t_t_id;
        }

        void solve() {
            g.maxflow();
        }

        partition_t get(vertex_t u) {
            return g.what_segment(u) == Graph<weight_t, weight_t, weight_t>::SOURCE ? s_id : t_id;
        }
    };

    class FlowBasedRefinementConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit FlowBasedRefinementConfiguration(const std::string &t_name) : ISerialRefinerConfiguration(t_name) {}

        u64 max_global_iteration = 1;
        u64 max_local_iteration  = 3;
        f64 alpha                = 2.0;
        f64 alpha_upper_bound    = 8.0;
        f64 alpha_modifier       = 2.0;
    };

    class FlowBasedRefinement final : public ISerialRefiner {
        vertex_t                 m_n         = 0;
        vertex_t                 m_m         = 0;
        partition_t              m_k         = 0;
        f64                      m_imbalance = 0.0;
        weight_t                 m_lmax      = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t>    m_distance;

        // active block scheduling
        u8         *active_this_round = nullptr;
        u8         *active_next_round = nullptr;
        PairWeight *pairs             = nullptr;
        size_t pairs_size = 0;

        // array for boundary vertices
        vertex_t *left_boundary = nullptr;
        size_t left_boundary_size = 0;

        vertex_t *right_boundary = nullptr;
        size_t right_boundary_size = 0;

        // array for regions
        vertex_t *left_region = nullptr;
        size_t left_region_size = 0;

        vertex_t *right_region = nullptr;
        size_t right_region_size = 0;

        u32 *is_left_region  = nullptr;
        u32 *is_right_region = nullptr;
        u32 is_region_mark = 0;

        u32 *bfs_level = nullptr;

        vertex_t *queue = nullptr;
        size_t queue_size = 0;

        u32 *seen = nullptr;
        u32 seen_mark = 0;

        // array for penalties
        weight_t *left_penalties  = nullptr;
        weight_t *right_penalties = nullptr;

        //Translation Table for mapping
        TranslationTable<vertex_t> translation_table;

        FlowNetwork flow_network;

        RandomEngine                           *random_engine    = nullptr;
        const FlowBasedRefinementConfiguration *config           = nullptr;
        StatisticCollector                     *m_stat_collector = nullptr;

    public:
        FlowBasedRefinement() = default;

        ~FlowBasedRefinement() override {
            free(active_this_round);
            free(active_next_round);

            free(left_boundary);
            free(right_boundary);

            free(left_region);
            free(right_region);
            free(is_left_region);
            free(is_right_region);
            free(queue);
            free(seen);

            free(left_penalties);
            free(right_penalties);
        }

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const f64 t_imbalance,
                        const weight_t t_lmax,
                        const std::vector<partition_t> &t_hierarchy,
                        const std::vector<weight_t> &t_distance,
                        RandomEngine &t_random_engine,
                        const ISerialRefinerConfiguration &i_config,
                        StatisticCollector &t_stat_collect) override {
            m_n         = t_n;
            m_m         = t_m;
            m_k         = t_k;
            m_lmax      = t_lmax;
            m_imbalance = t_imbalance;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;

            random_engine    = &t_random_engine;
            config           = dynamic_cast<const FlowBasedRefinementConfiguration *>(&i_config);
            m_stat_collector = &t_stat_collect;

            vertex_t    m_n_64     = round_up_64(m_n);
            partition_t m_k_64     = round_up_64(m_k);
            partition_t m_k_m_k_64 = round_up_64(m_k * m_k);

            // active block scheduling
            active_this_round = (u8 *) aligned_alloc(64, m_k_64 * sizeof(u8));
            active_next_round = (u8 *) aligned_alloc(64, m_k_64 * sizeof(u8));
            pairs             = (PairWeight *) aligned_alloc(64, m_k_m_k_64 * sizeof(PairWeight));
            pairs_size        = 0;

            left_boundary      = (vertex_t *) aligned_alloc(64, m_n_64 * sizeof(vertex_t));
            left_boundary_size = 0;

            right_boundary      = (vertex_t *) aligned_alloc(64, m_n_64 * sizeof(vertex_t));
            right_boundary_size = 0;

            left_region      = (vertex_t *) aligned_alloc(64, m_n_64 * sizeof(vertex_t));
            left_region_size = 0;

            right_region      = (vertex_t *) aligned_alloc(64, m_n_64 * sizeof(vertex_t));
            right_region_size = 0;

            bfs_level = (u32 *) aligned_alloc(64, m_n_64 * sizeof(u32));

            is_left_region = (u32 *) aligned_alloc(64, m_n_64 * sizeof(u32));
            std::fill_n(is_left_region, m_n_64, 0);
            is_right_region = (u32 *) aligned_alloc(64, m_n_64 * sizeof(u32));
            std::fill_n(is_right_region, m_n_64, 0);
            is_region_mark = 0;

            queue      = (vertex_t *) aligned_alloc(64, m_n_64 * sizeof(vertex_t));
            queue_size = 0;

            seen = (u32 *) aligned_alloc(64, m_n_64 * sizeof(u32));
            std::fill_n(seen, m_n_64, 0);
            seen_mark = 0;

            left_penalties  = (weight_t *) aligned_alloc(64, m_n_64 * sizeof(weight_t));
            right_penalties = (weight_t *) aligned_alloc(64, m_n_64 * sizeof(weight_t));

            translation_table.reserve(m_n_64, m_n_64);
        }

        void refine(const u64 level,
                    const u64 max_level,
                    const graph_t &g,
                    const d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph) override {
            std::fill_n(active_this_round, m_k, 1);
            std::fill_n(active_next_round, m_k, 0);

            for (u64 iteration = 0; iteration < config->max_global_iteration; ++iteration) {
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

                if (pairs_size == 0) { return; }

                for (size_t i = 0; i < pairs_size; ++i) {
                    partition_t left_id  = pairs[i].id1;
                    partition_t right_id = pairs[i].id2;
                    refine_blocks(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, left_id, right_id);
                }

                std::swap(active_this_round, active_next_round);
                std::fill_n(active_next_round, m_k, 0);

                HEAVYASSERT(assert_state_after_partitioning(g, p_manager, bv_manager, q_graph, m_k));
            }
        }

        void refine_blocks(const u64 level,
                           const u64 max_level,
                           const graph_t &g,
                           const d_oracle_t &d_oracle,
                           bv_manager_t &bv_manager,
                           p_manager_t &p_manager,
                           q_graph_t &q_graph,
                           partition_t left_id,
                           partition_t right_id) {
            ASSERT(left_id != right_id);

            f64 alpha             = config->alpha;
            f64 alpha_upper_bound = config->alpha_upper_bound;
            f64 alpha_modifier    = config->alpha_modifier;

            u64 max_local_iteration = config->max_local_iteration;
            u64 iteration           = 0;
            while (iteration < max_local_iteration) {
                iteration += 1;
                // get boundary vertices
                determine_boundary_vertices(g, bv_manager, p_manager, left_id, right_id);

                // calc max weight for each bfs
                weight_t lmax             = std::ceil((1.0 + m_imbalance * alpha) * ((f64) g.weight() / (f64) m_k));
                weight_t left_max_weight  = lmax - p_manager.get_bweight(right_id);
                weight_t right_max_weight = lmax - p_manager.get_bweight(left_id);

                // get both regions
                determine_regions(g, p_manager, left_id, left_max_weight, right_id, right_max_weight);

                if (left_region_size + right_region_size == 0) {
                    // if both regions are empty, increase their sizes
                    if (alpha == alpha_upper_bound) {
                        // the regions will not change next iteration, so just return
                        return;
                    }
                    // alpha = std::min(random_engine->get_f64(1.5, 2.0) * alpha, alpha_upper_bound);
                    alpha = std::min(alpha_modifier * alpha, alpha_upper_bound);
                    continue;
                }

                // determine penalties for all vertices
                determine_penalties(g, p_manager, d_oracle, left_id, right_id);

                // build a translation table from graph to flow network
                vertex_t    new_u = 0;
                for (size_t i     = 0; i < left_region_size; ++i) { translation_table.add(left_region[i], new_u++); }
                for (size_t i     = 0; i < right_region_size; ++i) { translation_table.add(right_region[i], new_u++); }

                // build flownetwork
                build_flow_network(g, d_oracle, left_id, right_id);

                // solve the flow network
                flow_network.solve();

                // check if it is a valid cut
                if (!valid_cut(g, p_manager, left_id, right_id)) {
                    // alpha = std::min(alpha / random_engine->get_f64(1.5, 2.0), alpha_upper_bound);
                    alpha = std::max(alpha / alpha_modifier, 1.0);
                    continue;
                }

                // check if the cut actually changes the partition
                if (!cut_changes_partition(left_id, right_id)) {
                    // cut is valid, but does not change anything
                    return;
                }

                // cut is valid and changes the partition, increase alpha
                // alpha = std::min(random_engine->get_f64(1.5, 2.0) * alpha, alpha_upper_bound);
                alpha = std::min(alpha * alpha_modifier, alpha_upper_bound);

                // make the changes
                change_boundary(g, bv_manager, p_manager, q_graph, left_id, right_id);

                active_next_round[left_id]  = 1;
                active_next_round[right_id] = 1;
            }
        }

        void determine_boundary_vertices(const graph_t &g,
                                         const bv_manager_t &bv_manager,
                                         const p_manager_t &p_manager,
                                         partition_t left_id,
                                         partition_t right_id) {
            left_boundary_size = 0;
            forall_bv_id_iu(bv_manager, left_id, i, u)
                {
                    forall_guiv(g, u, j, v)
                        {
                            partition_t v_id = p_manager[v];
                            if (v_id == right_id) {
                                left_boundary[left_boundary_size++] = u;
                                break;
                            }
                        }
                    endfor
                }
            endfor

            right_boundary_size = 0;
            forall_bv_id_iu(bv_manager, right_id, i, u)
                {
                    forall_guiv(g, u, j, v)
                        {
                            partition_t v_id = p_manager[v];
                            if (v_id == left_id) {
                                right_boundary[right_boundary_size++] = u;
                                break;
                            }
                        }
                    endfor
                }
            endfor
        }

        void determine_regions(const graph_t &g,
                               const p_manager_t &p_manager,
                               partition_t left_id,
                               weight_t left_max_weight,
                               partition_t right_id,
                               weight_t right_max_weight) {
            is_region_mark += 1;
            seen_mark += 2;
            // seen[u] == seen_mark     means u is processed
            // seen[u] == seen_mark - 1 means u is in the queue

            weight_t left_curr_weight = 0;

            queue_size = 0;
            for (size_t i = 0; i < left_boundary_size; ++i) {
                vertex_t u = left_boundary[i];
                queue[queue_size++] = u;
                seen[u]             = seen_mark - 1;
                bfs_level[u]        = 0;
            }

            size_t queue_idx = 0;
            left_region_size = 0;
            while (queue_idx < queue_size) {
                vertex_t u = queue[queue_idx++];
                if (seen[u] == seen_mark) { continue; }

                if (left_curr_weight + g.weight(u) <= left_max_weight) {
                    left_region[left_region_size++] = u;
                    is_left_region[u]               = is_region_mark;
                    left_curr_weight += g.weight(u);
                    forall_guiv(g, u, i, v)
                        {
                            bfs_level[v] = std::min(bfs_level[v], bfs_level[u] + 1);

                            partition_t v_id = p_manager[v];
                            if (v_id == left_id && seen[v] != seen_mark && seen[v] != seen_mark - 1) {
                                queue[queue_size++] = v;
                                seen[v]             = seen_mark - 1;
                                bfs_level[v]        = bfs_level[u] + 1;
                            }
                        }
                    endfor
                }
                seen[u]    = seen_mark;
            }

            weight_t right_curr_weight = 0;

            queue_size = 0;
            for (size_t i = 0; i < right_boundary_size; ++i) {
                vertex_t u = right_boundary[i];
                queue[queue_size++] = u;
                seen[u]             = seen_mark - 1;
                bfs_level[u]        = 0;
            }

            right_region_size = 0;
            queue_idx         = 0;
            while (queue_idx < queue_size) {
                vertex_t u = queue[queue_idx++];
                if (seen[u] == seen_mark) { continue; }

                if (right_curr_weight + g.weight(u) <= right_max_weight) {
                    right_region[right_region_size++] = u;
                    is_right_region[u]                = is_region_mark;
                    right_curr_weight += g.weight(u);
                    forall_guiv(g, u, i, v)
                        {
                            partition_t v_id = p_manager[v];
                            if (v_id == right_id && seen[v] != seen_mark && seen[v] != seen_mark - 1) {
                                queue[queue_size++] = v;
                                seen[v]             = seen_mark - 1;
                                bfs_level[v]        = bfs_level[u] + 1;
                            }
                        }
                    endfor
                }
                seen[u]    = seen_mark;
            }
        }

        void determine_penalties(const graph_t &g,
                                 const p_manager_t &p_manager,
                                 const d_oracle_t &d_oracle,
                                 partition_t left_id,
                                 partition_t right_id) {
            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u = left_region[j];
                left_penalties[u]  = 0;
                right_penalties[u] = 0;
                forall_guivw(g, u, i, v, w)
                    {
                        if (is_left_region[v] == is_region_mark || is_right_region[v] == is_region_mark) { continue; } // ignore neighbors that are in the region
                        partition_t v_id = p_manager[v];
                        left_penalties[u] += w * d_oracle.get(left_id, v_id);
                        right_penalties[u] += w * d_oracle.get(right_id, v_id);
                    }
                endfor
                left_penalties[u] *= 2;
                right_penalties[u] *= 2;
            }
            for (size_t j = 0; j < right_region_size; ++j) {
                vertex_t u         = right_region[j];
                left_penalties[u]  = 0;
                right_penalties[u] = 0;
                forall_guivw(g, u, i, v, w)
                    {
                        if (is_right_region[v] == is_region_mark || is_left_region[v] == is_region_mark) { continue; } // ignore neighbors that are in the region
                        partition_t v_id = p_manager[v];
                        left_penalties[u] += w * d_oracle.get(left_id, v_id);
                        right_penalties[u] += w * d_oracle.get(right_id, v_id);
                    }
                endfor
                left_penalties[u] *= 2;
                right_penalties[u] *= 2;
            }
        }

        void build_flow_network(const graph_t &g,
                                const d_oracle_t &d_oracle,
                                partition_t left_id,
                                partition_t right_id) {
            weight_t distance = d_oracle.get(left_id, right_id);

            // build flownetwork
            size_t n = left_region_size + right_region_size;
            flow_network.initialize(n);
            flow_network.set_s_t_vertex(left_id, right_id);

            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u = left_region[j];
                forall_guivw(g, u, i, v, w)
                    {
                        weight_t mult = 2;
                        if (is_left_region[v] != is_region_mark && is_right_region[v] != is_region_mark) { continue; } // vertex gets handled by penalties
                        if (is_left_region[v] == is_region_mark && bfs_level[u] < bfs_level[v]) { continue; } // only forward edges allowed
                        if (is_left_region[v] == is_region_mark && bfs_level[u] == bfs_level[v]) { mult = 1; } // the edge from v to u will also be added, therefore only one time the weight

                        vertex_t new_u = translation_table.get_n(u);
                        vertex_t new_v = translation_table.get_n(v);

                        ASSERT(new_u != new_v);

                        flow_network.add(new_u, new_v, mult * w * distance);
                    }
                endfor
            }
            for (size_t j = 0; j < right_region_size; ++j) {
                vertex_t u = right_region[j];
                forall_guivw(g, u, i, v, w)
                    {
                        weight_t mult = 2;
                        if (is_right_region[v] != is_region_mark) { continue; }  // vertex gets handled by penalties, or if v belongs to the left region, no edge is made
                        if (bfs_level[u] > bfs_level[v]) { continue; } // only forward edges allowed
                        if (bfs_level[u] == bfs_level[v]) { mult = 1; } // the edge from v to u will also be added, therefore only one time the weight
                        vertex_t new_u = translation_table.get_n(u);
                        vertex_t new_v = translation_table.get_n(v);

                        ASSERT(new_u != new_v);

                        flow_network.add(new_u, new_v, mult * w * distance);
                    }
                endfor
            }

            // add the penalties
            for (size_t i = 0; i < left_region_size; ++i) {
                vertex_t u             = left_region[i];
                vertex_t new_u         = translation_table.get_n(u);
                weight_t left_penalty  = left_penalties[u];
                weight_t right_penalty = right_penalties[u];
                if (left_penalty > 0) { flow_network.add_t_edge(new_u, left_penalty); }
                if (right_penalty > 0) { flow_network.add_s_edge(new_u, right_penalty); }
            }
            for (size_t i = 0; i < right_region_size; ++i) {
                vertex_t u             = right_region[i];
                vertex_t new_u         = translation_table.get_n(u);
                weight_t left_penalty  = left_penalties[u];
                weight_t right_penalty = right_penalties[u];
                if (left_penalty > 0) { flow_network.add_t_edge(new_u, left_penalty); }
                if (right_penalty > 0) { flow_network.add_s_edge(new_u, right_penalty); }
            }
        }

        bool valid_cut(const graph_t &g,
                       p_manager_t &p_manager,
                       partition_t left_id,
                       partition_t right_id) {
            weight_t    left_weight  = p_manager.get_bweight(left_id);
            weight_t    right_weight = p_manager.get_bweight(right_id);
            for (size_t j            = 0; j < left_region_size; ++j) {
                vertex_t u        = left_region[j];
                weight_t u_weight = g.weight(u);
                vertex_t new_u    = translation_table.get_n(u);
                if (flow_network.get(new_u) == right_id) {
                    left_weight -= u_weight;
                    right_weight += u_weight;
                }
            }
            for (size_t j            = 0; j < right_region_size; ++j) {
                vertex_t u        = right_region[j];
                weight_t u_weight = g.weight(u);
                vertex_t new_u    = translation_table.get_n(u);
                if (flow_network.get(new_u) == left_id) {
                    right_weight -= u_weight;
                    left_weight += u_weight;
                }
            }

            return left_weight <= m_lmax && right_weight <= m_lmax;
        }

        bool cut_changes_partition(partition_t left_id,
                                   partition_t right_id) {
            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u     = left_region[j];
                vertex_t new_u = translation_table.get_n(u);
                if (flow_network.get(new_u) == right_id) {
                    return true;
                }
            }
            for (size_t j = 0; j < right_region_size; ++j) {
                vertex_t u     = right_region[j];
                vertex_t new_u = translation_table.get_n(u);
                if (flow_network.get(new_u) == left_id) {
                    return true;
                }
            }

            return false;
        }

        void change_boundary(const graph_t &g,
                             bv_manager_t &bv_manager,
                             p_manager_t &p_manager,
                             q_graph_t &q_graph,
                             partition_t left_id,
                             partition_t right_id) {
            for (size_t j = 0; j < left_region_size; ++j) {
                vertex_t u     = left_region[j];
                vertex_t new_u = translation_table.get_n(u);
                if (flow_network.get(new_u) == right_id) {
                    if (bv_manager.is_boundary(u)) {
                        bv_manager.move(g, p_manager, u, left_id, right_id);
                    } else {
                        bv_manager.add_new(g, p_manager, u, right_id);
                    }

                    q_graph.move(g, p_manager, u, left_id, right_id);
                    p_manager.move(u, g.weight(u), left_id, right_id);
                }
            }
            for (size_t j = 0; j < right_region_size; ++j) {
                vertex_t u     = right_region[j];
                vertex_t new_u = translation_table.get_n(u);
                if (flow_network.get(new_u) == left_id) {
                    if (bv_manager.is_boundary(u)) {
                        bv_manager.move(g, p_manager, u, right_id, left_id);
                    } else {
                        bv_manager.add_new(g, p_manager, u, left_id);
                    }

                    q_graph.move(g, p_manager, u, right_id, left_id);
                    p_manager.move(u, g.weight(u), right_id, left_id);
                }
            }
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
