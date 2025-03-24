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
    class RandomEngine;

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
                    partition_t id_1 = pairs[i].id1;
                    partition_t id_2 = pairs[i].id2;
                    refine_blocks(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, id_1, id_2);
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
                           partition_t id_1,
                           partition_t id_2) {

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
