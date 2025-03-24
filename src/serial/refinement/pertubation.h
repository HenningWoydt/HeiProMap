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

#ifndef HEIPROMAP_PERTUBATION_H
#define HEIPROMAP_PERTUBATION_H

#include <algorithm>

#include "../../commons/utils.h"
#include "../datastructures/functions.h"
#include "../interfaces/ISerialRefiner.h"
#include "../rebalance/k_way_rebalancer.h"
#include "../utility/qap.h"

namespace HeiProMap {
    class Pertubation {
        vertex_t m_n    = 0;
        vertex_t m_m    = 0;
        partition_t m_k = 0;
        weight_t m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;

        u32* vertex_used  = nullptr;
        u32 vertex_marker = 0;

        u32* block_used  = nullptr;
        u32 block_marker = 0;

        Move* moves       = nullptr;
        size_t moves_size = 0;

        RandomEngine* random_engine          = nullptr;
        StatisticCollector* m_stat_collector = nullptr;

    public:
        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const weight_t t_lmax,
                        const std::vector<partition_t>& t_hierarchy,
                        const std::vector<weight_t>& t_distance,
                        RandomEngine& t_random_engine,
                        StatisticCollector& t_stat_collect) {
            m_n         = t_n;
            m_m         = t_m;
            m_k         = t_k;
            m_lmax      = t_lmax;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;

            vertex_t t_n_64 = round_up_64(t_n);
            vertex_t t_k_64 = round_up_64(t_k);

            vertex_used = (u32*)aligned_alloc(64, sizeof(u32) * t_n_64);
            std::fill_n(vertex_used, t_n_64, 0);
            vertex_marker = 0;

            block_used = (u32*)aligned_alloc(64, t_k_64 * sizeof(u32));
            std::fill_n(block_used, t_k_64, 0);
            block_marker = 0;

            random_engine    = &t_random_engine;
            m_stat_collector = &t_stat_collect;
        }

        void pertubate(const u64 level,
                       const u64 max_level,
                       const graph_t& g,
                       const d_oracle_t& d_oracle,
                       bv_manager_t& bv_manager,
                       p_manager_t& p_manager,
                       q_graph_t& q_graph) {
            size_t max_iteration   = 1;
            float move_probability = 0.01f;

            std::vector<partition_t> ids(m_k);
            std::iota(ids.begin(), ids.end(), 0);

            for (size_t iteration = 0; iteration < max_iteration; ++iteration) {
                // iterate over all blocks in random order
                std::shuffle(ids.begin(), ids.end(), random_engine->gen);

                for (partition_t u_id : ids) {
                    std::vector<vertex_t> curr_boundary;
                    forall_bv_id_iu(bv_manager, u_id, i, u)
                        {
                            curr_boundary.push_back(u);
                        }
                    endfor

                    for (vertex_t u : curr_boundary) {
                        if (!bv_manager.is_boundary(u)) { continue; }
                        weight_t u_weight = g.weight(u);

                        forall_guiv(g, u, j, v)
                            {
                                partition_t v_id = p_manager[v];
                                if (v_id == u_id) { continue; }
                                if (p_manager.get_bweight(v_id) + u_weight > m_lmax) { continue; }

                                if (random_engine->get_f32() < move_probability) {
                                    bv_manager.move(g, p_manager, u, u_id, v_id);
                                    q_graph.move(g, p_manager, u, u_id, v_id);
                                    p_manager.move(u, u_weight, u_id, v_id);
                                    break;
                                }
                            }
                        endfor
                    }
                }
            }
        }
    };
}

#endif //HEIPROMAP_PERTUBATION_H
