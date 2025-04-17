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

#ifndef HEIPROMAP_HIERARCHY_AWARE_ILP_REFINEMENT_H
#define HEIPROMAP_HIERARCHY_AWARE_ILP_REFINEMENT_H
/*

#include <gurobi_c++.h>

#include "../../commons/aligned_array.h"
#include "../../commons/random_engine.h"
#include "../../commons/statistic_collector.h"
#include "../../commons/utils.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/qap.h"

namespace HeiProMap {
    class HierarchyAwareILPRefinementConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit HierarchyAwareILPRefinementConfiguration(const std::string& t_name) : ISerialRefinerConfiguration(t_name) {}
        u64 max_iteration  = 1;
        u64 max_n_vertices = 100;
    };

    class HierarchyAwareILPRefinement final : public ISerialRefiner {
        vertex_t m_n    = 0;
        vertex_t m_m    = 0;
        partition_t m_k = 0;
        f64 m_imbalance = 0.0;
        weight_t m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;

        AlignedArray<vertex_t> vertices;
        size_t vertices_size = 0;
        AlignedArray<weight_t> penalties;
        AlignedArray<GRBVar> vars;

        RandomEngine* random_engine                            = nullptr;
        const HierarchyAwareILPRefinementConfiguration* config = nullptr;
        StatisticCollector* m_stat_collector                   = nullptr;

    public:
        HierarchyAwareILPRefinement() = default;

        ~HierarchyAwareILPRefinement() override = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const f64 t_imbalance,
                        const weight_t t_lmax,
                        const std::vector<partition_t>& t_hierarchy,
                        const std::vector<weight_t>& t_distance,
                        RandomEngine& t_random_engine,
                        const ISerialRefinerConfiguration& i_config,
                        StatisticCollector& t_stat_collect) override {
            vertex_t t_n_64        = round_up_64(t_n);
            partition_t t_k_64     = round_up_64(t_k);
            partition_t t_k_t_k_64 = round_up_64(t_k * t_k);

            m_n         = t_n;
            m_m         = t_m;
            m_k         = t_k;
            m_lmax      = t_lmax;
            m_imbalance = t_imbalance;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;

            random_engine    = &t_random_engine;
            config           = dynamic_cast<const HierarchyAwareILPRefinementConfiguration*>(&i_config);
            m_stat_collector = &t_stat_collect;

            vertices = AlignedArray<vertex_t>(t_n);

            penalties = AlignedArray<weight_t>(t_n * t_k);
            vars      = AlignedArray<GRBVar>(t_n * t_k);
        }

        void refine(const u64 level,
                    const u64 max_level,
                    const graph_t& g,
                    const d_oracle_t& d_oracle,
                    bv_manager_t& bv_manager,
                    p_manager_t& p_manager,
                    q_graph_t& q_graph) override {
            u64 iteration = 0;
            while (iteration < config->max_iteration) {
                iteration += 1;
                for (size_t i = 0; i < m_hierarchy.size() - 1; ++i) {
                    refine_layer(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, m_hierarchy.size() - 1 - i);
                }
            }
        }

        void refine_layer(const u64 level,
                          const u64 max_level,
                          const graph_t& g,
                          const d_oracle_t& d_oracle,
                          bv_manager_t& bv_manager,
                          p_manager_t& p_manager,
                          q_graph_t& q_graph,
                          size_t layer) {
            partition_t n_total_super_blocks = 1;
            for (size_t i = layer; i < m_hierarchy.size(); ++i) { n_total_super_blocks *= m_hierarchy[i]; }
            partition_t ids_per_super_block = m_k / n_total_super_blocks;

            partition_t n_local_super_blocks = m_hierarchy[layer];
            partition_t n_upper_blocks       = 1;
            for (size_t i = layer + 1; i < m_hierarchy.size(); ++i) { n_upper_blocks *= m_hierarchy[i]; }

            partition_t neighborhood_id_start = 0;
            partition_t neighborhood_id_end   = ids_per_super_block * n_local_super_blocks;

            for (size_t upper_block_id = 0; upper_block_id < n_upper_blocks; ++upper_block_id) {
                refine_neighborhood(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, neighborhood_id_start, neighborhood_id_end, n_local_super_blocks, ids_per_super_block);

                neighborhood_id_start = neighborhood_id_end;
                neighborhood_id_end += ids_per_super_block * n_local_super_blocks;
            }
        }

        void refine_neighborhood(const u64 level,
                                 const u64 max_level,
                                 const graph_t& g,
                                 const d_oracle_t& d_oracle,
                                 bv_manager_t& bv_manager,
                                 p_manager_t& p_manager,
                                 q_graph_t& q_graph,
                                 partition_t neighborhood_id_start,
                                 partition_t neighborhood_id_end,
                                 partition_t n_local_super_blocks,
                                 partition_t ids_per_super_block) {
            weight_t blocks_lmax = (weight_t)ids_per_super_block * m_lmax;
            std::vector<weight_t> blocks_weights(n_local_super_blocks, 0);

            for (partition_t i = 0; i < n_local_super_blocks; ++i) {
                for (partition_t j = 0; j < ids_per_super_block; ++j) {
                    partition_t id = neighborhood_id_start + i * ids_per_super_block + j;
                    blocks_weights[i] += p_manager.get_bweight(id);
                }
            }

            for (vertex_t u = 0; u < g.get_n(); ++u) {
                partition_t u_id = p_manager[u];
                if (!IN_NEIGHBORING_BLOCK(u_id, neighborhood_id_start, neighborhood_id_end)) { continue; }


            }
        }

        JSONString get_stats() override {
            std::string stats = "{ \n";
            stats.pop_back();
            stats.pop_back();
            stats += "\n}";

            JSONString json_stats;
            json_stats.s = stats;
            return json_stats;
        }
    };
}

 */
#endif //HEIPROMAP_HIERARCHY_AWARE_ILP_REFINEMENT_H
