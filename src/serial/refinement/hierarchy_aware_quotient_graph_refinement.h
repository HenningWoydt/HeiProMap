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

#ifndef HEIPROMAP_HIERARCHY_AWARE_QUOTIENT_GRAPH_REFINEMENT_H
#define HEIPROMAP_HIERARCHY_AWARE_QUOTIENT_GRAPH_REFINEMENT_H

#include <queue>

#include "../../commons/utils.h"
#include "../../commons/random_engine.h"
#include "../../commons/statistic_collector.h"
#include "../utility/functions.h"
#include "ISerialRefiner.h"
#include "../utility/qap.h"

namespace HeiProMap {
    class HierarchyAwareQuotientGraphRefinementConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit HierarchyAwareQuotientGraphRefinementConfiguration(const std::string &t_name) : ISerialRefinerConfiguration(t_name) {}
    };

    class HierarchyAwareQuotientGraphRefinement final : public ISerialRefiner {
    private:
        vertex_t                 m_n         = 0;
        vertex_t                 m_m         = 0;
        partition_t              m_k         = 0;
        f64                      m_imbalance = 0.0;
        weight_t                 m_lmax      = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t>    m_distance;

        AlignedArray<u32> vertex_used;
        u32               vertex_marker = 0;

        AlignedArray<u32> block_used;
        u32               block_marker = 0;

        std::priority_queue<KWayFMMove> heap;

        AlignedArray<Move> moves;
        size_t             moves_size = 0;

        RandomEngine                                             *random_engine    = nullptr;
        const HierarchyAwareQuotientGraphRefinementConfiguration *config           = nullptr;
        StatisticCollector                                       *m_stat_collector = nullptr;

    public:
        HierarchyAwareQuotientGraphRefinement() = default;

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
            m_imbalance = t_imbalance;
            m_lmax      = t_lmax;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;

            random_engine    = &t_random_engine;
            config           = dynamic_cast<const HierarchyAwareQuotientGraphRefinementConfiguration *>(&i_config);
            m_stat_collector = &t_stat_collect;

            vertex_used.initialize(m_n, 0);
            vertex_marker = 0;

            block_used.initialize(m_k, 0);
            block_marker = 0;

            moves.initialize(m_n);
            moves_size = 0;
        }

        void refine(const u64 level,
                    const u64 max_level,
                    graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph) override {

        }

        JSONString get_stats() override { return {}; };
    };
}

#endif //HEIPROMAP_HIERARCHY_AWARE_QUOTIENT_GRAPH_REFINEMENT_H
