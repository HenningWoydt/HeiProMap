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

#ifndef HEIPROMAP_LABEL_PROPAGATION_REFINEMENT_FARAJ20_H
#define HEIPROMAP_LABEL_PROPAGATION_REFINEMENT_FARAJ20_H

#include <random>

#include "../../definitions.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/qap.h"

namespace HeiProMap {
    class LabelPropagationFaraj20Configuration final : public ISerialRefinerConfiguration {
    public:
        explicit LabelPropagationFaraj20Configuration(const std::string &t_name) : ISerialRefinerConfiguration(t_name) {}
        u64 max_iteration = 25; // how many iterations to run the algorithm at most
    };

    /**
     * Executes label propagation refinement as described in
     * > Marcelo Fonseca Faraj, Alexander van der Grinten, Henning Meyerhenke, Jesper Larsson Träff, and Christian Schulz.
     * > High-quality Hierarchical Process Mapping.
     * > In 18th International Symposium on Experimental Algorithms, SEA 2020, June 16-18, 2020, Catania, Italy, volume 160 of LIPIcs, pages 4:1–4:15.
     */
     /*
    class LabelPropagationRefinementFaraj20 final : public ISerialRefiner {
        vertex_t m_n    = 0;
        vertex_t m_m    = 0;
        partition_t m_k = 0;
        weight_t m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;
        u64 m_seed = 0;

        u32* vertex_used  = nullptr;
        u32 vertex_marker = 0;

        u32* block_used  = nullptr;
        u32 block_marker = 0;

        RandomEngine* random_engine                        = nullptr;
        const LabelPropagationFaraj20Configuration* config = nullptr;
        StatisticCollector* m_stat_collector               = nullptr;

    public:
        LabelPropagationRefinementFaraj20() = default;

        ~LabelPropagationRefinementFaraj20() override {
            free(vertex_used);
            free(block_used);
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
            config           = dynamic_cast<const LabelPropagationFaraj20Configuration*>(&i_config);
            m_stat_collector = &t_stat_collect;

            vertex_t m_n_64 = round_up_64(m_n);
            vertex_used     = (u32*)aligned_alloc(64, m_n_64 * sizeof(u32));
            std::fill_n(vertex_used, m_n_64, vertex_marker);

            partition_t m_k_64 = round_up_64(m_k);
            block_used         = (u32*)aligned_alloc(64, m_k_64 * sizeof(u32));
            std::fill_n(block_used, m_k_64, block_marker);
        }

        void refine(const u64 level,
                    const graph_t& g,
                    const av_manager_t& av_manager,
                    const d_oracle_t& d_oracle,
                    bv_manager_t& bv_manager,
                    p_manager_t& p_manager,
                    q_graph_t& q_graph) override {
            bool move_occurred = true;
            for (u64 iteration = 0; iteration < config->max_iteration && move_occurred; ++iteration) {
                move_occurred = false;

                std::vector<vertex_t> curr_boundary;
                forall_bv_iu(bv_manager, i, u)
                    {
                        curr_boundary.push_back(u);
                    }
                endfor
                std::shuffle(curr_boundary.begin(), curr_boundary.end(), random_engine->gen);

                vertex_marker += 1;
                for (vertex_t u : curr_boundary) {
                    if (vertex_used[u] == vertex_marker) { continue; }
                    if (!bv_manager.is_boundary(u)) { continue; }

                    weight_t u_weight = g.get_weight(u);
                    partition_t u_id  = p_manager[u];

                    // make the move that reduces qap the most
                    partition_t best_id     = u_id;
                    weight_t best_id_weight = 0;
                    s64 best_qap_delta      = -1;
                    f32 counter             = 0;

                    block_marker += 1;
                    block_used[u_id] = block_marker;
                    forall_guiv(g, u, i, v)
                        {
                            partition_t v_id     = p_manager[v];
                            weight_t v_id_weight = p_manager.get_bweight(v_id);

                            if (block_used[v_id] != block_marker) {
                                if (v_id_weight + u_weight <= m_lmax) {
                                    s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);

                                    if (qap_delta > best_qap_delta || (qap_delta == best_qap_delta && v_id_weight < best_id_weight)) {
                                        best_id        = v_id;
                                        best_id_weight = v_id_weight;
                                        best_qap_delta = qap_delta;
                                        counter        = 1.0;
                                    } else if (qap_delta == best_qap_delta && qap_delta != -1) {
                                        counter += 1.0;
                                        // choose with probability 1/counter as it ensures uniform distribution
                                        if (random_engine->get_f32() < 1.0f / counter) {
                                            best_id = v_id;
                                        }
                                    }
                                }
                                block_used[v_id] = block_marker;
                            }
                        }
                    endfor

                    if (best_id != u_id) {
                        // choose if positive, if 0-gain choose 50% of the time
                        if (best_qap_delta > 0 || random_engine->get_f32() < 0.5) {
                            bv_manager.move(g, p_manager, u, u_id, best_id);
                            q_graph.move(g, p_manager, u, u_id, best_id);
                            p_manager.move(u, u_weight, u_id, best_id);
                            move_occurred = true;
                        }
                    }
                    vertex_used[u] = vertex_marker;
                }
            }
        }

        JSONString get_stats() override { return {}; };
    };
      */
}

#endif //HEIPROMAP_LABEL_PROPAGATION_REFINEMENT_FARAJ20_H
