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

#ifndef HEIPROMAP_DEEP_LIGHTNING_REFINEMENT_H
#define HEIPROMAP_DEEP_LIGHTNING_REFINEMENT_H

#include "ISerialDeepRefiner.h"
#include "../../../commons/definitions.h"
#include "../../../commons/JSON_utils.h"
#include "../../../commons/random_engine.h"
#include "../../../commons/statistic_collector.h"
#include "../../utility/qap.h"

namespace HeiProMap {
    class DeepLightningRefinementConfiguration final : public ISerialDeepRefinerConfiguration {
    public:
        explicit DeepLightningRefinementConfiguration(const std::string& t_name) : ISerialDeepRefinerConfiguration(t_name) {}
        u64 max_iteration = 100; // how many iterations to run the algorithm at most
        u64 max_recursion = 10;
    };

    class DeepLightningRefinement final : public ISerialDeepRefiner {
        vertex_t m_n    = 0;
        vertex_t m_m    = 0;
        partition_t m_k = 0;
        f64 m_imbalance = 0.0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;

        AlignedArray<u32> vertex_used;
        u32 vertex_marker = 0;

        AlignedArray<u32> block_used;
        u32 block_marker = 0;

        AlignedArray<vertex_t> curr_boundary;
        size_t curr_boundary_size = 0;

        AlignedArray<partition_t> blocks;
        AlignedArray<s64> blocks_qap_delta;
        size_t blocks_size = 0;

        RandomEngine* random_engine                        = nullptr;
        const DeepLightningRefinementConfiguration* config = nullptr;

        s64 current_qap_delta  = 0;
        size_t current_n_moves = 0;

    public:
        DeepLightningRefinement() = default;

        ~DeepLightningRefinement() override = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const f64 t_imbalance,
                        const std::vector<partition_t>& t_hierarchy,
                        const std::vector<weight_t>& t_distance,
                        RandomEngine& t_random_engine,
                        const ISerialDeepRefinerConfiguration& i_config) override {
            m_n         = t_n;
            m_m         = t_m;
            m_k         = t_k;
            m_imbalance = t_imbalance;
            m_hierarchy = t_hierarchy;
            m_distance  = t_distance;

            random_engine    = &t_random_engine;
            config           = dynamic_cast<const DeepLightningRefinementConfiguration*>(&i_config);

            vertex_used.initialize(m_n, 0);
            block_used.initialize(m_k, 0);

            curr_boundary.initialize(m_n);
            curr_boundary_size = 0;

            blocks.initialize(m_k);
            blocks_qap_delta.initialize(m_k);
            blocks_size = 0;
        }

        void refine(u64 level,
                    u64 max_level,
                    const graph_t& g,
                    deep_d_oracle_t& d_oracle,
                    deep_bv_manager_t& bv_manager,
                    deep_p_manager_t& p_manager,
                    deep_q_graph_t& q_graph) override {
            if (level + 3 < max_level) { return; }

            for (size_t i = 0; i < config->max_iteration; ++i) {
                std::vector<KWayFMMove> possible_moves;
                current_n_moves   = 0;
                current_qap_delta = 0;

                forall_bv_iu(bv_manager, i, u)
                    {
                        partition_t u_id  = p_manager[u];
                        weight_t u_weight = g.weight(u);

                        block_marker += 1;
                        forall_guiv(g, u, i, v)
                            {
                                partition_t v_id = p_manager[v];
                                if (u_id == v_id) { continue; }
                                if (block_used[v_id] == block_marker) { continue; }

                                if (p_manager.get_bweight(v_id) + u_weight > p_manager.get_lmax(v_id)) {
                                    // would overload
                                    s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                                    // if (qap_delta < 0) { continue; }

                                    possible_moves.push_back({u, u_id, v_id, qap_delta});
                                    block_used[v_id] = block_marker;
                                }
                            }
                        endfor
                    }
                endfor

                std::sort(possible_moves.begin(), possible_moves.end(), std::greater<KWayFMMove>());
                if (possible_moves.size() > 10) {
                    possible_moves.resize(10);
                }

                bool moves_found = false;
                for (auto& move : possible_moves) {
                    bv_manager.move(g, p_manager, move.u, move.u_id, move.to_move_id);
                    q_graph.move(g, p_manager, move.u, move.u_id, move.to_move_id);
                    p_manager.move(move.u, g.weight(move.u), move.u_id, move.to_move_id);

                    current_n_moves += 1;
                    current_qap_delta += move.qap_delta;

                    bool found = find_move(0, move.to_move_id, g, d_oracle, bv_manager, p_manager, q_graph);

                    if (found) {
                        moves_found = true;
                        break;
                    }

                    current_n_moves -= 1;
                    current_qap_delta -= move.qap_delta;

                    bv_manager.move(g, p_manager, move.u, move.to_move_id, move.u_id);
                    q_graph.move(g, p_manager, move.u, move.to_move_id, move.u_id);
                    p_manager.move(move.u, g.weight(move.u), move.to_move_id, move.u_id);
                }

                if (!moves_found) {
                    break;
                }
            }
        }

        bool find_move(size_t depth,
                       partition_t id,
                       const graph_t& g,
                       deep_d_oracle_t& d_oracle,
                       deep_bv_manager_t& bv_manager,
                       deep_p_manager_t& p_manager,
                       deep_q_graph_t& q_graph) {
            if (depth > config->max_recursion) {
                return false;
            }

            // from the overloaded block, either
            // positive move that balances block and not overloads non-used block -> stop
            // positive move that balances block and overloads non-used block -> continue
            // positive move that improves balances, but not overloads non-used block -> continue
            // if no move found than revert one step

            std::vector<KWayFMMove> possible_moves;
            KWayFMMove best_stop_move = {0, 0, std::numeric_limits<partition_t>::max(), -1};

            forall_bv_id_iu(bv_manager, id, i, u)
                {
                    weight_t u_weight = g.weight(u);

                    block_marker += 1;
                    forall_guiv(g, u, j, v)
                        {
                            partition_t v_id = p_manager[v];
                            if (id == v_id) { continue; }
                            if (block_used[v_id] == block_marker) { continue; }

                            s64 qap_delta = get_u_qap_delta(g, u, id, v_id, p_manager, d_oracle);

                            if (p_manager.get_bweight(id) - u_weight <= p_manager.get_lmax(id) && p_manager.get_bweight(v_id) + u_weight <= p_manager.get_lmax(v_id)) {
                                if (qap_delta > best_stop_move.qap_delta) {
                                    block_used[v_id] = block_marker;
                                    best_stop_move   = {u, id, v_id, qap_delta};
                                }
                            }

                            if (p_manager.get_bweight(id) - u_weight <= p_manager.get_lmax(id) && p_manager.get_bweight(v_id) + u_weight > p_manager.get_lmax(v_id)) {
                                block_used[v_id] = block_marker;
                                possible_moves.push_back({u, id, v_id, qap_delta});
                            }
                            if (p_manager.get_bweight(id) - u_weight > p_manager.get_lmax(id) && p_manager.get_bweight(v_id) + u_weight <= p_manager.get_lmax(v_id)) {
                                block_used[v_id] = block_marker;
                                possible_moves.push_back({u, id, v_id, qap_delta});
                            }
                        }
                    endfor
                }
            endfor

            std::sort(possible_moves.begin(), possible_moves.end(), std::greater<KWayFMMove>());

            if (possible_moves.size() > 3) {
                possible_moves.resize(3);
            }

            // std::cout << depth << " " << possible_moves.size() << " " << current_qap_delta << std::endl;

            for (auto& move : possible_moves) {
                bv_manager.move(g, p_manager, move.u, move.u_id, move.to_move_id);
                q_graph.move(g, p_manager, move.u, move.u_id, move.to_move_id);
                p_manager.move(move.u, g.weight(move.u), move.u_id, move.to_move_id);

                current_n_moves += 1;
                current_qap_delta += move.qap_delta;

                partition_t next_id = 0;
                // positive move that balances block and overloads another block -> continue
                if (p_manager.get_bweight(id) <= p_manager.get_lmax(id) && p_manager.get_bweight(move.to_move_id) > p_manager.get_lmax(move.to_move_id)) {
                    next_id = move.to_move_id;
                }

                // positive move that improves balances, but not overloads non-used block -> continue
                if (p_manager.get_bweight(id) > p_manager.get_lmax(id) && p_manager.get_bweight(move.to_move_id) <= p_manager.get_lmax(move.to_move_id)) {
                    next_id = id;
                }

                bool found = find_move(depth + 1, next_id, g, d_oracle, bv_manager, p_manager, q_graph);

                if (found) {
                    return true;
                }

                current_n_moves -= 1;
                current_qap_delta -= move.qap_delta;

                bv_manager.move(g, p_manager, move.u, move.to_move_id, move.u_id);
                q_graph.move(g, p_manager, move.u, move.to_move_id, move.u_id);
                p_manager.move(move.u, g.weight(move.u), move.to_move_id, move.u_id);
            }

            // positive move that balances block and not overloads non-used block -> stop
            if (best_stop_move.to_move_id != std::numeric_limits<partition_t>::max()) {
                if (best_stop_move.qap_delta + current_qap_delta < 0) { return false; }

                bv_manager.move(g, p_manager, best_stop_move.u, best_stop_move.u_id, best_stop_move.to_move_id);
                q_graph.move(g, p_manager, best_stop_move.u, best_stop_move.u_id, best_stop_move.to_move_id);
                p_manager.move(best_stop_move.u, g.weight(best_stop_move.u), best_stop_move.u_id, best_stop_move.to_move_id);

                current_n_moves += 1;
                current_qap_delta += best_stop_move.qap_delta;

                return true;
            }

            // if no move found than revert one step
            return false;
        }
    };
}

#endif //HEIPROMAP_DEEP_LIGHTNING_REFINEMENT_H
