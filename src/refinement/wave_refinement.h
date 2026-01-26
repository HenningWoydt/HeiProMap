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

#ifndef HEIPROMAP_WAVE_REFINEMENT_H
#define HEIPROMAP_WAVE_REFINEMENT_H

#include "ISerialRefiner.h"
#include "../definitions.h"
#include "../utility/JSON_utils.h"
#include "../utility/random_engine.h"
#include "../utility/qap.h"

namespace HeiProMap {
    class WaveRefinementConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit WaveRefinementConfiguration(const std::string &t_name) : ISerialRefinerConfiguration(t_name) {
        }

        u64 max_iteration = 25; // how many iterations to run the algorithm at most
    };

    class WaveRefinement final : public ISerialRefiner {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
        f64 m_imbalance = 0.0;
        weight_t m_lmax = 0;
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

        RandomEngine *random_engine = nullptr;
        const WaveRefinementConfiguration *config = nullptr;

        std::vector<KWayFMMove> global_moves;
        std::vector<KWayFMMove> moves;
        std::vector<KWayFMMove> possible_pos_moves;
        std::vector<KWayFMMove> possible_neg_moves;

    public:
        WaveRefinement() = default;

        ~WaveRefinement() override = default;

        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const f64 t_imbalance,
                        const weight_t t_lmax,
                        const std::vector<partition_t> &t_hierarchy,
                        const std::vector<weight_t> &t_distance,
                        RandomEngine &t_random_engine,
                        const ISerialRefinerConfiguration &i_config) override {
            m_n = t_n;
            m_m = t_m;
            m_k = t_k;
            m_imbalance = t_imbalance;
            m_lmax = t_lmax;
            m_hierarchy = t_hierarchy;
            m_distance = t_distance;

            random_engine = &t_random_engine;
            config = dynamic_cast<const WaveRefinementConfiguration *>(&i_config);

            vertex_used.initialize(m_n, 0);
            block_used.initialize(m_m, 0);

            curr_boundary.initialize(m_n);
            curr_boundary_size = 0;

            blocks.initialize(m_k);
            blocks_qap_delta.initialize(m_k);
            blocks_size = 0;
        }

        void refine([[maybe_unused]] const u64 level,
                    [[maybe_unused]] const u64 max_level,
                    graph_t &g,
                    d_oracle_t &d_oracle,
                    bv_manager_t &bv_manager,
                    p_manager_t &p_manager,
                    q_graph_t &q_graph) override {
            ScopedTimer _t("refinement", "WaveRefinement", "refine");

            u64 max_n_waves = 10000;
            u64 max_n_repetitions = 1;
            u64 max_n_add_nodes = 100;

            for (size_t wave_i = 0; wave_i < max_n_waves; ++wave_i) {
                partition_t u_id = random_engine->get_u32() % m_k;
                partition_t v_id = random_engine->get_u32() % m_k;

                if (u_id == v_id || !q_graph.has_edge(u_id, v_id)) {
                    wave_i -= 1;
                    continue;
                }

                s64 global_best_qap_delta = -1;
                global_moves.clear();

                for (size_t rep = 0; rep < max_n_repetitions; ++rep) {
                    s64 moves_qap_sum = 0;
                    s64 moves_best_qap_sum = 0;
                    size_t moves_best_idx = 0;
                    moves.clear();

                    for (size_t add_i = 0; add_i < max_n_add_nodes; ++add_i) {
                        s64 sum_qap_delta = 0;
                        size_t n_pos_moves = 0;

                        if (p_manager.get_bweight(v_id) == m_lmax) { continue; }

                        // collect all possible moves
                        possible_pos_moves.clear();
                        possible_neg_moves.clear();
                        forall_bv_id_iu(bv_manager, u_id, i, u) {
                                if (g.v_weights[u] + p_manager.get_bweight(v_id) > m_lmax) { continue; }
                                forall_guiv(g, u, j, v) {
                                        if (p_manager[v] != v_id) { continue; }

                                        s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);

                                        if (qap_delta > 0) {
                                            possible_pos_moves.push_back({u, u_id, v_id, qap_delta});
                                            sum_qap_delta += qap_delta;
                                            n_pos_moves += 1;
                                        } else {
                                            possible_neg_moves.push_back({u, u_id, v_id, qap_delta});
                                        }
                                        break;
                                    }
                                endfor
                            }
                        endfor

                        if (possible_pos_moves.size() + possible_neg_moves.size() == 0) { break; }

                        KWayFMMove m;
                        // randomly choose one move
                        if (n_pos_moves > 0) {
                            // choose a positive move
                            s64 threshold = random_engine->get_f32(0.0, 1.0) * sum_qap_delta;
                            s64 c = 0;
                            for (size_t i = 0; i < possible_pos_moves.size(); ++i) {
                                if (c <= threshold && c + possible_pos_moves[i].qap_delta > threshold) {
                                    m = possible_pos_moves[i];
                                    break;
                                }
                            }
                        } else {
                            // choose a negative move
                            s64 temp = -std::numeric_limits<s64>::max();
                            size_t idx = 0;
                            for (size_t i = 0; i < possible_neg_moves.size(); ++i) {
                                if (possible_neg_moves[i].qap_delta > temp) {
                                    temp = possible_neg_moves[i].qap_delta;
                                    idx = i;
                                }
                            }
                            m = possible_neg_moves[idx];
                        }

                        moves.push_back(m);
                        moves_qap_sum += m.qap_delta;
                        bv_manager.move(g, p_manager, m.u, u_id, v_id);
                        q_graph.move(g, p_manager, m.u, u_id, v_id);
                        p_manager.move(m.u, g.v_weights[m.u], u_id, v_id);

                        if (moves_qap_sum >= moves_best_qap_sum) {
                            moves_best_qap_sum = moves_qap_sum;
                            moves_best_idx = moves.size();
                        }
                    }

                    // revert all moves in partitioning manager
                    for (size_t i = 0; i < moves.size(); i++) {
                        vertex_t vertex = moves[moves.size() - 1 - i].u;
                        weight_t vertex_weight = g.v_weights[vertex];
                        partition_t vertex_id = p_manager[vertex];
                        partition_t move_id = u_id == vertex_id ? v_id : u_id;

                        bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                        q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                        p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                    }

                    if (moves_best_qap_sum > global_best_qap_delta) {
                        moves.resize(moves_best_idx);
                        global_moves = moves;
                        global_best_qap_delta = moves_best_qap_sum;
                    }
                }

                // make all moves to best index
                for (size_t i = 0; i < global_moves.size(); ++i) {
                    vertex_t vertex = global_moves[i].u;
                    weight_t vertex_weight = g.v_weights[vertex];
                    partition_t vertex_id = p_manager[vertex];
                    partition_t move_id = u_id == vertex_id ? v_id : u_id;

                    bv_manager.move(g, p_manager, vertex, vertex_id, move_id);
                    q_graph.move(g, p_manager, vertex, vertex_id, move_id);
                    p_manager.move(vertex, vertex_weight, vertex_id, move_id);
                }
            }
        }

        void refine_layer([[maybe_unused]] const u64 level,
                          [[maybe_unused]] const u64 max_level,
                          [[maybe_unused]] graph_t &g,
                          [[maybe_unused]] d_oracle_t &d_oracle,
                          [[maybe_unused]] bv_manager_t &bv_manager,
                          [[maybe_unused]] p_manager_t &p_manager,
                          [[maybe_unused]] q_graph_t &q_graph,
                          [[maybe_unused]] size_t layer) override {
        }
    };
}
#endif //HEIPROMAP_WAVE_REFINEMENT_H
