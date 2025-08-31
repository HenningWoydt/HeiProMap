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

#ifndef HEIPROMAP_DEEP_REBALANCER_H
#define HEIPROMAP_DEEP_REBALANCER_H

#include <cstdint>
#include <limits>
#include <vector>
#include <algorithm>
#include <chrono>
#include <iostream>

#include "../../serial/serial_definitions_1.h"
#include "../../serial/serial_definitions_2.h"
#include "../../serial/serial_definitions_3.h"
#include "../../serial/utility/qap.h"

namespace HeiProMap {
    class DeepRebalancer {
        u64 m_threads = 1;

        RandomEngine *random_engine = nullptr;

    public:
        void initialize(const u64 t_threads,
                        RandomEngine &t_random_engine) {
            m_threads = t_threads;

            random_engine = &t_random_engine;
        }

        void rebalance(const deep_graph_t &g,
                       deep_p_manager_t &p_manager,
                       deep_bv_manager_t &bv_manager,
                       deep_q_graph_t &q_graph,
                       deep_d_oracle_t &d_oracle,
                       partition_t k) {
            std::vector<vertex_t> boundary;

            while (true) {
                // collect all vertices
                boundary.clear();
                for (partition_t id = 0; id < k; ++id) {
                    if (p_manager.get_bweight(id) > p_manager.get_lmax(id)) {
                        forall_bv_id_iu(bv_manager, id, i, u) {
                                boundary.push_back(u);
                            }
                        endfor
                    }
                }

                // shuffle them
                std::shuffle(boundary.begin(), boundary.end(), random_engine->generator);

                bool made_move = false;
                for (vertex_t u: boundary) {
                    partition_t u_id = p_manager[u];
                    weight_t u_weight = g.weight(u);

                    if (p_manager.get_bweight(u_id) <= p_manager.get_lmax(u_id)) { continue; }

                    partition_t best_id = 0;
                    s64 best_qap_delta = -std::numeric_limits<s64>::max();
                    forall_guiv(g, u, j, v) {
                            partition_t v_id = p_manager[v];
                            if (p_manager.get_bweight(v_id) + u_weight > p_manager.get_lmax(v_id)) { continue; }

                            s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                            if (qap_delta > best_qap_delta) {
                                best_qap_delta = qap_delta;
                                best_id = v_id;
                            }
                        }
                    endfor

                    for (partition_t v_id = q_graph.lowest_level_neighborhood_start(u_id); v_id < q_graph.lowest_level_neighborhood_end(u_id); ++v_id) {
                        if (v_id == u_id) { continue; }
                        if (!p_manager.is_active(v_id)) { continue; }
                        if (p_manager.get_bweight(v_id) + u_weight > p_manager.get_lmax(v_id)) { continue; }

                        s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                        if (qap_delta > best_qap_delta) {
                            best_qap_delta = qap_delta;
                            best_id = v_id;
                        }
                    }


                    if (best_qap_delta != -std::numeric_limits<s64>::max()) {
                        bv_manager.move(g, p_manager, u, u_id, best_id);
                        q_graph.move(g, p_manager, u, u_id, best_id);
                        p_manager.move(u, u_weight, u_id, best_id);
                        made_move = true;
                    }
                }
                if (!made_move) {
                    break;
                }
            }
        }

        void rebalance_last_layer(const deep_graph_t &g,
                                  deep_p_manager_t &p_manager,
                                  deep_bv_manager_t &bv_manager,
                                  deep_q_graph_t &q_graph,
                                  deep_d_oracle_t &d_oracle,
                                  partition_t k) {
            rebalance(g, p_manager, bv_manager, q_graph, d_oracle, k);

            std::vector<vertex_t> boundary;

            bool global_move = true;
            while (global_move) {
                global_move = false;
                for (partition_t id = 0; id < k; ++id) {
                    if (p_manager.get_bweight(id) <= p_manager.get_lmax(id)) { continue; }

                    boundary.clear();
                    forall_bv_id_iu(bv_manager, id, i, u) {
                            boundary.push_back(u);
                        }
                    endfor

                    std::shuffle(boundary.begin(), boundary.end(), random_engine->generator);

                    for (vertex_t u: boundary) {
                        partition_t u_id = p_manager[u];
                        if (p_manager.get_bweight(u_id) <= p_manager.get_lmax(u_id)) { continue; }
                        partition_t best_id = 0;
                        s64 best_qap_delta = -std::numeric_limits<s64>::max();
                        weight_t u_weight = g.weight(u);
                        forall_guiv(g, u, j, v) {
                                partition_t v_id = p_manager[v];
                                if (v_id == id) { continue; }
                                if (p_manager.get_bweight(v_id) + u_weight > p_manager.get_lmax(v_id)) { continue; }

                                s64 qap_delta = get_u_qap_delta(g, u, id, v_id, p_manager, d_oracle);
                                if (qap_delta > best_qap_delta) {
                                    best_qap_delta = qap_delta;
                                    best_id = v_id;
                                }
                            }
                        endfor

                        for (partition_t v_id = q_graph.lowest_level_neighborhood_start(u_id); v_id < q_graph.lowest_level_neighborhood_end(u_id); ++v_id) {
                            if (v_id == u_id) { continue; }
                            if (v_id == id) { continue; }
                            if (p_manager.get_bweight(v_id) + u_weight > p_manager.get_lmax(v_id)) { continue; }

                            s64 qap_delta = get_u_qap_delta(g, u, id, v_id, p_manager, d_oracle);
                            if (qap_delta > best_qap_delta) {
                                best_qap_delta = qap_delta;
                                best_id = v_id;
                            }
                        }

                        if (best_qap_delta != -std::numeric_limits<s64>::max()) {
                            bv_manager.move(g, p_manager, u, u_id, best_id);
                            q_graph.move(g, p_manager, u, u_id, best_id);
                            p_manager.move(u, u_weight, u_id, best_id);
                            global_move = true;
                        }
                    }
                }
            }

            // desperate moves
            bool overloaded = true;
            while (overloaded) {
                overloaded = false;

                boundary.clear();
                for (partition_t id = 0; id < k; ++id) {
                    if (p_manager.get_bweight(id) > p_manager.get_lmax(id)) {
                        overloaded = true;
                        forall_bv_id_iu(bv_manager, id, i, u) {
                                boundary.push_back(u);
                            }
                        endfor
                    }
                }

                std::shuffle(boundary.begin(), boundary.end(), random_engine->generator);

                for (vertex_t u: boundary) {
                    partition_t u_id = p_manager[u];
                    if (p_manager.get_bweight(u_id) <= p_manager.get_lmax(u_id)) { continue; }

                    partition_t best_id = 0;
                    s64 best_qap_delta = -std::numeric_limits<s64>::max();
                    weight_t u_weight = g.weight(u);
                    size_t n_hits = 0;
                    for (size_t i = 0; i < k; ++i) {
                        partition_t v_id = ((u) * 1315423911ull ^ (u_id) * 2654435761ull ^ (i) * 97531ull) % k;
                        if (v_id == u_id) { continue; }
                        if (p_manager.get_bweight(v_id) + u_weight > p_manager.get_lmax(v_id)) { continue; }

                        n_hits += 1;
                        s64 qap_delta = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                        if (qap_delta > best_qap_delta) {
                            best_qap_delta = qap_delta;
                            best_id = v_id;
                        }
                        if (n_hits > 3) { break; }
                    }

                    if (best_qap_delta != -std::numeric_limits<s64>::max()) {
                        bv_manager.move(g, p_manager, u, u_id, best_id);
                        q_graph.move(g, p_manager, u, u_id, best_id);
                        p_manager.move(u, u_weight, u_id, best_id);
                    }
                }
            }
        }
    };
}

#endif //HEIPROMAP_DEEP_REBALANCER_H
