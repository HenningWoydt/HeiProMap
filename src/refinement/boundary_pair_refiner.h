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

#ifndef HEIPROMAP_BOUNDARY_PAIR_REFINER_H
#define HEIPROMAP_BOUNDARY_PAIR_REFINER_H

namespace HeiProMap {
    inline void boundary_pair_refiner(graph_t &g,
                                      d_oracle_t &d_oracle,
                                      bv_manager_t &bv_manager,
                                      p_manager_t &p_manager,
                                      q_graph_t &q_graph,
                                      block_conn_t &block_conn,
                                      f64 imbalance) {
        weight_t X = -std::numeric_limits<weight_t>::max();
        std::vector<weight_t> gains;
        std::vector<partition_t> target;
        //
        {
            ScopedTimer _t("refinement", "BoundaryPairRefiner", "get_deltas");

            gains.resize(g.n, X);
            target.resize(g.n, p_manager.k);

            vertex_t n_vertices = 0;
            vertex_t n_neg_vertices = 0;
            for (partition_t u_id = 0; u_id < p_manager.k; u_id++) {
                forall_bv_id_iu(bv_manager, u_id, i, u)
                    {
                        size_t n_conns = 0;
                        partition_t last_id = u_id;
                        forall_bc_ui_id(block_conn, u, j, id)
                            {
                                if (u_id == id) { continue; }
                                n_conns += 1;
                                last_id = id;
                            }
                        endfor

                        if (n_conns == 1) {
                            n_vertices += 1;

                            weight_t delta = get_u_qap_delta(g, u, u_id, last_id, p_manager, d_oracle, block_conn);
                            if (delta < 0) {
                                n_neg_vertices += 1;
                                gains[u] = delta;
                                target[u] = last_id;
                            }
                        }
                    }
                endfor
            }

            std::cout << "Found (" << n_vertices << ", " << n_neg_vertices << ") of " << bv_manager.size() << " with exactly one other connected block" << std::endl;
        }

        //
        {
            ScopedTimer _t("refinement", "BoundaryPairRefiner", "find_pairs");

            vertex_t n_pairs = 0;
            vertex_t n_pos_pairs = 0;
            weight_t max_pos = 0;
            std::vector<vertex_t> vertices(2);
            for (partition_t u_id = 0; u_id < p_manager.k; u_id++) {
                forall_bv_id_iu(bv_manager, u_id, i, u)
                    {
                        if (gains[u] == X) { continue; } // not interesting vertex

                        vertices[0] = u;
                        forall_guiv(g, u, j, v)
                            {
                                partition_t v_id = p_manager[v];
                                partition_t t_id = target[u];

                                if (v_id != u_id) { continue; } // do not move if not in same block
                                if (v >= u) { continue; } // do not count double
                                if (gains[v] == X) { continue; } // not interesting vertex
                                if (target[u] != target[v]) { continue; } // not the same target

                                vertices[1] = v;
                                n_pairs += 1;

                                weight_t qap_before = get_qap(g, vertices, p_manager, d_oracle);

                                p_manager.move(u, g.v_weights[u], u_id, t_id);
                                p_manager.move(v, g.v_weights[v], v_id, t_id);

                                weight_t qap_after = get_qap(g, vertices, p_manager, d_oracle);

                                p_manager.move(u, g.v_weights[u], t_id, u_id);
                                p_manager.move(v, g.v_weights[v], t_id, v_id);

                                weight_t delta = qap_before - qap_after;
                                if (delta >= 0) {
                                    n_pos_pairs += 1;
                                    max_pos = std::max(max_pos, delta);
                                }
                            }
                        endfor
                    }
                endfor
            }

            std::cout << "Found " << n_pairs << " pairs, found " << n_pos_pairs << " pos-pairs with " << max_pos << std::endl;
        }
    }
}


#endif //HEIPROMAP_BOUNDARY_PAIR_REFINER_H
