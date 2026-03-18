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

#ifndef HEIPROMAP_CYCLE_BASED_REFINEMENT_H
#define HEIPROMAP_CYCLE_BASED_REFINEMENT_H

namespace HeiProMap {
    struct Edge {
        vertex_t v;
        weight_t w;
    };

    struct CycleMove {
        vertex_t u;
        partition_t from_id;
        partition_t to_id;
        weight_t delta; // negative is good, positive is bad
    };

    class DirectedGraph {
    public:
        vertex_t n = 0;
        vertex_t m = 0;

        std::vector<std::vector<Edge> > edges;

        DirectedGraph() = default;

        DirectedGraph(vertex_t n) : n(n), m(0) {
            edges.resize(n);
        }

        void initialize(vertex_t t_n) {
            n = t_n;
            m = 0;
            edges.resize(n);
            for (vertex_t u = 0; u < n; ++u) {
                edges[u].clear();
            }
        }

        void add(vertex_t u, vertex_t v, weight_t w) {
            edges[u].push_back({v, w});
            m += 1;
        }
    };

    inline std::vector<vertex_t> find_negative_weight_cycle(DirectedGraph &dir_g) {
        std::vector<vertex_t> cycle;

        // Start all distances at 0 so we detect a negative cycle in any component.
        std::vector<weight_t> dist(dir_g.n, 0);
        std::vector<vertex_t> parent(dir_g.n, static_cast<vertex_t>(-1));

        vertex_t updated = static_cast<vertex_t>(-1);

        // Bellman-Ford: if something is still relaxed in the nth iteration,
        // then a negative cycle exists.
        for (vertex_t iter = 0; iter < dir_g.n; ++iter) {
            updated = static_cast<vertex_t>(-1);

            for (vertex_t u = 0; u < dir_g.n; ++u) {
                for (const Edge &e: dir_g.edges[u]) {
                    vertex_t v = e.v;
                    weight_t w = e.w;

                    if (dist[u] + w < dist[v]) {
                        dist[v] = dist[u] + w;
                        parent[v] = u;
                        updated = v;
                    }
                }
            }

            if (updated == static_cast<vertex_t>(-1)) {
                // No negative cycle.
                cycle.clear();
                return cycle;
            }
        }

        if (updated == static_cast<vertex_t>(-1)) {
            // No negative cycle.
            cycle.clear();
            return cycle;
        }

        // Move 'updated' inside the cycle.
        vertex_t x = updated;
        for (vertex_t i = 0; i < dir_g.n; ++i) {
            x = parent[x];
        }

        // Recover the cycle by following parents until we return to x.
        cycle.clear();
        vertex_t cur = x;

        cycle.push_back(cur);
        cur = parent[cur];
        while (cur != x) {
            cycle.push_back(cur);
            cur = parent[cur];
        }

        std::reverse(cycle.begin(), cycle.end());

        return cycle;
    }

    inline void cycle_refine(graph_t &g,
                             d_oracle_t &d_oracle,
                             bv_manager_t &bv_manager,
                             p_manager_t &p_manager,
                             q_graph_t &q_graph,
                             block_conn_t &block_conn,
                             f64 imbalance) {
        weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) p_manager.k));
        RandomEngine random_engine = RandomEngine(0);

        for (size_t iteration = 0; iteration < 10; ++iteration) {
            // get the current boundary
            std::vector<vertex_t> curr_boundary;
            //
            {
                ScopedTimer _t("refinement", "CycleRefinement", "get_boundary");

                for (partition_t id = 0; id < p_manager.k; ++id) {
                    forall_bv_id_iu(bv_manager, id, i, u)
                        {
                            curr_boundary.push_back(u);
                        }
                    }
                endfor
                fast_shuffle_unchecked(curr_boundary.data(), curr_boundary.data() + curr_boundary.size(), random_engine.generator);
            }

            // choose an independent set
            std::vector<vertex_t> independent_set;
            std::vector<u8> in_independent_set(g.n, 0);
            //
            {
                ScopedTimer _t("refinement", "CycleRefinement", "independent_set");

                for (auto u: curr_boundary) {
                    bool all_free = true;
                    forall_guiv(g, u, i, v)
                        {
                            all_free &= in_independent_set[v] == 0;
                        }
                    endfor
                    if (all_free) {
                        independent_set.push_back(u);
                        in_independent_set[u] = 1;
                    }
                }
            }

            // for each vertex in the independent set determine the gains of moving it to neighboring blocks
            size_t n_moves = 0;
            std::vector<std::vector<std::vector<CycleMove> > > moves;
            std::vector<CycleMove> flat_moves;
            std::vector<std::vector<std::vector<vertex_t> > > move_id;

            //
            {
                ScopedTimer _t("refinement", "CycleRefinement", "determine_moves");

                moves.resize(p_manager.k);
                move_id.resize(p_manager.k);
                for (partition_t id = 0; id < p_manager.k; ++id) {
                    moves[id].resize(p_manager.k);
                    move_id[id].resize(p_manager.k);
                }

                for (auto u: independent_set) {
                    partition_t u_id = p_manager[u];
                    forall_bc_ui_id(block_conn, u, i, id)
                        {
                            if (id == u_id) { continue; }

                            weight_t delta = get_u_qap_delta(g, u, u_id, id, p_manager, d_oracle, block_conn);

                            // if (delta <= -10) { continue; }

                            moves[u_id][id].push_back({u, u_id, id, -delta});
                            flat_moves.push_back({u, u_id, id, -delta});
                            move_id[u_id][id].push_back(n_moves);
                            n_moves += 1;
                        }
                    endfor
                }
            }

            // create the graph
            DirectedGraph dir_g;
            //
            {
                ScopedTimer _t("refinement", "CycleRefinement", "build_graph");

                dir_g.initialize(n_moves);

                // Insert edges only between buckets (a -> b) and (b -> c).
                for (partition_t a = 0; a < p_manager.k; ++a) {
                    for (partition_t b = 0; b < p_manager.k; ++b) {
                        if (a == b) { continue; }

                        for (partition_t c = 0; c < p_manager.k; ++c) {
                            if (c == b) { continue; }
                            if (a == c) { continue; }

                            for (size_t i = 0; i < moves[a][b].size(); ++i) {
                                CycleMove &move_1 = moves[a][b][i];
                                vertex_t id_1 = move_id[a][b][i];

                                for (size_t j = 0; j < moves[b][c].size(); ++j) {
                                    CycleMove &move_2 = moves[b][c][j];
                                    vertex_t id_2 = move_id[b][c][j];

                                    // don't connect moves that move the same vertex
                                    if (move_1.u == move_2.u) { continue; }

                                    // don't connect moves that would leave block b overloaded
                                    if (p_manager.get_bweight(b) + g.v_weights[move_1.u] - g.v_weights[move_2.u] > lmax) { continue; }

                                    // edge weight is the delta of the source move
                                    dir_g.add(id_1, id_2, move_1.delta);
                                }
                            }
                        }
                    }
                }
            }

            // find a negative cycle
            std::vector<vertex_t> cycle;
            //
            {
                ScopedTimer _t("refinement", "CycleRefinement", "find_cycle");

                /*
                auto chosen_cycles = select_best_disjoint_negative_2_cycles(g, p_manager, moves, move_id, lmax);

                if (!chosen_cycles.empty()) {
                    for (const auto &cy: chosen_cycles) {
                        const CycleMove &m1 = flat_moves[cy.id1];
                        const CycleMove &m2 = flat_moves[cy.id2];

                        for (const CycleMove *pmove: {&m1, &m2}) {
                            const CycleMove &move = *pmove;

                            vertex_t u = move.u;
                            weight_t u_w = g.v_weights[u];
                            partition_t u_id = p_manager[u];
                            partition_t new_id = move.to_id;

                            bv_manager.move(g, p_manager, u, u_id, new_id);
                            q_graph.move(g, p_manager, u, u_id, new_id);
                            block_conn.move(g, u, u_id, new_id);
                            p_manager.move(u, u_w, u_id, new_id);
                        }
                    }

                    continue;
                }
                */

                cycle = find_negative_weight_cycle(dir_g);
            }

            std::cout << curr_boundary.size() << " " << independent_set.size() << " " << dir_g.n << " " << dir_g.m << std::endl;
            std::cout << "Found: ";
            print(cycle);

            if (cycle.empty()) {
                // no cycle found
                return;
            }

            weight_t old_qap = get_qap(g, p_manager, d_oracle);
            std::cout << "A: " << old_qap << " " << p_manager.is_overloaded(lmax) << std::endl;

            // make all moves in the cycle
            {
                ScopedTimer _t("refinement", "CycleRefinement", "make_moves");

                for (vertex_t dir_u: cycle) {
                    CycleMove move = flat_moves[dir_u];

                    vertex_t u = move.u;
                    weight_t u_w = g.v_weights[u];
                    partition_t u_id = p_manager[u];
                    partition_t new_id = move.to_id;

                    bv_manager.move(g, p_manager, u, u_id, new_id);
                    q_graph.move(g, p_manager, u, u_id, new_id);
                    block_conn.move(g, u, u_id, new_id);
                    p_manager.move(u, u_w, u_id, new_id);
                }
            }

            weight_t new_qap = get_qap(g, p_manager, d_oracle);
            std::cout << "B: " << new_qap << " " << p_manager.is_overloaded(lmax) << " " << old_qap - new_qap << std::endl;
        }
    }
}


#endif //HEIPROMAP_CYCLE_BASED_REFINEMENT_H
