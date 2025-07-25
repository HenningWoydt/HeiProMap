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

#ifndef HEIPROMAP_SUITOR_ALGORITHM_H
#define HEIPROMAP_SUITOR_ALGORITHM_H

#include <vector>

#include "../../commons/definitions.h"
#include "../../commons/random_engine.h"
#include "../../serial/serial_definitions_1.h"
#include "../deep_definitions_1.h"
#include "../deep_definitions_2.h"
#include "../deep_definitions_3.h"

namespace HeiProMap {
    class SuitorMatcherConfiguration {
    public:
    };

    class SuitorMatcher {
        vertex_t m_n = 0;
        vertex_t m_m = 0;
        partition_t m_k = 0;
        weight_t m_l_max = 0;
        u64 m_threads = 1;

        const SuitorMatcherConfiguration *config = nullptr;
        RandomEngine *random_engine = nullptr;

    public:
        void initialize(const vertex_t t_n,
                        const vertex_t t_m,
                        const partition_t t_k,
                        const weight_t t_l_max,
                        const u64 t_threads,
                        RandomEngine &t_random_engine,
                        const SuitorMatcherConfiguration &i_config) {
            m_n = t_n;
            m_m = t_m;
            m_k = t_k;
            m_l_max = t_l_max;
            m_threads = t_threads;

            config = dynamic_cast<const SuitorMatcherConfiguration *>(&i_config);
            random_engine = &t_random_engine;
        }

        void match(const size_t level,
                   const deep_graph_t &g,
                   const deep_p_manager_t &p_manager,
                   Matching &matching) {
            // std::vector<vertex_t> partner(g.get_n(), g.get_n());
            std::vector<vertex_t> suitor(g.get_n(), g.get_n());
            std::vector<weight_t> ws(g.get_n(), 0);

            for (vertex_t u = 0; u < g.get_n(); u++) {
                vertex_t current = u;
                bool done = false;

                while (!done) {
                    vertex_t partner = suitor[current];
                    weight_t heaviest = ws[current];

                    forall_guivw(g, current, i, v, w) {
                            if (w > heaviest && w > ws[v]) {
                                partner = v;
                                heaviest = w;
                            }
                        }
                    endfor
                    done = true;

                    if (heaviest > 0) {
                        vertex_t y = suitor[partner];
                        suitor[partner] = current;
                        ws[partner] = heaviest;
                        if (y != g.get_n()) {
                            current = y;
                            done = false;
                        }
                    }
                }
            }

            for (vertex_t u = 0; u < g.get_n(); u++) {
                vertex_t u_partner = suitor[u];
                if (u < u_partner) { continue; }
                if (u_partner != g.get_n()) {
                    vertex_t v_partner = suitor[u_partner];
                    if (v_partner != g.get_n()) {
                        if (v_partner == u) {
                            // std::cout << "Add " << u << " to " << u_partner << std::endl;
                            matching.add(u, u_partner);
                        }
                    }
                }
            }
        }

        void match_rec(const size_t level,
                       const deep_graph_t &g,
                       const deep_p_manager_t &p_manager,
                       Matching &matching) {
            std::vector<vertex_t> partner(g.get_n(), g.get_n());
            std::vector<vertex_t> suitor(g.get_n(), g.get_n());

            std::vector<vertex_t> hits(g.get_n(), 0);

            for (vertex_t u = 0; u < g.get_n(); u++) {
                find_suitor_rec(0, u, g, partner, suitor);
            }

            for (vertex_t u = 0; u < g.get_n(); u++) {
                vertex_t u_partner = partner[u];
                if (u < u_partner) { continue; }
                if (u_partner != g.get_n()) {
                    vertex_t v_partner = partner[u_partner];
                    if (v_partner != g.get_n()) {
                        if (v_partner == u) {
                            // std::cout << "Add " << u << " to " << u_partner << std::endl;
                            matching.add(u, u_partner);

                            hits[u] += 1;
                            hits[u_partner] += 1;
                        }
                    }
                }
            }

            for (vertex_t u = 0; u < g.get_n(); u++) {
                if (hits[u] > 1) {
                    std::cout << "hit " << u << " " << hits[u] << "times" << std::endl;
                }
            }
        }

        struct Edge {
            vertex_t u;
            vertex_t v;
            weight_t w;

            Edge(vertex_t t_u, vertex_t t_v, weight_t t_w) {
                u = t_u;
                v = t_v;
                w = t_w;
            }

            bool operator<(const Edge &other) const {
                return w < other.w || (w == other.w && v < other.v);
            }
        };

        void find_suitor_rec(const u64 depth,
                             const vertex_t u,
                             const deep_graph_t &g,
                             std::vector<vertex_t> &partner,
                             std::vector<vertex_t> &suitor) {
            std::cout << depth << " ";
            for (size_t i = 0; i < depth; ++i) {
                std::cout << " ";
            }
            std::cout << "find_suitor " << u << " " << g.size(u) << std::endl;

            Edge best_edge(u, g.get_n(), 0);
            forall_guivw(g, u, i, v, w) {
                    Edge e(u, v, w);
                    Edge suitor_e(u, suitor[v], suitor[v] == g.get_n() ? 0 : g.weight(suitor[v]));

                    if (suitor_e < e) {
                        if (best_edge < e) {
                            best_edge = e;
                        }
                    }
                }
            endfor
            if (best_edge.v != g.get_n()) {
                partner[u] = best_edge.v;
            }

            if (partner[u] != g.get_n()) {
                vertex_t y = suitor[partner[u]];
                suitor[partner[u]] = u;

                if (y != g.get_n()) {
                    find_suitor_rec(depth + 1, y, g, partner, suitor);
                }
            }
        }
    };
}

#endif //HEIPROMAP_SUITOR_ALGORITHM_H
