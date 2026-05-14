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

#ifndef HEIPROMAP_ADAPTIVE_SOLVER_H
#define HEIPROMAP_ADAPTIVE_SOLVER_H

#include <iostream>
#include "solver.h"
#include "distance_oracle.h"
#include "../utility/algorithm_configuration.h"
#include "../partitioning/kaffpa_partitioner.h"
#include "../partitioning/kway_partitioner/kway_core.h"
#include "../utility/utils.h"
#include "../utility/qap.h"
#include "../refinement/flow_based_refinement.h"

namespace HeiProMap {

    class AdaptiveSolver {
    private:
        AlgorithmConfiguration m_ac;
        PartitionManager m_p_manager;

        void recursive_solve(graph_t& g, PartitionManager& p_manager, std::vector<partition_t> hierarchy,
                             std::vector<weight_t> distance, u64 current_level, u64 offset, const TranslationTable<vertex_t>& tt,
                             const weight_t total_weight) {
            partition_t k_of_subgraph = prod<partition_t>(hierarchy);
            f64 total_remaining_slack = ((1.0 + m_ac.imbalance) * (f64) k_of_subgraph * (f64) total_weight) / ((f64) m_ac.k * (f64) g.g_weight) - 1.0;
            total_remaining_slack = std::max(0.0, total_remaining_slack);

            if (current_level >= m_ac.hm_level) {
                AlgorithmConfiguration sub_ac = m_ac;
                sub_ac.hierarchy = hierarchy;
                sub_ac.distance = distance;
                sub_ac.k = k_of_subgraph;
                sub_ac.imbalance = total_remaining_slack;

                if (m_ac.get("--config") == "fast") {
                    sub_ac.set_fast();
                } else if (m_ac.get("--config") == "eco") {
                    sub_ac.set_eco();
                } else if (m_ac.get("--config") == "strong") {
                    sub_ac.set_strong();
                } else if (m_ac.get("--config") == "super-strong") {
                    sub_ac.set_super_strong();
                }
                
                std::vector<weight_t> v_weights(g.n);
                const weight_t* v_weights_ptr = g.v_weights.get_ptr();
                std::copy(v_weights_ptr, v_weights_ptr + g.n, v_weights.begin());
                Solver sub_solver(std::move(g), sub_ac);
                const PartitionManager& sub_p_manager = sub_solver.solve_subproblem();

                for (vertex_t u = 0; u < sub_p_manager.n; ++u) {
                    p_manager.set(tt.get_o(u), v_weights[u], offset + sub_p_manager[u]);
                }
                return;
            }

            partition_t k = hierarchy.back();
            hierarchy.pop_back();
            distance.pop_back();

            f64 per_level_epsilon = std::pow(1.0 + total_remaining_slack, 1.0 / (f64) (hierarchy.size() + 1)) - 1.0;
            per_level_epsilon = std::max(0.0, per_level_epsilon);

            AlignedArray<partition_t> partition;
            partition.initialize(g.n, 0);
            //
            {
                ScopedTimer _t("adaptive_solver", "adaptive_solver", "partition");
                if (m_ac.global_multisection_config.mode == GLOBAL_MULTISECTION_KAFFPA_STRONG) {
                    kaffpa_partition(g, k, per_level_epsilon, KAFFPA_PARTITION_STRONG, m_ac.seed, partition, m_ac.global_multisection_config.kappa);
                } else if (m_ac.global_multisection_config.mode == GLOBAL_MULTISECTION_KAFFPA_ECO) {
                    kaffpa_partition(g, k, per_level_epsilon, KAFFPA_PARTITION_ECO, m_ac.seed, partition, m_ac.global_multisection_config.kappa);
                } else {
                    kaffpa_partition(g, k, per_level_epsilon, KAFFPA_PARTITION_FAST, m_ac.seed, partition, m_ac.global_multisection_config.kappa);
                }
            }

            std::vector<vertex_t> new_ns(k, 0);
            std::vector<vertex_t> new_ms(k, 0);
            std::vector<weight_t> new_ws(k, 0);
            for (vertex_t u = 0; u < g.n; ++u) {
                partition_t u_id = partition[u];
                new_ns[u_id] += 1;
                new_ws[u_id] += g.v_weights[u];
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    vertex_t v = g.edges_v[i];
                    partition_t v_id = partition[v];
                    if (v_id == u_id) { new_ms[u_id] += 1; }
                }
            }

            partition_t k_per_subgraph = prod<partition_t>(hierarchy);

            for (partition_t i = 0; i < k; ++i) {
                ScopedTimer _t("adaptive_solver", "adaptive_solver", "extrcact_graph");
                graph_t sub_g(new_ns[i], new_ms[i], new_ws[i]);
                TranslationTable<vertex_t> sub_tt;
                sub_tt.reserve(new_ns[i], p_manager.n);

                std::vector<vertex_t> new_us(k, 0);
                for (vertex_t old_u = 0; old_u < g.n; ++old_u) {
                    if (partition[old_u] == i) {
                        sub_tt.add(tt.get_o(old_u), new_us[i]);
                        new_us[i] += 1;
                    }
                }

                std::vector<vertex_t> degrees(sub_g.n, 0);
                for (vertex_t old_u = 0; old_u < g.n; ++old_u) {
                    if (partition[old_u] == i) {
                        vertex_t new_u = sub_tt.get_n(tt.get_o(old_u));
                        for (size_t j = g.neighborhoods[old_u]; j < g.neighborhoods[old_u + 1]; ++j) {
                            vertex_t old_v = g.edges_v[j];
                            if (partition[old_v] == i) {
                                degrees[new_u]++;
                            }
                        }
                    }
                }
                
                sub_g.neighborhoods[0] = 0;
                for(vertex_t j = 0; j < sub_g.n; ++j) {
                    sub_g.neighborhoods[j+1] = sub_g.neighborhoods[j] + degrees[j];
                }

                std::vector<vertex_t> cursor(sub_g.n, 0);
                for (vertex_t old_u = 0; old_u < g.n; ++old_u) {
                    if (partition[old_u] == i) {
                        vertex_t new_u = sub_tt.get_n(tt.get_o(old_u));
                        sub_g.v_weights[new_u] = g.v_weights[old_u];

                        for (size_t j = g.neighborhoods[old_u]; j < g.neighborhoods[old_u + 1]; ++j) {
                            vertex_t old_v = g.edges_v[j];
                            if (partition[old_v] == i) {
                                vertex_t new_v = sub_tt.get_n(tt.get_o(old_v));
                                size_t pos = sub_g.neighborhoods[new_u] + cursor[new_u];
                                sub_g.edges_v[pos] = new_v;
                                sub_g.edges_w[pos] = g.edges_w[j];
                                cursor[new_u]++;
                            }
                        }
                    }
                }
                _t.stop();
                recursive_solve(sub_g, p_manager, hierarchy, distance, current_level + 1, offset + i * k_per_subgraph, sub_tt, total_weight);
            }
        }

    public:
        explicit AdaptiveSolver(const AlgorithmConfiguration& ac) : m_ac(ac) {}

        void solve() {
            const auto sp = std::chrono::high_resolution_clock::now();
            graph_t g(m_ac.graph_in);
            const weight_t total_weight = g.g_weight;
            m_p_manager.initialize(g.n, m_ac.k, g.g_weight);
            m_p_manager.reset_weights();
            TranslationTable<vertex_t> tt;
            tt.reserve(g.n, g.n);
            for(vertex_t u = 0; u < g.n; ++u) {
                tt.add(u, u);
            }

            recursive_solve(g, m_p_manager, m_ac.hierarchy, m_ac.distance, 0, 0, tt, total_weight);

            std::cout << m_ac.get("--config") << " " << m_ac.hm_level << std::endl;

            if (m_ac.hm_level > 0 && m_ac.get("--config") == "super-strong") {
                ScopedTimer _t_final_refine("final_refinement", "AdaptiveSolver", "refine");
                
                BoundaryVertexManager bv_manager;
                bv_manager.initialize(g.n, m_ac.k);
                QuotientGraph q_graph;
                q_graph.initialize(m_ac.k);
                BlockConn block_conn;
                block_conn.initialize(g.n, g.m, m_ac.k);
                DistanceOracle d_oracle;
                d_oracle.initialize(m_ac.hierarchy, m_ac.distance);

                // build data structures from scratch
                for (vertex_t u = 0; u < g.n; ++u) {
                    partition_t u_id = m_p_manager[u];
                    for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                        vertex_t v = g.edges_v[i];
                        partition_t v_id = m_p_manager[v];
                        if (u_id != v_id) {
                            bv_manager.add(u, u_id);
                            if (u < v) {
                                q_graph.add_edge(u_id, v_id, g.edges_w[i]);
                            }
                        }
                    }
                }
                block_conn.compute_from_scratch(g, m_p_manager);

                std::cout << "Before " << get_qap(g,  m_p_manager, d_oracle) << std::endl;

                FlowBasedRefinement final_refiner;
                final_refiner.initialize(g.n, g.m, m_ac.k, m_ac.threads, m_ac.seed, m_ac.flow_based_refinement_config);
                final_refiner.refine(g, d_oracle, bv_manager, m_p_manager, q_graph, block_conn, m_ac.imbalance, g.uniform_v_weights, g.uniform_e_weights);

                std::cout << "After " << get_qap(g,  m_p_manager, d_oracle) << std::endl;
            }

            std::vector<partition_t> p(g.n);
            for (vertex_t u = 0; u < g.n; ++u) { p[u] = m_p_manager[u]; }
            write_partition(p, m_ac.mapping_out);

            const auto ep = std::chrono::high_resolution_clock::now();
            f64 duration = get_seconds(sp, ep);

            DistanceOracle d_oracle;
            d_oracle.initialize(m_ac.hierarchy, m_ac.distance);
            weight_t qap = get_qap(g, m_p_manager, d_oracle);

            weight_t lmax = std::ceil((1.0 + m_ac.imbalance) * ((f64) g.g_weight / (f64) m_ac.k));

            std::cout << "Total time        : " << duration << std::endl;
            std::cout << "#Nodes            : " << g.n << std::endl;
            std::cout << "#Edges            : " << g.m << std::endl;
            std::cout << "k                 : " << m_ac.k << std::endl;
            std::cout << "Lmax              : " << lmax << std::endl;
            std::cout << "Final QAP         : " << qap << std::endl;
            std::cout << "max block w       : " << m_p_manager.max_weight() << std::endl;
            std::cout << "#empty partitions : " << (u32) m_p_manager.n_empty_blocks() << std::endl;
            std::cout << "#oload partitions : " << (u32) m_p_manager.n_oload_blocks(lmax) << std::endl;
            std::cout << "Sum oload weights : " << m_p_manager.sum_oload_weight(lmax) << std::endl;
        }
    };

}

#endif //HEIPROMAP_ADAPTIVE_SOLVER_H
