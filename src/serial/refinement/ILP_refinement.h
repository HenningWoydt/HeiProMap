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

#ifndef HEIPROMAP_ILP_REFINEMENT_H
#define HEIPROMAP_ILP_REFINEMENT_H

/*
#include <gurobi_c++.h>

#include "../../commons/aligned_array.h"
#include "../../commons/random_engine.h"
#include "../../commons/statistic_collector.h"
#include "../../commons/utils.h"
#include "../interfaces/ISerialRefiner.h"
#include "../utility/qap.h"

namespace HeiProMap {
    class ILPRefinementConfiguration final : public ISerialRefinerConfiguration {
    public:
        explicit ILPRefinementConfiguration(const std::string& t_name) : ISerialRefinerConfiguration(t_name) {}
        u64 max_iteration  = 1;
        u64 max_n_vertices = 100;
    };

    class ILPRefinement final : public ISerialRefiner {
        vertex_t m_n    = 0;
        vertex_t m_m    = 0;
        partition_t m_k = 0;
        f64 m_imbalance = 0.0;
        weight_t m_lmax = 0;
        std::vector<partition_t> m_hierarchy;
        std::vector<weight_t> m_distance;

        // active block scheduling
        u8* active_this_round = nullptr;
        u8* active_next_round = nullptr;
        PairWeight* pairs     = nullptr;
        size_t pairs_size     = 0;

        AlignedArray<vertex_t> vertices;
        size_t vertices_size = 0;
        AlignedArray<weight_t> u_id_penalties;
        AlignedArray<weight_t> v_id_penalties;
        AlignedArray<GRBVar> vars;

        RandomEngine* random_engine              = nullptr;
        const ILPRefinementConfiguration* config = nullptr;
        StatisticCollector* m_stat_collector     = nullptr;

    public:
        ILPRefinement() = default;

        ~ILPRefinement() override = default;

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
            config           = dynamic_cast<const ILPRefinementConfiguration*>(&i_config);
            m_stat_collector = &t_stat_collect;

            // active block scheduling
            active_this_round = (u8*)aligned_alloc(64, t_k_64 * sizeof(u8));
            active_next_round = (u8*)aligned_alloc(64, t_k_64 * sizeof(u8));
            pairs             = (PairWeight*)aligned_alloc(64, t_k_t_k_64 * sizeof(PairWeight));
            pairs_size        = 0;

            vertices = AlignedArray<vertex_t>(t_n);

            u_id_penalties = AlignedArray<weight_t>(t_n);
            v_id_penalties = AlignedArray<weight_t>(t_n);
            vars           = AlignedArray<GRBVar>(t_n);
        }

        void refine(const u64 level,
                    const u64 max_level,
                    const graph_t& g,
                    const d_oracle_t& d_oracle,
                    bv_manager_t& bv_manager,
                    p_manager_t& p_manager,
                    q_graph_t& q_graph) override {
            std::fill_n(active_this_round, m_k, 1);
            std::fill_n(active_next_round, m_k, 0);

            u64 iteration = 0;
            while (iteration < config->max_iteration) {
                METRICS(iteration_time.back().push_back(0.0);)
                METRICS(iteration_time_get_pairs.back().push_back(0.0);)
                METRICS(iteration_time_initialize.back().push_back(0.0);)
                METRICS(iteration_time_queue.back().push_back(0.0);)
                METRICS(iteration_time_moves.back().push_back(0.0);)
                METRICS(iteration_qap_delta.back().push_back(0);)

                METRICS_TIME(sp)

                iteration += 1;

                METRICS_TIME(sp_get_pairs)

                // determine all pairs in the quotient graph
                pairs_size = 0;
                for (partition_t u_id = 0; u_id < m_k; ++u_id) {
                    for (partition_t v_id = u_id + 1; v_id < m_k; ++v_id) {
                        if (q_graph.has_edge(u_id, v_id) && (active_this_round[u_id] || active_this_round[v_id])) {
                            pairs[pairs_size++] = {u_id, v_id, d_oracle.get(u_id, v_id)};
                        }
                    }
                }
                if (pairs_size == 0) { return; }
                std::sort(pairs, pairs + pairs_size, std::greater<>());

                METRICS_TIME(ep_get_pairs)
                METRICS(iteration_time_get_pairs.back().back() += get_seconds(sp_get_pairs, ep_get_pairs);)

                for (size_t j = 0; j < pairs_size; ++j) {
                    auto [u_id, v_id, distance] = pairs[j];
                    refine_blocks(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, u_id, v_id);
                }

                std::swap(active_this_round, active_next_round);
                std::fill_n(active_next_round, m_k, 0);
            }
        }

        void refine_blocks(const u64 level,
                           const u64 max_level,
                           const graph_t& g,
                           const d_oracle_t& d_oracle,
                           bv_manager_t& bv_manager,
                           p_manager_t& p_manager,
                           q_graph_t& q_graph,
                           partition_t u_id,
                           partition_t v_id) {
            refine_blocks_new(level, max_level, g, d_oracle, bv_manager, p_manager, q_graph, u_id, v_id);
            return;

            static f64 time_setup = 0.0;
            static f64 time_solve = 0.0;

            if (bv_manager.size(u_id) + bv_manager.size(v_id) > config->max_n_vertices) { return; }

            auto sp_setup = std::chrono::high_resolution_clock::now();

            // weight_t curr_qap = get_qap(g, p_manager, d_oracle);

            weight_t total_u_id_weight = 0;
            weight_t total_v_id_weight = 0;

            // ilp goes here
            GRBEnv env = GRBEnv(true);
            env.set("OutputFlag", "0"); // disable output
            env.start();
            GRBModel model = GRBModel(env);

            GRBQuadExpr penalty = 0;
            GRBLinExpr u_weight = 0;
            GRBLinExpr v_weight = 0;

            vertices_size = 0;
            forall_bv_id_iu(bv_manager, u_id, i, u)
                {
                    vertices[vertices_size++] = u;
                    vars[u]                   = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                }
            endfor
            forall_bv_id_iu(bv_manager, v_id, i, v)
                {
                    vertices[vertices_size++] = v;
                    vars[v]                   = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                }
            endfor

            forall_bv_id_iu(bv_manager, u_id, i, u)
                {
                    total_u_id_weight += g.weight(u);

                    weight_t penalty_u_id = 0;
                    weight_t penalty_v_id = 0;
                    forall_guivw(g, u, j, vv, w)
                        {
                            partition_t vv_id = p_manager[vv];
                            if (vv_id == v_id) {
                                if (vv > u) { continue; }
                                weight_t p = 2 * w * d_oracle.get(u_id, v_id);
                                penalty += p * (vars[u] + vars[vv] - 2 * vars[u] * vars[vv]);
                            } else {
                                penalty_u_id += 2 * w * d_oracle.get(vv_id, u_id);
                                penalty_v_id += 2 * w * d_oracle.get(vv_id, v_id);
                            }
                        }
                    endfor

                    penalty += penalty_u_id * (1 - vars[u]); // penalty if in u_id
                    penalty += penalty_v_id * vars[u]; // penalty if in v_id

                    u_weight += g.weight(u) * (1 - vars[u]);
                    v_weight += g.weight(u) * vars[u];
                }
            endfor

            forall_bv_id_iu(bv_manager, v_id, i, v)
                {
                    total_v_id_weight += g.weight(v);

                    weight_t penalty_u_id = 0;
                    weight_t penalty_v_id = 0;
                    forall_guivw(g, v, j, vv, w)
                        {
                            partition_t vv_id = p_manager[vv];
                            if (vv_id == u_id) {
                                if (vv > v) { continue; }
                                weight_t p = 2 * w * d_oracle.get(u_id, v_id);
                                penalty += p * (vars[v] + vars[vv] - 2 * vars[v] * vars[vv]);
                            } else {
                                penalty_u_id += 2 * w * d_oracle.get(vv_id, u_id);
                                penalty_v_id += 2 * w * d_oracle.get(vv_id, v_id);
                            }
                        }
                    endfor

                    penalty += penalty_u_id * (1 - vars[v]); // penalty if in u_id
                    penalty += penalty_v_id * vars[v]; // penalty if in v_id

                    u_weight += g.weight(v) * (1 - vars[v]);
                    v_weight += g.weight(v) * vars[v];
                }
            endfor

            model.addConstr(u_weight <= m_lmax - (p_manager.get_bweight(u_id) - total_u_id_weight));
            model.addConstr(v_weight <= m_lmax - (p_manager.get_bweight(v_id) - total_v_id_weight));

            auto ep_setup = std::chrono::high_resolution_clock::now();

            auto sp_solve = std::chrono::high_resolution_clock::now();
            model.setObjective(penalty, GRB_MINIMIZE);
            model.optimize();
            auto ep_solve = std::chrono::high_resolution_clock::now();

            u64 n_changes = 0;

            if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
                for (size_t i = 0; i < vertices_size; ++i) {
                    const vertex_t v  = vertices[i];
                    const GRBVar& var = vars[v];

                    int assignment     = round(var.get(GRB_DoubleAttr_X));
                    partition_t new_id = (assignment == 0) ? u_id : v_id;
                    partition_t old_id = p_manager[v];

                    if (old_id != new_id) {
                        n_changes++;
                        if (bv_manager.is_boundary(v)) {
                            bv_manager.move(g, p_manager, v, old_id, new_id);
                        } else {
                            bv_manager.add_new(g, p_manager, v, new_id);
                        }

                        q_graph.move(g, p_manager, v, old_id, new_id);
                        p_manager.move(v, g.weight(v), old_id, new_id);
                    }
                }
            }

            // time_setup += get_seconds(sp_setup, ep_setup);
            // time_solve += get_seconds(sp_solve, ep_solve);
            // std::cout << "Number of changes: " << n_changes << " " << curr_qap - get_qap(g, p_manager, d_oracle) << " " << time_setup << " " << time_solve << " " << max(p_manager.get_bweights()) << std::endl;
        }

        void refine_blocks_new(const u64 level,
                               const u64 max_level,
                               const graph_t& g,
                               const d_oracle_t& d_oracle,
                               bv_manager_t& bv_manager,
                               p_manager_t& p_manager,
                               q_graph_t& q_graph,
                               partition_t u_id,
                               partition_t v_id) {
            static f64 time_setup = 0.0;
            static f64 time_solve = 0.0;

            if (p_manager.size(u_id) + p_manager.size(v_id) > config->max_n_vertices) { return; }

            auto sp_setup = std::chrono::high_resolution_clock::now();

            weight_t curr_qap = get_qap(g, p_manager, d_oracle);

            // ilp goes here
            GRBEnv env = GRBEnv(true);
            env.set("OutputFlag", "0"); // disable output
            env.start();
            GRBModel model = GRBModel(env);

            GRBQuadExpr penalty = 0;
            GRBLinExpr u_weight = 0;
            GRBLinExpr v_weight = 0;

            vertices_size = 0;
            for (vertex_t u = 0; u < g.get_n(); ++u) {
                if (p_manager[u] == u_id || p_manager[u] == v_id) {
                    vertices[vertices_size++] = u;
                    vars[u]                   = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                }
            }

            if (vertices_size > config->max_n_vertices) { return; }

            for (size_t i = 0; i < vertices_size; ++i) {
                const vertex_t u    = vertices[i];
                partition_t curr_id = p_manager[u];

                weight_t penalty_u_id = 0;
                weight_t penalty_v_id = 0;

                forall_guivw(g, u, j, vv, w)
                    {
                        partition_t vv_id = p_manager[vv];
                        if (vv_id == u_id || vv_id == v_id) {
                            if (vv > u) { continue; }
                            weight_t p = 2 * w * d_oracle.get(u_id, v_id);
                            penalty += p * (vars[u] + vars[vv] - 2 * vars[u] * vars[vv]);
                        } else {
                            penalty_u_id += 2 * w * d_oracle.get(vv_id, u_id);
                            penalty_v_id += 2 * w * d_oracle.get(vv_id, v_id);
                        }
                    }
                endfor

                penalty += penalty_u_id * (1 - vars[u]); // penalty if in u_id
                penalty += penalty_v_id * vars[u]; // penalty if in v_id

                u_weight += g.weight(u) * (1 - vars[u]);
                v_weight += g.weight(u) * vars[u];
            }

            model.addConstr(u_weight <= m_lmax);
            model.addConstr(v_weight <= m_lmax);

            auto ep_setup = std::chrono::high_resolution_clock::now();

            auto sp_solve = std::chrono::high_resolution_clock::now();
            model.setObjective(penalty, GRB_MINIMIZE);
            model.optimize();
            auto ep_solve = std::chrono::high_resolution_clock::now();

            u64 n_changes = 0;

            if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
                for (size_t i = 0; i < vertices_size; ++i) {
                    const vertex_t v  = vertices[i];
                    const GRBVar& var = vars[v];

                    int assignment     = round(var.get(GRB_DoubleAttr_X));
                    partition_t new_id = (assignment == 0) ? u_id : v_id;
                    partition_t old_id = p_manager[v];

                    if (old_id != new_id) {
                        n_changes++;
                        if (bv_manager.is_boundary(v)) {
                            bv_manager.move(g, p_manager, v, old_id, new_id);
                        } else {
                            bv_manager.add_new(g, p_manager, v, new_id);
                        }

                        q_graph.move(g, p_manager, v, old_id, new_id);
                        p_manager.move(v, g.weight(v), old_id, new_id);
                    }
                }
            }

            time_setup += get_seconds(sp_setup, ep_setup);
            time_solve += get_seconds(sp_solve, ep_solve);
            std::cout << "Number of changes: " << n_changes << " " << curr_qap - get_qap(g, p_manager, d_oracle) << " " << time_setup << " " << time_solve << " " << max(p_manager.get_bweights()) << std::endl;
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

#endif //HEIPROMAP_ILP_REFINEMENT_H
