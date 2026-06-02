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

#ifndef HEIPROMAP_GLOBAL_MULTISECTION_H
#define HEIPROMAP_GLOBAL_MULTISECTION_H

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>

#include "../definitions.h"
#include "../utility/aligned_array.h"
#include "../utility/utils.h"
#include "../datastructures/partition_manager.h"
#include "../datastructures/boundary_vertex_manger.h"
#include "../datastructures/quotient_graph.h"
#include "../datastructures/block_conn.h"
#include "../datastructures/distance_oracle.h"
#include "kaffpa_partitioner.h"
#include "greedy_partitioner.h"
#include "recursive_bisection.h"
#include "../coarsening/heavy_edge_matching.h"
#include "../refinement/flow_based_refinement.h"
#include "../refinement/label_propagation_refinement.h"
#include "../refinement/quotient_graph_refinement.h"
#include "../utility/qap.h"
#include "kway_partitioner/kway_core.h"
#include "../utility/translation_table.h"

namespace HeiProMap {
    enum GlobalMultisectionMode {
        GLOBAL_MULTISECTION_UNDEFINED,
        GLOBAL_MULTISECTION_KAFFPA_STRONG,
        GLOBAL_MULTISECTION_KAFFPA_ECO,
        GLOBAL_MULTISECTION_KAFFPA_FAST,
        GLOBAL_MULTISECTION_METIS_KWAY,
        GLOBAL_MULTISECTION_HEIPA_FAST,
        GLOBAL_MULTISECTION_HEIPA_ECO,
        GLOBAL_MULTISECTION_HEIPA_STRONG,
        GLOBAL_MULTISECTION_HEIPA_SUPER_STRONG,
        GLOBAL_MULTISECTION_GREEDY,
        GLOBAL_MULTISECTION_RECURSIVE_BISECTION,
        GLOBAL_MULTISECTION_GGG,
        GLOBAL_MULTISECTION_HYBRID,
    };

    inline GlobalMultisectionMode string_to_global_multisection_mode(const std::string &str) {
        if (str == "UNDEFINED") return GLOBAL_MULTISECTION_UNDEFINED;
        if (str == "kaffpa-strong") return GLOBAL_MULTISECTION_KAFFPA_STRONG;
        if (str == "kaffpa-eco") return GLOBAL_MULTISECTION_KAFFPA_ECO;
        if (str == "kaffpa-fast") return GLOBAL_MULTISECTION_KAFFPA_FAST;
        if (str == "metis-kway") return GLOBAL_MULTISECTION_METIS_KWAY;
        if (str == "heipa-fast") return GLOBAL_MULTISECTION_HEIPA_FAST;
        if (str == "heipa-eco") return GLOBAL_MULTISECTION_HEIPA_ECO;
        if (str == "heipa-strong") return GLOBAL_MULTISECTION_HEIPA_STRONG;
        if (str == "heipa-super-strong") return GLOBAL_MULTISECTION_HEIPA_SUPER_STRONG;
        if (str == "greedy") return GLOBAL_MULTISECTION_GREEDY;
        if (str == "recursive-bisection") return GLOBAL_MULTISECTION_RECURSIVE_BISECTION;
        if (str == "ggg") return GLOBAL_MULTISECTION_GGG;
        if (str == "hybrid") return GLOBAL_MULTISECTION_HYBRID;
        return GLOBAL_MULTISECTION_UNDEFINED;
    }

    inline std::string global_multisection_mode_to_string(GlobalMultisectionMode mode) {
        switch (mode) {
            case GLOBAL_MULTISECTION_UNDEFINED:
                return "UNDEFINED";
            case GLOBAL_MULTISECTION_KAFFPA_STRONG:
                return "kaffpa-strong";
            case GLOBAL_MULTISECTION_KAFFPA_ECO:
                return "kaffpa-eco";
            case GLOBAL_MULTISECTION_KAFFPA_FAST:
                return "kaffpa-fast";
            case GLOBAL_MULTISECTION_METIS_KWAY:
                return "metis-kway";
            case GLOBAL_MULTISECTION_HEIPA_FAST:
                return "heipa-fast";
            case GLOBAL_MULTISECTION_HEIPA_ECO:
                return "heipa-eco";
            case GLOBAL_MULTISECTION_HEIPA_STRONG:
                return "heipa-strong";
            case GLOBAL_MULTISECTION_HEIPA_SUPER_STRONG:
                return "heipa-super-strong";
            case GLOBAL_MULTISECTION_GREEDY:
                return "greedy";
            case GLOBAL_MULTISECTION_RECURSIVE_BISECTION:
                return "recursive-bisection";
            case GLOBAL_MULTISECTION_GGG:
                return "ggg";
            case GLOBAL_MULTISECTION_HYBRID:
                return "hybrid";
            default:
                return "UNDEFINED";
        }
    }

    class GlobalMultisectionConfiguration {
    public:
        std::string mode_string;
        GlobalMultisectionMode mode; // Which mode to use STRONG, ECO, FAST
        u64 kappa = 1;
        bool refine = false;
        u64 v_cycles = 0;
        u64 v_cycle_depth = 1;
        bool collect_dataset = false;
        std::string data_dir = "data";
        LabelPropagationConfiguration label_propagation_config = LabelPropagationConfiguration("Label Propagation");
        QuotientGraphRefinementConfiguration quotient_graph_refinement_config = QuotientGraphRefinementConfiguration("Quotient Graph");
        FlowBasedRefinementConfiguration flow_based_refinement_config = FlowBasedRefinementConfiguration("Flow Based");
    };

    // Forward declaration of the HeiPa wrapper to resolve circular dependency
    void heipa_multisection_partition_wrapper(graph_t &g, partition_t k, f64 imb, u64 seed, AlignedArray<partition_t> &partition, GlobalMultisectionMode mode, u64 kappa);

    struct Item {
        graph_t *g = nullptr;
        partition_t k = 0;
        f64 imb = 0.0;
        u64 seed = 0;
        TranslationTable<vertex_t> *tt = nullptr;
        std::vector<partition_t> *identifier = nullptr;
        bool free = true;
    };

    class GlobalMultisectionPartitioner {
    public:
        static void partition(graph_t &g,
                              p_manager_t &p_manager,
                              const std::vector<partition_t> &hierarchy,
                              [[maybe_unused]] const std::vector<weight_t> &distance,
                              const f64 imbalance,
                              const GlobalMultisectionConfiguration &i_config,
                              u64 seed) {
            HEIPROMAP_PROFILE_SCOPE("partition", "GlobalMultisectionPartitioner", "partition");
            GlobalMultisectionConfiguration config = i_config;

            RandomEngine rnd_engine(seed);

            AlignedArray<partition_t> partition;
            partition.initialize(g.n);

            partition_t l = (partition_t) hierarchy.size();

            std::vector<partition_t> index_vec = {1};
            for (partition_t i = 0; i < l - 1; ++i) { index_vec.push_back(index_vec[i] * hierarchy[i]); }

            std::vector<partition_t> k_rem_vec(l);
            u64 p = 1;

            for (partition_t i = 0; i < l; ++i) {
                k_rem_vec[i] = p * hierarchy[i];
                p *= hierarchy[i];
            }

            const f64 global_imbalance = imbalance;
            const weight_t global_g_weight = g.g_weight;
            const partition_t global_k = prod<partition_t>(hierarchy);
            const vertex_t global_n = g.n;

            // create the first graph
            Item first_graph;
            first_graph.free = false;
            first_graph.g = &g;
            first_graph.identifier = new std::vector<partition_t>();
            first_graph.tt = new TranslationTable<vertex_t>();
            first_graph.tt->reserve(g.n, global_n);

            // initialize the translation table of the first graph
            vertex_t temp_new_u = 0;
            for (vertex_t old_u = 0; old_u < g.n; ++old_u) {
                {
                    first_graph.tt->add(old_u, temp_new_u);
                    temp_new_u += 1;
                }
            }

            // fill in other information
            first_graph.k = hierarchy.back();
            first_graph.imb = determine_adaptive_imbalance(global_imbalance, global_g_weight, global_k, first_graph.g->g_weight, k_rem_vec[l - 1], l);
            first_graph.seed = rnd_engine.get_s32();

            // initialize stack;
            std::vector<Item> stack = {first_graph};

            // process the stack
            while (!stack.empty()) {
                Item item = stack.back(); // process first item
                stack.pop_back();         // remove top item

                if (config.mode == GLOBAL_MULTISECTION_KAFFPA_STRONG) {
                    kaffpa_partition(*item.g, item.k, item.imb, KAFFPA_PARTITION_STRONG, item.seed, partition, config.kappa, config.collect_dataset, config.data_dir);
                } else if (config.mode == GLOBAL_MULTISECTION_KAFFPA_ECO) {
                    kaffpa_partition(*item.g, item.k, item.imb, KAFFPA_PARTITION_ECO, item.seed, partition, config.kappa, config.collect_dataset, config.data_dir);
                } else if (config.mode == GLOBAL_MULTISECTION_KAFFPA_FAST) {
                    kaffpa_partition(*item.g, item.k, item.imb, KAFFPA_PARTITION_FAST, item.seed, partition, config.kappa, config.collect_dataset, config.data_dir);
                } else if (config.mode == GLOBAL_MULTISECTION_METIS_KWAY) {
                    kway_partition(*item.g, item.k, item.imb, item.seed, partition, config.kappa);
                } else if (config.mode >= GLOBAL_MULTISECTION_HEIPA_FAST && config.mode <= GLOBAL_MULTISECTION_HEIPA_SUPER_STRONG) {
                    heipa_multisection_partition_wrapper(*item.g, item.k, item.imb, item.seed, partition, config.mode, config.kappa);
                } else {
                    std::cerr << "Mode " << config.mode << " not implemented" << std::endl;
                    abort();
                }

                if (config.refine) {
                    u64 cycles = std::max((u64)1, config.v_cycles);
                    for (u64 v = 0; v < cycles; ++v) {
                        refine_partition(*item.g, item.k, item.imb, item.seed + v, partition, config);
                    }
                }

                if (item.identifier->size() == l - 1) {
                    // insert solution
                    u64 offset = 0;
                    for (partition_t i = 0; i < l - 1; ++i) { offset += item.identifier->operator[](i) * index_vec[index_vec.size() - 1 - i]; }
                    for (vertex_t u = 0; u < item.g->n; ++u) { p_manager.set(item.tt->get_o(u), item.g->v_weights[u], offset + partition[u]); }
                } else {
                    HEIPROMAP_PROFILE_SCOPE("partition", "GlobalMultisectionPartitioner", "split");
                    // split problem
                    std::vector<vertex_t> new_ns(item.k, 0);
                    std::vector<vertex_t> new_ms(item.k, 0);
                    std::vector<weight_t> new_ws(item.k, 0);

                    for (vertex_t u = 0; u < item.g->n; ++u) {
                        partition_t id = partition[u];
                        new_ns[id] += 1;
                        new_ws[id] += item.g->v_weights[u];
                        for (size_t i = item.g->neighborhoods[u]; i < item.g->neighborhoods[u + 1]; ++i) {
                            vertex_t v = item.g->edges_v[i];
                            if (partition[v] == id) { new_ms[id] += 1; }
                        }
                    }

                    for (partition_t i = 0; i < item.k; ++i) {
                        Item next_item;
                        next_item.free = true;
                        next_item.g = new graph_t(new_ns[i], new_ms[i], new_ws[i]);
                        next_item.tt = new TranslationTable<vertex_t>();
                        next_item.tt->reserve(new_ns[i], global_n);

                        std::vector<vertex_t> new_us(item.k, 0);
                        for (vertex_t old_u = 0; old_u < item.g->n; ++old_u) {
                            if (partition[old_u] == i) {
                                next_item.tt->add(item.tt->get_o(old_u), new_us[i]);
                                new_us[i] += 1;
                            }
                        }

                        std::vector<vertex_t> degrees(next_item.g->n, 0);
                        for (vertex_t old_u = 0; old_u < item.g->n; ++old_u) {
                            if (partition[old_u] == i) {
                                vertex_t new_u = next_item.tt->get_n(item.tt->get_o(old_u));
                                for (size_t j = item.g->neighborhoods[old_u]; j < item.g->neighborhoods[old_u + 1]; ++j) {
                                    vertex_t old_v = item.g->edges_v[j];
                                    if (partition[old_v] == i) {
                                        degrees[new_u]++;
                                    }
                                }
                            }
                        }

                        next_item.g->neighborhoods[0] = 0;
                        for (vertex_t j = 0; j < next_item.g->n; ++j) {
                            next_item.g->neighborhoods[j + 1] = next_item.g->neighborhoods[j] + degrees[j];
                        }

                        std::vector<vertex_t> cursor(next_item.g->n, 0);
                        for (vertex_t old_u = 0; old_u < item.g->n; ++old_u) {
                            if (partition[old_u] == i) {
                                vertex_t new_u = next_item.tt->get_n(item.tt->get_o(old_u));
                                next_item.g->v_weights[new_u] = item.g->v_weights[old_u];

                                for (size_t j = item.g->neighborhoods[old_u]; j < item.g->neighborhoods[old_u + 1]; ++j) {
                                    vertex_t old_v = item.g->edges_v[j];
                                    if (partition[old_v] == i) {
                                        vertex_t new_v = next_item.tt->get_n(item.tt->get_o(old_v));
                                        size_t pos = next_item.g->neighborhoods[new_u] + cursor[new_u];
                                        next_item.g->edges_v[pos] = new_v;
                                        next_item.g->edges_w[pos] = item.g->edges_w[j];
                                        cursor[new_u]++;
                                    }
                                }
                            }
                        }

                        next_item.identifier = new std::vector<partition_t>(*item.identifier);
                        next_item.identifier->push_back(i);
                        next_item.k = hierarchy[hierarchy.size() - 1 - next_item.identifier->size()];
                        next_item.imb = determine_adaptive_imbalance(global_imbalance, global_g_weight, global_k, next_item.g->g_weight, k_rem_vec[l - 1 - next_item.identifier->size()], l - (partition_t) next_item.identifier->size());
                        next_item.seed = rnd_engine.get_s32();

                        stack.push_back(next_item);
                    }
                }

                // free
                delete item.identifier;
                delete item.tt;
                if (item.free) {
                    delete item.g;
                }
            }
        }

        static f64 determine_adaptive_imbalance(const f64 global_imbalance,
                                                const weight_t global_g_weight,
                                                const partition_t global_k,
                                                const weight_t local_g_weight,
                                                const partition_t local_k,
                                                const partition_t depth) {
            f64 local_imbalance = (((1.0 + global_imbalance) * (f64) local_k * (f64) global_g_weight) / (f64) (global_k * local_g_weight));
            local_imbalance = std::pow(local_imbalance, (f64) 1 / (f64) depth) - 1.0;
            return local_imbalance;
        }

    private:
        static void run_core_refinement(graph_t &g,
                                        partition_t k,
                                        f64 imbalance,
                                        u64 seed,
                                        PartitionManager &pm,
                                        const GlobalMultisectionConfiguration &config) {
            DistanceOracle d_oracle;
            d_oracle.initialize({k}, {1});

            BoundaryVertexManager bv_manager;
            bv_manager.initialize(g.n, k);

            QuotientGraph q_graph;
            q_graph.initialize(k);

            BlockConn block_conn;
            block_conn.initialize(g.n, 2 * g.m + g.n, k);
            block_conn.reset_build();

            for (vertex_t u = 0; u < g.n; ++u) {
                block_conn.begin_vertex(g, u);
                partition_t u_id = pm[u];
                for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                    vertex_t v = g.edges_v[i];
                    weight_t w = g.edges_w[i];
                    partition_t v_id = pm[v];
                    block_conn.add_connection(u, v_id, w);
                    if (u_id != v_id) {
                        bv_manager.add(u, u_id);
                        if (u < v) q_graph.add_edge(u_id, v_id, w);
                    }
                }
            }

            AlignedArray<weight_t> lmax_constraints;
            lmax_constraints.initialize(k);
            weight_t lmax = std::ceil((1.0 + imbalance) * ((f64) g.g_weight / (f64) k));
            for (partition_t i = 0; i < k; ++i) {
                lmax_constraints[i] = lmax;
            }

            if (config.label_propagation_config.enabled) {
                LabelPropagationRefinement lp_refine;
                lp_refine.initialize(g.n, g.m, k, 1, seed, config.label_propagation_config);
                lp_refine.refine(g, d_oracle, bv_manager, pm, q_graph, block_conn, lmax_constraints, g.uniform_v_weights, g.uniform_e_weights);
            }

            if (config.quotient_graph_refinement_config.enabled) {
                QuotientGraphRefinement qg_refine;
                qg_refine.initialize(g.n, g.m, k, 1, seed, config.quotient_graph_refinement_config);
                qg_refine.refine(g, d_oracle, bv_manager, pm, q_graph, block_conn, lmax_constraints, g.uniform_v_weights, g.uniform_e_weights);
            }

            if (config.flow_based_refinement_config.enabled) {
                FlowBasedRefinement flow_refine;
                flow_refine.initialize(g.n, g.m, k, 1, seed, config.flow_based_refinement_config);
                flow_refine.refine(g, d_oracle, bv_manager, pm, q_graph, block_conn, lmax_constraints, g.uniform_v_weights, g.uniform_e_weights);
            }
        }

        static void refine_partition(graph_t &g,
                                     partition_t k,
                                     f64 imbalance,
                                     u64 seed,
                                     AlignedArray<partition_t> &partition,
                                     const GlobalMultisectionConfiguration &config,
                                     u64 depth = 0) {
            HEIPROMAP_PROFILE_SCOPE("refinement", "GlobalMultisectionPartitioner", "refine_partition");

            PartitionManager pm;
            pm.initialize(g.n, k, g.g_weight);
            pm.reset_weights();
            for (vertex_t u = 0; u < g.n; ++u) {
                pm.set(u, g.v_weights[u], partition[u]);
            }

            // 1. Initial Refinement
            run_core_refinement(g, k, imbalance, seed, pm, config);

            // 2. V-Cycles
            if (depth < config.v_cycle_depth) {
                Mapping mapping;
                HeavyEdgeMatching hem;
                HeavyEdgeMatchingConfiguration hem_config;
                hem.match(g, pm, mapping, imbalance, seed + depth * 100, hem_config);

                if (g.n > 2 * k) {
                    graph_t coarse_g;
                    if (g.uniform_v_weights && g.uniform_e_weights) {
                        coarse_g.initialize<true, true>(g, mapping);
                    } else if (g.uniform_v_weights) {
                        coarse_g.initialize<true, false>(g, mapping);
                    } else if (g.uniform_e_weights) {
                        coarse_g.initialize<false, true>(g, mapping);
                    } else {
                        coarse_g.initialize<false, false>(g, mapping);
                    }

                    AlignedArray<partition_t> coarse_partition;
                    coarse_partition.initialize(coarse_g.n);
                    for (vertex_t u = 0; u < g.n; ++u) {
                        coarse_partition[mapping.get(u)] = pm[u];
                    }

                    refine_partition(coarse_g, k, imbalance, seed + depth * 100, coarse_partition, config, depth + 1);

                    // Project back
                    for (vertex_t u = 0; u < g.n; ++u) {
                        partition_t old_id = pm[u];
                        partition_t new_id = coarse_partition[mapping.get(u)];
                        if (old_id != new_id) {
                            pm.move(u, g.uniform_v_weights ? 1 : g.v_weights[u], old_id, new_id);
                        }
                    }

                    // Refine Fine again
                    run_core_refinement(g, k, imbalance, seed + depth * 100 + 10, pm, config);
                }
            }

            for (vertex_t u = 0; u < g.n; ++u) {
                partition[u] = pm[u];
            }
        }
    };
}

#endif //HEIPROMAP_GLOBAL_MULTISECTION_H
