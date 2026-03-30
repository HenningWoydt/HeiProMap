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

#include "../definitions.h"
#include "../utility/utils.h"
#include "kaffpa_partitioner.h"
#include "metis_partitioner.h"
#include "mtkahypar_partition.h"
#include "../refinement/flow_based_refinement.h"
#include "../refinement/multi_try_fm_refinement.h"
#include "../utility/qap.h"
#include "kway_partitioner/kway_core.h"

namespace HeiProMap {
    enum GlobalMultisectionMode {
        GLOBAL_MULTISECTION_UNDEFINED,
        GLOBAL_MULTISECTION_KAFFPA_STRONG,
        GLOBAL_MULTISECTION_KAFFPA_ECO,
        GLOBAL_MULTISECTION_KAFFPA_FAST,
        GLOBAL_MULTISECTION_METIS_RECURSIVE,
        GLOBAL_MULTISECTION_METIS_KWAY,
        GLOBAL_MULTISECTION_MTKAHYPAR_DEFAULT,
        GLOBAL_MULTISECTION_MTKAHYPAR_QUALITY,
        GLOBAL_MULTISECTION_MTKAHYPAR_HIGHEST_QUALITY,
    };

    inline GlobalMultisectionMode string_to_global_multisection_mode(const std::string &str) {
        if (str == "UNDEFINED") return GLOBAL_MULTISECTION_UNDEFINED;
        if (str == "kaffpa-strong") return GLOBAL_MULTISECTION_KAFFPA_STRONG;
        if (str == "kaffpa-eco") return GLOBAL_MULTISECTION_KAFFPA_ECO;
        if (str == "kaffpa-fast") return GLOBAL_MULTISECTION_KAFFPA_FAST;
        if (str == "metis-recursive") return GLOBAL_MULTISECTION_METIS_RECURSIVE;
        if (str == "metis-kway") return GLOBAL_MULTISECTION_METIS_KWAY;
        if (str == "mtkahypar-default") return GLOBAL_MULTISECTION_MTKAHYPAR_DEFAULT;
        if (str == "mtkahypar-quality") return GLOBAL_MULTISECTION_MTKAHYPAR_QUALITY;
        if (str == "mtkahypar-highestquality") return GLOBAL_MULTISECTION_MTKAHYPAR_HIGHEST_QUALITY;
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
            case GLOBAL_MULTISECTION_METIS_RECURSIVE:
                return "metis-recursive";
            case GLOBAL_MULTISECTION_METIS_KWAY:
                return "metis-kway";
            case GLOBAL_MULTISECTION_MTKAHYPAR_DEFAULT:
                return "mtkahypar-default";
            case GLOBAL_MULTISECTION_MTKAHYPAR_QUALITY:
                return "mtkahypar-quality";
            case GLOBAL_MULTISECTION_MTKAHYPAR_HIGHEST_QUALITY:
                return "mtkahypar-highestquality";
            default:
                return "UNDEFINED";
        }
    }

    class GlobalMultisectionConfiguration {
    public:
        std::string mode_string;
        GlobalMultisectionMode mode; // Which mode to use STRONG, ECO, FAST
        u64 kappa = 1;
    };

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
            GlobalMultisectionConfiguration config = *dynamic_cast<const GlobalMultisectionConfiguration *>(&i_config);

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

            // create the first graph
            Item first_graph;
            first_graph.free = false;
            first_graph.g = &g;
            first_graph.identifier = new std::vector<partition_t>();
            first_graph.tt = new TranslationTable<vertex_t>();
            first_graph.tt->reserve(g.n, g.n);

            // initialize the translation table of the first graph
            vertex_t temp_new_u = 0;
            forall_gu(g, old_u)
                {
                    first_graph.tt->add(old_u, temp_new_u);
                    temp_new_u += 1;
                }
            endfor

            // fill in other information
            first_graph.k = hierarchy.back();
            first_graph.imb = determine_adaptive_imbalance(global_imbalance, global_g_weight, global_k, first_graph.g->g_weight, k_rem_vec[l - 1], l);
            first_graph.seed = rnd_engine.get_s32();

            // initialize stack;
            std::vector<Item> stack = {first_graph};

            // process the stack
            while (!stack.empty()) {
                Item item = stack.back(); // process first item
                stack.pop_back(); // remove top item

                if (config.mode == GLOBAL_MULTISECTION_KAFFPA_STRONG) {
                    kaffpa_partition(*item.g, item.k, item.imb, KAFFPA_PARTITION_STRONG, item.seed, partition, config.kappa);
                } else if (config.mode == GLOBAL_MULTISECTION_KAFFPA_ECO) {
                    kaffpa_partition(*item.g, item.k, item.imb, KAFFPA_PARTITION_ECO, item.seed, partition, config.kappa);
                } else if (config.mode == GLOBAL_MULTISECTION_KAFFPA_FAST) {
                    kaffpa_partition(*item.g, item.k, item.imb, KAFFPA_PARTITION_FAST, item.seed, partition, config.kappa);
                } else if (config.mode == GLOBAL_MULTISECTION_METIS_RECURSIVE) {
                    // metis_partition(*item.g, item.k, item.imb, METIS_PARTITION_RECURSIVE, item.seed, partition);
                } else if (config.mode == GLOBAL_MULTISECTION_METIS_KWAY) {
                    kway_partition(*item.g, item.k, item.imb, item.seed, partition, config.kappa);
                    // metis_partition(*item.g, item.k, item.imb, METIS_PARTITION_KWAY, item.seed, partition);
                } else if (config.mode == GLOBAL_MULTISECTION_MTKAHYPAR_DEFAULT) {
                    // mtkahypar_partition(*item.g, item.k, item.imb, MTKAHYPAR_DEFAULT, item.seed, partition);
                } else if (config.mode == GLOBAL_MULTISECTION_MTKAHYPAR_QUALITY) {
                    // mtkahypar_partition(*item.g, item.k, item.imb, MTKAHYPAR_QUALITY, item.seed, partition);
                } else if (config.mode == GLOBAL_MULTISECTION_MTKAHYPAR_HIGHEST_QUALITY) {
                    // mtkahypar_partition(*item.g, item.k, item.imb, MTKAHYPAR_HIGHEST_QUALITY, item.seed, partition);
                } else {
                    std::cerr << "Mode " << config.mode << " not implemented" << std::endl;
                    abort();
                }

                {
                    // refine
                    PartitionManager local_p_manager;
                    BoundaryVertexManager local_boundary_manager;
                    QuotientGraph local_quotient_graph;
                    BlockConn local_block_conn;
                    DistanceOracle local_distance_oracle;

                    local_p_manager.initialize(item.g->n, item.k, item.g->g_weight);
                    local_p_manager.n_vertices[0] = 0;
                    local_p_manager.bweights[0] = 0;
                    local_boundary_manager.initialize(item.g->n, item.k);
                    local_quotient_graph.initialize(item.k);
                    local_block_conn.initialize(item.g->n, item.g->m, item.k);
                    local_distance_oracle.initialize({item.k}, {1});

                    forall_gu((*item.g), u)
                        {
                            const partition_t u_id = partition[u];
                            const weight_t u_w = item.g->v_weights[u];
                            local_p_manager.set(u, u_w, u_id);

                            forall_guivw((*item.g), u, i, v, w)
                                {
                                    const partition_t v_id = partition[v];

                                    if (u_id != v_id) {
                                        local_boundary_manager.add(u, u_id); // boundary vertex
                                        if (u < v) {
                                            local_quotient_graph.add_edge(u_id, v_id, w); // quotient graph
                                        }
                                    }
                                }
                            endfor
                        }
                    endfor
                    local_block_conn.compute_from_scratch((*item.g), local_p_manager);

                    weight_t edge_cut = get_edge_cut((*item.g), local_p_manager);
                    weight_t diff = 1;
                    while (diff > 0) {
                        MultiTryFmRefinementConfiguration multi_try_config("refinements");
                        multi_try_config.max_iteration = 10;
                        multi_try_config.min_n_steps = 10;
                        MultiTryFMRefinement multi_try_fm_refinement;
                        multi_try_fm_refinement.initialize(item.g->n, item.g->m, item.k, 1, seed, multi_try_config);
                        multi_try_fm_refinement.refine((*item.g), local_distance_oracle, local_boundary_manager, local_p_manager, local_quotient_graph, local_block_conn, item.imb);

                        FlowBasedRefinementConfiguration flow_config("refinements");
                        flow_config.alpha = 1.0;
                        flow_config.max_local_iteration = 10;
                        flow_config.max_global_iteration = 2;
                        FlowBasedRefinement flow_refinement;
                        flow_refinement.initialize(item.g->n, item.g->m, item.k, 1, seed, flow_config);
                        flow_refinement.refine((*item.g), local_distance_oracle, local_boundary_manager, local_p_manager, local_quotient_graph, local_block_conn, item.imb);

                        weight_t temp_edge_cut = get_edge_cut((*item.g), local_p_manager);
                        diff = edge_cut - temp_edge_cut;
                        edge_cut = temp_edge_cut;
                    }
                }

                if (item.identifier->size() == l - 1) {
                    // insert solution
                    u64 offset = 0;
                    for (partition_t i = 0; i < l - 1; ++i) { offset += item.identifier->operator[](i) * index_vec[index_vec.size() - 1 - i]; }
                    for (vertex_t u = 0; u < item.g->n; ++u) { p_manager.set(item.tt->get_o(u), item.g->v_weights[u], offset + partition[u]); }
                } else {
                    // create the subgraphs and place them in the next stack

                    // collect the number of vertices and edges for each new subgraph
                    std::vector<vertex_t> new_ns(item.k, 0);
                    std::vector<vertex_t> new_ms(item.k, 0);
                    std::vector<weight_t> new_ws(item.k, 0);
                    for (vertex_t u = 0; u < item.g->n; ++u) {
                        partition_t u_id = partition[u];
                        new_ns[u_id] += 1;
                        new_ws[u_id] += item.g->v_weights[u];
                        for (size_t i = item.g->neighborhoods[u]; i < item.g->neighborhoods[u + 1]; ++i) {
                            vertex_t v = item.g->edges_v[i];
                            partition_t v_id = partition[v];
                            if (v_id == u_id) { new_ms[u_id] += 1; }
                        }
                    }

                    // create the new subgraphs on the stack
                    for (partition_t i = 0; i < item.k; ++i) {
                        stack.emplace_back();
                        Item &new_item = stack.back();
                        new_item.g = new CSRGraph(new_ns[i], new_ms[i], new_ws[i]);
                        new_item.identifier = new std::vector<partition_t>(*item.identifier);
                        new_item.identifier->push_back(i);
                        new_item.tt = new TranslationTable<vertex_t>();
                        new_item.tt->reserve(new_ns[i], g.n);
                    }

                    // fill the translation tables
                    std::vector<vertex_t> new_us(item.k, 0);
                    for (vertex_t old_u = 0; old_u < item.g->n; ++old_u) {
                        partition_t u_id = partition[old_u];
                        size_t idx = stack.size() - (item.k - u_id);

                        stack[idx].tt->add(item.tt->get_o(old_u), new_us[u_id]);
                        new_us[u_id] += 1;
                    }

                    // create the graphs
                    for (vertex_t old_u = 0; old_u < item.g->n; ++old_u) {
                        partition_t u_id = partition[old_u];
                        size_t idx = stack.size() - (item.k - u_id);

                        graph_t &sub_g = *stack[idx].g;
                        TranslationTable<vertex_t> &sub_tt = *stack[idx].tt;

                        vertex_t new_u = sub_tt.get_n(item.tt->get_o(old_u)); // vertex in new graph

                        // set the weight
                        sub_g.v_weights[new_u] = item.g->v_weights[old_u];

                        // set the edges
                        sub_g.neighborhoods[new_u + 1] = sub_g.neighborhoods[new_u];
                        for (size_t i = item.g->neighborhoods[old_u]; i < item.g->neighborhoods[old_u + 1]; ++i) {
                            vertex_t old_v = item.g->edges_v[i];

                            if (partition[old_v] == u_id) {
                                // add the edge
                                vertex_t new_v = sub_tt.get_n(item.tt->get_o(old_v)); // vertex in new graph

                                sub_g.edges_v[sub_g.neighborhoods[new_u + 1]] = new_v;
                                sub_g.edges_w[sub_g.neighborhoods[new_u + 1]] = item.g->edges_w[i];
                                sub_g.neighborhoods[new_u + 1] += 1;
                            }
                        }
                    }

                    // fill in other information
                    for (partition_t i = 0; i < item.k; ++i) {
                        size_t idx = stack.size() - 1 - i;
                        stack[idx].k = hierarchy[l - 1 - stack[idx].identifier->size()];
                        stack[idx].imb = determine_adaptive_imbalance(global_imbalance, global_g_weight, global_k, stack[idx].g->g_weight, k_rem_vec[l - 1 - stack[idx].identifier->size()], l - stack[idx].identifier->size());
                        stack[idx].seed = rnd_engine.get_s32();
                    }
                }
                if (item.free) {
                    delete item.g;
                }
                delete item.identifier;
                delete item.tt;
            }
        }

        static f64 determine_adaptive_imbalance(const f64 global_imbalance,
                                                const u64 global_g_weight,
                                                const u64 global_k,
                                                const u64 local_g_weight,
                                                const u64 local_k_rem,
                                                const u64 depth) {
            f64 local_imbalance = (1.0 + global_imbalance) * ((f64) (local_k_rem * global_g_weight) / (f64) (global_k * local_g_weight));
            local_imbalance = std::pow(local_imbalance, (f64) 1 / (f64) depth) - 1.0;
            return local_imbalance;
        }
    };
}

#endif //HEIPROMAP_GLOBAL_MULTISECTION_H
