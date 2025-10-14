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

#include "../../commons/definitions.h"
#include "../../commons/utils.h"
#include "../../../extern/local/kahip/include/kaHIP_interface.h"

namespace HeiProMap {
    enum GlobalMultisectionMode {
        GLOBAL_MULTISECTION_UNDEFINED,
        GLOBAL_MULTISECTION_STRONG,
        GLOBAL_MULTISECTION_ECO,
        GLOBAL_MULTISECTION_FAST,
    };

    inline GlobalMultisectionMode string_to_global_multisection_mode(const std::string &str) {
        if (str == "UNDEFINED") return GLOBAL_MULTISECTION_UNDEFINED;
        if (str == "strong") return GLOBAL_MULTISECTION_STRONG;
        if (str == "eco") return GLOBAL_MULTISECTION_ECO;
        if (str == "fast") return GLOBAL_MULTISECTION_FAST;
        return GLOBAL_MULTISECTION_UNDEFINED;
    }

    inline std::string global_multisection_mode_to_string(GlobalMultisectionMode mode) {
        switch (mode) {
            case GLOBAL_MULTISECTION_UNDEFINED:
                return "UNDEFINED";
            case GLOBAL_MULTISECTION_STRONG:
                return "strong";
            case GLOBAL_MULTISECTION_ECO:
                return "eco";
            case GLOBAL_MULTISECTION_FAST:
                return "fast";
            default:
                return "UNDEFINED";
        }
    }

    class GlobalMultisectionConfiguration {
    public:
        std::string            mode_string;
        GlobalMultisectionMode mode; // Which mode to use STRONG, ECO, FAST
    };

    struct Item {
        // variables needed for KaFFPa
        int    n               = 0;
        int    *vwgt           = nullptr;
        int    *xadj           = nullptr;
        int    *adjcwgt        = nullptr;
        int    *adjncy         = nullptr;
        int    nparts          = 0;
        double imbalance       = 0;
        bool   suppress_output = true;
        int    seed            = 0;
        int    mode            = 0;
        int    edge_cut_temp   = 0;
        int    *part_temp      = nullptr;
        int    edge_cut        = 0;
        int    *part           = nullptr;

        // variables needed for multisection
        int                      total_weight = 0;
        TranslationTable<int>    tt;
        std::vector<partition_t> identifier;

        ~Item() {
            free(vwgt);
            free(xadj);
            free(adjcwgt);
            free(adjncy);
            free(part_temp);
            free(part);
        }
    };

    class GlobalMultisectionPartitioner {
    private:
        std::vector<Item *> free_items;

    public:
        ~GlobalMultisectionPartitioner()  {
            for (auto &item: free_items) {
                delete item;
            }
        }

        void partition(const graph_t &g,
                       p_manager_t &p_manager,
                       const std::vector<partition_t> &hierarchy,
                       const std::vector<weight_t> &distance,
                       const f64 imbalance,
                       RandomEngine &t_random_engine,
                       const GlobalMultisectionConfiguration &i_config)  {
            GlobalMultisectionConfiguration config = *dynamic_cast<const GlobalMultisectionConfiguration *>(&i_config);

            int mode;

            if (config.mode == GLOBAL_MULTISECTION_STRONG) {
                mode = STRONG;
            } else if (config.mode == GLOBAL_MULTISECTION_ECO) {
                mode = ECO;
            } else if (config.mode == GLOBAL_MULTISECTION_FAST) {
                mode = FAST;
            } else {
                std::cout << "ERROR: GlobalMultisectionPartitioner: unknown mode " << config.mode << " " << config.mode_string << std::endl;
                exit(EXIT_FAILURE);
            }

            int n = (int) g.get_n();
            int m = (int) g.get_m();
            int l = (int) hierarchy.size();

            std::vector<partition_t> index_vec = {1};
            for (int                 i         = 0; i < l - 1; ++i) { index_vec.push_back(index_vec[i] * hierarchy[i]); }

            std::vector<partition_t> k_rem_vec(l);
            u64                      p = 1;

            for (int i = 0; i < l; ++i) {
                k_rem_vec[i] = p * hierarchy[i];
                p *= hierarchy[i];
            }

            const f64         global_imbalance = imbalance;
            const weight_t    global_g_weight  = g.weight();
            const partition_t global_k         = prod<partition_t>(hierarchy);

            // create the first graph
            Item *first_graph = get_available_item(n, m, l);
            first_graph->tt.reserve(n, g.get_n());

            // initialize the translation table of the first graph
            vertex_t new_u = 0;
            forall_gu(g, old_u)
                {
                    first_graph->tt.add(old_u, new_u);
                    new_u += 1;
                }
            endfor

            // write the graph
            first_graph->n = (int) new_u;
            first_graph->xadj[0] = 0;
            first_graph->total_weight = 0;
            forall_gu(g, old_u)
                {
                    new_u = first_graph->tt.get_n(old_u);

                    first_graph->vwgt[new_u] = (int) g.weight(old_u);
                    first_graph->total_weight += (int) g.weight(old_u);
                    first_graph->xadj[new_u + 1] = first_graph->xadj[new_u];

                    forall_guivw(g, old_u, i, old_v, w)
                        {
                            vertex_t new_v = first_graph->tt.get_n(old_v);

                            first_graph->adjcwgt[first_graph->xadj[new_u + 1]] = (int) w;
                            first_graph->adjncy[first_graph->xadj[new_u + 1]]  = (int) new_v;
                            first_graph->xadj[new_u + 1] += 1;
                        }
                    endfor
                }
            endfor

            // fill in other information
            first_graph->nparts          = (int) hierarchy.back();
            first_graph->imbalance       = determine_adaptive_imbalance(global_imbalance, global_g_weight, global_k, first_graph->total_weight, k_rem_vec[l - 1], l);
            first_graph->suppress_output = true;
            first_graph->seed            = t_random_engine.get_s32();
            first_graph->mode            = mode;

            // initialize stack;
            std::vector<Item *> stack = {first_graph};

            // process the stack
            while (!stack.empty()) {
                Item *item = stack.back(); // process first item
                stack.pop_back(); // remove top item

                // TIME_POINT(sp_kaffpa);
                item->edge_cut = std::numeric_limits<int>::max();
                for (int i = 0; i < 1; ++i) {
                    kaffpa(&item->n, item->vwgt, item->xadj, item->adjcwgt, item->adjncy, &item->nparts, &item->imbalance, item->suppress_output, item->seed + i, item->mode, &item->edge_cut_temp, item->part_temp);
                    if (item->edge_cut_temp < item->edge_cut) {
                        item->edge_cut = item->edge_cut_temp;
                        std::swap(item->part_temp, item->part);
                    }
                }

                if ((int) item->identifier.size() == l - 1) {
                    // insert solution
                    u64           offset = 0;
                    for (int      i      = 0; i < l - 1; ++i) { offset += item->identifier[i] * index_vec[index_vec.size() - 1 - i]; }
                    for (vertex_t u      = 0; u < (vertex_t) item->n; ++u) { p_manager.set(item->tt.get_o(u), item->vwgt[u], offset + item->part[u]); }
                } else {
                    // create the subgraphs and place them in the next stack

                    // collect the number of vertices and edges for each new subgraph
                    std::vector<int> new_n(item->nparts, 0);
                    for (int         u = 0; u < item->n; ++u) {
                        int p_id = item->part[u];
                        new_n[p_id] += 1; // increase number of vertices
                    }

                    // create the new subgraphs on the stack
                    for (int i = 0; i < item->nparts; ++i) {
                        Item *new_item = get_available_item(n, m, l);
                        new_item->identifier = item->identifier;
                        new_item->identifier.push_back(i);
                        new_item->tt.reserve(new_n[i], g.get_n());
                        stack.push_back(new_item);
                    }

                    // fill the translation tables
                    std::vector<int> new_us(item->nparts, 0);

                    for (int old_u = 0; old_u < item->n; ++old_u) {
                        int p_id = item->part[old_u];
                        int idx  = (int) stack.size() - (item->nparts - p_id);

                        stack[idx]->tt.add(item->tt.get_o(old_u), new_us[p_id]);
                        new_us[p_id] += 1;
                    }

                    // finalize the translation tables
                    for (int i = 0; i < item->nparts; ++i) {
                        stack[stack.size() - 1 - i]->n            = (int) new_n[new_n.size() - 1 - i];
                        stack[stack.size() - 1 - i]->total_weight = 0;
                        stack[stack.size() - 1 - i]->xadj[0] = 0;
                    }

                    // create the graphs
                    for (int old_u = 0; old_u < item->n; ++old_u) {
                        int p_id = item->part[old_u];
                        int idx  = (int) stack.size() - (item->nparts - p_id);

                        int new_u = stack[idx]->tt.get_n(item->tt.get_o(old_u)); // vertex in new graph

                        // set the weight
                        stack[idx]->vwgt[new_u] = item->vwgt[old_u];
                        stack[idx]->total_weight += item->vwgt[old_u];
                        stack[idx]->xadj[new_u + 1] = stack[idx]->xadj[new_u];

                        // set the edges
                        for (int i = item->xadj[old_u]; i < item->xadj[old_u + 1]; ++i) {
                            int old_v = item->adjncy[i];

                            if (item->part[old_v] == p_id) {
                                // add the edge
                                int new_v = stack[idx]->tt.get_n(item->tt.get_o(old_v)); // vertex in new graph

                                stack[idx]->adjncy[stack[idx]->xadj[new_u + 1]]  = new_v;
                                stack[idx]->adjcwgt[stack[idx]->xadj[new_u + 1]] = item->adjcwgt[i];
                                stack[idx]->xadj[new_u + 1] += 1;
                            }
                        }
                    }

                    // fill in other information
                    for (int i = 0; i < item->nparts; ++i) {
                        int idx = (int) stack.size() - 1 - i;
                        stack[idx]->nparts          = (int) hierarchy[l - 1 - stack[idx]->identifier.size()];
                        stack[idx]->imbalance       = determine_adaptive_imbalance(global_imbalance, global_g_weight, global_k, stack[idx]->total_weight, k_rem_vec[l - 1 - stack[idx]->identifier.size()], l - stack[idx]->identifier.size());
                        stack[idx]->suppress_output = true;
                        stack[idx]->seed            = t_random_engine.get_s32();
                        stack[idx]->mode            = mode;
                    }
                }
                free_items.push_back(item);
            }
        }

        /**
         * Returns a pointer to an item that can hold at least n vertices with m edges and an identifier of size l.
         *
         * @param n
         * @param m
         * @param l
         * @return
         */
        Item *get_available_item(int n, int m, int l) {
            if (free_items.empty()) {
                int n_64 = round_up_64(n + 1);
                int m_64 = round_up_64(m);

                // if no items are available, then create a new one
                Item *new_item = new Item();
                new_item->vwgt      = (int *) aligned_alloc(64, n_64 * sizeof(int));
                new_item->xadj      = (int *) aligned_alloc(64, n_64 * sizeof(int));
                new_item->adjcwgt   = (int *) aligned_alloc(64, m_64 * sizeof(int));
                new_item->adjncy    = (int *) aligned_alloc(64, m_64 * sizeof(int));
                new_item->part_temp = (int *) aligned_alloc(64, n_64 * sizeof(int));
                new_item->part      = (int *) aligned_alloc(64, n_64 * sizeof(int));

                return new_item;
            } else {
                // if an item exists, then just use the last one
                Item *last_item = free_items.back();
                free_items.pop_back();
                return last_item;
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
