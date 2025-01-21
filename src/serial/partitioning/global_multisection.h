#ifndef HEIDELBERGPROCESSMAPPING_GLOBAL_MULTISECTION_H
#define HEIDELBERGPROCESSMAPPING_GLOBAL_MULTISECTION_H

#include "../../definitions.h"
#include "../../macros.h"
#include "../utility/utils.h"
#include "../datastructures/translation_table.h"
#include "../interfaces/ISerialPartitioner.h"
#include "interface/kaHIP_interface.h"

namespace HeiProMap {
    class GlobalMultisectionPartitioner final : public ISerialPartitioner {
    public:
        template<typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialPartitionManager>
        void partition(TSerialGraph &g,
                       TSerialActiveVertexManager &av_manager,
                       TSerialPartitionManager &p_manager,
                       const std::vector<partition_t> &hierarchy,
                       [[maybe_unused]] const std::vector<weight_t> &distance,
                       const f64 imbalance) {
            // references for better code readability
            const size_t l = hierarchy.size();

            std::vector<partition_t> index_vec = {1};
            for (size_t              i         = 0; i < l - 1; ++i) { index_vec.push_back(index_vec[i] * hierarchy[i]); }

            std::vector<partition_t> k_rem_vec(l);
            u64                      p = 1;
            for (size_t              i = 0; i < l; ++i) {
                k_rem_vec[i] = p * hierarchy[i];
                p *= hierarchy[i];
            }

            const f64         global_imbalance = imbalance;
            const weight_t    global_g_weight  = g.get_weight();
            const partition_t global_k         = prod<partition_t>(hierarchy);

            // initialize stack;
            std::vector<Item<TSerialGraph, TSerialActiveVertexManager>> stack      = {{new std::vector<partition_t>(), new KaFFPaGraph(g, av_manager), true}};
            int                                                         *partition = (int *) malloc(av_manager.get_n_active() * sizeof(int));

            // process the stack
            while (!stack.empty()) {
                Item item = stack.back(); // process first item
                stack.pop_back(); // remove top item

                // load item to process
                KaFFPaGraph<TSerialGraph, TSerialActiveVertexManager> &kaffpa_g   = (*item.g);
                TranslationTable                                      &kaffpa_tt  = (*item.g).tt;
                std::vector<partition_t>                              &identifier = (*item.identifier);

                // get depth info
                size_t      depth           = l - 1 - identifier.size();
                partition_t local_k         = hierarchy[depth];
                partition_t local_k_rem     = k_rem_vec[depth];
                f64         local_imbalance = determine_adaptive_imbalance(global_imbalance, global_g_weight, global_k, kaffpa_g.total_v_weight, local_k_rem, depth + 1);

                // partition the subgraph
                int kaffpa_k        = (int) local_k;
                int kaffpa_edge_cut = 0;
                int mode            = STRONG; // FAST, ECO, STRONG
                kaffpa(&kaffpa_g.n, kaffpa_g.v_weights, kaffpa_g.adj_ptr, kaffpa_g.e_weights, kaffpa_g.adj, &kaffpa_k, &local_imbalance, true, 0, mode, &kaffpa_edge_cut, partition);

                if (depth == 0) {
                    // insert solution
                    u64           offset = 0;
                    for (u64      i      = 0; i < identifier.size(); ++i) { offset += identifier[i] * index_vec[index_vec.size() - 1 - i]; }
                    for (vertex_t u      = 0; u < (vertex_t) kaffpa_g.n; ++u) { p_manager.set(kaffpa_tt.get_o(u), kaffpa_g.v_weights[u], offset + partition[u]); }
                } else {
                    // create the subgraphs and place them in the next stack

                    // collect the number of vertices and edges for each new subgraph
                    std::vector<vertex_t> new_n(local_k, 0);
                    std::vector<vertex_t> new_m(local_k, 0);
                    for (int              u = 0; u < kaffpa_g.n; ++u) {
                        ASSERT(0 <= partition[u] && partition[u] < (int) local_k);
                        new_n[partition[u]] += 1; // increase number of vertices
                        for (int i = kaffpa_g.adj_ptr[u]; i < kaffpa_g.adj_ptr[u + 1]; ++i) {
                            int v = kaffpa_g.adj[i];
                            if (partition[u] == partition[v]) {
                                new_m[partition[u]] += 1; // increase number of edges
                            }
                        }
                    }

                    // create the new subgraphs on the stack
                    for (partition_t i = 0; i < local_k; ++i) {
                        stack.emplace_back(new std::vector<partition_t>(identifier), new KaFFPaGraph<TSerialGraph, TSerialActiveVertexManager>(new_n[i], new_m[i]), true);
                        stack.back().identifier->push_back(i);
                    }

                    // fill the translation tables
                    std::vector<vertex_t> new_us(local_k, 0);
                    for (int              u = 0; u < kaffpa_g.n; ++u) {
                        partition_t p_id = partition[u];
                        size_t      idx  = stack.size() - (local_k - p_id);

                        stack[idx].g->tt.add(kaffpa_tt.get_o(u), new_us[p_id]);
                        new_us[p_id] += 1;
                    }

                    // create the graphs
                    for (int u = 0; u < kaffpa_g.n; ++u) {
                        partition_t p_id = partition[u];
                        size_t      idx  = stack.size() - (local_k - p_id);

                        int sub_u = stack[idx].g->tt.get_n(kaffpa_tt.get_o(u)); // vertex in new graph

                        // set the weight
                        stack[idx].g->set_weight(sub_u, kaffpa_g.v_weights[u]);

                        // set the edges
                        for (int i = kaffpa_g.adj_ptr[u]; i < kaffpa_g.adj_ptr[u + 1]; ++i) {
                            int v = kaffpa_g.adj[i];

                            if (partition[u] == partition[v]) {
                                // add the edge
                                int sub_v    = stack[idx].g->tt.get_n(kaffpa_tt.get_o(v)); // vertex in new graph
                                int curr_end = stack[idx].g->adj_ptr[sub_u + 1];
                                stack[idx].g->adj[curr_end]       = sub_v;
                                stack[idx].g->e_weights[curr_end] = kaffpa_g.e_weights[i];
                                stack[idx].g->adj_ptr[sub_u + 1] += 1;
                            }
                        }
                        if (sub_u + 2 < stack[idx].g->n + 1) {
                            stack[idx].g->adj_ptr[sub_u + 2] = stack[idx].g->adj_ptr[sub_u + 1];
                        }
                    }
                }
            }
            free(partition);
        }

    private:
        template<typename TSerialGraph, typename TSerialActiveVertexManager>
        class KaFFPaGraph {
        public:
            int n = 0;
            int m = 0;

            int *v_weights     = nullptr;
            int total_v_weight = 0;

            int *adj_ptr   = nullptr;
            int *adj       = nullptr;
            int *e_weights = nullptr;
            int last_u     = 0;

            TranslationTable tt;

            KaFFPaGraph(TSerialGraph &g, TSerialActiveVertexManager &av_manager) {
                // remap all active vertices to [0, ..., n-1]
                vertex_t            new_u = 0;
                for (const vertex_t old_u: av_manager) {
                    tt.add(old_u, new_u);
                    new_u += 1;
                }

                // set n and m
                n = (int) new_u;
                m = (int) g.get_m();

                // allocate enough space
                v_weights = (int *) malloc(n * sizeof(int));
                adj_ptr   = (int *) malloc((n + 1) * sizeof(int));
                adj       = (int *) malloc(m * sizeof(int));
                e_weights = (int *) malloc(m * sizeof(int));

                // fill in v_weights
                for (new_u = 0; new_u < (vertex_t) n; ++new_u) {
                    const int w = (int) g.get_weight(tt.get_o(new_u));
                    set_weight((int) new_u, w);
                }

                // fill in adj
                adj_ptr[0] = 0;
                adj_ptr[1] = 0;
                for (new_u = 0; new_u < (vertex_t) n; ++new_u) {
                    vertex_t old_u = tt.get_o(new_u);
                    for (int i     = 0; i < (int) g.size(old_u); ++i) {
                        ASSERT(adj_ptr[new_u] + i < m);
                        const vertex_t old_v = g.neighbor(old_u, i);
                        const int      new_v = (int) tt.get_n(old_v);
                        const int      w     = (int) g.get_weight(old_u, i);
                        add_edge((int) new_u, new_v, w);
                    }
                }
            }

            KaFFPaGraph(const int t_n, const int t_m) {
                n = t_n;
                m = t_m;

                v_weights = (int *) malloc(n * sizeof(int));
                adj_ptr   = (int *) malloc((n + 1) * sizeof(int));
                adj       = (int *) malloc(m * sizeof(int));
                e_weights = (int *) malloc(m * sizeof(int));

                adj_ptr[0] = 0;
                adj_ptr[1] = 0;
            }

            void set_weight(const int u, const int weight) {
                ASSERT(u >= 0 && u < n && weight >= 0);
                v_weights[u] = weight;
                total_v_weight += weight;
            }

            void add_edge(const int u, const int v, const int weight) {
                ASSERT(u >= 0 && u < n && v >= 0 && weight >= 0);
                ASSERT(last_u == u || last_u == u - 1);

                adj[adj_ptr[u + 1]]       = v;
                e_weights[adj_ptr[u + 1]] = weight;
                adj_ptr[u + 1] += 1;

                if (u + 2 < n + 1) {
                    adj_ptr[u + 2] = adj_ptr[u + 1];
                }
                last_u = u;
            }

            ~KaFFPaGraph() {
                free(v_weights);
                free(adj_ptr);
                free(adj);
                free(e_weights);
            }
        };

        template<typename TSerialGraph, typename TSerialActiveVertexManager>
        struct Item {
            std::vector<partition_t>                              *identifier;
            KaFFPaGraph<TSerialGraph, TSerialActiveVertexManager> *g;
            bool                                                  to_delete;

            Item(std::vector<partition_t> *t_identifier, KaFFPaGraph<TSerialGraph, TSerialActiveVertexManager> *t_g, bool t_to_delete) {
                identifier = t_identifier;
                g          = t_g;
                to_delete  = t_to_delete;
            }
        };

        f64 determine_adaptive_imbalance(const f64 global_imbalance,
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

#endif //HEIDELBERGPROCESSMAPPING_GLOBAL_MULTISECTION_H
