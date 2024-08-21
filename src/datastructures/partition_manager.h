#ifndef SERIALPROCESSMAPPING_PARTITION_MANAGER_H
#define SERIALPROCESSMAPPING_PARTITION_MANAGER_H

#include "../utility/definitions.h"
#include "../utility/macros.h"
#include "../utility/utils.h"
#include "graph.h"
#include "../utility/qap.h"
#include "iterators/active_vertex_iterator.h"
#include "distance_oracle.h"
#include "boundary_vertex_manager.h"
#include "quotient_graph.h"
#include "iterators/boundary_vertex_iterator.h"

namespace SPM {

    class PartitionManager {
    private:
        Graph *p_g = nullptr;
        u64 k = 0;

        // actual partition
        std::vector<partition_t> partition;

        // boundary vertices
        BoundaryVertexManager bvm;

        // partition weights
        std::vector<u64> pweights;

        // quotient graph
        QuotientGraph qg;

    public:
        PartitionManager() : bvm(partition) {}

        void initialize(Graph *t_g, u64 t_k) {
            p_g = t_g;
            k = t_k;

            // actual partition
            partition.resize(p_g->get_n());

            // boundary vertices
            bvm.initialize(t_g, t_k);

            // partition weights
            pweights.resize(k);

            // quotient graph
            qg.initialize(t_k);
        }

        u64 get_k() const {
            return k;
        }

        vertex_t &operator[](vertex_t u) {
            Graph &g = *p_g;
            ASSERT(u < g.get_n());

            return partition[u];
        }

        const vertex_t &operator[](vertex_t u) const {
            Graph &g = *p_g;
            ASSERT(u < g.get_n());

            return partition[u];
        }

        std::vector<vertex_t> get_partition() {
            ASSERT(check_pweights());
            return partition;
        }

        void set_partition(const std::vector<partition_t> &t_partition) {
            partition = t_partition;
        }

        std::vector<u64> get_pweights() {
            ASSERT(check_pweights());
            return pweights;
        }

        bool is_boundary_vertex(vertex_t u) const {
            return bvm.is_boundary_vertex(u);
        }

        BoundaryVertexManager &get_bvm() {
            return bvm;
        }

        QuotientGraph &get_qg() {
            return qg;
        }

        u64 get_pweight(vertex_t u_id) const {
            ASSERT(u_id < k);

            return pweights[u_id];
        }

        void init_after_partition() {
            Graph &g = *p_g;

            // Determine boundary vertices
            for (ActiveVertexIterator avi(g); avi.not_end(); avi.next()) {
                vertex_t u = avi.get();
                partition_t u_id = partition[u];
                bvm.insert(u, u_id);
            }
            ASSERT(bvm.holds_all_boundary_vertices());

            // Determine partition weights
            std::fill(pweights.begin(), pweights.end(), 0);
            for (ActiveVertexIterator avi(g); avi.not_end(); avi.next()) {
                vertex_t u = avi.get();
                pweights[partition[u]] += g.get_vertex_weight(u);
            }
            HEAVYASSERT(check_pweights());

            // Determine quotient graph
            for (BoundaryVertexIterator bvi(bvm); bvi.not_end(); bvi.next()) {
                vertex_t u = bvi.get();
                partition_t u_id = partition[u];
                for (EdgeW &e: g[u]) {
                    partition_t v_id = partition[e.v];
                    if (u_id != v_id) {
                        // different partitions so u is a boundary vertex
                        qg.add_edge(u_id, v_id, e.w);
                    }
                }
            }

            ASSERT(bvm.holds_all_boundary_vertices());
        }

        void uncontract_edge(vertex_t u, vertex_t v) {
            Graph &g = *p_g;

            // put v in same partition as u
            partition[v] = partition[u];

            // uncontract in graph and boundary vertices
            bvm.pre_uncontract(u);
            g.uncontract_edge(u, v);
            bvm.after_uncontract(u, v);

            // Note: Weights dont change
        }

        void move(vertex_t u, partition_t new_u_id) {
            Graph &g = *p_g;
            ASSERT(u < g.get_n());
            ASSERT(new_u_id < k);
            ASSERT(g.get_vertex_state(u) == 1);

            partition_t old_u_id = partition[u];
            if (old_u_id == new_u_id) {
                // do nothing
                return;
            }

            // update weights
            weight_t u_weight = g.get_vertex_weight(u);
            ASSERT(pweights[old_u_id] >= u_weight);
            pweights[old_u_id] -= u_weight;
            pweights[new_u_id] += u_weight;

            // update boundary vertices
            bvm.move(u, new_u_id);

            // update quotient graph
            for (EdgeW &e: (*p_g)[u]) {
                partition_t v_id = partition[e.v];
                if (old_u_id != v_id) {
                    qg.remove_edge(old_u_id, v_id, e.w);
                }
                if (new_u_id != v_id) {
                    qg.add_edge(new_u_id, v_id, e.w);
                }
            }

            // move u
            partition[u] = new_u_id;

            ASSERT(bvm.holds_all_boundary_vertices());
        }

    private:
        u64 get_bf_pweight(vertex_t u_id) const {
            Graph &g = *p_g;
            u64 w = 0;
            for (ActiveVertexIterator avi(g); avi.not_end(); avi.next()) {
                vertex_t v = avi.get();
                vertex_t v_id = partition[v];
                weight_t v_w = g.get_vertex_weight(v);
                if (v_id == u_id) {
                    w += v_w;
                }
            }
            return w;
        }

        bool check_pweights() const {
            for (vertex_t i = 0; i < k; ++i) {
                if (pweights[i] != get_bf_pweight(i)) {
                    return false;
                }
            }
            return true;
        }
    };

}

#endif //SERIALPROCESSMAPPING_PARTITION_MANAGER_H
