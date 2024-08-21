#ifndef SERIALPROCESSMAPPING_BOUNDARY_VERTEX_MANAGER_H
#define SERIALPROCESSMAPPING_BOUNDARY_VERTEX_MANAGER_H

#include "../utility/definitions.h"
#include "../utility/macros.h"
#include "../utility/utils.h"
#include "graph.h"
#include "../utility/qap.h"
#include "iterators/active_vertex_iterator.h"
#include "distance_oracle.h"

namespace SPM {

    /**
     * A class that allows for quick access to the current boundary vertices,
     * either for the complete graph (BoundaryVertexIterator) or for individual
     * blocks (BlockBoundaryVertexIterator).
     *
     * HOW TO USE:
     *  - Use BoundaryVertexIterator to iterate across all boundary vertices.
     *  - Use BlockBoundaryVertexIterator to iterate across all boundary
     *    vertices of one block.
     *  - All modifications to the boundary happen instantly! If you move/remove
     *    vertices during an iteration newly discovered boundary vertices will
     *    occur at the end of the iteration and removed vertices will not appear.
     *    The same vertex can appear multiple times, but never twice in the same
     *    block.
     *  - move(u, b) moves vertex u to block b.
     *  - remove(u) removes vertex u from the boundary.
     *  - tidy_up() removes duplicates and unnecessary entries for all blocks.
     *  - tidy_up(u) removes duplicates and unnecessary entries for block b.
     *
     *
     * DEVELOPER NOTES:
     *  - BoundaryVertexIterator and BlockBoundaryVertexIterator are additionally
     *    responsible for cleaning the vectors of vertices that are wrongly placed.
     */
    class BoundaryVertexManager {
    private:
        Graph *p_g = nullptr;
        const std::vector<partition_t> &partition; // current partition
        u64 k = 0;

        // Current number of edges to other boundary vertices
        std::vector<u64> n_boundary_edges;

        // Current boundary for each block. Can be inconsistent!
        std::vector<std::vector<vertex_t>> boundaries;

    public:
        explicit BoundaryVertexManager(std::vector<partition_t> &t_partition) : partition(t_partition) {};

        /**
         * Initializes the object.
         *
         * @param t_g The graph.
         * @param t_k The number of blocks.
         */
        void initialize(Graph *t_g, u64 t_k) {
            p_g = t_g;
            k = t_k;

            n_boundary_edges.resize(p_g->get_n(), 0);
            boundaries.resize(k);
        }

        /**
         * Inserts the vertex u to the block b. Vertex u has to be in no block.
         * If vertex u is currently in another block, then use move(u, b)!
         *
         * @param u The vertex u.
         * @param u_id The block b.
         */
        void insert(vertex_t u, partition_t u_id){
            ASSERT(u < p_g->get_n());
            ASSERT(u_id < k);

            // add connections to other boundary vertices
            bool is_boundary = false;
            for(EdgeW &e : (*p_g)[u]){
                partition_t v_id = partition[e.v];
                if(v_id != u_id){
                    // u and v are boundary vertices in different blocks
                    n_boundary_edges[u] += 1;
                    is_boundary = true;
                }
            }

            if(is_boundary){
                // put u in its respective boundary queue
                boundaries[u_id].emplace_back(u);
            }
        }

        /**
         * Moves the vertex u to the block b. Vertex u has to be in another
         * block. If vertex u is currently not in another block, then use
         * insert(u, b)!
         *
         * @param u The vertex u.
         * @param new_b The block b.
         */
        void move(vertex_t u, partition_t new_b){
            ASSERT(u < p_g->get_n());
            ASSERT(new_b < k);

            partition_t old_b = partition[u];
            if(old_b == new_b){
                // u is already in b, do nothing
                return;
            }

            // put u in b
            emplace_if_not_exists(new_b, u);

            // new boundary vertices could be discovered and other could be removed
            for(EdgeW &e : (*p_g)[u]){
                partition_t v_id = partition[e.v];
                if (v_id == new_b){
                    // u was moved to the same block as v, both loose 1 boundary edge
                    ASSERT(n_boundary_edges[u] > 0);
                    ASSERT(n_boundary_edges[e.v] > 0);
                    n_boundary_edges[u] -= 1;
                    n_boundary_edges[e.v] -= 1;
                } else if (v_id == old_b){
                    // u was moved to a different block as v, both gain 1 boundary edge
                    n_boundary_edges[u] += 1;
                    n_boundary_edges[e.v] += 1;
                    emplace_if_not_exists(v_id, e.v);
                }
                // else, v and b are in different blocks and still connected, nothing changes
            }
        }

        /**
         * Use this function before uncontraction of u and v in the graph.
         *
         * @param u The vertex u.
         */
        void pre_uncontract(vertex_t u){
            partition_t u_id = partition[u];
            if(n_boundary_edges[u] == 0){
                // u is not a boundary edge, so nothing to do
                return;
            }

            // decrement all boundary edges from u
            for(EdgeW &e: (*p_g)[u]){
                partition_t v_id = partition[e.v];
                if(v_id != u_id){
                    // is a boundary edge, decrement
                    n_boundary_edges[u] -= 1;
                    n_boundary_edges[e.v] -= 1;
                    // NOTE: after_contract() will be called shortly, this will
                    // reconnect u and e.v, if they are still neighbors
                }
            }
        }

        /**
         * Use this function after uncontraction of u and v in the graph.
         *
         * @param u The vertex u.
         * @param v The vertex v.
         */
        void after_uncontract(vertex_t u, vertex_t v){
            partition_t u_id = partition[u];
            partition_t v_id = partition[v];
            ASSERT(u_id == v_id);

            // increment boundary edges from u
            for(EdgeW &e: (*p_g)[u]){
                partition_t ev_id = partition[e.v];
                if(ev_id != u_id){
                    // is a boundary edge, increment
                    n_boundary_edges[u] += 1;
                    n_boundary_edges[e.v] += 1;
                }
            }

            // increment boundary edges from v
            for(EdgeW &e: (*p_g)[v]){
                partition_t ev_id = partition[e.v];
                if(ev_id != v_id){
                    // is a boundary edge, increment
                    n_boundary_edges[v] += 1;
                    n_boundary_edges[e.v] += 1;
                }
            }

            // check if v becomes a boundary vertex
            if(n_boundary_edges[v] > 0){
                boundaries[v_id].emplace_back(v);
            }
        }

        /**
         * Determines whether u is a boundary vertex.
         *
         * @param u The vertex u.
         * @return True if u is a boundary vertex, false else.
         */
        bool is_boundary_vertex(vertex_t u) const {
            return n_boundary_edges[u] > 0;
        }

        /**
         * Returns the block of vertex u.
         *
         * @param u The vertex u.
         * @return The block and BLOCK_TOMBSTONE if u is not a boundary vertex.
         */
        partition_t get_b(vertex_t u){
            return partition[u];
        }

        /**
         * Removes the vertex u from the boundary.
         *
         * @param u The vertex u.
         */
         /*
        void remove(vertex_t u){
            ASSERT(u < p_g->get_n());
            ids[u] = BLOCK_TOMBSTONE; // mark as not in any block
        }
          */

        /**
         * Removes wrongly placed and duplicate entries from block b.
         *
         * @param b The block b.
         */
        void tidy_up(partition_t b){
            // remove unnecessary entries
            for(size_t i = 0; i < boundaries[b].size(); ++i){
                if(partition[i] != b){
                    boundaries[b] = boundaries.back();
                    boundaries.pop_back();
                    i -= 1;
                    continue;
                }
            }

            // remove duplicate entries
            std::sort(boundaries[b].begin(), boundaries[b].end());
            boundaries[b].erase(std::unique(boundaries[b].begin(), boundaries[b].end()), boundaries[b].end());
        }

        /**
         * Removes wrongly placed and duplicate entries from all blocks.
         */
        void tidy_up(){
            for(partition_t b = 0; b < k; ++b){
                tidy_up(b);
            }
        }

        const std::vector<partition_t> &get_partition(){
            return partition;
        }

        const std::vector<u64> &get_n_boundary_edges(){
            return n_boundary_edges;
        }

        std::vector<std::vector<vertex_t>> &get_boundaries(){
            return boundaries;
        }

        partition_t get_k(){
            return k;
        }

        bool holds_all_boundary_vertices(){
            for (ActiveVertexIterator avi((*p_g)); avi.not_end(); avi.next()) {
                vertex_t u = avi.get();
                if(n_boundary_edges[u] != n_boundary_edges_via_graph(u)){
                    return false;
                }
            }
            return true;
        }

    private:
        /**
         * Determines whether a vertex is currently a boundary vertex.
         * Determines via own structure.
         *
         * @param u The vertex.
         * @return True if it is a boundary vertex, false else.
         */
        bool is_boundary_vertex_via_self(vertex_t u){
            if(n_boundary_edges[u] == 0){
                // u must have boundary edges
                return false;
            }
            std::cout << u << " " << partition[u] << " " << n_boundary_edges[u] << std::endl;
            std::cout << to_string(boundaries[partition[u]]) << std::endl;
            ASSERT(std::find(boundaries[partition[u]].begin(), boundaries[partition[u]].end(), u) != boundaries[partition[u]].end());
            return true;
        }

        bool is_boundary_vertex_via_graph(vertex_t u){
            partition_t u_id = partition[u];
            for(EdgeW &e : (*p_g)[u]){
                partition_t v_id = partition[e.v];
                if(u_id != v_id){
                    return true;
                }
            }
            return false;
        }

        u64 n_boundary_edges_via_graph(vertex_t u){
            u64 count = 0;
            partition_t u_id = partition[u];
            for(EdgeW &e : (*p_g)[u]){
                partition_t v_id = partition[e.v];
                if(u_id != v_id){
                    count += 1;
                }
            }
            return count;
        }

        void emplace_if_not_exists(partition_t b, vertex_t u){
            for(size_t i = 0; i< boundaries[b].size(); ++i){
                if(boundaries[b][i] == u){
                    return;
                }
                if(partition[boundaries[b][i]] != b){
                    boundaries[b][i] = boundaries[b].back();
                    boundaries[b].pop_back();
                    i -= 1;
                    continue;
                }
            }
            boundaries[b].emplace_back(u);
        }
    };

}

#endif //SERIALPROCESSMAPPING_BOUNDARY_VERTEX_MANAGER_H
