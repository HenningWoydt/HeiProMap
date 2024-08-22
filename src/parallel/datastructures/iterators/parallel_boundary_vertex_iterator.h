#ifndef SERIALPROCESSMAPPING_PARALLEL_BOUNDARY_VERTEX_ITERATOR_H
#define SERIALPROCESSMAPPING_PARALLEL_BOUNDARY_VERTEX_ITERATOR_H

#include <vector>
#include <fstream>
#include <numeric>

#include "../../../utility/definitions.h"
#include "../../../datastructures/boundary_vertex_manager.h"

namespace HeiProMap {

    class ParallelBoundaryVertexIterator {
    private:
        u64 n_threads;

        partition_t k;

        size_t b_idx;
        std::vector<size_t> indices;

        const std::vector<partition_t> &partition; // ref to vector in BoundaryVertexManager
        const std::vector<u64> &n_boundary_edges;
        std::vector<std::vector<vertex_t>> &boundaries; // ref to vector in BoundaryVertexManager

    public:
        explicit ParallelBoundaryVertexIterator(BoundaryVertexManager &bvm, u64 t_n_threads) : partition(bvm.get_partition()), n_boundary_edges(bvm.get_n_boundary_edges()), boundaries(bvm.get_boundaries()) {
            n_threads = t_n_threads;
            k = bvm.get_k();
            b_idx = 0;
            indices.resize(k, 0);
        }

        vertex_t get() {
            return boundaries[b_idx][indices[b_idx]];
        }

        void next() {
            indices[b_idx] += 1;
        }

        bool not_end() {
            for(partition_t i = 0; i < k; ++i){
                while(indices[b_idx] < boundaries[b_idx].size()){
                    if(n_boundary_edges[boundaries[b_idx][indices[b_idx]]] == 0){
                        boundaries[b_idx][indices[b_idx]] = boundaries[b_idx].back();
                        boundaries[b_idx].pop_back();
                        continue;
                    }
                    return true;
                }

                b_idx = (b_idx + 1) % k;
            }
            return false;
        }
    };
}

#endif //SERIALPROCESSMAPPING_PARALLEL_BOUNDARY_VERTEX_ITERATOR_H
