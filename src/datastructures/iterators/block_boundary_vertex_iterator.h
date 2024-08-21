#ifndef SERIALPROCESSMAPPING_BLOCK_BOUNDARY_VERTEX_ITERATOR_H
#define SERIALPROCESSMAPPING_BLOCK_BOUNDARY_VERTEX_ITERATOR_H

#include <vector>
#include <fstream>
#include <numeric>

#include "../../utility/definitions.h"
#include "../../utility/utils.h"
#include "../../utility/macros.h"
#include "../graph.h"
#include "../partition_manager.h"

namespace SPM {

    class BlockBoundaryVertexIterator {
    private:
        partition_t k;
        size_t idx;

        const std::vector<partition_t> &partition; // ref to vector in BoundaryVertexManager
        const std::vector<u64> &n_boundary_edges;
        std::vector<vertex_t> &boundary; // ref to vector in BoundaryVertexManager
    public:
        BlockBoundaryVertexIterator(BoundaryVertexManager &bvm, partition_t b) : partition(bvm.get_partition()), n_boundary_edges(bvm.get_n_boundary_edges()), boundary(bvm.get_boundaries()[b]) {
            k = bvm.get_k();
            idx = 0;
        }

        vertex_t get() {
            return boundary[idx];
        }

        void next() {
            idx += 1;
        }

        bool not_end() {
            while (idx < boundary.size()) {
                if (n_boundary_edges[boundary[idx]] == 0) {
                    boundary[idx] = boundary.back();
                    boundary.pop_back();
                    continue;
                }
                return true;
            }
            return false;
        }

    };

}

#endif //SERIALPROCESSMAPPING_BLOCK_BOUNDARY_VERTEX_ITERATOR_H
