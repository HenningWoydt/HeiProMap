#include <gtest/gtest.h>
#include <vector>
#include "src/partitioning/recursive_bisection.h"
#include "src/datastructures/csr_graph.h"
#include "src/datastructures/partition_manager.h"

namespace HeiProMap {

TEST(RecursiveBisectionTest, GridGraphRefinement) {
    // Create a 10x10 grid graph (100 vertices, ~180 edges)
    const vertex_t dim = 10;
    const vertex_t n = dim * dim;
    std::vector<weight_t> v_weights(n, 1);
    std::vector<size_t> neighborhoods;
    std::vector<vertex_t> edges_v;
    std::vector<weight_t> edges_w;

    neighborhoods.push_back(0);
    for (vertex_t i = 0; i < dim; ++i) {
        for (vertex_t j = 0; j < dim; ++j) {
            vertex_t u = i * dim + j;
            if (i > 0) { edges_v.push_back((i - 1) * dim + j); edges_w.push_back(1); }
            if (i < dim - 1) { edges_v.push_back((i + 1) * dim + j); edges_w.push_back(1); }
            if (j > 0) { edges_v.push_back(i * dim + (j - 1)); edges_w.push_back(1); }
            if (j < dim - 1) { edges_v.push_back(i * dim + (j + 1)); edges_w.push_back(1); }
            neighborhoods.push_back(edges_v.size());
        }
    }

    CSRGraph g;
    g.initialize(n, edges_v.size(), v_weights.data(), neighborhoods.data(), edges_w.data(), edges_v.data());

    const partition_t k = 8;
    PartitionManager pm;
    pm.initialize(n, k, n);

    RecursiveBisectionPartitioner partitioner;
    RecursiveBisectionConfiguration config;
    config.use_full_refine = true;
    config.kappa = 2;
    config.method = BisectionMethod::HYBRID;
    
    // LP refinement
    config.lp_config.enabled = true;
    config.lp_config.max_iteration = 5;

    partitioner.partition(g, pm, k, 42, 0.03, config);

    // Check balance
    weight_t max_w = 0;
    for (partition_t i = 0; i < k; ++i) {
        max_w = std::max(max_w, pm.get_bweight(i));
    }
    
    f64 expected_max = std::ceil((1.03) * (f64)n / (f64)k);
    EXPECT_LE(max_w, (weight_t)expected_max);

    // Check if all vertices are assigned
    size_t total_assigned = 0;
    for (partition_t i = 0; i < k; ++i) {
        total_assigned += pm.size(i);
    }
    EXPECT_EQ(total_assigned, n);
}

} // namespace HeiProMap
