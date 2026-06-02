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

#include <gtest/gtest.h>
#include <vector>

#include "src/datastructures/csr_graph.h"
#include "src/datastructures/subgraph_extractor.h"
#include "src/utility/translation_table.h"

namespace HeiProMap {

class SubgraphExtractorTest : public ::testing::Test {
protected:
    CSRGraph g;
    std::vector<u8> side;

    void SetUp() override {
        // Create a simple graph: 0-1, 1-2, 2-3 (Path P4)
        // Vertex weights: 1, 2, 3, 4
        // Edge weights: 10, 20, 30
        vertex_t n = 4;
        vertex_t m = 6;
        weight_t v_weights[] = {1, 2, 3, 4};
        size_t neighborhoods[] = {0, 1, 3, 5, 6};
        vertex_t edges_v[] = {1, 0, 2, 1, 3, 2};
        weight_t edges_w[] = {10, 10, 20, 20, 30, 30};

        g.initialize(n, m, v_weights, neighborhoods, edges_w, edges_v);
        
        // Split: {0, 1} vs {2, 3}
        side = {0, 0, 1, 1};
    }
};

TEST_F(SubgraphExtractorTest, ExtractSide0) {
    CSRGraph sub_g;
    TranslationTable<vertex_t> tt;
    SubgraphExtractor::extract(g, side, 0, sub_g, tt);

    // Vertices 0 and 1 are in sub_g
    EXPECT_EQ(sub_g.n, 2);
    EXPECT_EQ(sub_g.g_weight, 1 + 2);
    
    // Check translation
    // Mapping could be 0->0, 1->1 (based on original order)
    vertex_t sub_0 = tt.get_n(0);
    vertex_t sub_1 = tt.get_n(1);
    EXPECT_NE(sub_0, sub_1);
    
    EXPECT_EQ(tt.get_o(sub_0), 0);
    EXPECT_EQ(tt.get_o(sub_1), 1);
    
    EXPECT_EQ(sub_g.v_weights[sub_0], 1);
    EXPECT_EQ(sub_g.v_weights[sub_1], 2);

    // Edges in sub_g: only 0-1 remains
    // 0 has neighbor 1, 1 has neighbor 0
    EXPECT_EQ(sub_g.m, 2);
    EXPECT_EQ(sub_g.neighborhoods[1] - sub_g.neighborhoods[0], 1);
    EXPECT_EQ(sub_g.neighborhoods[2] - sub_g.neighborhoods[1], 1);
    
    EXPECT_EQ(sub_g.edges_v[sub_g.neighborhoods[sub_0]], sub_1);
    EXPECT_EQ(sub_g.edges_w[sub_g.neighborhoods[sub_0]], 10);
    
    EXPECT_EQ(sub_g.edges_v[sub_g.neighborhoods[sub_1]], sub_0);
    EXPECT_EQ(sub_g.edges_w[sub_g.neighborhoods[sub_1]], 10);
}

TEST_F(SubgraphExtractorTest, ExtractSide1) {
    CSRGraph sub_g;
    TranslationTable<vertex_t> tt;
    SubgraphExtractor::extract(g, side, 1, sub_g, tt);

    // Vertices 2 and 3 are in sub_g
    EXPECT_EQ(sub_g.n, 2);
    EXPECT_EQ(sub_g.g_weight, 3 + 4);
    
    vertex_t sub_2 = tt.get_n(2);
    vertex_t sub_3 = tt.get_n(3);
    
    EXPECT_EQ(tt.get_o(sub_2), 2);
    EXPECT_EQ(tt.get_o(sub_3), 3);

    // Edges in sub_g: only 2-3 remains
    EXPECT_EQ(sub_g.m, 2);
    EXPECT_EQ(sub_g.edges_v[sub_g.neighborhoods[sub_2]], sub_3);
    EXPECT_EQ(sub_g.edges_w[sub_g.neighborhoods[sub_2]], 30);
}

TEST_F(SubgraphExtractorTest, ExtractEmptySide) {
    CSRGraph sub_g;
    TranslationTable<vertex_t> tt;
    std::vector<u8> all_zero_side = {0, 0, 0, 0};
    SubgraphExtractor::extract(g, all_zero_side, 1, sub_g, tt);

    EXPECT_EQ(sub_g.n, 0);
    EXPECT_EQ(sub_g.m, 0);
    EXPECT_EQ(sub_g.g_weight, 0);
}

TEST_F(SubgraphExtractorTest, DisconnectedSubgraph) {
    // 0-1, 2-3. Split {0, 2} vs {1, 3}
    std::vector<u8> checkerboard_side = {0, 1, 0, 1};
    CSRGraph sub_g;
    TranslationTable<vertex_t> tt;
    SubgraphExtractor::extract(g, checkerboard_side, 0, sub_g, tt);

    // Subgraph should have vertices 0 and 2, but NO edges because 0 was connected to 1 and 2 was connected to 1 and 3.
    EXPECT_EQ(sub_g.n, 2);
    EXPECT_EQ(sub_g.m, 0);
    EXPECT_EQ(sub_g.v_weights[tt.get_n(0)], 1);
    EXPECT_EQ(sub_g.v_weights[tt.get_n(2)], 3);
}

} // namespace HeiProMap
