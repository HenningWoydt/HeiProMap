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

#include "test_utils.h"

namespace HeiProMap {
    void graphs_are_equal(const ISerialGraph& g1, const ISerialGraph& g2) {
        // compare number of vertices
        EXPECT_EQ(g1.get_n(), g2.get_n());

        // compare number of edges
        EXPECT_EQ(g1.get_m(), g2.get_m());

        // compare graph weight
        EXPECT_EQ(g1.get_weight(), g2.get_weight());

        // compare vertex weights
        for (vertex_t u = 0; u < g1.get_n(); u++) {
            EXPECT_EQ(g1.get_weight(u), g2.get_weight(u));
        }

        // compare neighborhood sizes
        for (vertex_t u = 0; u < g1.get_n(); u++) {
            EXPECT_EQ(g1.size(u), g2.size(u));
        }

        // compare neighborhoods
        std::vector<EdgeVW> g1_neighborhood;
        std::vector<EdgeVW> g2_neighborhood;
        for (vertex_t u = 0; u < g1.get_n(); u++) {
            g1_neighborhood.clear();
            g2_neighborhood.clear();
            for (size_t i = 0; i < g1.size(u); i++) {
                vertex_t v = g1.neighbor(u, i);
                weight_t w = g1.get_weight(u, i);
                g1_neighborhood.emplace_back(v, w);
            }
            for (size_t i = 0; i < g2.size(u); i++) {
                vertex_t v = g2.neighbor(u, i);
                weight_t w = g2.get_weight(u, i);
                g2_neighborhood.emplace_back(v, w);
            }
            std::sort(g1_neighborhood.begin(), g1_neighborhood.end());
            std::sort(g2_neighborhood.begin(), g2_neighborhood.end());
            EXPECT_EQ(g1_neighborhood, g2_neighborhood);
        }
    }
}
