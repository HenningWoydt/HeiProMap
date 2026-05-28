#include <gtest/gtest.h>
#include "src/utility/hungarian.h"

namespace HeiProMap {

TEST(HungarianTest, BasicAssignment) {
    std::vector<std::vector<weight_t>> cost_matrix = {
        {9, 2, 7, 8},
        {6, 4, 3, 7},
        {5, 8, 1, 8},
        {7, 6, 9, 4}
    };
    
    auto result = solve_hungarian(cost_matrix);
    std::vector<int> expected = {1, 0, 2, 3}; // 2+6+1+4 = 13
    
    EXPECT_EQ(result, expected);
    
    weight_t total_cost = 0;
    for(size_t i = 0; i < result.size(); ++i) {
        total_cost += cost_matrix[i][result[i]];
    }
    EXPECT_EQ(total_cost, 13);
}

TEST(HungarianTest, IdentityMatrix) {
    std::vector<std::vector<weight_t>> cost_matrix = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };
    
    auto result = solve_hungarian(cost_matrix);
    std::vector<int> expected = {1, 0, 2}; // 0+0+1 = 1
                                           // or 1,2,0 -> 0+0+0 = 0
                                           // or 2,0,1 -> 0+0+0 = 0
    
    weight_t total_cost = 0;
    for(size_t i = 0; i < result.size(); ++i) {
        total_cost += cost_matrix[i][result[i]];
    }
    EXPECT_EQ(total_cost, 0);
}

} // namespace HeiProMap
