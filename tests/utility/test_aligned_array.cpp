#include <gtest/gtest.h>
#include "src/utility/aligned_array.h"

namespace HeiProMap {

TEST(AlignedArrayTest, BasicFunctionality) {
    AlignedArray<int> arr;
    arr.initialize(100, 42);
    
    EXPECT_GE(arr.size(), 100);
    for (size_t i = 0; i < 100; ++i) {
        EXPECT_EQ(arr[i], 42);
    }
    
    arr[0] = 10;
    EXPECT_EQ(arr[0], 10);
}

TEST(AlignedArrayTest, CopyAndMove) {
    AlignedArray<int> arr1;
    arr1.initialize(10, 1);
    
    // Copy
    AlignedArray<int> arr2 = arr1;
    EXPECT_EQ(arr2.size(), arr1.size());
    EXPECT_EQ(arr2[0], 1);
    
    arr2[0] = 2;
    EXPECT_EQ(arr1[0], 1); // Ensure it's a deep copy
    
    // Move
    AlignedArray<int> arr3 = std::move(arr1);
    EXPECT_EQ(arr3[0], 1);
    EXPECT_EQ(arr1.get_ptr(), nullptr);
}

} // namespace HeiProMap
