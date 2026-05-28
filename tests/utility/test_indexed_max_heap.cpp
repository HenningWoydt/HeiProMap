#include <gtest/gtest.h>
#include "src/utility/indexed_max_heap.h"

namespace HeiProMap {

TEST(IndexedMaxHeapTest, BasicFunctionality) {
    IndexedMaxHeap<int> heap;
    heap.initialize(10);
    
    heap.push(0, 10);
    heap.push(1, 20);
    heap.push(2, 15);
    
    EXPECT_EQ(heap.size(), 3);
    EXPECT_EQ(heap.top_key(), 1);
    EXPECT_EQ(heap.top(), 20);
    
    heap.pop();
    EXPECT_EQ(heap.top_key(), 2);
    EXPECT_EQ(heap.top(), 15);
    
    heap.update(0, 30);
    EXPECT_EQ(heap.top_key(), 0);
    EXPECT_EQ(heap.top(), 30);
    
    heap.clear();
    EXPECT_TRUE(heap.empty());
}

} // namespace HeiProMap
