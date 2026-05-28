#include <gtest/gtest.h>
#include "src/utility/indexed_update_heap.h"

namespace HeiProMap {

TEST(IndexedUpdateHeapTest, PushAndPop) {
    IndexedUpdateHeap heap;
    heap.initialize(10);

    heap.push(0, 1, 10);
    heap.push(1, 1, 20);
    heap.push(2, 1, 15);

    EXPECT_EQ(heap.size(), 3);
    EXPECT_EQ(heap.top_u(), 1);
    EXPECT_EQ(heap.top_qap_delta(), 20);

    heap.pop();
    EXPECT_EQ(heap.top_u(), 2);
    EXPECT_EQ(heap.top_qap_delta(), 15);
}

TEST(IndexedUpdateHeapTest, Update) {
    IndexedUpdateHeap heap;
    heap.initialize(10);

    heap.push(0, 1, 10);
    heap.push(1, 1, 5);

    EXPECT_EQ(heap.top_u(), 0);

    heap.update(1, 1, 30);
    EXPECT_EQ(heap.top_u(), 1);
    EXPECT_EQ(heap.top_qap_delta(), 30);
}

TEST(IndexedUpdateHeapTest, PushUpdate) {
    IndexedUpdateHeap heap;
    heap.initialize(10);
    
    heap.push_update(0, 1, 10);
    EXPECT_TRUE(heap.entry_exists(0));
    EXPECT_EQ(heap.top_u(), 0);
    
    heap.push_update(0, 1, 5); // Update existing
    EXPECT_EQ(heap.top_qap_delta(), 5);
    
    heap.push_update(1, 1, 20); // Push new
    EXPECT_EQ(heap.top_u(), 1);
}

TEST(IndexedUpdateHeapTest, Clear) {
    IndexedUpdateHeap heap;
    heap.initialize(10);

    heap.push(0, 1, 10);
    heap.push(1, 1, 20);
    
    EXPECT_FALSE(heap.empty());
    
    heap.clear();
    EXPECT_TRUE(heap.empty());
    
    // After clear, old entries should not exist
    EXPECT_FALSE(heap.entry_exists(0));
    
    heap.push(0, 1, 5);
    EXPECT_TRUE(heap.entry_exists(0));
}

} // namespace HeiProMap
