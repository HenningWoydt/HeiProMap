#include <gtest/gtest.h>
#include "src/utility/memory_stack.h"

namespace HeiProMap {

TEST(MemoryStackTest, BasicAllocation) {
    MemoryStack stack;
    stack.ensure(1024);
    
    int* p1 = (int*)stack.get_memory(sizeof(int));
    *p1 = 42;
    EXPECT_EQ(*p1, 42);
    
    double* p2 = (double*)stack.get_memory(sizeof(double));
    *p2 = 3.14;
    EXPECT_EQ(*p2, 3.14);
}

TEST(MemoryStackTest, ClearAndReuse) {
    MemoryStack stack;
    stack.ensure(128);
    
    int* p1 = (int*)stack.get_memory(sizeof(int));
    *p1 = 100;
    
    stack.clear();
    
    int* p2 = (int*)stack.get_memory(sizeof(int));
    // p2 should be at the same memory location as p1
    *p2 = 200;
    
    EXPECT_EQ(*p1, 200);
}

TEST(MemoryStackTest, Reallocation) {
    MemoryStack stack;
    stack.ensure(64);
    
    void* p1 = stack.get_memory(32);
    ASSERT_NE(p1, nullptr);
    
    // This should trigger a reallocation
    stack.ensure(128);
    
    void* p2 = stack.get_memory(32);
    ASSERT_NE(p2, nullptr);
}

} // namespace HeiProMap
