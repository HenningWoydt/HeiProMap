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

#ifndef HEIPROMAP_MACROS_H
#define HEIPROMAP_MACROS_H

#include <iostream>
#include <string>

#if defined(__GNUC__) || defined(__clang__)
#define HEIPROMAP_RESTRICT __restrict__
#define HEIPROMAP_ASSUME_ALIGNED(ptr, alignment) ((decltype(ptr))__builtin_assume_aligned((ptr), (alignment)))
#elif defined(_MSC_VER)
#define HEIPROMAP_RESTRICT __restrict
#define HEIPROMAP_ASSUME_ALIGNED(ptr, alignment) (ptr)
#else
#define HEIPROMAP_RESTRICT
#define HEIPROMAP_ASSUME_ALIGNED(ptr, alignment) (ptr)
#endif

namespace HeiProMap {
    #ifndef ENABLE_PROFILER
    #define ENABLE_PROFILER true
    #endif

    #ifndef ASSERT_ENABLED
    #define ASSERT_ENABLED false
    #endif

    #ifndef HEAVYASSERT_ENABLED
    #define HEAVYASSERT_ENABLED false
    #endif

    #if (ASSERT_ENABLED)
    // Use ASSERT for quick operations like O(1) operations, for other Asserts use HEAVYASSERT
    #define ASSERT(condition) if(!(condition)) {std::cerr << "Error in file " << __FILE__ << " in function " << __FUNCTION__ << " at line " << __LINE__ << "!" << std::endl; abort(); } ((void)0)
    #else
    #define ASSERT(condition) ((void)0)
    #endif


    #if (HEAVYASSERT_ENABLED)
    // Use HEAVYASSERT for expensive operations like O(n), O(n^2) operations, for faster Asserts use ASSERT
    #define HEAVYASSERT(condition) if(!(condition)) {std::cerr << "Error in file " << __FILE__ << " in function " << __FUNCTION__ << " at line " << __LINE__ << "!" << std::endl; abort(); } ((void)0)
    #else
    #define HEAVYASSERT(condition) ((void)0)
    #endif
}

#endif //HEIPROMAP_MACROS_H
