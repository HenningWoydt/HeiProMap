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

namespace HeiProMap {
#ifndef ASSERT_ENABLED
#define ASSERT_ENABLED true
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

#ifdef __GNUC__  // (GCC or Clang)
#define ASSUME_ALIGNED(dtype, ptr, alignment) (dtype) __builtin_assume_aligned((ptr), (alignment))
#elif defined(_MSC_VER)  // MSVC
#include <intrin.h>
    // MSVC does not have a direct equivalent to __builtin_assume_aligned,
    // but we can at least use __assume to tell the compiler something
    // about the pointer. This is not exactly the same, though.
    inline void *ASSUME_ALIGNED(void *ptr, size_t /*alignment*/) {
        __assume(ptr != nullptr);
        return ptr;
    }
#else
    // Fallback: do nothing
#define ASSUME_ALIGNED(dtype, ptr, alignment) (ptr)
#endif


    // Macro to iterate over the neighborhood of vertex u of a graph
#define forall_guivw(g, u, i, v, w) for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) { const vertex_t v = g.edges_v[i]; const weight_t w = g.edges_w[i];
#define forall_guiv(g, u, i, v)  for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) { const vertex_t v = g.edges_v[i];
#define forall_gu(g, u) for (vertex_t u = 0; u < g.n; ++u) {
#define endfor }

    // Macro to iterate over all boundary vertices
#define forall_bv_id_iu(bv_manager, id, i, u) for (size_t i = 0; i < bv_manager.size(id); ++i) { const vertex_t u = bv_manager.get(id, i);

    // Macro to iterate over the connections of a vertex
#define forall_bc_ui_id(block_conn, u, i, id) for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) { const partition_t id = block_conn.get_id(i);
#define forall_bc_ui_id_idw(block_conn, u, i, id, id_w) for (size_t i = block_conn.start(u); i < block_conn.end(u); ++i) { const partition_t id = block_conn.get_id(i); const weight_t id_w = block_conn.get_w(i);
}

#endif //HEIPROMAP_MACROS_H
