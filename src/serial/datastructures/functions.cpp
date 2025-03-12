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

#include "functions.h"

namespace HeiProMap {
    bool is_boundary(const graph_t& g,
                     const p_manager_t& p_manager,
                     const vertex_t u) {
        partition_t u_id = p_manager[u];

#pragma GCC unroll 4
        forall_guiv(g, u, i, v)
            {
                partition_t v_id = p_manager[v];
                if (u_id != v_id) {
                    return true;
                }
            }
        endfor
        return false;
    }

    bool is_connected_to(const graph_t& g,
                         const p_manager_t& p_manager,
                         const vertex_t u,
                         const partition_t id) {

#pragma GCC unroll 4
        forall_guiv(g, u, i, v)
            {
                partition_t v_id = p_manager[v];
                if (id == v_id) {
                    return true;
                }
            }
        endfor
        return false;
    }
}
