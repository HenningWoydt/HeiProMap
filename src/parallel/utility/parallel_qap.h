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

#ifndef HEIPROMAP_PARALLEL_QAP_H
#define HEIPROMAP_PARALLEL_QAP_H

#include "../parallel_definitions_1.h"
#include "../../definitions.h"

namespace HeiProMap {
    weight_t get_qap(p_graph_t &g,
                     p_av_manager_t &av_manager,
                     p_p_manager_t &p_manager,
                     p_d_oracle_t &d_oracle) {
        weight_t qap = 0;

        for (vertex_t u: av_manager) {
            ASSERT(av_manager.is_active(u));

            partition_t u_id = p_manager[u];

            for (size_t i = 0; i < g.size(u); ++i) {
                vertex_t    v    = g.neighbor(u, i);
                weight_t    ew   = g.get_weight(u, i);
                partition_t v_id = p_manager[v];
                weight_t    d    = d_oracle.get(u_id, v_id);
                qap += (d * ew);
            }
        }

        return qap;
    }
}

#endif //HEIPROMAP_PARALLEL_QAP_H
