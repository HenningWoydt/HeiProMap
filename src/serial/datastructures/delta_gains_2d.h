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

#ifndef HEIPROMAP_DELTA_GAINS_2D_H
#define HEIPROMAP_DELTA_GAINS_2D_H

#include <vector>

#include "../../definitions.h"
#include "../interfaces/ISerialBoundaryVertexManager.h"
#include "../interfaces/ISerialDistanceOracle.h"
#include "../interfaces/ISerialGraph.h"
#include "../interfaces/ISerialPartitionManager.h"

namespace HeiProMap {
    class DeltaGains2D {
    private:
        vertex_t n;
        partition_t k;
        partition_t id1, id2;

        std::vector<bool> active;

        std::vector<s64> gains;

    public:
        DeltaGains2D(vertex_t t_n, partition_t t_k) {
            n = t_n;
            k = t_k;
            active.resize(t_n);
            gains.resize(t_n * 2);
        }

        template <typename TSerialGraph, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle>
        void initialize(TSerialGraph& g,
                        TSerialBoundaryVertexManager& bv_manager,
                        TSerialPartitionManager& p_manager,
                        TSerialDistanceOracle& d_oracle,
                        partition_t u_id,
                        partition_t v_id) {
            static_assert(std::is_base_of_v<ISerialGraph, TSerialGraph>, "TSerialGraph must inherit from ISerialGraph");
            static_assert(std::is_base_of_v<ISerialBoundaryVertexManager, TSerialBoundaryVertexManager>, "TSerialBoundaryVertexManager must inherit from ISerialBoundaryVertexManager");
            static_assert(std::is_base_of_v<ISerialPartitionManager, TSerialPartitionManager>, "TSerialPartitionManager must inherit from ISerialPartitionManager");
            static_assert(std::is_base_of_v<ISerialDistanceOracle, TSerialDistanceOracle>, "TSerialDistanceOracle must inherit from ISerialDistanceOracle");

            // set ids
            id1 = u_id;
            id2 = v_id;

            // reset everything
            std::fill(active.begin(), active.end(), false);
            std::fill(gains.begin(), gains.end(), 0);

            // calculate the gains from bottom up, for all boundary vertices
            for (vertex_t u : bv_manager[u_id]) {
                gains[2 * u + 1] = get_u_qap_delta(g, u, u_id, v_id, p_manager, d_oracle);
                active[u]        = true;
            }
            for (vertex_t v : bv_manager[v_id]) {
                gains[2 * v] = get_u_qap_delta(g, v, v_id, u_id, p_manager, d_oracle);
                active[v]    = true;
            }
        }

        bool is_up_to_date(vertex_t u) {
            return active[u];
        }

        s64 get_gain(vertex_t u, partition_t new_u_id) {
            if (new_u_id == id1) {
                return gains[2 * u];
            }
            return gains[2 * u + 1];
        }

        void move(){}

        template <typename TSerialGraph, typename TSerialBoundaryVertexManager, typename TSerialPartitionManager, typename TSerialDistanceOracle>
        void update(vertex_t u,
                    partition_t u_id,
                    TSerialGraph& g,
                    TSerialBoundaryVertexManager& bv_manager,
                    TSerialPartitionManager& p_manager,
                    TSerialDistanceOracle& d_oracle) {

        }
    };
}

#endif //HEIPROMAP_DELTA_GAINS_2D_H
