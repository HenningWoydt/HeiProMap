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

#ifndef HEIPROMAP_DEFINITIONS_H
#define HEIPROMAP_DEFINITIONS_H

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>

namespace HeiProMap {
    class CSRGraph;
    typedef CSRGraph graph_t;

#ifdef USE_GRAPH_DISTANCE_ORACLE
    class GraphDistanceOracle;
    typedef GraphDistanceOracle d_oracle_t;
#else
    class DistanceOracle;
    typedef DistanceOracle d_oracle_t;
#endif

    class PartitionManager;
    typedef PartitionManager p_manager_t;
    class BoundaryVertexManager;
    typedef BoundaryVertexManager bv_manager_t;
    
    class LargeQuotientGraph;
typedef LargeQuotientGraph q_graph_t;
    class BlockConn;
    typedef BlockConn block_conn_t;

    typedef int8_t s8;
    typedef int16_t s16;
    typedef int32_t s32;
    typedef int64_t s64;

    typedef uint8_t u8;
    typedef uint16_t u16;
    typedef uint32_t u32;
    typedef uint64_t u64;

    typedef float f32;
    typedef double f64;

    typedef u64 vertex_t;
    typedef s64 weight_t;
    typedef u64 partition_t;

    constexpr size_t HEAP_TOMBSTONE = std::numeric_limits<size_t>::max();
    constexpr partition_t NO_ID = std::numeric_limits<partition_t>::max();

    struct CommandLineOption {
        std::string large_key;
        std::string small_key;
        std::string description;
        std::string default_val;
        std::string input;
        bool is_set;
    };

    /**
     * Struct that holds a vertex and a weight.
     */
    struct EdgeVW {
        vertex_t v;
        weight_t w;

        EdgeVW() = default;

        EdgeVW(const vertex_t v, const weight_t w) {
            this->v = v;
            this->w = w;
        }

        bool operator<(const EdgeVW &y) const {
            return v < y.v;
        }

        bool operator==(const EdgeVW &y) const {
            return v == y.v && w == y.w;
        }
    };

    /**
     * Struct that holds a vertex and a weight.
     */
    struct QGraphUV {
        partition_t u;
        partition_t v;

        QGraphUV() = default;

        QGraphUV(partition_t u, partition_t v) {
            this->u = u;
            this->v = v;
        }
    };

    /**
     * Struct that holds two vertices.
     */
    struct EdgeUV {
        vertex_t u;
        vertex_t v;

        EdgeUV(const vertex_t u, const vertex_t v) {
            this->u = u;
            this->v = v;
        }

        bool operator==(const EdgeUV &e) const {
            return (u == e.u && v == e.v);
        }

        bool operator<(const EdgeUV &e) const {
            return u < e.u;
        }
    };

    enum struct EdgeRatingFunction {
        WEIGHT,
        EXPANSION,
        EXPANSIONSTAR,
        EXPANSIONSTARSTAR,
        INNEROUTER
    };

    class EdgeUVW {
    public:
        vertex_t u;
        vertex_t v;
        f32 w;
        u32 tiebreaker = 0;

    public:
        EdgeUVW() = default;

        EdgeUVW(const vertex_t u, const vertex_t v, const f32 w) : u(u), v(v), w(w), tiebreaker(0) {
        }
        
        EdgeUVW(const vertex_t u, const vertex_t v, const f32 w, const u32 tb) : u(u), v(v), w(w), tiebreaker(tb) {
        }

        bool operator<(const EdgeUVW &e) const {
            if (w != e.w) return w < e.w;
            return tiebreaker < e.tiebreaker;
        }

        bool operator>(const EdgeUVW &e) const {
            if (w != e.w) return w > e.w;
            return tiebreaker > e.tiebreaker;
        }

        bool operator==(const EdgeUVW &e) const {
            bool a = u == e.u && v == e.v && w == e.w;
            bool b = v == e.u && u == e.v && w == e.w;
            return a || b;
        }
    };

    class EdgeUVWeight {
    public:
        vertex_t u;
        vertex_t v;
        weight_t w;

    public:
        EdgeUVWeight() = default;

        EdgeUVWeight(const vertex_t u, const vertex_t v, const weight_t w) : u(u), v(v), w(w) {
        }

        bool operator<(const EdgeUVWeight &e) const {
            return w < e.w;
        }

        bool operator>(const EdgeUVWeight &e) const {
            return w > e.w;
        }

        bool operator==(const EdgeUVWeight &e) const {
            bool a = u == e.u && v == e.v && w == e.w;
            bool b = v == e.u && u == e.v && w == e.w;
            return a || b;
        }
    };

    class Swap {
    public:
        vertex_t u;
        weight_t qap_delta;

    public:
        Swap() = default;

        Swap(const vertex_t t_u, const weight_t qap_delta) : u(t_u), qap_delta(qap_delta) {
        }

        bool operator>(const Swap &m) const {
            return qap_delta > m.qap_delta;
        }

        bool operator<(const Swap &m) const {
            return qap_delta < m.qap_delta;
        }

        bool operator<=(const Swap &m) const {
            return qap_delta <= m.qap_delta;
        }
    };

    class KWayFMMove {
    public:
        vertex_t u;
        partition_t u_id;
        partition_t to_move_id;
        weight_t qap_delta;

        KWayFMMove() = default;

        KWayFMMove(const vertex_t t_u, const partition_t t_u_id, const partition_t t_to_move, const weight_t t_qap_delta) {
            u = t_u;
            u_id = t_u_id;
            to_move_id = t_to_move;
            qap_delta = t_qap_delta;
        }

        bool operator>(const KWayFMMove &m) const {
            return qap_delta > m.qap_delta;
        }

        bool operator>=(const KWayFMMove &m) const {
            return qap_delta >= m.qap_delta;
        }

        bool operator<(const KWayFMMove &m) const {
            return qap_delta < m.qap_delta;
        }

        bool operator<=(const KWayFMMove &m) const {
            return qap_delta <= m.qap_delta;
        }
    };

    class SmallKWayFMMove {
    public:
        partition_t to_move_id;
        weight_t qap_delta;

        SmallKWayFMMove() = default;

        SmallKWayFMMove(const partition_t t_to_move, const weight_t t_qap_delta) {
            to_move_id = t_to_move;
            qap_delta = t_qap_delta;
        }

        bool operator>(const SmallKWayFMMove &m) const {
            return qap_delta > m.qap_delta;
        }

        bool operator>=(const SmallKWayFMMove &m) const {
            return qap_delta >= m.qap_delta;
        }

        bool operator<(const SmallKWayFMMove &m) const {
            return qap_delta < m.qap_delta;
        }

        bool operator<=(const SmallKWayFMMove &m) const {
            return qap_delta <= m.qap_delta;
        }
    };

    class Move {
    public:
        vertex_t u;
        partition_t u_id;
        partition_t to_move_id;

        Move() = default;

        Move(const vertex_t t_u, const vertex_t t_u_id, const partition_t t_to_move) {
            u = t_u;
            u_id = t_u_id;
            to_move_id = t_to_move;
        }
    };
}

#endif //HEIPROMAP_DEFINITIONS_H
