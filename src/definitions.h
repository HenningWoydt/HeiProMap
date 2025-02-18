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

namespace HeiProMap {
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

    typedef u32 vertex_t;
    typedef s32 weight_t;
    typedef u32 partition_t;

    constexpr size_t HEAP_TOMBSTONE = std::numeric_limits<size_t>::max();

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

        bool operator<(const EdgeVW& y) const {
            return v < y.v;
        }

        bool operator==(const EdgeVW& y) const {
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

        bool operator==(const EdgeUV& e) const {
            return (u == e.u && v == e.v);
        }

        bool operator<(const EdgeUV& e) const {
            return u < e.u;
        }
    };

    class EdgeUVW {
    public:
        vertex_t u;
        vertex_t v;
        f64 w;

    public:
        EdgeUVW() = default;

        EdgeUVW(const vertex_t u, const vertex_t v, const f64 w) : u(u), v(v), w(w) {}

        bool operator<(const EdgeUVW& e) const {
            return w < e.w;
        }

        bool operator>(const EdgeUVW& e) const {
            return w > e.w;
        }
    };

    class EdgeF64VW {
    public:
        vertex_t v;
        f64 w;

    public:
        EdgeF64VW() = default;

        EdgeF64VW(const vertex_t v, const f64 w) : v(v), w(w) {}

        bool operator<(const EdgeF64VW& e) const {
            return w < e.w;
        }

        bool operator>(const EdgeF64VW& e) const {
            return w < e.w;
        }
    };

    class MovePQ {
    public:
        vertex_t p_id;
        s64 qap_delta;

    public:
        MovePQ() = default;

        MovePQ(const vertex_t p_id, const s64 qap_delta) : p_id(p_id), qap_delta(qap_delta) {}

        bool operator>(const MovePQ& m) const {
            return qap_delta > m.qap_delta;
        }

        bool operator<(const MovePQ& m) const {
            return qap_delta < m.qap_delta;
        }
    };

    class Move {
    public:
        vertex_t u;
        vertex_t p_id;
        s64 qap_delta;

    public:
        Move() = default;

        Move(const vertex_t t_u, const vertex_t p_id, const s64 qap_delta) : u(t_u), p_id(p_id), qap_delta(qap_delta) {}

        bool operator>(const Move& m) const {
            return qap_delta > m.qap_delta;
        }

        bool operator<(const Move& m) const {
            return qap_delta < m.qap_delta;
        }
    };

    class Swap {
    public:
        vertex_t u;
        s64 qap_delta;

    public:
        Swap() = default;

        Swap(const vertex_t t_u, const s64 qap_delta) : u(t_u), qap_delta(qap_delta) {}

        bool operator>(const Swap& m) const {
            return qap_delta > m.qap_delta;
        }

        bool operator<(const Swap& m) const {
            return qap_delta < m.qap_delta;
        }

        bool operator<=(const Swap& m) const {
            return qap_delta <= m.qap_delta;
        }
    };

    // partitioning algorithms
    enum PartitioningAlgorithms {
        GREEDY,
        KAFFPA_STRONG,
        KAFFPA_ECO,
        KAFFPA_FAST
    };
}

#endif //HEIPROMAP_DEFINITIONS_H
