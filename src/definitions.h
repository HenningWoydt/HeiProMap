#ifndef MT_RECPROMAP_DEFINITIONS_H
#define MT_RECPROMAP_DEFINITIONS_H

#include <vector>
#include <cstdint>
#include <cstdlib>
#include <algorithm>
#include <limits>

namespace HeiProMap {
    typedef int8_t s8;
    typedef int16_t s16;
    typedef int32_t s32;
    typedef int64_t s64;

    typedef u_int8_t u8;
    typedef u_int16_t u16;
    typedef u_int32_t u32;
    typedef u_int64_t u64;

    typedef float f32;
    typedef double f64;

    typedef u32 vertex_t;
    typedef s32 weight_t;
    typedef u32 partition_t;

    const size_t HEAP_TOMBSTONE = std::numeric_limits<size_t>::max();

    /**
     * Struct that holds a vertex and a weight.
     */
    struct EdgeVW {
        vertex_t v;
        weight_t w;

        EdgeVW() = default;

        EdgeVW(vertex_t v, weight_t w) {
            this->v = v;
            this->w = w;
        }

        bool operator<(const EdgeVW &y) const {
            return v < y.v;
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

        EdgeUV(vertex_t u, vertex_t v) {
            this->u = u;
            this->v = v;
        }

        bool operator==(const EdgeUV &e) const {
            return (u == e.u && v == e.v);
        }
    };

    class EdgeUVW {
    public:
        vertex_t u;
        vertex_t v;
        f64 w;

    public:
        EdgeUVW() = default;

        EdgeUVW(vertex_t u, vertex_t v, vertex_t w) : u(u), v(v), w(w) {}

        bool operator<(const EdgeUVW &e) const {
            return w < e.w;
        }

        bool operator>(const EdgeUVW &e) const {
            return w < e.w;
        }
    };

    class MovePQ {
    public:
        vertex_t p_id;
        s64 qap_delta;

    public:
        MovePQ() = default;

        MovePQ(vertex_t p_id, s64 qap_delta) : p_id(p_id), qap_delta(qap_delta) {}

        bool operator>(const MovePQ &m) const {
            return qap_delta > m.qap_delta;
        }

        bool operator<(const MovePQ &m) const {
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

        Move(vertex_t t_u, vertex_t p_id, s64 qap_delta) : u(t_u), p_id(p_id), qap_delta(qap_delta) {}

        bool operator>(const Move &m) const {
            return qap_delta > m.qap_delta;
        }

        bool operator<(const Move &m) const {
            return qap_delta < m.qap_delta;
        }
    };

    class Swap {
    public:
        vertex_t u;
        s64 qap_delta;

    public:
        Swap() = default;

        Swap(vertex_t t_u, s64 qap_delta) : u(t_u), qap_delta(qap_delta) {}

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

    // partitioning algorithms
    enum PartitioningAlgorithms {
        GREEDY,
        KAFFPA_STRONG,
        KAFFPA_ECO,
        KAFFPA_FAST
    };

}

#endif //MT_RECPROMAP_DEFINITIONS_H
