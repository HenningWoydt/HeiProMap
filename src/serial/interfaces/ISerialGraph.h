#ifndef HEIDELBERGPROCESSMAPPING_ISERIALGRAPH_H
#define HEIDELBERGPROCESSMAPPING_ISERIALGRAPH_H

namespace HeiProMap {
    class ISerialGraph {
    public:
        virtual ~ISerialGraph() = default;
        virtual vertex_t get_n() const = 0;
        virtual vertex_t get_m() const = 0;
        virtual weight_t get_weight() const = 0;
        virtual weight_t get_weight(vertex_t u) const = 0;
        virtual size_t size(vertex_t u) const = 0;
        virtual vertex_t neighbor(vertex_t u, size_t idx) const = 0;
        virtual weight_t get_weight(vertex_t u, size_t idx) const = 0;
        virtual bool edge_exists(vertex_t u, vertex_t v) const = 0;
    };
}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALGRAPH_H
