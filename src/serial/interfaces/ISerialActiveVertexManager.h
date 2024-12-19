#ifndef HEIDELBERGPROCESSMAPPING_ISERIALACTIVEVERTEXMANAGER_H
#define HEIDELBERGPROCESSMAPPING_ISERIALACTIVEVERTEXMANAGER_H

#include <vector>

namespace HeiProMap {
    class ISerialActiveVertexManager {
    public:
        virtual ~ISerialActiveVertexManager() = default;
        virtual void initialize(size_t n) = 0;
        virtual vertex_t get_n_active() const = 0;
        virtual void activate_vertex(vertex_t u) = 0;
        virtual void disable_vertex(vertex_t u) = 0;
        virtual bool is_active(vertex_t u) const = 0;
        virtual bool is_disabled(vertex_t u) const = 0;
        virtual bool get_state(vertex_t u) const = 0;
        virtual void contract(const std::vector<EdgeUV>& matches) = 0;
        virtual void uncontract(const std::vector<EdgeUV>& matches) = 0;
        //virtual void reset_iterator() = 0;
        //virtual vertex_t get() const = 0;
        //virtual void next() = 0;
        //virtual bool available() = 0;
    };
}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALACTIVEVERTEXMANAGER_H
