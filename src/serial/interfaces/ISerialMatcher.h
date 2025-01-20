#ifndef HEIDELBERGPROCESSMAPPING_ISERIALMATCHER_H
#define HEIDELBERGPROCESSMAPPING_ISERIALMATCHER_H

#include "ISerialActiveVertexManager.h"
#include "ISerialGraph.h"

namespace HeiProMap {
    class ISerialMatcher {
    public:
        virtual ~ISerialMatcher() = default;

        virtual void initialize(size_t n, weight_t l_max) = 0;

        template<typename TSerialGraph, typename TSerialActiveVertexManager>
        void match(const TSerialGraph &g,
                   TSerialActiveVertexManager &av_manager,
                   std::vector<EdgeUV> &matches) {}
    };
}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALMATCHER_H
