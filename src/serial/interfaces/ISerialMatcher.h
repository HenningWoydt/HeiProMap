#ifndef HEIDELBERGPROCESSMAPPING_ISERIALMATCHER_H
#define HEIDELBERGPROCESSMAPPING_ISERIALMATCHER_H

#include "ISerialActiveVertexManager.h"
#include "ISerialGraph.h"

namespace HeiProMap {
    // template <typename TSerialGraph, typename TSerialActiveVertexManager>
    class ISerialMatcher {
        // static_assert(std::is_base_of_v<ISerialGraph, TSerialGraph>, "TSerialGraph must inherit from ISerialGraph");
        // static_assert(std::is_base_of_v<ISerialActiveVertexManager, TSerialActiveVertexManager>, "TSerialActiveVertexManager must inherit from ISerialActiveVertexManager");

    public:
        virtual ~ISerialMatcher() = default;
        virtual void initialize(size_t n, weight_t l_max) = 0;

        template <typename TSerialGraph, typename TSerialActiveVertexManager>
        void match(const TSerialGraph& g,
                   TSerialActiveVertexManager& av_manager,
                   std::vector<EdgeUV>& matches) {}
    };
}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALMATCHER_H
