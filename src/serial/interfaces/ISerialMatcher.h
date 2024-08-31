#ifndef HEIDELBERGPROCESSMAPPING_ISERIALMATCHER_H
#define HEIDELBERGPROCESSMAPPING_ISERIALMATCHER_H

#include "../../interfaces/IMatcher.h"
#include "ISerialGraph.h"
#include "ISerialActiveVertexManager.h"

namespace HeiProMap {

    template<typename TSerialGraph, typename TSerialActiveVertexManager>
    class ISerialMatcher : IMatcher {
        static_assert(std::is_base_of<ISerialGraph, TSerialGraph>::value, "TSerialGraph must inherit from ISerialGraph");
        static_assert(std::is_base_of<ISerialActiveVertexManager<TSerialGraph>, TSerialActiveVertexManager>::value, "TSerialActiveVertexManager must inherit from ISerialActiveVertexManager");
    public:
        // initialize
        virtual void initialize(TSerialGraph *t_p_g,
                                TSerialActiveVertexManager *t_p_av_manager) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALMATCHER_H
