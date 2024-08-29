#ifndef HEIDELBERGPROCESSMAPPING_ISERIALMATCHER_H
#define HEIDELBERGPROCESSMAPPING_ISERIALMATCHER_H

#include "../../interfaces/IMatcher.h"
#include "ISerialGraph.h"
#include "ISerialActiveVertexManager.h"

namespace HeiProMap {

    class ISerialMatcher : IMatcher {
    public:
        // initialize
        virtual void initialize(ISerialGraph *t_p_g,
                                ISerialActiveVertexManager *t_p_av_manager) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALMATCHER_H
