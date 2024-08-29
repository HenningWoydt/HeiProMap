#ifndef HEIDELBERGPROCESSMAPPING_IPARALLELMATCHER_H
#define HEIDELBERGPROCESSMAPPING_IPARALLELMATCHER_H

#include "../../interfaces/IMatcher.h"
#include "IParallelGraph.h"
#include "IParallelActiveVertexManager.h"

namespace HeiProMap {

    class IParallelMatcher : IMatcher {
    public:
        // initialize
        virtual void initialize(IParallelGraph *t_p_g,
                                IParallelActiveVertexManager *t_p_av_manager,
                                u64 n_threads) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IPARALLELMATCHER_H
