#ifndef HEIDELBERGPROCESSMAPPING_IPARALLELMATCHER_H
#define HEIDELBERGPROCESSMAPPING_IPARALLELMATCHER_H

#include "../../interfaces/IMatcher.h"
#include "IParallelGraph.h"
#include "IParallelActiveVertexManager.h"

namespace HeiProMap {

    template<typename TParallelGraph, typename TParallelActiveVertexManager>
    class IParallelMatcher : IMatcher {
        static_assert(std::is_base_of<IParallelGraph, TParallelGraph>::value, "TParallelGraph must inherit from IParallelGraph");
        static_assert(std::is_base_of<IParallelActiveVertexManager<TParallelGraph>, TParallelActiveVertexManager>::value, "TParallelActiveVertexManager must inherit from IParallelActiveVertexManager");
    public:
        // initialize
        virtual void initialize(TParallelGraph *t_p_g,
                                TParallelActiveVertexManager *t_p_av_manager,
                                u64 n_threads) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IPARALLELMATCHER_H
