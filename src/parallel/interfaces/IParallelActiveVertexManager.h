#ifndef HEIDELBERGPROCESSMAPPING_IPARALLELACTIVEVERTEXMANAGER_H
#define HEIDELBERGPROCESSMAPPING_IPARALLELACTIVEVERTEXMANAGER_H

#include "IParallelGraph.h"
#include "../../interfaces/IActiveVertexManager.h"

namespace HeiProMap {

    template<typename TParallelGraph>
    class IParallelActiveVertexManager : public IActiveVertexManager {
        static_assert(std::is_base_of<IParallelGraph, TParallelGraph>::value, "TParallelGraph must inherit from IParallelGraph");
    public:
        // initialize
        virtual void initialize(TParallelGraph *t_p_g,
                                u64 n_threads) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IPARALLELACTIVEVERTEXMANAGER_H
