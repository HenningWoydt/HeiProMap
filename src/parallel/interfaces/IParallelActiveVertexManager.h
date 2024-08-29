#ifndef HEIDELBERGPROCESSMAPPING_IPARALLELACTIVEVERTEXMANAGER_H
#define HEIDELBERGPROCESSMAPPING_IPARALLELACTIVEVERTEXMANAGER_H

#include "IParallelGraph.h"
#include "../../interfaces/IActiveVertexManager.h"

namespace HeiProMap {

    class IParallelActiveVertexManager : public IActiveVertexManager {
    public:
        // initialize
        virtual void initialize(IParallelGraph *t_p_g,
                                u64 n_threads) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IPARALLELACTIVEVERTEXMANAGER_H
