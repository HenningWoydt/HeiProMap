#ifndef HEIDELBERGPROCESSMAPPING_ISERIALACTIVEVERTEXMANAGER_H
#define HEIDELBERGPROCESSMAPPING_ISERIALACTIVEVERTEXMANAGER_H

#include "../../interfaces/IActiveVertexManager.h"
#include "ISerialGraph.h"

namespace HeiProMap {

    class ISerialActiveVertexManager : public IActiveVertexManager {
    public:
        // initialize
        virtual void initialize(ISerialGraph *t_p_g) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALACTIVEVERTEXMANAGER_H
