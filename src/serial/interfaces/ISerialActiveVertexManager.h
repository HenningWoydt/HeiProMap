#ifndef HEIDELBERGPROCESSMAPPING_ISERIALACTIVEVERTEXMANAGER_H
#define HEIDELBERGPROCESSMAPPING_ISERIALACTIVEVERTEXMANAGER_H

#include "../../interfaces/IActiveVertexManager.h"
#include "ISerialGraph.h"

namespace HeiProMap {

    template<typename TSerialGraph>
    class ISerialActiveVertexManager : public IActiveVertexManager {
        static_assert(std::is_base_of<ISerialGraph, TSerialGraph>::value, "TSerialGraph must inherit from ISerialGraph");
    public:
        // initialize
        virtual void initialize(TSerialGraph *t_p_g) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALACTIVEVERTEXMANAGER_H
