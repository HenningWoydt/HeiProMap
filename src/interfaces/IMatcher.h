#ifndef HEIDELBERGPROCESSMAPPING_IMATCHER_H
#define HEIDELBERGPROCESSMAPPING_IMATCHER_H

#include <vector>

#include "../definitions.h"

namespace HeiProMap {

    class IMatcher {
    public:
        virtual ~IMatcher() = default;
        // matching
        virtual void match(std::vector<EdgeUV> &t_matches) = 0;
    };

}

#endif //HEIDELBERGPROCESSMAPPING_IMATCHER_H
