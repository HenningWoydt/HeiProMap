//
// Created by henning on 9/29/26.
//

#ifndef HEIPROMAP_ACTIVE_BLOCK_SCHEDULING_H
#define HEIPROMAP_ACTIVE_BLOCK_SCHEDULING_H

#include "../utility/aligned_array.h"

namespace HeiProMap {
    class ActiveBlockScheduling {
    public:
        partition_t m_k = 0;
        // active block scheduling
        AlignedArray<u8> active_this_round;
        AlignedArray<u8> active_next_round;

        void initialize(const partition_t t_k) {
            m_k = t_k;
            active_this_round.initialize(t_k);
            active_next_round.initialize(t_k);
        }

        void reset(const partition_t t_k) {
            active_this_round.initialize(t_k, 1);
            active_next_round.initialize(t_k, 0);
        }

        void next_round() {
            std::swap(active_this_round, active_next_round);
            active_next_round.initialize(m_k, 0);
        }

        void activate(partition_t id1, partition_t id2) {
            active_next_round[id1] = 1;
            active_this_round[id2] = 1;
        }
    };
}

#endif //HEIPROMAP_ACTIVE_BLOCK_SCHEDULING_H
