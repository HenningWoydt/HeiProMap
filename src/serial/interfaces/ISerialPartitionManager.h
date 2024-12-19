#ifndef HEIDELBERGPROCESSMAPPING_ISERIALPARTITIONMANAGER_H
#define HEIDELBERGPROCESSMAPPING_ISERIALPARTITIONMANAGER_H

#include "../../definitions.h"

namespace HeiProMap {
    class ISerialPartitionManager {
    public:
        virtual ~ISerialPartitionManager() = default;
        virtual void initialize(vertex_t n, partition_t k) = 0;
        virtual const partition_t& operator[](vertex_t u) const = 0;
        virtual void set(vertex_t u, weight_t w, partition_t id) = 0;
        virtual void move(vertex_t u, weight_t w, partition_t old_id, partition_t new_id) = 0;
        virtual weight_t get_bweight(partition_t id) const = 0;
        virtual std::vector<weight_t> get_bweights() const = 0;
        virtual void uncontract(vertex_t u, vertex_t v) = 0;
        virtual void uncontract(const std::vector<EdgeUV>& matches) = 0;
    };
}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALPARTITIONMANAGER_H
