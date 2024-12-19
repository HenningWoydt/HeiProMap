#ifndef HEIDELBERGPROCESSMAPPING_PARTITION_MANAGER_H
#define HEIDELBERGPROCESSMAPPING_PARTITION_MANAGER_H

#include "../../definitions.h"
#include "../interfaces/ISerialPartitionManager.h"

namespace HeiProMap {
    class PartitionManager final : public ISerialPartitionManager {
        vertex_t m_n    = 0;
        partition_t m_k = 0;

        std::vector<partition_t> partition;
        std::vector<weight_t> bweights;

    public:
        void initialize(const vertex_t t_n,
                        const partition_t t_k) override {
            m_n = t_n;
            m_k = t_k;

            partition.resize(m_n);
            bweights.resize(m_k, 0);
        }

        // read
        const partition_t& operator[](const vertex_t u) const override { return partition[u]; }

        // write
        void set(const vertex_t u, const weight_t w, const partition_t id) override {
            bweights[id] += w;
            partition[u] = id;
        }

        void move(const vertex_t u, const weight_t w, const partition_t old_id, const partition_t new_id) override {
            bweights[old_id] -= w;
            bweights[new_id] += w;
            partition[u] = new_id;
        }

        weight_t get_bweight(const partition_t id) const override { return bweights[id]; }
        std::vector<weight_t> get_bweights() const override { return bweights; }
        void uncontract(const vertex_t u, const vertex_t v) override { partition[v] = partition[u]; }

        void uncontract(const std::vector<EdgeUV>& matches) override {
            for (const auto [u, v] : matches) {
                partition[v] = partition[u];
            }
        }
    };
}

#endif //HEIDELBERGPROCESSMAPPING_PARTITION_MANAGER_H
