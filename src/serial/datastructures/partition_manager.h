#ifndef SERIALPROCESSMAPPING_PARTITION_MANAGER_H
#define SERIALPROCESSMAPPING_PARTITION_MANAGER_H

#include "../../definitions.h"
#include "../../macros.h"
#include "../utility/utils.h"
#include "../utility/qap.h"
#include "../interfaces/ISerialActiveVertexManager.h"
#include "../interfaces/ISerialPartitionManager.h"

namespace HeiProMap {

    template<typename TSerialGraph, typename TSerialActiveVertexManager>
    class PartitionManager : public ISerialPartitionManager<TSerialGraph, TSerialActiveVertexManager> {
    private:
        TSerialGraph *m_p_g = nullptr;
        TSerialActiveVertexManager *m_p_av_manager = nullptr;
        partition_t m_k = 0;

        // actual partition
        std::vector<partition_t> partition;

        // partition weights
        std::vector<weight_t> bweights;

    public:
        void initialize(TSerialGraph *t_p_g,
                        TSerialActiveVertexManager *t_p_av_manager,
                        partition_t t_k) final {
            ASSERT(t_p_g != nullptr);
            ASSERT(t_p_av_manager != nullptr);

            m_p_g = t_p_g;
            m_p_av_manager = t_p_av_manager;
            m_k = t_k;

            // actual partition
            partition.resize(m_p_g->get_n());

            // partition weights
            bweights.resize(m_k, 0);
        }

        // read
        const partition_t &operator[](vertex_t u) const final { return partition[u]; }

        // write
        void set(vertex_t u, partition_t id) final {
            ASSERT(m_p_g != nullptr);
            bweights[id] += m_p_g->get_weight(u);
            partition[u] = id;
        }

        void move(vertex_t u, partition_t old_id, partition_t new_id) final {
            ASSERT(m_p_g != nullptr);
            bweights[old_id] -= m_p_g->get_weight(u);
            bweights[new_id] += m_p_g->get_weight(u);
            partition[u] = new_id;
        }

        // weights
        weight_t get_bweight(partition_t id) const final { return bweights[id]; }

        std::vector<weight_t> get_bweights() const final { return bweights; }

        // uncoarsing
        void uncontract(vertex_t u, vertex_t v) final { partition[v] = partition[u]; }
    };

}

#endif //SERIALPROCESSMAPPING_PARTITION_MANAGER_H
