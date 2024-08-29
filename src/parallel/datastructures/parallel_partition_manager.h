#ifndef HEIDELBERGPROCESSMAPPING_PARALLEL_PARTITION_MANAGER_H
#define HEIDELBERGPROCESSMAPPING_PARALLEL_PARTITION_MANAGER_H

#include <atomic>

#include "../../definitions.h"
#include "../../macros.h"
#include "../../serial/utility/utils.h"
#include "../../serial/utility/qap.h"
#include "../../interfaces/IPartitionManager.h"
#include "../interfaces/IParallelActiveVertexManager.h"
#include "../interfaces/IParallelPartitionManager.h"

namespace HeiProMap {

    class ParallelPartitionManager : public IParallelPartitionManager {
    private:
        IParallelGraph *m_p_g = nullptr;
        IParallelActiveVertexManager *m_p_av_manager = nullptr;
        partition_t m_k = 0;

        // actual partition
        std::vector<partition_t> partition;

        // partition weights
        std::atomic<weight_t> *m_p_bweights = nullptr;

    public:

        ~ParallelPartitionManager() { free(m_p_bweights); }

        void initialize(IParallelGraph *t_p_g,
                        IParallelActiveVertexManager *t_p_av_manager,
                        partition_t t_k,
                        u64 n_threads) final {
            ASSERT(t_p_g != nullptr);
            ASSERT(t_p_av_manager != nullptr);

            m_p_g = t_p_g;
            m_p_av_manager = t_p_av_manager;
            m_k = t_k;

            // actual partition
            partition.resize(m_p_g->get_n());

            // partition weights
            m_p_bweights = (std::atomic<weight_t> *) malloc(m_k * sizeof(std::atomic<weight_t>));
            for(size_t i = 0; i < m_k; ++i){ m_p_bweights[i] = 0; }
        }

        // read
        const partition_t &operator[](vertex_t u) const final { return partition[u]; }

        // write
        void set(vertex_t u, partition_t id) final {
            ASSERT(m_p_g != nullptr);
            m_p_bweights[id] += m_p_g->get_weight(u);
            partition[u] = id;
        }

        void move(vertex_t u, partition_t old_id, partition_t new_id) final {
            ASSERT(m_p_g != nullptr);
            m_p_bweights[old_id] -= m_p_g->get_weight(u);
            m_p_bweights[new_id] += m_p_g->get_weight(u);
            partition[u] = new_id;
        }

        // weights
        weight_t get_bweight(partition_t id) const final { return m_p_bweights[id]; }

        std::vector<weight_t> get_bweights() const final {
            std::vector<weight_t> temp(m_k);
            for (size_t i = 0; i < m_k; ++i) { temp[i] = m_p_bweights[i]; }
            return temp;
        }

        // uncoarsing
        void uncontract(vertex_t u, vertex_t v) final { partition[v] = partition[u]; }
    };

}

#endif //HEIDELBERGPROCESSMAPPING_PARALLEL_PARTITION_MANAGER_H
