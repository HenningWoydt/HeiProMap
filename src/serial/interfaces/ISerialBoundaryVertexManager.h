#ifndef HEIDELBERGPROCESSMAPPING_ISERIALBOUNDARYVERTEXMANAGER_H
#define HEIDELBERGPROCESSMAPPING_ISERIALBOUNDARYVERTEXMANAGER_H

#include "ISerialGraph.h"
#include "ISerialPartitionManager.h"

namespace HeiProMap {
    class ISerialBoundaryVertexManager {

    public:
        virtual ~ISerialBoundaryVertexManager() = default;
        virtual void initialize(vertex_t t_n, partition_t t_k) = 0;

        template <typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialPartitionManager>
        void uncontract(std::vector<EdgeUV>& matches, TSerialGraph& new_g, TSerialGraph& old_g, TSerialActiveVertexManager& av_manager, TSerialPartitionManager& p_manager) {}

        template <typename TSerialGraph, typename TSerialPartitionManager>
        void move(TSerialGraph& g, TSerialPartitionManager& p_manager, vertex_t u, partition_t old_id, partition_t new_id) {}

        virtual vertex_t get_n_boundary() const = 0;
        virtual bool is_boundary(vertex_t u) const = 0;
        virtual void add(vertex_t u, partition_t id) = 0;
    };
}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALBOUNDARYVERTEXMANAGER_H
