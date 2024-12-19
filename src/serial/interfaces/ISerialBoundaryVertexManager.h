#ifndef HEIDELBERGPROCESSMAPPING_ISERIALBOUNDARYVERTEXMANAGER_H
#define HEIDELBERGPROCESSMAPPING_ISERIALBOUNDARYVERTEXMANAGER_H

#include "ISerialGraph.h"
#include "ISerialPartitionManager.h"

namespace HeiProMap {
    template <typename TSerialGraph, typename TSerialActiveVertexManager, typename TSerialPartitionManager>
    class ISerialBoundaryVertexManager {
        static_assert(std::is_base_of_v<ISerialGraph, TSerialGraph>, "TSerialGraph must inherit from ISerialGraph");
        static_assert(std::is_base_of_v<ISerialActiveVertexManager, TSerialActiveVertexManager>, "TSerialActiveVertexManager must inherit from ISerialActiveVertexManager");
        static_assert(std::is_base_of_v<ISerialPartitionManager, TSerialPartitionManager>, "TSerialPartitionManager must inherit from ISerialPartitionManager");

    public:
        virtual ~ISerialBoundaryVertexManager() = default;
        virtual void initialize(vertex_t t_n, partition_t t_k) = 0;
        virtual void uncontract(std::vector<EdgeUV>& matches, TSerialGraph& new_g, TSerialGraph& old_g, TSerialActiveVertexManager& av_manager, TSerialPartitionManager& p_manager) = 0;
        virtual void move(TSerialGraph &g, TSerialPartitionManager &p_manager, vertex_t u, partition_t old_id, partition_t new_id) = 0;
        virtual vertex_t get_n_boundary() const = 0;
        virtual bool is_boundary(vertex_t u) const = 0;
        virtual void add(vertex_t u, partition_t id) = 0;
    };
}

#endif //HEIDELBERGPROCESSMAPPING_ISERIALBOUNDARYVERTEXMANAGER_H
