#ifndef HEIPROMAP_GRID_DISTANCE_ORACLE_H
#define HEIPROMAP_GRID_DISTANCE_ORACLE_H

#include <vector>
#include <cmath>

#include "../utility/aligned_array.h"
#include "../definitions.h"
#include "../utility/macros.h"
#include "../utility/utils.h"
#include "../utility/profiler.h"
#include "csr_graph.h"

namespace HeiProMap {
    class GridDistanceOracle {
        partition_t m_k = 0;
        const CSRGraph* m_graph = nullptr;

    public:
        ~GridDistanceOracle() {}

        void initialize(const CSRGraph& topology_graph, size_t max_threads, bool use_grid = false) {
            HEIPROMAP_PROFILE_SCOPE("misc", "GridDistanceOracle", "initialize");
            m_graph = &topology_graph;
            m_k = topology_graph.n;
        }

        weight_t get(partition_t u_id, partition_t v_id) const {
            return get(u_id, v_id, 0);
        }

        weight_t get(partition_t u_id, partition_t v_id, size_t thread_id) const {
            if (u_id == v_id) return 0;
            
            if (!m_graph->coords.empty()) {
                double dx = std::abs(m_graph->coords[u_id].first - m_graph->coords[v_id].first);
                double dy = std::abs(m_graph->coords[u_id].second - m_graph->coords[v_id].second);
                return static_cast<weight_t>(std::round(dx + dy));
            }
            
            return 1; // Fallback
        }

        partition_t get_h(partition_t u_id, partition_t v_id) const {
            return 0; 
        }

        partition_t get_k() const { return m_k; }

        bool last_level_pair(partition_t u_id, partition_t v_id) const {
            return false;
        }
    };
}

#endif //HEIPROMAP_GRID_DISTANCE_ORACLE_H
