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

        size_t m_W = 0;

        void initialize(const CSRGraph& topology_graph, size_t max_threads, bool use_grid = false, size_t grid_W = 0) {
            HEIPROMAP_PROFILE_SCOPE("misc", "GridDistanceOracle", "initialize");
            m_graph = &topology_graph;
            m_k = topology_graph.n;
            
            m_W = grid_W;
            if (m_W == 0 && m_k > 0) {
                for (size_t i = 1; i < m_k; ++i) {
                    if (topology_graph.neighborhoods[i+1] - topology_graph.neighborhoods[i] == 2) {
                        m_W = i + 1;
                        break;
                    }
                }
                if (m_W == 0) m_W = 1;
            }
        }

        weight_t get(partition_t u_id, partition_t v_id) const {
            return get(u_id, v_id, 0);
        }

        weight_t get(partition_t u_id, partition_t v_id, size_t thread_id) const {
            if (u_id == v_id) return 0;
            if (m_W > 0) {
                partition_t ux = u_id % m_W;
                partition_t uy = u_id / m_W;
                partition_t vx = v_id % m_W;
                partition_t vy = v_id / m_W;
                partition_t dx = (ux > vx) ? (ux - vx) : (vx - ux);
                partition_t dy = (uy > vy) ? (uy - vy) : (vy - uy);
                return static_cast<weight_t>(dx + dy);
            }
            return 1;
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
