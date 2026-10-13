#ifndef HEIPROMAP_GRID_DISTANCE_ORACLE_H
#define HEIPROMAP_GRID_DISTANCE_ORACLE_H

#include <vector>
#include <cmath>

#include "../utility/aligned_array.h"
#include "../utility/mapping.h"
#include "../definitions.h"
#include "../utility/macros.h"
#include "../utility/utils.h"
#include "../utility/profiler.h"
#include "../datastructures/csr_graph.h"

namespace HeiProMap {
    class GridDistanceOracle {
        struct Level {
            CSRGraph topology;
            size_t W = 1;
        };

        std::vector<Level> m_levels;
        std::vector<Mapping> m_mappings; // m_mappings[i]: level[i] (fine) -> level[i+1] (coarse)

        static size_t infer_W(const CSRGraph &g) {
            for (size_t i = 1; i < g.n; ++i) {
                if (g.neighborhoods[i + 1] - g.neighborhoods[i] == 2) {
                    return i + 1;
                }
            }
            return 1;
        }

    public:
        void initialize(CSRGraph topology, size_t /*max_threads*/) {
            HEIPROMAP_PROFILE_SCOPE("misc", "GridDistanceOracle", "initialize");
            m_levels.clear();
            m_mappings.clear();

            Level lvl;
            lvl.topology = std::move(topology);
            lvl.topology.infer_grid_coords();
            lvl.W = infer_W(lvl.topology);
            m_levels.push_back(std::move(lvl));
        }

        // Current number of partitions (= topology nodes at current level).
        partition_t current_k() const { return m_levels.back().topology.n; }

        // k at the next finer level (valid only when depth() > 0).
        partition_t fine_k() const { return m_levels[m_levels.size() - 2].topology.n; }

        const CSRGraph &current_topology() const { return m_levels.back().topology; }

        // Number of coarsenings applied so far.
        size_t depth() const { return m_levels.size() - 1; }

        bool can_coarsen() const { return m_levels.back().topology.n > 64; }

        // Push one coarser level onto the stack.
        void coarsen() {
            HEIPROMAP_PROFILE_SCOPE("misc", "GridDistanceOracle", "coarsen");
            const Level &cur = m_levels.back();
            size_t W_old = cur.W;
            size_t H_old = (cur.topology.n + W_old - 1) / W_old;
            size_t W_new = (W_old + 1) / 2;
            size_t H_new = (H_old + 1) / 2;
            size_t coarse_n = W_new * H_new;

            Mapping mapping;
            mapping.initialize(cur.topology.n);
            mapping.set_coarse_n(coarse_n);
            for (vertex_t u = 0; u < cur.topology.n; ++u) {
                vertex_t x_new = (u % W_old) / 2;
                vertex_t y_new = (u / W_old) / 2;
                mapping.set(u, y_new * W_new + x_new);
            }

            Level coarse_lvl;
            coarse_lvl.topology.contract_topology(cur.topology, mapping, CSRGraph::ContractionMode::MAX);
            coarse_lvl.W = W_new;

            m_mappings.push_back(std::move(mapping));
            m_levels.push_back(std::move(coarse_lvl));
        }

        // The mapping from the current finest level to the next coarser level.
        // Valid only when depth() > 0. Use before calling uncoarsen().
        const Mapping &top_mapping() const { return m_mappings.back(); }

        // Pop the current coarsest level, returning to the previous finer one.
        void uncoarsen() {
            HEIPROMAP_PROFILE_SCOPE("misc", "GridDistanceOracle", "uncoarsen");
            m_levels.pop_back();
            m_mappings.pop_back();
        }

        weight_t get(partition_t u_id, partition_t v_id) const {
            return get(u_id, v_id, 0);
        }

        weight_t get(partition_t u_id, partition_t v_id, size_t /*thread_id*/) const {
            if (u_id == v_id) return 0;
            size_t W = m_levels.back().W;
            if (W > 0) {
                partition_t dx = (u_id % W > v_id % W) ? (u_id % W - v_id % W) : (v_id % W - u_id % W);
                partition_t dy = (u_id / W > v_id / W) ? (u_id / W - v_id / W) : (v_id / W - u_id / W);
                return static_cast<weight_t>(dx + dy);
            }
            return 1;
        }

        partition_t get_h(partition_t /*u_id*/, partition_t /*v_id*/) const { return 0; }
        partition_t get_k() const { return current_k(); }
        bool last_level_pair(partition_t /*u_id*/, partition_t /*v_id*/) const { return false; }
    };
}

#endif //HEIPROMAP_GRID_DISTANCE_ORACLE_H
