#include <cmath>
#ifndef HEIPROMAP_GRAPH_DISTANCE_ORACLE_H
#define HEIPROMAP_GRAPH_DISTANCE_ORACLE_H

#include <queue>
#include <vector>
#include <limits>
#include <iostream>
#include <iomanip>

#include "../utility/aligned_array.h"
#include "../definitions.h"
#include "../utility/macros.h"
#include "../utility/utils.h"
#include "../utility/profiler.h"
#include "csr_graph.h"

namespace HeiProMap {
    class GraphDistanceOracle {
        partition_t m_k = 0;
        const CSRGraph* m_graph = nullptr;
        bool m_use_grid = false;

        static constexpr size_t NUM_CACHED_ROWS = 1024*128;

        template <class T, class Container = std::vector<T>, class Compare = std::less<typename Container::value_type>>
        class clearable_pq : public std::priority_queue<T, Container, Compare> {
        public:
            void clear() { this->c.clear(); }
        };

        struct CachedRow {
            partition_t u;
            uint32_t last_used;
            std::vector<weight_t> dist;
        };

        struct ThreadState {
            std::vector<CachedRow> rows;
            std::vector<int32_t> row_locator;
            uint32_t tick;
            clearable_pq<std::pair<weight_t, partition_t>,
                         std::vector<std::pair<weight_t, partition_t>>,
                         std::greater<>> pq;
                         
            uint64_t stat_queries = 0;
            uint64_t stat_hits = 0;
            uint64_t stat_misses = 0;
            uint64_t stat_dijkstra_visited = 0;
        };
        mutable std::vector<ThreadState> m_thread_states;

    public:
        ~GraphDistanceOracle() {
            print_stats();
        }

        void print_stats() const {
            uint64_t total_queries = 0;
            uint64_t total_hits = 0;
            uint64_t total_misses = 0;
            uint64_t total_visited = 0;
            for (const auto& state : m_thread_states) {
                total_queries += state.stat_queries;
                total_hits += state.stat_hits;
                total_misses += state.stat_misses;
                total_visited += state.stat_dijkstra_visited;
            }
            if (total_queries > 0) {
                std::cout << "[DistanceOracle] Stats for k=" << m_k 
                          << " | Queries: " << total_queries 
                          << " | Hits: " << total_hits 
                          << " (" << std::fixed << std::setprecision(2) << (100.0 * total_hits / total_queries) << "%)"
                          << " | Misses: " << total_misses 
                          << " | Dijkstra nodes visited: " << total_visited 
                          << std::endl;
            }
        }

        void initialize(const CSRGraph& topology_graph, size_t max_threads, bool use_grid = false) {
            HEIPROMAP_PROFILE_SCOPE("misc", "GraphDistanceOracle", "initialize");
            print_stats();
            m_graph = &topology_graph;
            m_k = topology_graph.n;
            m_use_grid = use_grid;
            
            m_thread_states.resize(max_threads);
            for (auto& state : m_thread_states) {
                if (!m_use_grid) {
                    // Limit total cache memory per thread to ~256MB to avoid OOM on huge graphs
                    size_t max_mem_per_thread = 256ULL * 1024 * 1024;
                    size_t row_size = m_k * sizeof(weight_t);
                    size_t max_rows_by_mem = std::max<size_t>(1, max_mem_per_thread / row_size);
                    
                    size_t num_rows = std::min<size_t>({NUM_CACHED_ROWS, m_k, max_rows_by_mem});
                    state.rows.resize(num_rows);
                    for (size_t i = 0; i < num_rows; ++i) {
                        state.rows[i].u = -1;
                        state.rows[i].last_used = 0;
                        state.rows[i].dist.resize(m_k);
                    }
                    state.row_locator.assign(m_k, -1);
                } else {
                    state.rows.clear();
                    state.row_locator.clear();
                }
                state.tick = 0;
                
                state.stat_queries = 0;
                state.stat_hits = 0;
                state.stat_misses = 0;
                state.stat_dijkstra_visited = 0;
            }
        }

        weight_t get(partition_t u_id, partition_t v_id) const {
            size_t thread_id = omp_in_parallel() ? omp_get_thread_num() : 0;
            return get(u_id, v_id, thread_id);
        }

        weight_t get(partition_t u_id, partition_t v_id, size_t thread_id) const {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);
            ASSERT(thread_id < m_thread_states.size());
            
            if (u_id == v_id) return 0;
            if (m_use_grid && !m_graph->coords.empty()) {
                double dx = std::abs(m_graph->coords[u_id].first - m_graph->coords[v_id].first);
                double dy = std::abs(m_graph->coords[u_id].second - m_graph->coords[v_id].second);
                return static_cast<weight_t>(std::round(dx + dy));
            }

            ThreadState& state = m_thread_states[thread_id];
            state.tick++;
            state.stat_queries++;

            int32_t idx_u = state.row_locator[u_id];
            if (idx_u != -1) {
                state.rows[idx_u].last_used = state.tick;
                state.stat_hits++;
                return state.rows[idx_u].dist[v_id];
            }
            int32_t idx_v = state.row_locator[v_id];
            if (idx_v != -1) {
                state.rows[idx_v].last_used = state.tick;
                state.stat_hits++;
                return state.rows[idx_v].dist[u_id];
            }
            
            state.stat_misses++;

            // Not found. Evict LRU row
            size_t lru_idx = 0;
            uint32_t min_tick = state.tick + 1; // Start higher than possible
            for (size_t i = 0; i < state.rows.size(); ++i) {
                if (state.rows[i].u == (partition_t)-1) {
                    lru_idx = i;
                    break;
                }
                if (state.rows[i].last_used < min_tick) {
                    min_tick = state.rows[i].last_used;
                    lru_idx = i;
                }
            }

            CachedRow& row = state.rows[lru_idx];
            if (row.u != (partition_t)-1) {
                state.row_locator[row.u] = -1; // clear old mapping
            }
            
            row.u = u_id;
            row.last_used = state.tick;
            state.row_locator[u_id] = lru_idx;
            
            // Fast reset
            std::fill(row.dist.begin(), row.dist.end(), std::numeric_limits<weight_t>::max() / 2);

            state.pq.push({0, u_id});
            row.dist[u_id] = 0;

            while (!state.pq.empty()) {
                auto [d, u] = state.pq.top();
                state.pq.pop();
                state.stat_dijkstra_visited++;

                if (d > row.dist[u]) continue;

                for (size_t i = m_graph->neighborhoods[u]; i < m_graph->neighborhoods[u + 1]; ++i) {
                    partition_t v = m_graph->edges_v[i];
                    weight_t w = m_graph->edges_w[i];

                    if (row.dist[u] + w < row.dist[v]) {
                        row.dist[v] = row.dist[u] + w;
                        state.pq.push({row.dist[v], v});
                    }
                }
            }

            state.pq.clear();
            return row.dist[v_id];
        }

        partition_t get_h(partition_t u_id, partition_t v_id) const {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);
            return 0; 
        }

        partition_t get_k() const { return m_k; }

        bool last_level_pair(partition_t u_id, partition_t v_id) const {
            return false;
        }
    };
}

#endif //HEIPROMAP_GRAPH_DISTANCE_ORACLE_H
