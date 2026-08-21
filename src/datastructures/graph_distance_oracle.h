/*******************************************************************************
 * MIT License
 *
 * This file is part of HeiProMap.
 *
 * Copyright (C) 2025 Henning Woydt <henning.woydt@informatik.uni-heidelberg.de>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#ifndef HEIPROMAP_GRAPH_DISTANCE_ORACLE_H
#define HEIPROMAP_GRAPH_DISTANCE_ORACLE_H

#include <queue>
#include <vector>
#include <limits>

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

        struct CacheEntry {
            partition_t u;
            partition_t v;
            weight_t dist;
        };
        struct CacheSet {
            CacheEntry entries[2];
            u8 lru;
        };
        static constexpr size_t NUM_CACHED_ROWS = 1024*32;

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

        void initialize(const CSRGraph& topology_graph, size_t max_threads) {
            HEIPROMAP_PROFILE_SCOPE("misc", "GraphDistanceOracle", "initialize");
            print_stats();
            m_graph = &topology_graph;
            m_k = topology_graph.n;
            
            m_thread_states.resize(max_threads);
            for (auto& state : m_thread_states) {
                size_t num_rows = std::min<size_t>(1024*32, m_k);
                state.rows.resize(num_rows);
                for (size_t i = 0; i < num_rows; ++i) {
                    state.rows[i].u = -1;
                    state.rows[i].dist.resize(m_k);
                }
                state.row_locator.assign(m_k, -1);
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
            uint32_t min_tick = state.rows[0].last_used;
            for (size_t i = 1; i < state.rows.size(); ++i) {
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
