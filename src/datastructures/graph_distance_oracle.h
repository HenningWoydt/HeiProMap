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
        // 524288 sets * 2 ways = 1,048,576 total entries
        static constexpr size_t NUM_SETS = 524288; // Power of 2

        template <class T, class Container = std::vector<T>, class Compare = std::less<typename Container::value_type>>
        class clearable_pq : public std::priority_queue<T, Container, Compare> {
        public:
            void clear() { this->c.clear(); }
        };

        struct ThreadState {
            std::vector<weight_t> dist;
            std::vector<partition_t> touched;
            clearable_pq<std::pair<weight_t, partition_t>,
                         std::vector<std::pair<weight_t, partition_t>>,
                         std::greater<>> pq;
            std::vector<CacheSet> cache;
        };
        mutable std::vector<ThreadState> m_thread_states;

    public:
        void initialize(const CSRGraph& topology_graph, size_t max_threads) {
            HEIPROMAP_PROFILE_SCOPE("misc", "GraphDistanceOracle", "initialize");
            m_graph = &topology_graph;
            m_k = topology_graph.n;
            
            m_thread_states.resize(max_threads);
            for (auto& state : m_thread_states) {
                state.dist.resize(m_k, std::numeric_limits<weight_t>::max() / 2);
                CacheSet empty_set = { {{ (partition_t)-1, (partition_t)-1, 0 }, { (partition_t)-1, (partition_t)-1, 0 }}, 0 };
                state.cache.resize(NUM_SETS, empty_set);
            }
        }

        weight_t get(partition_t u_id, partition_t v_id) const {
            size_t thread_id = omp_in_parallel() ? omp_get_thread_num() : 0;
            return get(u_id, v_id, thread_id);
        }

        inline uint64_t mix_hash(uint64_t key) const {
            key ^= key >> 33;
            key *= 0xff51afd7ed558ccdULL;
            key ^= key >> 33;
            key *= 0xc4ceb9fe1a85ec53ULL;
            key ^= key >> 33;
            return key;
        }

        weight_t get(partition_t u_id, partition_t v_id, size_t thread_id) const {
            ASSERT(u_id < m_k);
            ASSERT(v_id < m_k);
            ASSERT(thread_id < m_thread_states.size());
            
            if (u_id == v_id) return 0;

            partition_t min_u = std::min(u_id, v_id);
            partition_t max_v = std::max(u_id, v_id);
            uint64_t key = ((uint64_t)min_u << 32) | (uint64_t)max_v;
            size_t idx = mix_hash(key) & (NUM_SETS - 1);

            ThreadState& state = m_thread_states[thread_id];
            
            if (state.cache[idx].entries[0].u == min_u && state.cache[idx].entries[0].v == max_v) {
                state.cache[idx].lru = 1;
                return state.cache[idx].entries[0].dist;
            }
            if (state.cache[idx].entries[1].u == min_u && state.cache[idx].entries[1].v == max_v) {
                state.cache[idx].lru = 0;
                return state.cache[idx].entries[1].dist;
            }

            state.pq.push({0, u_id});
            state.dist[u_id] = 0;
            state.touched.push_back(u_id);

            weight_t result = std::numeric_limits<weight_t>::max() / 2;

            while (!state.pq.empty()) {
                auto [d, u] = state.pq.top();
                state.pq.pop();

                if (u == v_id) {
                    result = d;
                    break;
                }

                if (d > state.dist[u]) continue;

                for (size_t i = m_graph->neighborhoods[u]; i < m_graph->neighborhoods[u + 1]; ++i) {
                    partition_t v = m_graph->edges_v[i];
                    weight_t w = m_graph->edges_w[i];

                    if (state.dist[u] + w < state.dist[v]) {
                        if (state.dist[v] == std::numeric_limits<weight_t>::max() / 2) {
                            state.touched.push_back(v);
                        }
                        state.dist[v] = state.dist[u] + w;
                        state.pq.push({state.dist[v], v});
                    }
                }
            }

            // Cleanup for the next query
            for (partition_t node : state.touched) {
                state.dist[node] = std::numeric_limits<weight_t>::max() / 2;
            }
            state.touched.clear();
            state.pq.clear(); // Doesn't reallocate!
            
            u8 lru_idx = state.cache[idx].lru;
            state.cache[idx].entries[lru_idx].u = min_u;
            state.cache[idx].entries[lru_idx].v = max_v;
            state.cache[idx].entries[lru_idx].dist = result;
            state.cache[idx].lru = 1 - lru_idx;

            return result;
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
