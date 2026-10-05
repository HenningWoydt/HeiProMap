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

#ifndef HEIPROMAP_DISTANCE_1_MATCHING_H
#define HEIPROMAP_DISTANCE_1_MATCHING_H

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

#include "../definitions.h"
#include "../datastructures/partition_manager.h"
#include "../datastructures/bit_vector_array.h"
#include "../utility/aligned_array.h"
#include "../utility/profiler.h"

namespace HeiProMap {
    template<bool LARGE_K>
    class Distance1Matching {
        partition_t m_k = 0;
        u32 m_matched_epoch = 1;
        AlignedArray<u32> m_vertex_matched_epoch;

        // Used edge tracking
        std::vector<size_t> m_used_edge_indices;
        AlignedArray<u8> m_used_this_round;
        std::vector<std::vector<partition_t> > m_used_neighbors;
        AlignedArray<partition_t> m_touched_vertices;
        size_t m_num_touched = 0;
        u32 m_edge_round_epoch = 1;

        // Candidate queues
        AlignedArray<partition_t> m_candidates;
        size_t m_num_candidates = 0;
        AlignedArray<partition_t> m_next_candidates;
        size_t m_num_next_candidates = 0;
        AlignedArray<u32> m_in_candidates_epoch;
        u32 m_candidate_epoch = 1;
        bool m_candidates_initialized = false;

        // Static batches
        BitVectorArray m_matched_batches;
        std::vector<u64> m_forbidden_buffer;

    public:
        Distance1Matching() = default;

        size_t heap_bytes() const {
            size_t bytes = m_vertex_matched_epoch.heap_bytes() + m_used_this_round.heap_bytes()
                         + m_touched_vertices.heap_bytes() + m_candidates.heap_bytes()
                         + m_next_candidates.heap_bytes() + m_in_candidates_epoch.heap_bytes()
                         + m_matched_batches.heap_bytes()
                         + m_used_edge_indices.capacity() * sizeof(size_t)
                         + m_forbidden_buffer.capacity() * sizeof(u64);
            for (partition_t id = 0; id < m_k; ++id) {
                bytes += m_used_neighbors[id].capacity() * sizeof(partition_t);
            }
            return bytes;
        }

        void initialize(const partition_t k) {
            m_k = k;

            m_vertex_matched_epoch.initialize(m_k, 0);
            m_matched_epoch = 1;

            m_edge_round_epoch = 1;
            m_matched_batches.initialize(m_k, 512);
            m_forbidden_buffer.resize(m_matched_batches.words_per_element());
            m_used_edge_indices.clear();
            m_candidates.initialize(m_k);
            m_num_candidates = 0;
            m_next_candidates.initialize(m_k);
            m_num_next_candidates = 0;
            m_in_candidates_epoch.initialize(m_k, 0);
            m_candidate_epoch = 1;
            m_candidates_initialized = false;

            if constexpr (!LARGE_K) {
                m_used_this_round.initialize(m_k * m_k, 0);
            } else {
                m_used_neighbors.clear();
                m_used_neighbors.resize(m_k);
                m_touched_vertices.initialize(m_k);
                m_num_touched = 0;
            }
        }

        void reset_used_edges() {
            m_num_candidates = 0;
            m_candidates_initialized = false;
            m_edge_round_epoch++;
            if (m_edge_round_epoch == 0) {
                m_edge_round_epoch = 1;
            }

            if constexpr (!LARGE_K) {
                for (const size_t eidx: m_used_edge_indices) {
                    if (eidx < m_used_this_round.size()) {
                        m_used_this_round[eidx] = 0;
                    }
                }
                m_used_edge_indices.clear();
            } else {
                for (size_t i = 0; i < m_num_touched; ++i) {
                    m_used_neighbors[m_touched_vertices[i]].clear();
                }
                m_num_touched = 0;
            }
        }

        template<typename QGraphT>
        void init_candidates(const QGraphT &q_graph,
                             const PartitionManager &p_manager,
                             const AlignedArray<u8> &active_this_round) {
            HEIPROMAP_PROFILE_SCOPE("d1_matching", "Distance1Matching", "init_candidates");

            if (!m_candidates_initialized) {
                m_num_candidates = 0;

                for (partition_t i = 0; i < p_manager.n_active; ++i) {
                    partition_t u = p_manager.active_ids[i];

                    if (q_graph.degree(u) == 0) { continue; }

                    bool has_valid_edge = false;
                    q_graph.for_each_neighbor(u, [&](const partition_t v, const weight_t) {
                        if (has_valid_edge) return;
                        if (v > u && (active_this_round[u] != 0 || active_this_round[v] != 0)) {
                            has_valid_edge = true;
                        }
                    });
                    if (has_valid_edge) {
                        m_candidates[m_num_candidates++] = u;
                    }
                }

                m_candidates_initialized = true;
            }
        }

        template<typename QGraphT>
        bool find_matching(const QGraphT &q_graph,
                           const PartitionManager &p_manager,
                           const AlignedArray<u8> &active_this_round,
                           std::vector<std::pair<partition_t, partition_t> > &matching) {
            matching.clear();

            init_candidates(q_graph, p_manager, active_this_round);

            if (m_num_candidates == 0) { return false; }

            if (m_matched_epoch >= std::numeric_limits<u32>::max() - 1) {
                m_vertex_matched_epoch.fill(0);
                m_matched_epoch = 1;
            }
            const u32 cur_epoch = m_matched_epoch++;

            m_candidate_epoch++;
            if (m_candidate_epoch == 0) {
                m_in_candidates_epoch.fill(0);
                m_candidate_epoch = 1;
            }
            const u32 cur_candidate_epoch = m_candidate_epoch;

            m_num_next_candidates = 0;
            auto add_to_next = [&](const partition_t x) {
                if (m_in_candidates_epoch[x] != cur_candidate_epoch) {
                    m_in_candidates_epoch[x] = cur_candidate_epoch;
                    m_next_candidates[m_num_next_candidates++] = x;
                }
            };

            for (size_t i = 0; i < m_num_candidates; ++i) {
                const partition_t u_id = m_candidates[i];
                if (m_vertex_matched_epoch[u_id] >= cur_epoch) {
                    add_to_next(u_id);
                    continue;
                }

                bool matched_u = false;
                size_t unused_valid_edges = 0;

                if constexpr (LARGE_K) {
                    const auto &adj_u = q_graph.neighbors(u_id);
                    for (const auto &edge: adj_u) {
                        const partition_t v_id = edge.target;
                        if (v_id <= u_id) continue;
                        if (active_this_round[u_id] == 0 && active_this_round[v_id] == 0) continue;
                        if (edge.used_epoch == m_edge_round_epoch) continue;

                        unused_valid_edges++;

                        if (m_vertex_matched_epoch[v_id] >= cur_epoch) continue;

                        matching.emplace_back(u_id, v_id);
                        edge.used_epoch = m_edge_round_epoch;
                        m_vertex_matched_epoch[u_id] = cur_epoch;
                        m_vertex_matched_epoch[v_id] = cur_epoch;
                        matched_u = true;
                        break;
                    }
                } else {
                    q_graph.for_each_neighbor(u_id, [&](const partition_t v_id, const weight_t) {
                        if (matched_u) return;
                        if (v_id <= u_id) return;
                        if (active_this_round[u_id] == 0 && active_this_round[v_id] == 0) return;
                        if (is_edge_used(q_graph, u_id, v_id)) return;

                        unused_valid_edges++;

                        if (m_vertex_matched_epoch[v_id] >= cur_epoch) return;

                        matching.emplace_back(u_id, v_id);
                        mark_edge_used(q_graph, u_id, v_id);
                        m_vertex_matched_epoch[u_id] = cur_epoch;
                        m_vertex_matched_epoch[v_id] = cur_epoch;
                        matched_u = true;
                    });
                }

                if (matched_u) {
                    add_to_next(u_id);
                    add_to_next(matching.back().second);
                } else if (unused_valid_edges > 0) {
                    add_to_next(u_id);
                }
            }

            swap(m_candidates, m_next_candidates);
            std::swap(m_num_candidates, m_num_next_candidates);
            return !matching.empty();
        }

        template<typename QGraphT>
        void compute_static_matchings(const QGraphT &q_graph,
                                      const PartitionManager &p_manager,
                                      const AlignedArray<u8> &active_this_round,
                                      std::vector<std::vector<std::pair<partition_t, partition_t> > > &matchings,
                                      const size_t max_rounds = 512,
                                      const u64 min_threshold = 10) {
            HEIPROMAP_PROFILE_SCOPE("d1_matching", "Distance1Matching", "compute_static_matchings");
            matchings.clear();
            if (m_k == 0) return;

            const size_t rounds = std::min(max_rounds, static_cast<size_t>(512));
            if (m_matched_batches.num_bits() != rounds || m_matched_batches.num_elements() != m_k) {
                m_matched_batches.initialize(m_k, rounds);
                m_forbidden_buffer.resize(m_matched_batches.words_per_element());
            } else {
                m_matched_batches.clear();
            }

            matchings.resize(rounds);
            const size_t words = m_matched_batches.words_per_element();
            u64 *forbidden = m_forbidden_buffer.data();

            const auto process_vertex = [&](const partition_t u) {
                if constexpr (LARGE_K) {
                    const auto &adj_u = q_graph.neighbors(u);
                    const u64 *matched_u = m_matched_batches.get_words(u);

                    for (const auto &edge: adj_u) {
                        const partition_t v = edge.target;
                        if (v <= u) continue;
                        if (active_this_round[u] == 0 && active_this_round[v] == 0) continue;

                        const u64 *matched_v = m_matched_batches.get_words(v);

                        for (size_t w = 0; w < words; ++w) {
                            forbidden[w] = matched_u[w] | matched_v[w];
                        }

                        const size_t b = m_matched_batches.find_first_zero_bit(forbidden);
                        if (b < rounds) {
                            matchings[b].emplace_back(u, v);
                            m_matched_batches.set_bit(u, b);
                            m_matched_batches.set_bit(v, b);
                        }
                    }
                } else {
                    q_graph.for_each_neighbor(u, [&](const partition_t v, const weight_t) {
                        if (v <= u) return;
                        if (active_this_round[u] == 0 && active_this_round[v] == 0) return;

                        const u64 *matched_u = m_matched_batches.get_words(u);
                        const u64 *matched_v = m_matched_batches.get_words(v);

                        for (size_t w = 0; w < words; ++w) {
                            forbidden[w] = matched_u[w] | matched_v[w];
                        }

                        const size_t b = m_matched_batches.find_first_zero_bit(forbidden);
                        if (b < rounds) {
                            matchings[b].emplace_back(u, v);
                            m_matched_batches.set_bit(u, b);
                            m_matched_batches.set_bit(v, b);
                        }
                    });
                }
            };

            if (p_manager.n_active > 0) {
                for (partition_t i = 0; i < p_manager.n_active; ++i) {
                    process_vertex(p_manager.active_ids[i]);
                }
            } else {
                for (partition_t u = 0; u < m_k; ++u) {
                    if (p_manager.is_active(u)) {
                        process_vertex(u);
                    }
                }
            }

            size_t valid_count = 0;
            for (size_t b = 0; b < rounds; ++b) {
                if (matchings[b].size() >= min_threshold) {
                    if (valid_count != b) {
                        matchings[valid_count] = std::move(matchings[b]);
                    }
                    valid_count++;
                }
            }
            matchings.resize(valid_count);
        }

    private:
        template<typename QGraphT>
        bool is_edge_used(const QGraphT &q_graph, const partition_t u, const partition_t v) const {
            if constexpr (!LARGE_K) {
                const size_t eidx = q_graph.edge_index(u, v);
                return eidx < m_used_this_round.size() && m_used_this_round[eidx] == 1;
            } else {
                const auto &used = m_used_neighbors[u];
                return std::find(used.begin(), used.end(), v) != used.end();
            }
        }

        template<typename QGraphT>
        void mark_edge_used(const QGraphT &q_graph, const partition_t u, const partition_t v) {
            if constexpr (!LARGE_K) {
                const size_t eidx = q_graph.edge_index(u, v);
                if (eidx < m_used_this_round.size()) {
                    m_used_this_round[eidx] = 1;
                    m_used_edge_indices.push_back(eidx);
                }
            } else {
                if (m_used_neighbors[u].empty()) {
                    m_touched_vertices[m_num_touched++] = u;
                }
                m_used_neighbors[u].push_back(v);
            }
        }
    };
}

#endif //HEIPROMAP_DISTANCE_1_MATCHING_H
