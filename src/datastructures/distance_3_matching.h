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

#ifndef HEIPROMAP_DISTANCE_3_MATCHING_H
#define HEIPROMAP_DISTANCE_3_MATCHING_H

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
    /**
     * Distance3Matching computes a distance-3 matching on a quotient graph.
     * All endpoints of matched edges are at distance >= 3 from each other in the quotient graph.
     */
    template<bool LARGE_K>
    class Distance3Matching {
        partition_t m_k = 0;
        u32 m_frozen_epoch = 1;
        AlignedArray<u32> m_vertex_frozen_epoch;

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

        // valid checks
        AlignedArray<u32> m_valid_frozen;
        u32 m_valid_epoch = 1;

        // static batches
        BitVectorArray m_matched_batches;
        BitVectorArray m_dist1_batches;
        BitVectorArray m_dist2_batches;
        std::vector<u64> m_forbidden_buffer;

    public:
        Distance3Matching() = default;

        void initialize(const partition_t k) {
            m_k = k;

            m_vertex_frozen_epoch.initialize(m_k, 0);
            m_frozen_epoch = 1;

            m_edge_round_epoch = 1;
            m_valid_frozen.initialize(m_k, 0);
            m_valid_epoch = 1;
            m_matched_batches.initialize(m_k, 512);
            m_dist1_batches.initialize(m_k, 512);
            m_dist2_batches.initialize(m_k, 512);
            m_forbidden_buffer.resize((512 + 63) / 64);
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
            HEIPROMAP_PROFILE_SCOPE("d3_matching", "Distance3Matching", "init_candidates");

            if (!m_candidates_initialized) {
                m_num_candidates = 0;

                for (partition_t i = 0; i < p_manager.n_active; ++i) {
                    partition_t u = p_manager.active_ids[i];

                    if (q_graph.degree(u) == 0) { continue ; }

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

            if (m_frozen_epoch >= std::numeric_limits<u32>::max() - 2) {
                m_vertex_frozen_epoch.fill(0);
                m_frozen_epoch = 1;
            }
            const u32 epoch_dist2 = m_frozen_epoch;
            const u32 epoch_dist1 = m_frozen_epoch + 1;
            m_frozen_epoch += 2;

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

            auto freeze_from = [&](const partition_t x) {
                if constexpr (LARGE_K) {
                    const auto &adj_x = q_graph.neighbors(x);
                    for (const auto &e1: adj_x) {
                        const partition_t n1 = e1.target;
                        if (m_vertex_frozen_epoch[n1] >= epoch_dist1) {
                            continue;
                        }
                        m_vertex_frozen_epoch[n1] = epoch_dist1;

                        const auto &adj_n1 = q_graph.neighbors(n1);
                        for (const auto &e2: adj_n1) {
                            const partition_t n2 = e2.target;
                            if (m_vertex_frozen_epoch[n2] < epoch_dist2) {
                                m_vertex_frozen_epoch[n2] = epoch_dist2;
                            }
                        }
                    }
                } else {
                    q_graph.for_each_neighbor(x, [&](const partition_t n1, const weight_t) {
                        if (m_vertex_frozen_epoch[n1] >= epoch_dist1) {
                            return;
                        }
                        m_vertex_frozen_epoch[n1] = epoch_dist1;

                        q_graph.for_each_neighbor(n1, [&](const partition_t n2, const weight_t) {
                            if (m_vertex_frozen_epoch[n2] < epoch_dist2) {
                                m_vertex_frozen_epoch[n2] = epoch_dist2;
                            }
                        });
                    });
                }
            };

            {
                HEIPROMAP_PROFILE_SCOPE("d3_matching", "Distance3Matching", "candidate_loop");
                for (size_t i = 0; i < m_num_candidates; ++i) {
                    const partition_t u_id = m_candidates[i];
                    if (m_vertex_frozen_epoch[u_id] >= epoch_dist2) {
                        add_to_next(u_id);
                        continue;
                    }

                    bool matched_u = false;
                    size_t unused_valid_edges = 0;

                    if constexpr (LARGE_K) {
                        const auto &adj_u = q_graph.neighbors(u_id);
                        for (const auto &edge: adj_u) {
                            const partition_t v_id = edge.target;
                            if (v_id <= u_id) {
                                continue;
                            }
                            if (active_this_round[u_id] == 0 && active_this_round[v_id] == 0) {
                                continue;
                            }

                            if (edge.used_epoch == m_edge_round_epoch) {
                                continue;
                            }

                            unused_valid_edges++;

                            if (m_vertex_frozen_epoch[v_id] >= epoch_dist2) {
                                continue;
                            }

                            matching.emplace_back(u_id, v_id);
                            edge.used_epoch = m_edge_round_epoch;

                            m_vertex_frozen_epoch[u_id] = epoch_dist1;
                            m_vertex_frozen_epoch[v_id] = epoch_dist1;
                            freeze_from(u_id);
                            freeze_from(v_id);
                            matched_u = true;
                            break;
                        }
                    } else {
                        q_graph.for_each_neighbor(u_id, [&](const partition_t v_id, const weight_t) {
                            if (matched_u) {
                                return;
                            }
                            if (v_id <= u_id) {
                                return;
                            }
                            if (active_this_round[u_id] == 0 && active_this_round[v_id] == 0) {
                                return;
                            }

                            if (is_edge_used(q_graph, u_id, v_id)) {
                                return;
                            }

                            unused_valid_edges++;

                            if (m_vertex_frozen_epoch[v_id] >= epoch_dist2) {
                                return;
                            }

                            matching.emplace_back(u_id, v_id);
                            mark_edge_used(q_graph, u_id, v_id);

                            m_vertex_frozen_epoch[u_id] = epoch_dist1;
                            m_vertex_frozen_epoch[v_id] = epoch_dist1;
                            freeze_from(u_id);
                            freeze_from(v_id);
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
            }

            swap(m_candidates, m_next_candidates);
            std::swap(m_num_candidates, m_num_next_candidates);
            return !matching.empty();
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

        template<typename QGraphT>
        void freeze_validation(const QGraphT &q_graph, const partition_t x, const u32 cur_dist1, const u32 cur_dist2) {
            m_valid_frozen[x] = cur_dist1;
            if constexpr (LARGE_K) {
                const auto &adj_x = q_graph.neighbors(x);
                for (const auto &e1: adj_x) {
                    const partition_t n1 = e1.target;
                    if (m_valid_frozen[n1] >= cur_dist1) {
                        continue;
                    }
                    m_valid_frozen[n1] = cur_dist1;

                    const auto &adj_n1 = q_graph.neighbors(n1);
                    for (const auto &e2: adj_n1) {
                        const partition_t n2 = e2.target;
                        if (m_valid_frozen[n2] < cur_dist2) {
                            m_valid_frozen[n2] = cur_dist2;
                        }
                    }
                }
            } else {
                q_graph.for_each_neighbor(x, [&](const partition_t n1, const weight_t) {
                    if (m_valid_frozen[n1] >= cur_dist1) {
                        return;
                    }
                    m_valid_frozen[n1] = cur_dist1;

                    q_graph.for_each_neighbor(n1, [&](const partition_t n2, const weight_t) {
                        if (m_valid_frozen[n2] < cur_dist2) {
                            m_valid_frozen[n2] = cur_dist2;
                        }
                    });
                });
            }
        }

    public:
        template<typename QGraphT>
        void compute_static_matchings(const QGraphT &q_graph,
                                      const PartitionManager &p_manager,
                                      const AlignedArray<u8> &active_this_round,
                                      std::vector<std::vector<std::pair<partition_t, partition_t> > > &matchings,
                                      const size_t max_rounds = 512,
                                      const u64 min_threshold = 10) {
            HEIPROMAP_PROFILE_SCOPE("d3_matching", "Distance3Matching", "compute_static_matchings");
            matchings.clear();
            if (m_k == 0) return;

            const size_t rounds = std::min(max_rounds, static_cast<size_t>(512));
            if (m_matched_batches.num_bits() != rounds || m_matched_batches.num_elements() != m_k) {
                m_matched_batches.initialize(m_k, rounds);
                m_dist1_batches.initialize(m_k, rounds);
                m_dist2_batches.initialize(m_k, rounds);
                m_forbidden_buffer.resize(m_matched_batches.words_per_element());
            } else {
                m_matched_batches.clear();
                m_dist1_batches.clear();
                m_dist2_batches.clear();
            }

            matchings.resize(rounds);
            const size_t words = m_matched_batches.words_per_element();
            u64 *forbidden = m_forbidden_buffer.data();

            const auto freeze_endpoint = [&](const partition_t x, const size_t b) {
                m_matched_batches.set_bit(x, b);
                m_dist1_batches.set_bit(x, b);
                m_dist2_batches.set_bit(x, b);

                if constexpr (LARGE_K) {
                    for (const auto &e1: q_graph.neighbors(x)) {
                        const partition_t n1 = e1.target;
                        if (!m_dist1_batches.test_bit(n1, b)) {
                            m_dist1_batches.set_bit(n1, b);
                            m_dist2_batches.set_bit(n1, b);

                            for (const auto &e2: q_graph.neighbors(n1)) {
                                m_dist2_batches.set_bit(e2.target, b);
                            }
                        }
                    }
                } else {
                    q_graph.for_each_neighbor(x, [&](const partition_t n1, const weight_t) {
                        if (!m_dist1_batches.test_bit(n1, b)) {
                            m_dist1_batches.set_bit(n1, b);
                            m_dist2_batches.set_bit(n1, b);

                            q_graph.for_each_neighbor(n1, [&](const partition_t n2, const weight_t) {
                                m_dist2_batches.set_bit(n2, b);
                            });
                        }
                    });
                }
            };

            const auto process_vertex = [&](const partition_t u) {
                if constexpr (LARGE_K) {
                    const auto &adj_u = q_graph.neighbors(u);
                    const u64 *matched_u = m_matched_batches.get_words(u);
                    const u64 *d2_u = m_dist2_batches.get_words(u);

                    for (const auto &edge: adj_u) {
                        const partition_t v = edge.target;
                        if (v <= u) continue;
                        if (active_this_round[u] == 0 && active_this_round[v] == 0) continue;

                        const u64 *matched_v = m_matched_batches.get_words(v);
                        const u64 *d2_v = m_dist2_batches.get_words(v);

                        for (size_t w = 0; w < words; ++w) {
                            forbidden[w] = matched_u[w] | matched_v[w] | d2_u[w] | d2_v[w];
                        }

                        const size_t b = m_matched_batches.find_first_zero_bit(forbidden);
                        if (b < rounds) {
                            matchings[b].emplace_back(u, v);
                            freeze_endpoint(u, b);
                            freeze_endpoint(v, b);
                        }
                    }
                } else {
                    q_graph.for_each_neighbor(u, [&](const partition_t v, const weight_t) {
                        if (v <= u) return;
                        if (active_this_round[u] == 0 && active_this_round[v] == 0) return;

                        const u64 *matched_u = m_matched_batches.get_words(u);
                        const u64 *d2_u = m_dist2_batches.get_words(u);
                        const u64 *matched_v = m_matched_batches.get_words(v);
                        const u64 *d2_v = m_dist2_batches.get_words(v);

                        for (size_t w = 0; w < words; ++w) {
                            forbidden[w] = matched_u[w] | matched_v[w] | d2_u[w] | d2_v[w];
                        }

                        const size_t b = m_matched_batches.find_first_zero_bit(forbidden);
                        if (b < rounds) {
                            matchings[b].emplace_back(u, v);
                            freeze_endpoint(u, b);
                            freeze_endpoint(v, b);
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

        template<typename QGraphT>
        void filter_valid_matching(const QGraphT &q_graph,
                                   const std::vector<std::pair<partition_t, partition_t> > &input_matching,
                                   std::vector<std::pair<partition_t, partition_t> > &filtered_matching) {
            filtered_matching.clear();
            if (input_matching.empty()) {
                return;
            }

            if (m_valid_epoch >= std::numeric_limits<u32>::max() - 2) {
                m_valid_frozen.fill(0);
                m_valid_epoch = 1;
            }
            const u32 epoch_dist2 = m_valid_epoch;
            const u32 epoch_dist1 = m_valid_epoch + 1;
            m_valid_epoch += 2;

            for (const auto &pair: input_matching) {
                const partition_t u = pair.first;
                const partition_t v = pair.second;

                if (!q_graph.has_edge(u, v)) {
                    continue;
                }

                if (m_valid_frozen[u] >= epoch_dist2 || m_valid_frozen[v] >= epoch_dist2) {
                    continue;
                }

                filtered_matching.push_back(pair);
                freeze_validation(q_graph, u, epoch_dist1, epoch_dist2);
                freeze_validation(q_graph, v, epoch_dist1, epoch_dist2);
            }
        }

        /**
         * Incremental filter: two-pass approach.
         * Pass 1: Accept clean pairs (both endpoints not dirty) without any freeze.
         * Pass 2: Validate dirty pairs against clean endpoints + each other.
         * dirty_blocks should include the 2-hop neighborhood of structurally changed blocks.
         */
        template<typename QGraphT>
        void filter_valid_matching_incremental(const QGraphT &q_graph,
                                               const std::vector<std::pair<partition_t, partition_t> > &input_matching,
                                               const AlignedArray<u8> &dirty_blocks,
                                               std::vector<std::pair<partition_t, partition_t> > &filtered_matching) {
            filtered_matching.clear();
            if (input_matching.empty()) {
                return;
            }

            if (m_valid_epoch >= std::numeric_limits<u32>::max() - 3) {
                m_valid_frozen.fill(0);
                m_valid_epoch = 1;
            }
            const u32 epoch_dist2 = m_valid_epoch;
            const u32 epoch_dist1 = m_valid_epoch + 1;
            const u32 epoch_clean = m_valid_epoch + 2;
            m_valid_epoch += 3;

            // Pass 1: Accept clean pairs, mark endpoints at epoch_clean
            for (const auto &pair: input_matching) {
                const partition_t u = pair.first;
                const partition_t v = pair.second;
                if (!dirty_blocks[u] && !dirty_blocks[v]) {
                    filtered_matching.push_back(pair);
                    m_valid_frozen[u] = epoch_clean;
                    m_valid_frozen[v] = epoch_clean;
                }
            }

            // Pass 2: Validate dirty pairs
            for (const auto &pair: input_matching) {
                const partition_t u = pair.first;
                const partition_t v = pair.second;
                if (!dirty_blocks[u] && !dirty_blocks[v]) continue;

                if (!q_graph.has_edge(u, v)) continue;

                // Check against already-accepted dirty pairs (epoch_dist2 marks)
                if (m_valid_frozen[u] >= epoch_dist2 || m_valid_frozen[v] >= epoch_dist2) {
                    continue;
                }

                // Check distance-3 against clean pairs: 2-hop BFS for clean endpoints
                if (has_clean_conflict(q_graph, u, epoch_clean) ||
                    has_clean_conflict(q_graph, v, epoch_clean)) {
                    continue;
                }

                filtered_matching.push_back(pair);
                freeze_validation(q_graph, u, epoch_dist1, epoch_dist2);
                freeze_validation(q_graph, v, epoch_dist1, epoch_dist2);
            }
        }

    private:
        template<typename QGraphT>
        bool has_clean_conflict(const QGraphT &q_graph, const partition_t x, const u32 epoch_clean) const {
            if (m_valid_frozen[x] >= epoch_clean) return true;
            if constexpr (LARGE_K) {
                for (const auto &e1: q_graph.neighbors(x)) {
                    if (m_valid_frozen[e1.target] >= epoch_clean) return true;
                    for (const auto &e2: q_graph.neighbors(e1.target)) {
                        if (m_valid_frozen[e2.target] >= epoch_clean) return true;
                    }
                }
            } else {
                bool found = false;
                q_graph.for_each_neighbor(x, [&](const partition_t n1, const weight_t) {
                    if (found) return;
                    if (m_valid_frozen[n1] >= epoch_clean) {
                        found = true;
                        return;
                    }
                    q_graph.for_each_neighbor(n1, [&](const partition_t n2, const weight_t) {
                        if (found) return;
                        if (m_valid_frozen[n2] >= epoch_clean) { found = true; }
                    });
                });
                return found;
            }
            return false;
        }
    };
}

#endif //HEIPROMAP_DISTANCE_3_MATCHING_H
