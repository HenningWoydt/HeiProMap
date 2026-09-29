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
#include <unordered_set>
#include <utility>
#include <vector>

#include "../definitions.h"
#include "../utility/aligned_array.h"

namespace HeiProMap {
    /**
     * Distance3Matching computes a distance-3 matching on a quotient graph.
     * All endpoints of matched edges are at distance >= 3 from each other in the quotient graph.
     */
    template<bool LARGE_K>
    class Distance3Matching {
        partition_t m_k = 0;
        u32 m_frozen_epoch = 1;
        std::vector<u32> m_vertex_frozen_epoch;
        std::vector<size_t> m_used_edge_indices;
        AlignedArray<u8> m_used_this_round;
        std::vector<std::vector<partition_t> > m_used_neighbors;
        std::vector<partition_t> m_touched_vertices;
        std::vector<partition_t> m_candidates;
        std::vector<partition_t> m_next_candidates;
        std::vector<u32> m_in_candidates_epoch;
        u32 m_candidate_epoch = 1;
        bool m_candidates_initialized = false;

    public:
        Distance3Matching() = default;

        void initialize(const partition_t k) {
            m_k = k;
            m_vertex_frozen_epoch.assign(m_k, 0);
            m_frozen_epoch = 1;
            m_used_edge_indices.clear();
            m_candidates.clear();
            m_next_candidates.clear();
            m_in_candidates_epoch.assign(m_k, 0);
            m_candidate_epoch = 1;
            m_candidates_initialized = false;

            if constexpr (!LARGE_K) {
                m_used_this_round.initialize(static_cast<size_t>(m_k) * static_cast<size_t>(m_k));
            } else {
                m_used_neighbors.clear();
                m_used_neighbors.resize(m_k);
                m_touched_vertices.clear();
            }
        }

        void reset_used_edges() {
            m_candidates.clear();
            m_candidates_initialized = false;

            if constexpr (!LARGE_K) {
                for (const size_t eidx : m_used_edge_indices) {
                    if (eidx < m_used_this_round.size()) {
                        m_used_this_round[eidx] = 0;
                    }
                }
                m_used_edge_indices.clear();
            } else {
                for (const partition_t u : m_touched_vertices) {
                    m_used_neighbors[u].clear();
                }
                m_touched_vertices.clear();
            }
        }

        template<typename QGraphT>
        bool find_matching(const QGraphT &q_graph,
                           const AlignedArray<u8> &active_this_round,
                           std::vector<std::pair<partition_t, partition_t> > &matching) {
            matching.clear();

            if (!m_candidates_initialized) {
                m_candidates.clear();
                for (partition_t u = 0; u < m_k; ++u) {
                    if (q_graph.degree(u) == 0) {
                        continue;
                    }
                    bool has_valid_edge = false;
                    q_graph.for_each_neighbor(u, [&](const partition_t v, const weight_t) {
                        if (has_valid_edge) return;
                        if (v > u && (active_this_round[u] != 0 || active_this_round[v] != 0)) {
                            has_valid_edge = true;
                        }
                    });
                    if (has_valid_edge) {
                        m_candidates.push_back(u);
                    }
                }
                m_candidates_initialized = true;
            }

            if (m_candidates.empty()) {
                return false;
            }

            if (m_frozen_epoch >= std::numeric_limits<u32>::max() - 2) {
                std::fill(m_vertex_frozen_epoch.begin(), m_vertex_frozen_epoch.end(), 0);
                m_frozen_epoch = 1;
            }
            const u32 epoch_dist2 = m_frozen_epoch;
            const u32 epoch_dist1 = m_frozen_epoch + 1;
            m_frozen_epoch += 2;

            m_candidate_epoch++;
            if (m_candidate_epoch == 0) {
                std::fill(m_in_candidates_epoch.begin(), m_in_candidates_epoch.end(), 0);
                m_candidate_epoch = 1;
            }
            const u32 cur_candidate_epoch = m_candidate_epoch;

            m_next_candidates.clear();
            auto add_to_next = [&](const partition_t x) {
                if (m_in_candidates_epoch[x] != cur_candidate_epoch) {
                    m_in_candidates_epoch[x] = cur_candidate_epoch;
                    m_next_candidates.push_back(x);
                }
            };

            auto freeze_from = [&](const partition_t x) {
                q_graph.for_each_neighbor(x, [&](const partition_t n1, const weight_t) {
                    if (m_vertex_frozen_epoch[n1] == epoch_dist1) {
                        return;
                    }
                    m_vertex_frozen_epoch[n1] = epoch_dist1;

                    q_graph.for_each_neighbor(n1, [&](const partition_t n2, const weight_t) {
                        if (m_vertex_frozen_epoch[n2] < epoch_dist2) {
                            m_vertex_frozen_epoch[n2] = epoch_dist2;
                        }
                    });
                });
            };

            for (const partition_t u_id : m_candidates) {
                if (m_vertex_frozen_epoch[u_id] >= epoch_dist2) {
                    add_to_next(u_id);
                    continue;
                }

                bool matched_u = false;
                size_t unused_valid_edges = 0;

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

                    if constexpr (!LARGE_K) {
                        const size_t eidx = q_graph.edge_index(u_id, v_id);
                        if (eidx < m_used_this_round.size() && m_used_this_round[eidx] == 1) {
                            return;
                        }
                    } else {
                        const auto &used = m_used_neighbors[u_id];
                        if (std::find(used.begin(), used.end(), v_id) != used.end()) {
                            return;
                        }
                    }

                    unused_valid_edges++;

                    if (m_vertex_frozen_epoch[v_id] >= epoch_dist2) {
                        return;
                    }

                    matching.emplace_back(u_id, v_id);
                    if constexpr (!LARGE_K) {
                        const size_t eidx = q_graph.edge_index(u_id, v_id);
                        if (eidx < m_used_this_round.size()) {
                            m_used_this_round[eidx] = 1;
                            m_used_edge_indices.push_back(eidx);
                        }
                    } else {
                        if (m_used_neighbors[u_id].empty()) {
                            m_touched_vertices.push_back(u_id);
                        }
                        m_used_neighbors[u_id].push_back(v_id);
                    }

                    m_vertex_frozen_epoch[u_id] = epoch_dist1;
                    m_vertex_frozen_epoch[v_id] = epoch_dist1;
                    freeze_from(u_id);
                    freeze_from(v_id);
                    matched_u = true;
                });

                if (matched_u) {
                    add_to_next(u_id);
                    add_to_next(matching.back().second);
                } else if (unused_valid_edges > 0) {
                    add_to_next(u_id);
                }
            }

            std::swap(m_candidates, m_next_candidates);
            return !matching.empty();
        }
    };
}

#endif //HEIPROMAP_DISTANCE_3_MATCHING_H
