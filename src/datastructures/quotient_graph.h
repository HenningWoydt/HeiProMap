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

#ifndef HEIPROMAP_QUOTIENT_GRAPH_H
#define HEIPROMAP_QUOTIENT_GRAPH_H

#include <unordered_set>
#include <random>

#include "../definitions.h"
#include "../utility/aligned_array.h"

namespace HeiProMap {
    class QuotientGraph {
        partition_t m_k = 0;

        AlignedArray<weight_t> m_adj_mtx;

    public:
        void initialize(const partition_t t_k) {
            ScopedTimer _t("io", "QuotientGraph", "initialize");

            m_k = t_k;

            size_t size = (size_t) m_k * (size_t) m_k;
            m_adj_mtx.initialize(size, 0);
        }

        void add_edge(const partition_t u_id, const partition_t v_id, const weight_t w) {
            partition_t min = std::min(u_id, v_id);
            partition_t max = std::max(u_id, v_id);
            m_adj_mtx[min * m_k + max] += w;
        }

        void remove_edge(const partition_t u_id, const partition_t v_id, const weight_t w) {
            partition_t min = std::min(u_id, v_id);
            partition_t max = std::max(u_id, v_id);
            m_adj_mtx[min * m_k + max] -= w;
        }

        bool has_edge(const partition_t u_id, const partition_t v_id) const {
            partition_t min = std::min(u_id, v_id);
            partition_t max = std::max(u_id, v_id);
            return m_adj_mtx[min * m_k + max] > 0;
        }

        weight_t get_weight(const partition_t u_id, const partition_t v_id) const {
            partition_t min = std::min(u_id, v_id);
            partition_t max = std::max(u_id, v_id);
            return m_adj_mtx[min * m_k + max];
        }

        /**
         * O(max_deg)
         *
         * @param g
         * @param p_manager
         * @param u
         * @param old_id
         * @param new_id
         */
        void move(const graph_t &g,
                  const p_manager_t &p_manager,
                  const vertex_t u,
                  const partition_t old_id,
                  const partition_t new_id) {
            ASSERT(new_id < m_k);
            ASSERT(new_id != old_id);

            forall_guivw(g, u, i, v, w)
                {
                    partition_t v_id = p_manager[v];

                    // remove old edge, if existed
                    if (old_id != v_id) {
                        remove_edge(old_id, v_id, w);
                    }
                    // add new edge, if has to exist
                    if (new_id != v_id) {
                        add_edge(new_id, v_id, w);
                    }
                }
            endfor
        }

        void copy_from(const QuotientGraph &q) {
            for (u64 i = 0; i < m_k * m_k; i++) {
                m_adj_mtx[i] = q.m_adj_mtx[i];
            }
        }

        bool find_distance_3_matching(AlignedArray<u8> &active_this_round,
                                      AlignedArray<u8> &used_edges_this_round,
                                      std::vector<std::pair<partition_t, partition_t> > &matching) {
            matching.clear();

            std::vector<u8> vertex_frozen(m_k, 0);

            auto freeze_distance_2 = [&](partition_t x) {
                vertex_frozen[x] = 1;

                for (partition_t n1 = 0; n1 < m_k; ++n1) {
                    if (n1 == x || !has_edge(x, n1)) {
                        continue;
                    }

                    vertex_frozen[n1] = 1;

                    for (partition_t n2 = 0; n2 < m_k; ++n2) {
                        if (n2 == n1 || !has_edge(n1, n2)) {
                            continue;
                        }
                        vertex_frozen[n2] = 1;
                    }
                }
            };

            for (partition_t u_id = 0; u_id < m_k; ++u_id) {
                if (vertex_frozen[u_id] == 1) { continue; }

                for (partition_t v_id = u_id + 1; v_id < m_k; ++v_id) {
                    if (!has_edge(u_id, v_id)) { continue; }
                    if (vertex_frozen[v_id] == 1) { continue; }
                    if (active_this_round[u_id] == 0 && active_this_round[v_id] == 0) { continue; }

                    size_t eidx = edge_index(u_id, v_id);
                    if (used_edges_this_round[eidx] == 1) { continue; }

                    matching.emplace_back(u_id, v_id);
                    used_edges_this_round[eidx] = 1;

                    freeze_distance_2(u_id);
                    freeze_distance_2(v_id);
                    break;
                }
            }

            return !matching.empty();
        }

        size_t edge_index(const partition_t u_id, const partition_t v_id) const {
            partition_t min_id = std::min(u_id, v_id);
            partition_t max_id = std::max(u_id, v_id);
            return (size_t) min_id * (size_t) m_k + (size_t) max_id;
        }

        bool is_valid_distance_3_matching(const std::vector<std::pair<partition_t, partition_t> > &matching) const {
            std::vector<u8> is_in_matching(m_k, 0);

            for (const auto &[u_id, v_id]: matching) {
                if (u_id >= m_k || v_id >= m_k) {
                    std::cout << "Not Valid 3 Distance Matching" << std::endl;
                    return false;
                }
                if (u_id == v_id) {
                    std::cout << "Not Valid 3 Distance Matching" << std::endl;
                    return false;
                }
                if (!has_edge(u_id, v_id)) {
                    std::cout << "Not Valid 3 Distance Matching" << std::endl;
                    return false;
                }
                if (is_in_matching[u_id] == 1 || is_in_matching[v_id] == 1) {
                    std::cout << "Not Valid 3 Distance Matching" << std::endl;
                    return false;
                }

                is_in_matching[u_id] = 1;
                is_in_matching[v_id] = 1;
            }

            for (const auto &[u_id, v_id]: matching) {
                // Check all vertices within distance <= 2 from u_id
                for (partition_t id = 0; id < m_k; ++id) {
                    if (id != u_id && id != v_id && has_edge(u_id, id) && is_in_matching[id] == 1) {
                        std::cout << "Not Valid 3 Distance Matching" << std::endl;
                        return false;
                    }

                    if (!has_edge(u_id, id)) {
                        continue;
                    }

                    for (partition_t id2 = 0; id2 < m_k; ++id2) {
                        if (id2 != u_id && id2 != v_id && has_edge(id, id2) && is_in_matching[id2] == 1) {
                            std::cout << "Not Valid 3 Distance Matching" << std::endl;
                            return false;
                        }
                    }
                }

                // Check all vertices within distance <= 2 from v_id
                for (partition_t id = 0; id < m_k; ++id) {
                    if (id != u_id && id != v_id && has_edge(v_id, id) && is_in_matching[id] == 1) {
                        std::cout << "Not Valid 3 Distance Matching" << std::endl;
                        return false;
                    }

                    if (!has_edge(v_id, id)) {
                        continue;
                    }

                    for (partition_t id2 = 0; id2 < m_k; ++id2) {
                        if (id2 != u_id && id2 != v_id && has_edge(id, id2) && is_in_matching[id2] == 1) {
                            std::cout << "Not Valid 3 Distance Matching" << std::endl;
                            return false;
                        }
                    }
                }
            }

            return true;
        }

        void write_as_metis(const std::string &file_name) {
            std::ofstream out(file_name);
            if (!out) {
                throw std::runtime_error("Could not open file for writing: " + file_name);
            }

            // Count undirected edges (ignore self-loops).
            size_t num_edges = 0;
            for (partition_t u = 0; u < m_k; ++u) {
                for (partition_t v = u + 1; v < m_k; ++v) {
                    if (get_weight(u, v) > 0) {
                        ++num_edges;
                    }
                }
            }

            // METIS header:
            // <num vertices> <num edges> <fmt>
            // fmt = 1 => edge weights present
            out << m_k << " " << num_edges << " 1\n";

            // For each vertex, write adjacency list:
            // neighbor_id weight neighbor_id weight ...
            // METIS uses 1-based vertex numbering.
            for (partition_t u = 0; u < m_k; ++u) {
                bool first = true;

                for (partition_t v = 0; v < m_k; ++v) {
                    if (u == v) {
                        continue; // skip self-loops
                    }

                    weight_t w = get_weight(u, v);
                    if (w > 0) {
                        if (!first) {
                            out << " ";
                        }
                        out << (v + 1) << " " << w;
                        first = false;
                    }
                }

                out << "\n";
            }
        }

        std::vector<std::vector<partition_t> > get_rnd_cycles(u64 max_n_cycles,
                                                              u64 n_samples_per_start,
                                                              u64 min_cycle_length,
                                                              u64 max_cycle_length) {
            std::vector<std::vector<partition_t> > cycles;

            if (m_k == 0 || max_n_cycles == 0 || n_samples_per_start == 0 || max_cycle_length < 3) {
                return cycles;
            }

            std::mt19937_64 rng(std::random_device{}());

            auto canonicalize_cycle = [](std::vector<partition_t> &cyc) {
                if (cyc.empty()) {
                    return;
                }

                const size_t n = cyc.size();

                // find smallest element
                size_t min_pos = 0;
                for (size_t i = 1; i < n; ++i) {
                    if (cyc[i] < cyc[min_pos]) {
                        min_pos = i;
                    }
                }

                // build rotated version
                std::vector<partition_t> rot(n);
                for (size_t i = 0; i < n; ++i) {
                    rot[i] = cyc[(min_pos + i) % n];
                }

                // build reversed version with same first element
                std::vector<partition_t> rev(n);
                rev[0] = rot[0];
                for (size_t i = 1; i < n; ++i) {
                    rev[i] = rot[n - i];
                }

                if (rev < rot) {
                    cyc.swap(rev);
                } else {
                    cyc.swap(rot);
                }
            };

            auto same_cycle = [](const std::vector<partition_t> &a,
                                 const std::vector<partition_t> &b) -> bool {
                return a == b;
            };

            // timestamp-based visited structure:
            // mark[v] == cur_mark  <=>  v is in current path
            std::vector<u64> mark(m_k, 0);
            u64 cur_mark = 1;

            std::vector<partition_t> path;
            path.reserve(std::min<u64>(m_k, max_cycle_length));

            std::vector<partition_t> neighbors;
            neighbors.reserve(m_k);

            for (u64 sample = 0; sample < n_samples_per_start; ++sample) {
                for (partition_t start = 0; start < m_k; ++start) {
                    // handle rare timestamp overflow
                    std::fill(mark.begin(), mark.end(), 0);

                    path.clear();
                    path.push_back(start);
                    mark[start] = cur_mark;

                    partition_t cur = start;

                    for (u64 depth = 1; depth < max_cycle_length; ++depth) {
                        neighbors.clear();

                        for (partition_t v = 0; v < m_k; ++v) {
                            if (v == cur) {
                                continue;
                            }
                            if (!has_edge(cur, v)) {
                                continue;
                            }

                            if (v == start) {
                                if (path.size() >= min_cycle_length) {
                                    neighbors.push_back(v);
                                }
                            } else if (mark[v] != cur_mark) {
                                neighbors.push_back(v);
                            }
                        }

                        if (neighbors.empty()) {
                            break;
                        }

                        std::uniform_int_distribution<size_t> dist(0, neighbors.size() - 1);
                        partition_t nxt = neighbors[dist(rng)];

                        if (nxt == start) {
                            std::vector<partition_t> cyc(path);
                            canonicalize_cycle(cyc);

                            bool exists = false;
                            for (const auto &existing: cycles) {
                                if (same_cycle(existing, cyc)) {
                                    exists = true;
                                    break;
                                }
                            }

                            if (!exists) {
                                cycles.push_back(std::move(cyc));
                            }
                            break;
                        }

                        path.push_back(nxt);
                        mark[nxt] = cur_mark;
                        cur = nxt;
                    }

                    ++cur_mark;
                }
            }

            if (cycles.size() > max_n_cycles) {
                cycles.resize(max_n_cycles);
            }

            return cycles;
        }

        bool cycle_exists(const std::vector<partition_t> &cycle) {
            for (size_t i = 0; i < cycle.size(); ++i) {
                partition_t u = cycle[i];
                partition_t v = cycle[(i + 1) % cycle.size()];

                if (!has_edge(u, v)) {
                    return false;
                }
            }

            return true;
        }
    };
}

#endif //HEIPROMAP_QUOTIENT_GRAPH_H
