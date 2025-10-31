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

#ifndef HEIPROMAP_CSR_GRAPH_H
#define HEIPROMAP_CSR_GRAPH_H

#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include "../utility/aligned_array.h"
#include "../definitions.h"
#include "../utility/matching.h"
#include "../utility/mapping.h"
#include "../utility/utils.h"
#include "../utility/profiler.h"

namespace HeiProMap {
    static inline uint64_t splitmix64(uint64_t x) {
        x += 0x9E3779B97F4A7C15ull;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
        return x ^ (x >> 31);
    }

    // Hash a vertex id (optionally salted by map_u)
    static inline uint64_t hash_vertex(vertex_t v, vertex_t salt = 0) {
        uint64_t x = static_cast<uint64_t>(v) ^ (static_cast<uint64_t>(salt) * 0x9E3779B97F4A7C15ull);
        return splitmix64(x);
    }

    /**
    * Standard undirected Graph that can hold vertex and edge weights.
    */
    class CSRGraph {
        vertex_t m_n = 0;
        vertex_t m_m = 0;

        weight_t m_graph_weight = 0;

        AlignedArray<weight_t> m_v_weights;
        AlignedArray<size_t> m_neighborhoods;
        AlignedArray<vertex_t> m_edges_v;
        AlignedArray<weight_t> m_edges_w;

    public:
        CSRGraph() = default;

        explicit CSRGraph(const std::string &file_path) {
            ScopedTimer _t_allocate("io", "CSRGraph", "allocate");
            if (!file_exists(file_path)) {
                std::cerr << "File " << file_path << " does not exist!" << std::endl;
                exit(EXIT_FAILURE);
            }

            // mmap the whole file
            MMap mm = mmap_file_ro(file_path);
            char *p = mm.data;
            const char *end = mm.data + mm.size; // (tiny fix: don't do -1)

            _t_allocate.stop();
            ScopedTimer _t_read_header("io", "CSRGraph", "read_header");

            // skip comment lines
            while (*p == '%') {
                while (*p != '\n') { ++p; }
                ++p;
            }

            // skip whitespace
            while (*p == ' ') { ++p; }

            // read number of vertices
            m_n = 0;
            while (*p != ' ' && *p != '\n') {
                m_n = m_n * 10 + (vertex_t) (*p - '0');
                ++p;
            }

            // skip whitespace
            while (*p == ' ') { ++p; }

            // read number of edges
            m_m = 0;
            while (*p != ' ' && *p != '\n') {
                m_m = m_m * 10 + (vertex_t) (*p - '0');
                ++p;
            }
            m_m *= 2;

            // search end of line or fmt
            std::string fmt = "000";
            bool has_v_weights = false;
            bool has_e_weights = false;
            while (*p == ' ') { ++p; }
            if (*p != '\n') {
                // found fmt
                fmt[0] = *p;
                ++p;
                if (*p != '\n') {
                    // found fmt
                    fmt[1] = *p;
                    ++p;
                    if (*p != '\n') {
                        // found fmt
                        fmt[2] = *p;
                        ++p;
                    }
                }
                // skip whitespaces
                while (*p == ' ') { ++p; }
            }
            m_graph_weight = 0;
            m_v_weights.initialize(m_n);
            m_neighborhoods.initialize(m_n + 1);
            m_neighborhoods[0] = 0;
            m_edges_v.initialize(m_m);
            m_edges_w.initialize(m_m);
            has_v_weights = fmt[1] == '1';
            has_e_weights = fmt[2] == '1';

            _t_read_header.stop();
            ScopedTimer _t_read_edges("io", "CSRGraph", "read_edges");

            ++p;
            vertex_t u = 0;
            size_t curr_m = 0;
            while (p < end) {
                // skip comment lines
                while (*p == '%') {
                    while (*p != '\n') { ++p; }
                    ++p;
                }

                // skip whitespaces
                while (*p == ' ') { ++p; }

                // read in vertex weight
                weight_t vw = 1;
                if (has_v_weights) {
                    vw = 0;
                    while (*p != ' ' && *p != '\n') {
                        vw = vw * 10 + (weight_t) (*p - '0');
                        ++p;
                    }

                    // skip whitespaces
                    while (*p == ' ') { ++p; }
                }
                m_v_weights[u] = vw;
                m_graph_weight += vw;

                // read in edges
                while (*p != '\n') {
                    vertex_t v = 0;
                    weight_t w = 1;

                    while (*p != ' ' && *p != '\n') {
                        v = v * 10 + (vertex_t) (*p - '0');
                        ++p;
                    }

                    // skip whitespaces
                    while (*p == ' ') { ++p; }

                    if (has_e_weights) {
                        w = 0;
                        while (*p != ' ' && *p != '\n') {
                            w = w * 10 + (weight_t) (*p - '0');
                            ++p;
                        }

                        // skip whitespaces
                        while (*p == ' ') { ++p; }
                    }

                    m_edges_v[curr_m] = v - 1;
                    m_edges_w[curr_m] = w;
                    ++curr_m;
                }
                m_neighborhoods[u + 1] = curr_m;
                ++u;
                ++p;
            }

            if (curr_m != m_m) {
                std::cerr << "Number of expected edges " << m_m << " not equal to number edges " << curr_m << " found!\n";
                munmap_file(mm);
                exit(EXIT_FAILURE);
            }

            _t_read_edges.stop();
            // done with the file
            munmap_file(mm);
        }

        void initialize(const CSRGraph &g,
                        const Mapping &mapping) {
            ScopedTimer _t_allocate("contraction", "CSRGraph", "allocate");

            m_n = mapping.get_coarse_n();
            m_graph_weight = g.m_graph_weight;
            m_v_weights.initialize(m_n, 0);

            _t_allocate.stop();
            ScopedTimer _t_overest_sizes("contraction", "CSRGraph", "overest_sizes");

            AlignedArray<vertex_t> overest_sizes;
            overest_sizes.initialize(m_n, 0);

            // overestimate neighborhood sizes and collect weights
            for (vertex_t u = 0; u < mapping.get_old_n(); ++u) {
                vertex_t map_u = mapping.get_map_u(u);
                overest_sizes[map_u] += g.size(u);
                m_v_weights[map_u] += g.weight(u);
            }

            _t_overest_sizes.stop();
            ScopedTimer _t_prefix_sum("contraction", "CSRGraph", "prefix_sum");

            // prefix-sum on overestimated neighborhood
            AlignedArray<vertex_t> overest_neighborhood;
            overest_neighborhood.initialize(m_n + 1);
            overest_neighborhood[0] = 0;
            for (vertex_t map_u = 0; map_u < m_n; ++map_u) {
                overest_neighborhood[map_u + 1] = overest_neighborhood[map_u] + overest_sizes[map_u];
            }

            _t_prefix_sum.stop();
            ScopedTimer _t_insert_edges("contraction", "CSRGraph", "insert_edges");

            // insert edges in overestimated array
            m_m = 0;
            AlignedArray<vertex_t> edges_v;
            AlignedArray<weight_t> edges_w;
            edges_v.initialize(overest_neighborhood[m_n], g.get_n());
            edges_w.initialize(overest_neighborhood[m_n]);

            AlignedArray<vertex_t> sizes;
            sizes.initialize(m_n, 0);

            for (vertex_t u = 0; u < mapping.get_old_n(); ++u) {
                vertex_t map_u = mapping.get_map_u(u);

                forall_guivw(g, u, i, v, w) {
                    vertex_t map_v = mapping.get_map_u(v);
                    if (map_u == map_v) { continue; }

                    vertex_t beg = overest_neighborhood[map_u];
                    vertex_t end = overest_neighborhood[map_u + 1];
                    vertex_t len = end - beg;
                    if (len == 0) { continue; }

                    // insert map_v
                    vertex_t j = beg + (hash_vertex(map_v) % len);
                    while (true) {
                        if (j == end) { j = beg; }
                        if (edges_v[j] == map_v) {
                            edges_w[j] += w;
                            break;
                        }
                        if (edges_v[j] == g.get_n()) {
                            edges_v[j] = map_v;
                            edges_w[j] = w;
                            m_m += 1;
                            sizes[map_u] += 1;
                            break;
                        }
                        j += 1;
                    }
                }
                endfor
            }

            _t_insert_edges.stop();
            ScopedTimer _t_real_neighborhood("contraction", "CSRGraph", "real_neighborhood");

            // insert edges in real array
            m_neighborhoods.initialize(m_n + 1);
            m_neighborhoods[0] = 0;
            for (vertex_t map_u = 0; map_u < m_n; ++map_u) {
                m_neighborhoods[map_u + 1] = m_neighborhoods[map_u] + sizes[map_u];
            }

            _t_real_neighborhood.stop();
            ScopedTimer _t_copy_edges("contraction", "CSRGraph", "copy_edges");

            m_edges_v.initialize(m_m);
            m_edges_w.initialize(m_m);

            for (vertex_t map_u = 0; map_u < mapping.get_coarse_n(); ++map_u) {
                size_t cursor = m_neighborhoods[map_u];
                for (vertex_t j = overest_neighborhood[map_u]; j < overest_neighborhood[map_u + 1]; ++j) {
                    vertex_t map_v = edges_v[j];
                    weight_t map_w = edges_w[j];

                    if (map_v != g.get_n()) {
                        m_edges_v[cursor] = map_v;
                        m_edges_w[cursor] = map_w;
                        cursor += 1;
                    }
                }
            }

            _t_copy_edges.stop();
        }

        // Move constructor
        CSRGraph(CSRGraph &&other) noexcept {
            m_n = other.m_n;
            m_m = other.m_m;

            m_graph_weight = other.m_graph_weight;

            std::swap(m_v_weights, other.m_v_weights);
            std::swap(m_neighborhoods, other.m_neighborhoods);
            std::swap(m_edges_v, other.m_edges_v);
            std::swap(m_edges_w, other.m_edges_w);
        }

        // Optionally disable copying.
        CSRGraph(const CSRGraph &) = delete;

        CSRGraph &operator=(const CSRGraph &) = delete;

        vertex_t get_n() const { return m_n; }

        vertex_t get_m() const { return m_m; }

        weight_t weight() const { return m_graph_weight; }

        weight_t weight(const vertex_t u) const { return m_v_weights[u]; }

        size_t size(const vertex_t u) const { return m_neighborhoods[u + 1] - m_neighborhoods[u]; }

        vertex_t neighbor(const vertex_t u, const size_t idx) const { return m_edges_v[m_neighborhoods[u] + idx]; }

        weight_t weight(const vertex_t u, const size_t idx) const { return m_edges_w[m_neighborhoods[u] + idx]; }

        void write_graph(std::string file_path) {
            std::ofstream file(file_path);

            file << m_n << " " << m_m / 2 << " 011" << std::endl;
            for (size_t i = 0; i < m_n; ++i) {
                file << m_v_weights[i] << " ";
                for (size_t j = m_neighborhoods[i]; j < m_neighborhoods[i + 1]; ++j) {
                    file << m_edges_v[j] + 1 << " " << m_edges_w[j] << " ";
                }
                file << std::endl;
            }
        }
    };
}

#endif //HEIPROMAP_CSR_GRAPH_H
