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

#include <omp.h>

#include "../utility/aligned_array.h"
#include "../definitions.h"
#include "../utility/matching.h"
#include "../utility/mapping.h"
#include "../utility/utils.h"
#include "../utility/profiler.h"

namespace HeiProMap {
    /**
    * Standard undirected Graph that can hold vertex and edge weights.
    */
    class CSRGraph {
    public:
        vertex_t n = 0;
        vertex_t m = 0;
        weight_t g_weight = 0;
        bool uniform_v_weights = true;
        bool uniform_e_weights = true;

        AlignedArray<weight_t> v_weights;
        AlignedArray<size_t> neighborhoods;
        AlignedArray<vertex_t> edges_v;
        AlignedArray<weight_t> edges_w;

        CSRGraph() = default;

        explicit CSRGraph(const std::string &file_path) {
            HEIPROMAP_PROFILE_SCOPE("io", "CSRGraph", "allocate");
            if (!file_exists(file_path)) {
                std::cerr << "File " << file_path << " does not exist!" << std::endl;
                exit(EXIT_FAILURE);
            }

            // mmap the whole file
            MMap mm = mmap_file_ro(file_path);
            char *p = mm.data;
            const char *end = mm.data + mm.size; // (tiny fix: don't do -1)

            HEIPROMAP_PROFILE_SCOPE("io", "CSRGraph", "read_header");

            // skip comment lines
            while (*p == '%') {
                while (*p != '\n') { ++p; }
                ++p;
            }

            // skip whitespace
            while (*p == ' ') { ++p; }

            // read number of vertices
            n = 0;
            while (*p != ' ' && *p != '\n') {
                n = n * 10 + (vertex_t) (*p - '0');
                ++p;
            }

            // skip whitespace
            while (*p == ' ') { ++p; }

            // read number of edges
            m = 0;
            while (*p != ' ' && *p != '\n') {
                m = m * 10 + (vertex_t) (*p - '0');
                ++p;
            }
            m *= 2;

            // search end of line or fmt
            std::string fmt = "000";
            bool has_v_weights = false;
            bool has_e_weights = false;
            while (*p == ' ') { ++p; }
            if (*p != '\n') {
                // found fmt
                std::string tmp_fmt;
                while (*p != ' ' && *p != '\n') {
                    tmp_fmt += *p;
                    ++p;
                }
                if (tmp_fmt.size() == 1) {
                    fmt[2] = tmp_fmt[0];
                } else if (tmp_fmt.size() == 2) {
                    fmt[1] = tmp_fmt[0];
                    fmt[2] = tmp_fmt[1];
                } else if (tmp_fmt.size() == 3) {
                    fmt[0] = tmp_fmt[0];
                    fmt[1] = tmp_fmt[1];
                    fmt[2] = tmp_fmt[2];
                }
                // skip whitespaces
                while (*p == ' ') { ++p; }
            }
            g_weight = 0;
            v_weights.initialize(n);
            neighborhoods.initialize(n + 1);
            neighborhoods[0] = 0;
            edges_v.initialize(m);
            edges_w.initialize(m);
            has_v_weights = fmt[1] == '1';
            has_e_weights = fmt[2] == '1';
            uniform_v_weights = !has_v_weights;
            uniform_e_weights = !has_e_weights;

            HEIPROMAP_PROFILE_SCOPE("io", "CSRGraph", "read_edges");

            ++p;
            if (has_v_weights && has_e_weights) {
                read_edges<true, true>(p, end);
            } else if (has_v_weights && !has_e_weights) {
                read_edges<true, false>(p, end);
            } else if (!has_v_weights && has_e_weights) {
                read_edges<false, true>(p, end);
            } else {
                read_edges<false, false>(p, end);
            }

            // done with the file
            munmap_file(mm);
        }

        explicit CSRGraph(vertex_t t_n, vertex_t t_m, weight_t t_g_weight) {
            n = t_n;
            m = t_m;
            g_weight = t_g_weight;
            uniform_v_weights = false;
            uniform_e_weights = false;

            v_weights.initialize(n, 0);
            neighborhoods.initialize(n + 1);
            neighborhoods[0] = 0;
            edges_v.initialize(m);
            edges_w.initialize(m);
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void initialize(const CSRGraph &g,
                        const Mapping &mapping) {
            AlignedArray<vertex_t> n_mapped;
            AlignedArray<vertex_t> n_mapped_prefix;
            AlignedArray<vertex_t> cursor;
            AlignedArray<vertex_t> mapped_vertices;

            struct SeenEntry {
                u32 epoch;
                size_t idx;
            };
            AlignedArray<SeenEntry> seen_idx; // one cache miss instead of two
            u32 epoch = 0;
            // allocate
            {
                HEIPROMAP_PROFILE_SCOPE("contraction", "CSRGraph", "allocate");

                n = mapping.get_coarse_n();
                g_weight = g.g_weight;
                uniform_v_weights = t_uniform_v_weights;
                uniform_e_weights = t_uniform_e_weights;
                v_weights.initialize(n, 0);
                neighborhoods.initialize(n + 1);
                neighborhoods[0] = 0;
                edges_v.initialize(g.m);
                edges_w.initialize(g.m);

                n_mapped.initialize(n, 0);
                n_mapped_prefix.initialize(n + 1);
                n_mapped_prefix[0] = 0;
                cursor.initialize(n + 1);
                mapped_vertices.initialize(g.n);

                seen_idx.initialize(n, {0, 0});
            }
            // count how many vertices are mapped to each new vertex
            {
                HEIPROMAP_PROFILE_SCOPE("contraction", "CSRGraph", "n_mapped");

                // count and collect weights
                for (vertex_t u = 0; u < mapping.get_old_n(); ++u) {
                    vertex_t map_u = mapping.get(u);
                    n_mapped[map_u] += 1;
                    if constexpr (t_uniform_v_weights) {
                        v_weights[map_u] += 1;
                    } else {
                        v_weights[map_u] += g.v_weights[u];
                    }
                }
            }
            // prefix sum on n_mapped
            {
                HEIPROMAP_PROFILE_SCOPE("contraction", "CSRGraph", "prefix_sum");

                for (vertex_t map_u = 0; map_u < n; ++map_u) {
                    n_mapped_prefix[map_u + 1] = n_mapped_prefix[map_u] + n_mapped[map_u];
                }
            }
            // copy so we have a cursor
            {
                HEIPROMAP_PROFILE_SCOPE("contraction", "CSRGraph", "copy_cursor");

                for (vertex_t map_u = 0; map_u <= n; ++map_u) {
                    cursor[map_u] = n_mapped_prefix[map_u];
                }
            }
            // insert mapped vertices
            {
                HEIPROMAP_PROFILE_SCOPE("contraction", "CSRGraph", "mapped_vertices");

                for (vertex_t u = 0; u < g.n; ++u) {
                    vertex_t map_u = mapping.get(u);
                    mapped_vertices[cursor[map_u]] = u;
                    cursor[map_u] += 1;
                }
            }
            // insert edges in real array
            {
                HEIPROMAP_PROFILE_SCOPE("contraction", "CSRGraph", "real_neighborhood");

                m = 0;
                for (vertex_t map_u = 0; map_u < n; ++map_u) {
                    epoch += 1;
                    for (u64 i = n_mapped_prefix[map_u]; i < n_mapped_prefix[map_u + 1]; ++i) {
                        vertex_t u = mapped_vertices[i];

                        if (i + 1 < n_mapped_prefix[map_u + 1])
                            __builtin_prefetch(&g.neighborhoods[mapped_vertices[i + 1]], 0, 1);


                        for (size_t j = g.neighborhoods[u]; j < g.neighborhoods[u + 1]; ++j) {
                            const vertex_t v = g.edges_v[j];

                            vertex_t map_v = mapping.get(v);
                            if (map_u == map_v) { continue; }
                            weight_t w = t_uniform_e_weights ? 1 : g.edges_w[j];

                            if (seen_idx[map_v].epoch == epoch) {
                                size_t k = seen_idx[map_v].idx;
                                edges_w[k] += w;
                            } else {
                                seen_idx[map_v].epoch = epoch;
                                seen_idx[map_v].idx = m;
                                edges_v[m] = map_v;
                                edges_w[m] = w;
                                m += 1;
                            }
                        }
                    }
                    neighborhoods[map_u + 1] = m;
                }
            }
        }

        template<bool t_uniform_v_weights, bool t_uniform_e_weights>
        void parallel_initialize(const CSRGraph &g,
                                 const Mapping &mapping,
                                 const u64 threads) {
            AlignedArray<vertex_t> overest_sizes;
            AlignedArray<vertex_t> overest_neighborhood;
            AlignedArray<vertex_t> m_per_thread;
            AlignedArray<vertex_t> temp_edges_v;
            AlignedArray<weight_t> temp_edges_w;
            AlignedArray<vertex_t> sizes;
            //
            {
                HEIPROMAP_PROFILE_SCOPE("contraction", "CSRGraph", "allocate");

                n = mapping.get_coarse_n();
                g_weight = g.g_weight;
                uniform_v_weights = t_uniform_v_weights;
                uniform_e_weights = t_uniform_e_weights;
                v_weights.initialize(n, 0);
                overest_sizes.initialize(n, 0);
                overest_neighborhood.initialize(n + 1);
            }
            //
            {
                HEIPROMAP_PROFILE_SCOPE("contraction", "CSRGraph", "overest_sizes");

                // overestimate neighborhood sizes and collect weights
                for (vertex_t u = 0; u < mapping.get_old_n(); ++u) {
                    vertex_t map_u = mapping.get(u);
                    overest_sizes[map_u] += g.deg(u);
                    v_weights[map_u] += t_uniform_v_weights ? 1 : g.v_weights[u];
                }
            }
            //
            {
                HEIPROMAP_PROFILE_SCOPE("contraction", "CSRGraph", "prefix_sum");
                overest_neighborhood[0] = 0;
                for (vertex_t map_u = 0; map_u < n; ++map_u) {
                    overest_neighborhood[map_u + 1] = overest_neighborhood[map_u] + overest_sizes[map_u];
                }
            }
            //
            {
                HEIPROMAP_PROFILE_SCOPE("contraction", "CSRGraph", "allocate_insert_edges");

                m_per_thread.initialize(threads, 0);
                temp_edges_v.initialize(overest_neighborhood[n], g.n);
                temp_edges_w.initialize(overest_neighborhood[n]);
                sizes.initialize(n, 0);
            }
            //
            {
                HEIPROMAP_PROFILE_SCOPE("contraction", "CSRGraph", "insert_edges");

                // insert edges in overestimated array
                m = 0;

                // Partition the *coarse* id space [0, m_n) into disjoint slices.
                auto slice_begin = [&](u64 t) -> vertex_t { return t * n / threads; };
                auto slice_end = [&](u64 t) -> vertex_t { return (t + 1) * n / threads; };

                #pragma omp parallel num_threads(threads)
                {
                    const u64 tid = (u64) omp_get_thread_num();
                    const vertex_t mu_beg = slice_begin(tid);
                    const vertex_t mu_end = tid == threads - 1 ? n : slice_end(tid);

                    size_t local_m = 0;
                    for (vertex_t u = 0; u < mapping.get_old_n(); ++u) {
                        vertex_t map_u = mapping.get(u);
                        if (map_u < mu_beg || map_u >= mu_end) continue; // not my bucket range

                        for (size_t i = g.neighborhoods[u]; i < g.neighborhoods[u + 1]; ++i) {
                            const vertex_t v = g.edges_v[i];
                            vertex_t map_v = mapping.get(v);
                            if (map_u == map_v) { continue; }
                            weight_t w = t_uniform_e_weights ? 1 : g.edges_w[i];

                            vertex_t beg = overest_neighborhood[map_u];
                            vertex_t end = overest_neighborhood[map_u + 1];
                            vertex_t len = end - beg;
                            if (len == 0) { continue; }

                            // insert map_v
                            vertex_t j = beg + (map_v % len);
                            while (true) {
                                if (j == end) { j = beg; }
                                if (temp_edges_v[j] == map_v) {
                                    temp_edges_w[j] += w;
                                    break;
                                }
                                if (temp_edges_v[j] == g.n) {
                                    temp_edges_v[j] = map_v;
                                    temp_edges_w[j] = w;
                                    local_m += 1;
                                    sizes[map_u] += 1;
                                    break;
                                }
                                j += 1;
                            }
                        }
                    }
                    m_per_thread[tid] = local_m;
                }

                m = 0;
                for (size_t i = 0; i < threads; ++i) { m += m_per_thread[i]; }
            }
            //
            {
                HEIPROMAP_PROFILE_SCOPE("contraction", "CSRGraph", "real_neighborhood");

                // insert edges in real array
                neighborhoods.initialize(n + 1);
                neighborhoods[0] = 0;
                for (vertex_t map_u = 0; map_u < n; ++map_u) {
                    neighborhoods[map_u + 1] = neighborhoods[map_u] + sizes[map_u];
                }
            }
            //
            {
                HEIPROMAP_PROFILE_SCOPE("contraction", "CSRGraph", "copy_edges");

                edges_v.initialize(m);
                edges_w.initialize(m);

                #pragma omp parallel for num_threads(threads)
                for (vertex_t map_u = 0; map_u < mapping.get_coarse_n(); ++map_u) {
                    size_t cursor = neighborhoods[map_u];
                    for (vertex_t j = overest_neighborhood[map_u]; j < overest_neighborhood[map_u + 1]; ++j) {
                        vertex_t map_v = temp_edges_v[j];
                        weight_t map_w = temp_edges_w[j];

                        if (map_v != g.n) {
                            edges_v[cursor] = map_v;
                            edges_w[cursor] = map_w;
                            cursor += 1;
                        }
                    }
                }
            }
        }

        // Move constructor
        CSRGraph(CSRGraph &&other) noexcept {
            n = other.n;
            m = other.m;
            g_weight = other.g_weight;
            uniform_v_weights = other.uniform_v_weights;
            uniform_e_weights = other.uniform_e_weights;

            v_weights = std::move(other.v_weights);
            neighborhoods = std::move(other.neighborhoods);
            edges_v = std::move(other.edges_v);
            edges_w = std::move(other.edges_w);
        }

        CSRGraph(const CSRGraph &other) {
            n = other.n;
            m = other.m;
            g_weight = other.g_weight;
            uniform_v_weights = other.uniform_v_weights;
            uniform_e_weights = other.uniform_e_weights;

            v_weights = other.v_weights;
            neighborhoods = other.neighborhoods;
            edges_v = other.edges_v;
            edges_w = other.edges_w;
        }

        CSRGraph &operator=(const CSRGraph &other) {
            if (this != &other) {
                n = other.n;
                m = other.m;
                g_weight = other.g_weight;
                uniform_v_weights = other.uniform_v_weights;
                uniform_e_weights = other.uniform_e_weights;

                v_weights = other.v_weights;
                neighborhoods = other.neighborhoods;
                edges_v = other.edges_v;
                edges_w = other.edges_w;
            }
            return *this;
        }

        size_t deg(const vertex_t u) const { return neighborhoods[u + 1] - neighborhoods[u]; }

        void write_graph(const std::string &file_path) const {
            std::ofstream file(file_path);

            file << n << " " << m / 2 << " 011" << std::endl;
            for (size_t i = 0; i < n; ++i) {
                file << v_weights[i] << " ";
                for (size_t j = neighborhoods[i]; j < neighborhoods[i + 1]; ++j) {
                    file << edges_v[j] + 1 << " " << edges_w[j] << " ";
                }
                file << std::endl;
            }
        }

        void print_degree_distribution() const {
            std::unordered_map<size_t, size_t> degree_count;

            size_t max_deg = 0;

            // compute degrees
            for (vertex_t u = 0; u < n; ++u) {
                size_t deg = neighborhoods[u + 1] - neighborhoods[u];
                degree_count[deg]++;
                if (deg > max_deg) {
                    max_deg = deg;
                }
            }

            // move to vector for sorted output
            std::vector<std::pair<size_t, size_t> > dist;
            dist.reserve(degree_count.size());

            for (auto &p: degree_count) {
                dist.emplace_back(p); // (degree, count)
            }

            std::sort(dist.begin(), dist.end());

            // print
            std::cout << "#degree count\n";
            for (auto &[deg, count]: dist) {
                std::cout << deg << " " << count << "\n";
            }
        }

    private:
        template<bool has_v_weights, bool has_e_weights>
        void read_edges(char *p, const char *end) {
            vertex_t u = 0;
            size_t curr_m = 0;
            while (p < end) {
                while (*p == '%') {
                    while (*p != '\n') { ++p; }
                    ++p;
                }

                while (*p == ' ') { ++p; }

                weight_t vw = 1;
                if constexpr (has_v_weights) {
                    vw = 0;
                    while (*p != ' ' && *p != '\n') {
                        vw = vw * 10 + (weight_t) (*p - '0');
                        ++p;
                    }
                    while (*p == ' ') { ++p; }
                    if (vw != 1) { uniform_v_weights = false; }
                }
                v_weights[u] = vw;
                g_weight += vw;

                while (*p != '\n' && p < end) {
                    vertex_t v = 0;
                    while (*p != ' ' && *p != '\n') {
                        v = v * 10 + (vertex_t) (*p - '0');
                        ++p;
                    }
                    while (*p == ' ') { ++p; }

                    weight_t w = 1;
                    if constexpr (has_e_weights) {
                        w = 0;
                        while (*p != ' ' && *p != '\n') {
                            w = w * 10 + (weight_t) (*p - '0');
                            ++p;
                        }
                        while (*p == ' ') { ++p; }
                        if (w != 1) { uniform_e_weights = false; }
                    }

                    edges_v[curr_m] = v - 1;
                    edges_w[curr_m] = w;
                    ++curr_m;
                }
                neighborhoods[u + 1] = curr_m;
                ++u;
                ++p;
            }

            if (curr_m != m) {
                std::cerr << "Number of expected edges " << m << " not equal to number edges " << curr_m << " found!\n";
                exit(EXIT_FAILURE);
            }
        }
    };
}

#endif //HEIPROMAP_CSR_GRAPH_H
