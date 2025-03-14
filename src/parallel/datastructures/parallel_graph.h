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

#ifndef HEIPROMAP_PARALLEL_GRAPH_H
#define HEIPROMAP_PARALLEL_GRAPH_H

#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <omp.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include "../../commons/definitions.h"
#include "../../commons/macros.h"
#include "../interfaces/IParallelGraph.h"

namespace HeiProMap {

    /**
    * Standard undirected Graph that can hold vertex and edge weights.
    */
    class ParallelGraph final : public IParallelGraph {

    private:
        vertex_t m_n = 0; // original number of vertices
        vertex_t m_m = 0; // original number of edges

        u64 m_n_threads = 1;

        // vertex weight
        std::vector<weight_t> m_v_weights;

        weight_t m_g_weight = 0;

        // adjacency and edge weights
        std::vector<std::vector<EdgeVW>> m_adj;

    public:
        ParallelGraph(const std::string &graph_in, u64 t_threads) {
            m_n_threads = t_threads;

            // Open the file
            int fd = open(graph_in.c_str(), O_RDONLY);
            if (fd == -1) {
                std::cerr << "File " << graph_in << " does not exist!" << std::endl;
                exit(EXIT_FAILURE);
            }

            // Get the file size
            struct stat fileInfo;
            if (fstat(fd, &fileInfo) == -1) {
                std::cerr << "File " << graph_in << " Could not get file size!" << std::endl;
                close(fd);
                exit(EXIT_FAILURE);
            }
            size_t file_size = fileInfo.st_size;

            // Memory-map the file
            char *file_arr = static_cast<char *>(mmap(nullptr, file_size, PROT_READ, MAP_PRIVATE, fd, 0));
            if (file_arr == MAP_FAILED) {
                std::cerr << "File " << graph_in << " Could not map the file!" << std::endl;
                close(fd);
                exit(EXIT_FAILURE);
            }

            size_t i = 0;
            while (file_arr[i] == '%') {
                // skip line, since it is a comment
                move_while_not(file_arr, i, '\n', file_size);
                ++i;
            }

            // skip leading white spaces
            move_while(file_arr, i, ' ', file_size);

            // read n
            m_n = 0;
            while (file_arr[i] != ' ') { m_n = m_n * 10 + (file_arr[i++] - '0'); }

            // skip whitespaces
            move_while(file_arr, i, ' ', file_size);

            // read m
            m_m = 0;
            while (file_arr[i] != ' ' && file_arr[i] != '\n') { m_m = m_m * 10 + (file_arr[i++] - '0'); }

            // skip whitespaces
            move_while(file_arr, i, ' ', file_size);

            // read fmt, since its optional special code
            // char fmt_0 = '0';
            char fmt_1 = '0';
            char fmt_2 = '0';
            while (file_arr[i] != '\n') {
                if (file_arr[i] != ' ') {
                    // found one fmt number
                    fmt_2 = file_arr[i++];

                    if (i < file_size && file_arr[i] != ' ' && file_arr[i] != '\n') {
                        fmt_1 = fmt_2;
                        fmt_2 = file_arr[i++];

                        if (i < file_size && file_arr[i] != ' ' && file_arr[i] != '\n') {
                            // fmt_0 = fmt_1;
                            fmt_1 = fmt_2;
                            fmt_2 = file_arr[i++];
                        }
                    }
                    break;
                }
                i += 1;
            }
            // now only ' ' expected until '\n'
            move_while(file_arr, i, ' ', file_size);
            ++i; // now on the next line

            m_v_weights.resize(m_n);
            m_adj.resize(m_n);
#pragma omp parallel default(none) firstprivate(i, file_size, file_arr, fmt_1, fmt_2, m_n_threads) num_threads(m_n_threads)
            {
                weight_t            local_g_weight = 0;
                std::vector<EdgeVW> edges;
                edges.reserve(100);

                size_t   thread_id  = omp_get_thread_num();
                vertex_t base_range = floor((f64) m_n / (f64) m_n_threads);
                vertex_t rem        = m_n % m_n_threads;

                vertex_t start_u;
                vertex_t end_u;
                if (thread_id < rem) {
                    start_u = thread_id * (base_range + 1);
                    end_u   = start_u + base_range + 1;
                } else {
                    start_u = rem * (base_range + 1) + (thread_id - rem) * base_range;
                    end_u   = start_u + base_range;
                }

                vertex_t u = 0;
                while (true) {
                    if (file_arr[i] == '%') {
                        // this line is a comment, ignore it
                        move_while_not(file_arr, i, '\n', file_size);
                        ++i;
                        continue;
                    }

                    if (u < start_u) {
                        // this line should not be read by this thread
                        move_while_not(file_arr, i, '\n', file_size);
                        ++i;
                        u += 1;
                        continue;
                    }

                    if (u >= end_u) {
                        // this thread has read everything
                        break;
                    }

                    // this line contains vertex information
                    move_while(file_arr, i, ' ', file_size); // skip leading whitespaces

                    weight_t u_w = 1;
                    if (fmt_1 == '1') {
                        u_w = 0;
                        // read in vertex weight
                        while (i < file_size && file_arr[i] != ' ' && file_arr[i] != '\n') { u_w = u_w * 10 + (file_arr[i++] - '0'); }
                        move_while(file_arr, i, ' ', file_size); // move to next number
                    }
                    m_v_weights[u] = u_w;
                    local_g_weight += u_w;

                    while (i < file_size && file_arr[i] != '\n') {
                        // read in the edges
                        vertex_t v = 0;
                        while (i < file_size && file_arr[i] != ' ' && file_arr[i] != '\n') { v = v * 10 + (file_arr[i++] - '0'); }
                        move_while(file_arr, i, ' ', file_size); // move to next number

                        weight_t w = 1;
                        if (fmt_2 == '1') {
                            w = 0;
                            while (i < file_size && file_arr[i] != ' ' && file_arr[i] != '\n') { w = w * 10 + (file_arr[i++] - '0'); }
                            move_while(file_arr, i, ' ', file_size); // move to next number
                        }

                        edges.emplace_back(v - 1, w);
                    }
                    m_adj[u].reserve(edges.size());
                    for (auto &e: edges) {
                        m_adj[u].emplace_back(e.v, e.w);
                    }
                    std::sort(m_adj[u].begin(), m_adj[u].end());
                    edges.clear();

                    ++i;
                    u += 1;

                    if (u == m_n) {
                        break;
                    }
                }

#pragma omp atomic
                m_g_weight += local_g_weight;
            }

            // Clean up
            munmap(file_arr, file_size);
            close(fd);
        };

        ParallelGraph(const ParallelGraph &g,
                      const EdgeUV *matches,
                      size_t &matches_size,
                      u64 threads) {
            matches = ASSUME_ALIGNED(EdgeUV *, matches, 64);

            vertex_t m_n_64 = round_up_64(g.get_n() + 1);
            vertex_t m_m_64 = round_up_64(g.get_m() + 1);

            m_n = g.get_n();
            m_m = 0;

            m_g_weight = g.m_g_weight;
            m_v_weights.resize(m_n);
            std::copy(g.m_v_weights.begin(), g.m_v_weights.end(), m_v_weights.begin());

            m_adj.resize(m_n);

            // define the state of each vertex
            constexpr u8 NOT_MATCHED    = 0;
            constexpr u8 FIRST_MATCHED  = 1;
            constexpr u8 SECOND_MATCHED = 2;

            u8 *vertex_state = (u8 *) aligned_alloc(64, m_n_64 * sizeof(u8));
            std::fill_n(vertex_state, m_n_64, NOT_MATCHED);

            vertex_t *vertex_neighbor = (vertex_t *) aligned_alloc(64, m_n_64 * sizeof(vertex_t));

            // check the matching
            for (size_t i = 0; i < matches_size; ++i) {
                const auto [u, v] = matches[i];

                vertex_state[u]    = FIRST_MATCHED;
                vertex_state[v]    = SECOND_MATCHED;
                vertex_neighbor[u] = v;
                vertex_neighbor[v] = u;

                m_v_weights[v] = 0;
                m_v_weights[u] = g.get_weight(u) + g.get_weight(v);
            }

            for (vertex_t u = 0; u < m_n; ++u) {
                if (vertex_state[u] == NOT_MATCHED) {
                    // copy it to the next graph
                    for (size_t i = 0; i < g.size(u); ++i) {
                        vertex_t vv = g.neighbor(u, i);
                        weight_t ww = g.get_weight(u, i);

                        // if the vv vertex is matched, then make an edge to the neighbor vertex
                        vv = vertex_state[vv] == SECOND_MATCHED ? vertex_neighbor[vv] : vv;

                        // if the edge is present, then add the weight, else expand it
                        bool found = false;
                        for (auto &[v, w]: m_adj[u]) {
                            if (v == vv) {
                                w += ww;
                                found = true;
                            }
                        }
                        if (!found) {
                            m_adj[u].emplace_back(vv, ww);
                        }
                    }
                } else if (vertex_state[u] == FIRST_MATCHED) {
                    // the vertex gets all neighbors of u and v
                    vertex_t v = vertex_neighbor[u];

                    for (size_t i = 0; i < g.size(u); ++i) {
                        vertex_t vv = g.neighbor(u, i);
                        weight_t ww = g.get_weight(u, i);
                        // do not add edge to matched vertex
                        if (vv == v) { continue; }

                        // if the vv vertex is matched, then make an edge to the neighbor vertex
                        vv = vertex_state[vv] == SECOND_MATCHED ? vertex_neighbor[vv] : vv;

                        // if the edge is present, then add the weight, else expand it
                        bool found = false;
                        for (auto &[vvv, www]: m_adj[u]) {
                            if (vvv == vv) {
                                www += ww;
                                found = true;
                            }
                        }
                        if (!found) {
                            m_adj[u].emplace_back(vv, ww);
                        }
                    }
                    for (size_t i = 0; i < g.size(v); ++i) {
                        vertex_t vv = g.neighbor(v, i);
                        weight_t ww = g.get_weight(v, i);
                        // do not add edge to matched vertex
                        if (vv == u) { continue; }

                        // if the vv vertex is matched, then make an edge to the neighbor vertex
                        vv = vertex_state[vv] == SECOND_MATCHED ? vertex_neighbor[vv] : vv;

                        // if the edge is present, then add the weight, else expand it
                        bool found = false;
                        for (auto &[vvv, www]: m_adj[u]) {
                            if (vvv == vv) {
                                www += ww;
                                found = true;
                            }
                        }
                        if (!found) {
                            m_adj[u].emplace_back(vv, ww);
                        }
                    }
                }
            }

            free(vertex_state);
            free(vertex_neighbor);

            m_m = 0;
            for (auto &vec: m_adj) { m_m += vec.size(); }

        }

        vertex_t get_n() const override { return m_n; }

        vertex_t get_m() const override { return m_m; }

        weight_t get_weight() const override { return m_g_weight; }

        weight_t get_weight(vertex_t u) const override { return m_v_weights[u]; }

        size_t size(vertex_t u) const override { return m_adj[u].size(); }

        vertex_t neighbor(vertex_t u, size_t idx) const override { return m_adj[u][idx].v; }

        weight_t get_weight(vertex_t u, size_t idx) const override { return m_adj[u][idx].w; }

        bool edge_exists(vertex_t u, vertex_t v) const override { return std::any_of(m_adj[u].begin(), m_adj[u].end(), [&](const EdgeVW &e) { return e.v == v; }); }
    };

}

#endif //HEIPROMAP_PARALLEL_GRAPH_H
