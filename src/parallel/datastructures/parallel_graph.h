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

#include <iostream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <omp.h>
#include <cstring>

#include "../../interfaces/IGraph.h"
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
        ParallelGraph() = default;

        ~ParallelGraph() final {
#pragma omp parallel for schedule(static) default(none) num_threads(m_n_threads)
            for(size_t i = 0; i < m_n; ++i){
                std::vector<EdgeVW>().swap(m_adj[i]);
            }
        }

        // initialization
        void initialize(const std::string &graph_in, u64 n_threads) final {
            m_n_threads = n_threads;

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
#pragma omp parallel default(none) firstprivate(i, file_size, file_arr, fmt_1, fmt_2, n_threads) num_threads(n_threads)
            {
                weight_t local_g_weight = 0;
                std::vector<EdgeVW> edges;
                edges.reserve(100);

                size_t thread_id = omp_get_thread_num();
                vertex_t base_range = floor((f64) m_n / (f64) n_threads);
                vertex_t rem = m_n % n_threads;

                vertex_t start_u;
                vertex_t end_u;
                if (thread_id < rem) {
                    start_u = thread_id * (base_range + 1);
                    end_u = start_u + base_range + 1;
                } else {
                    start_u = rem * (base_range + 1) + (thread_id - rem) * base_range;
                    end_u = start_u + base_range;
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

                    if(u == m_n){
                        break;
                    }
                }

#pragma omp atomic
                m_g_weight += local_g_weight;
            }

            // Clean up
            munmap(file_arr, file_size);
            close(fd);
        }

        // graph properties
        vertex_t get_n() const final { return m_n; }

        vertex_t get_m() const final { return m_m; }

        weight_t get_weight() const final { return m_g_weight; }

        // vertex weights
        weight_t get_weight(vertex_t u) const final { return m_v_weights[u]; }

        size_t size(vertex_t u) const final { return m_adj[u].size(); }

        vertex_t neighbor(vertex_t u, size_t idx) const final { return m_adj[u][idx].v; }
        weight_t get_weight(vertex_t u, size_t idx) const final { return m_adj[u][idx].w; }

        // edge manipulation
        bool edge_exists(vertex_t u, vertex_t v) const final { return std::any_of(m_adj[u].begin(), m_adj[u].end(), [&](const EdgeVW &e) { return e.v == v; }); }

        // coarsing and uncoarsing
        void contract(vertex_t u, vertex_t v) final {
            // add weight of v to u
            m_v_weights[u] += m_v_weights[v];  // thread safe, only this thread modifies u

            // remove v from all its neighbors
            for (const EdgeVW &e: m_adj[v]) {
                remove_edge(e.v, v); // remove v from e.v  // not thread safe, multiple threads could remove from e.v
            }

            // connect neighbors of v to u, but not u
            for (const EdgeVW &e: m_adj[v]) {
                if (u != e.v) {
                    add_edge_with_weight(u, e.v, e.w);  // not thread safe, remove_edge(e.v, v) could modify u
                    add_edge_with_weight(e.v, u, e.w);  // not thread safe, remove_edge(e.v, v) could modify e.v
                }
            }
        }

        void uncontract(vertex_t u, vertex_t v) final {
            // remove neighbors of v from u
            for (const EdgeVW &e: m_adj[v]) {
                if (u != e.v) {
                    remove_edge_with_weight(u, e.v, e.w);
                    remove_edge_with_weight(e.v, u, e.w);
                }
            }

            // connect v to all its neighbors
            for (const EdgeVW &e: m_adj[v]) {
                add_edge_with_weight_guaranteed(e.v, v, e.w); // add v to e.v
            }

            // subtract weight of v from u
            m_v_weights[u] -= m_v_weights[v];
        }

    private:
        // coarsing
        void remove_edge(vertex_t u, vertex_t v) {
            size_t lower_idx = own_lower_bound_guaranteed(m_adj[u], v);
            m_adj[u].erase(m_adj[u].cbegin() + lower_idx);
        }

        void add_edge_with_weight(vertex_t u, vertex_t v, weight_t weight = 1) {
            size_t lower_idx = own_lower_bound_not_guaranteed(m_adj[u], v);
            if (lower_idx != m_adj[u].size() && m_adj[u][lower_idx].v == v) {
                m_adj[u][lower_idx].w += weight;
            } else {
                m_adj[u].insert(m_adj[u].begin() + lower_idx, {v, weight});
            }
        }

        // uncoarsing
        void add_edge_with_weight_guaranteed(vertex_t u, vertex_t v, weight_t weight = 1) {
            size_t lower_idx = own_lower_bound_not_guaranteed(m_adj[u], v);
            m_adj[u].insert(m_adj[u].begin() + lower_idx, {v, weight});
        }

        void remove_edge_with_weight(vertex_t u, vertex_t v, weight_t weight = 1) {
            size_t lower_idx = own_lower_bound_guaranteed(m_adj[u], v);
            m_adj[u][lower_idx].w -= weight;
            if (m_adj[u][lower_idx].w == 0) {
                m_adj[u].erase(m_adj[u].begin() + lower_idx);
            }
        }

    };

}

#endif //HEIPROMAP_PARALLEL_GRAPH_H
