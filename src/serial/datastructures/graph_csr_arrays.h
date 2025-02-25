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

#ifndef HEIPROMAP_GRAPH_CSR_ARRAYS_H
#define HEIPROMAP_GRAPH_CSR_ARRAYS_H

#include <fcntl.h>
#include <iostream>
#include <unistd.h>
#include <vector>
#include <sys/mman.h>
#include <sys/stat.h>

#include "../../definitions.h"
#include "../interfaces/ISerialGraph.h"
#include "../utility/utils.h"

namespace HeiProMap {
    /**
    * Standard undirected Graph that can hold vertex and edge weights.
    */
    class GraphCSRArrays final : public ISerialGraph {
        vertex_t m_n = 0;
        vertex_t m_m = 0;

        weight_t m_vertex_weights = 0;
        weight_t m_edge_weights   = 0;

        weight_t* m_v_weights   = nullptr;
        size_t* m_neighborhoods = nullptr;
        EdgeVW* m_edges         = nullptr;
        size_t m_curr_m         = 0;

    public:
        explicit GraphCSRArrays(const std::string& graph_in) {
            // Open the file
            int fd = open(graph_in.c_str(), O_RDONLY);
            if (fd == -1) {
                std::cerr << "File " << graph_in << " does not exist!" << std::endl;
                exit(EXIT_FAILURE);
            }

            // Get the file size
            struct stat fileInfo{};
            if (fstat(fd, &fileInfo) == -1) {
                std::cerr << "File " << graph_in << " Could not get file size!" << std::endl;
                close(fd);
                exit(EXIT_FAILURE);
            }
            size_t file_size = fileInfo.st_size;

            // Memory-map the file
            char* file_arr = static_cast<char*>(mmap(nullptr, file_size, PROT_READ, MAP_PRIVATE, fd, 0));
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
            m_m *= 2;

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

            vertex_t m_n_64    = round_up_64(m_n + 1);
            vertex_t m_m_64    = round_up_64(m_m + 1);
            m_neighborhoods    = (size_t*)aligned_alloc(64, m_n_64 * sizeof(size_t));
            m_neighborhoods[0] = 0;
            m_edges            = (EdgeVW*)aligned_alloc(64, m_m_64 * sizeof(EdgeVW));

            vertex_t u = 0;
            if (fmt_1 == '0' && fmt_2 == '0') {
                m_v_weights = (weight_t*)aligned_alloc(64, m_n_64 * sizeof(weight_t));
                std::fill_n(m_v_weights, m_n_64, 1);

                m_vertex_weights = (weight_t)m_n;
                while (true) {
                    if (file_arr[i] == '%') {
                        // this line is a comment, ignore it
                        while (file_arr[i] != '\n') { ++i; }
                        ++i;
                        continue;
                    }
                    // this line contains vertex information

                    while (file_arr[i] == ' ') { ++i; }

                    while (file_arr[i] != '\n') {
                        // read in the edges
                        vertex_t v = 0;
                        while (file_arr[i] != ' ' && file_arr[i] != '\n') { v = v * 10 + (file_arr[i++] - '0'); }
                        while (file_arr[i] == ' ') { ++i; }

                        m_edges[m_curr_m++] = {v - 1, 1};
                    }

                    ++i;
                    m_neighborhoods[u + 1] = m_curr_m;
                    u += 1;

                    if (u + 32 >= m_n) {
                        break;
                    }
                }
                while (true) {
                    if (file_arr[i] == '%') {
                        // this line is a comment, ignore it
                        while (i < file_size && file_arr[i] != '\n') { ++i; }
                        // move_while_not(file_arr, i, '\n', file_size);
                        ++i;
                        continue;
                    }
                    // this line contains vertex information

                    while (i < file_size && file_arr[i] == ' ') { ++i; }

                    while (i < file_size && file_arr[i] != '\n') {
                        // read in the edges
                        vertex_t v = 0;
                        for (; i < file_size && file_arr[i] != ' ' && file_arr[i] != '\n'; ++i) { v = v * 10 + (file_arr[i] - '0'); }
                        while (i < file_size && file_arr[i] == ' ') { ++i; }

                        m_edges[m_curr_m++] = {v - 1, 1};
                    }

                    ++i;
                    m_neighborhoods[u + 1] = m_curr_m;
                    u += 1;

                    if (u == m_n) {
                        break;
                    }
                }
            } else {
                m_v_weights = (weight_t*)aligned_alloc(64, m_n_64 * sizeof(weight_t));
                std::fill_n(m_v_weights, m_n_64, 0);

                while (true) {
                    if (file_arr[i] == '%') {
                        // this line is a comment, ignore it
                        move_while_not(file_arr, i, '\n', file_size);
                        ++i;
                        continue;
                    }
                    // this line contains vertex information

                    move_while(file_arr, i, ' ', file_size); // skip leading whitespaces

                    weight_t u_w = 1;
                    if (fmt_1 == '1') {
                        // read in vertex weight
                        u_w = 0;
                        while (i < file_size && file_arr[i] != ' ' && file_arr[i] != '\n') { u_w = u_w * 10 + (file_arr[i++] - '0'); }
                        move_while(file_arr, i, ' ', file_size); // move to the next number
                    }
                    m_v_weights[u] = u_w;
                    m_vertex_weights += u_w;

                    while (i < file_size && file_arr[i] != '\n') {
                        // read in the edges
                        vertex_t v = 0;
                        while (i < file_size && file_arr[i] != ' ' && file_arr[i] != '\n') { v = v * 10 + (file_arr[i++] - '0'); }
                        move_while(file_arr, i, ' ', file_size); // move to the next number

                        weight_t w = 1;
                        if (fmt_2 == '1') {
                            w = 0;
                            while (i < file_size && file_arr[i] != ' ' && file_arr[i] != '\n') { w = w * 10 + (file_arr[i++] - '0'); }
                            move_while(file_arr, i, ' ', file_size); // move to the next number
                        }

                        m_edges[m_curr_m++] = {v - 1, w};
                    }

                    ++i;
                    m_neighborhoods[u + 1] = m_curr_m;
                    u += 1;

                    if (u == m_n) {
                        break;
                    }
                }
            }

            // Clean up
            munmap(file_arr, file_size);
            close(fd);
        }

        GraphCSRArrays(const GraphCSRArrays& g, const std::vector<EdgeUV>& matching) {
            vertex_t m_n_64 = round_up_64(g.get_n() + 1);
            vertex_t m_m_64 = round_up_64(g.get_m() + 1);

            m_n                = g.get_n();
            m_m                = 0;
            m_curr_m           = 0;
            m_v_weights        = (weight_t*)aligned_alloc(64, m_n_64 * sizeof(weight_t));
            m_neighborhoods    = (size_t*)aligned_alloc(64, m_n_64 * sizeof(size_t));
            m_neighborhoods[0] = 0;
            m_edges            = (EdgeVW*)aligned_alloc(64, m_m_64 * sizeof(size_t));
            m_vertex_weights   = g.m_vertex_weights;

            // define the state of each vertex
            constexpr u8 NOT_MATCHED    = 0;
            constexpr u8 FIRST_MATCHED  = 1;
            constexpr u8 SECOND_MATCHED = 2;
            std::vector<u8> vertex_state(m_n, NOT_MATCHED);
            std::vector<vertex_t> vertex_neighbor(m_n);

            // check the matching
            for (const auto& [u, v] : matching) {
                vertex_state[u]    = FIRST_MATCHED;
                vertex_state[v]    = SECOND_MATCHED;
                vertex_neighbor[u] = v;
                vertex_neighbor[v] = u;
            }

            for (vertex_t u = 0; u < m_n; ++u) {
                if (vertex_state[u] == NOT_MATCHED) {
                    // copy it to the next graph
                    m_v_weights[u] = g.get_weight(u);
                    for (size_t i = 0; i < g.size(u); ++i) {
                        vertex_t vv = g.neighbor(u, i);
                        weight_t ww = g.get_weight(u, i);

                        // if the vv vertex is matched, then make an edge to the neighbor vertex
                        vv                    = vertex_state[vv] == SECOND_MATCHED ? vertex_neighbor[vv] : vv;
                        m_edges[m_curr_m].v   = vv;
                        m_edges[m_curr_m++].w = ww;

                        // if the edge is present, then add the weight, else expand it
                        for (size_t j = m_neighborhoods[u]; j < m_curr_m - 1; ++j) {
                            if (m_edges[j].v == vv) {
                                m_edges[j].w += ww;
                                m_curr_m -= 1;
                                break;
                            }
                        }
                    }
                } else if (vertex_state[u] == FIRST_MATCHED) {
                    // the vertex gets all neighbors of v
                    vertex_t v = vertex_neighbor[u];

                    m_v_weights[v] = 0;
                    m_v_weights[u] = g.get_weight(u) + g.get_weight(v);

                    for (size_t i = 0; i < g.size(u); ++i) {
                        vertex_t vv = g.neighbor(u, i);
                        weight_t ww = g.get_weight(u, i);

                        // do not add edge to matched vertex
                        if (vv == v) { continue; }

                        // if the vv vertex is matched, then make an edge to the neighbor vertex
                        vv = vertex_state[vv] == SECOND_MATCHED ? vertex_neighbor[vv] : vv;

                        m_edges[m_curr_m].v   = vv;
                        m_edges[m_curr_m++].w = ww;

                        // if the edge is present, then add the weight, else expand it
                        for (size_t j = m_neighborhoods[u]; j < m_curr_m - 1; ++j) {
                            if (m_edges[j].v == vv) {
                                m_edges[j].w += ww;
                                m_curr_m -= 1;
                                break;
                            }
                        }
                    }
                    for (size_t i = 0; i < g.size(v); ++i) {
                        vertex_t vv = g.neighbor(v, i);
                        weight_t ww = g.get_weight(v, i);

                        // do not add edge to matched vertex
                        if (vv == u) { continue; }

                        // if the vv vertex is matched, then make an edge to the neighbor vertex
                        vv = vertex_state[vv] == SECOND_MATCHED ? vertex_neighbor[vv] : vv;

                        m_edges[m_curr_m].v   = vv;
                        m_edges[m_curr_m++].w = ww;

                        // if the edge is present, then add the weight, else expand it
                        for (size_t j = m_neighborhoods[u]; j < m_curr_m - 1; ++j) {
                            if (m_edges[j].v == vv) {
                                m_edges[j].w += ww;
                                m_curr_m -= 1;
                                break;
                            }
                        }
                    }
                }
                m_neighborhoods[u + 1] = m_curr_m;
            }
            m_m = m_curr_m;
        }

        // Move constructor
        GraphCSRArrays(GraphCSRArrays&& other) noexcept {
            m_n = other.m_n;
            m_m = other.m_m;

            m_vertex_weights = other.m_vertex_weights;
            m_edge_weights   = other.m_edge_weights;

            m_v_weights     = other.m_v_weights;
            m_neighborhoods = other.m_neighborhoods;
            m_edges         = other.m_edges;
            m_curr_m        = other.m_curr_m;

            other.m_v_weights     = nullptr;
            other.m_neighborhoods = nullptr;
            other.m_edges         = nullptr;
        }

        // Optionally disable copying.
        GraphCSRArrays(const GraphCSRArrays&)            = delete;
        GraphCSRArrays& operator=(const GraphCSRArrays&) = delete;

        ~GraphCSRArrays() override {
            free(m_v_weights);
            free(m_neighborhoods);
            free(m_edges);
        }

        vertex_t get_n() const override { return m_n; }
        vertex_t get_m() const override { return m_m; }
        weight_t get_weight() const override { return m_vertex_weights; }
        weight_t get_weight(const vertex_t u) const override { return m_v_weights[u]; }
        size_t size(const vertex_t u) const override { return m_neighborhoods[u + 1] - m_neighborhoods[u]; }
        vertex_t neighbor(const vertex_t u, const size_t idx) const override { return m_edges[m_neighborhoods[u] + idx].v; }
        weight_t get_weight(const vertex_t u, const size_t idx) const override { return m_edges[m_neighborhoods[u] + idx].w; }

        // edge manipulation
        bool edge_exists(const vertex_t u, const vertex_t v) const override {
            for (size_t i = m_neighborhoods[u]; i < m_neighborhoods[u + 1]; ++i) {
                if (m_edges[i].v == v) {
                    return true;
                }
            }
            return false;
        }

        class NeighborhoodIterator {
            vertex_t m_u;
            size_t* m_neighborhoods;
            EdgeVW* m_edges;

        public:
            NeighborhoodIterator(vertex_t u,
                                 size_t* neighborhoods,
                                 EdgeVW* edges) : m_u(u), m_neighborhoods(neighborhoods), m_edges(edges) {}

            class Iterator {
                EdgeVW* m_edges;
                size_t m_idx;

            public:
                // Constructor
                Iterator(EdgeVW* edges, size_t idx) {
                    m_edges = edges;
                    m_idx   = idx;
                }

                // Dereference operator
                EdgeVW operator*() const {
                    return m_edges[m_idx];
                }

                // Pre-increment operator
                Iterator& operator++() {
                    m_idx++;
                    return *this;
                }

                bool operator!=(const Iterator& other) const {
                    return m_idx != other.m_idx;
                }
            };

            Iterator begin() const { return {m_edges, m_neighborhoods[m_u]}; }
            Iterator end() const { return {m_edges, m_neighborhoods[m_u + 1]}; }
        };

        NeighborhoodIterator operator[](const vertex_t u) {
            return {u, m_neighborhoods, m_edges};
        }
    };
}

#endif //HEIPROMAP_GRAPH_CSR_ARRAYS_H
