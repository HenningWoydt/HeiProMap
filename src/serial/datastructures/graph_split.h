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

#ifndef HEIPROMAP_GRAPH_SPLIT_H
#define HEIPROMAP_GRAPH_SPLIT_H

#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include "../../commons/definitions.h"
#include "../../commons/utils.h"
#include "../coarsening/matching.h"
#include "../interfaces/ISerialGraph.h"

namespace HeiProMap {
    /**
    * Standard undirected Graph that can hold vertex and edge weights.
    */
    class GraphSplit final : public ISerialGraph {
        vertex_t m_n = 0;
        vertex_t m_m = 0;

        weight_t m_graph_weight = 0;

        weight_t* m_v_weights   = nullptr;
        size_t* m_neighborhoods = nullptr;
        vertex_t* m_edges_v     = nullptr;
        weight_t* m_edges_w     = nullptr;

    public:
        GraphSplit() = default;

        explicit GraphSplit(const std::string& graph_in) {
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
            m_edges_v          = (vertex_t*)aligned_alloc(64, m_m_64 * sizeof(vertex_t));
            m_edges_w          = (weight_t*)aligned_alloc(64, m_m_64 * sizeof(weight_t));

            size_t curr_m = 0;
            vertex_t u    = 0;
            if (fmt_1 == '0' && fmt_2 == '0') {
                m_v_weights = (weight_t*)aligned_alloc(64, m_n_64 * sizeof(weight_t));
                std::fill_n(m_v_weights, m_n_64, 1);

                m_graph_weight = (weight_t)m_n;
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

                        m_edges_v[curr_m] = v - 1;
                        m_edges_w[curr_m] = 1;
                        curr_m += 1;
                    }

                    ++i;
                    m_neighborhoods[u + 1] = curr_m;
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

                        m_edges_v[curr_m] = v - 1;
                        m_edges_w[curr_m] = 1;
                        curr_m += 1;
                    }

                    ++i;
                    m_neighborhoods[u + 1] = curr_m;
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
                    m_graph_weight += u_w;

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

                        m_edges_v[curr_m] = v - 1;
                        m_edges_w[curr_m] = w;
                        curr_m += 1;
                    }

                    ++i;
                    m_neighborhoods[u + 1] = curr_m;
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

        void initialize(const GraphSplit& g,
                        Matching& matching) {
            matching.set_translation();

            m_n = matching.get_n_coarse_nodes();

            const vertex_t m_n_64 = round_up_64(m_n + 1);
            const vertex_t m_m_64 = round_up_64(g.get_m() + 1);

            m_graph_weight = g.m_graph_weight;
            m_v_weights    = (weight_t*)aligned_alloc(64, m_n_64 * sizeof(weight_t));

            m_neighborhoods = (size_t*)aligned_alloc(64, m_n_64 * sizeof(size_t));
            m_edges_v       = (vertex_t*)aligned_alloc(64, m_m_64 * sizeof(vertex_t));
            m_edges_w       = (weight_t*)aligned_alloc(64, m_m_64 * sizeof(weight_t));

            struct IdxMark {
                vertex_t idx;
                u32 mark;
            };

            IdxMark* idx_mark = (IdxMark*)aligned_alloc(64, m_n_64 * sizeof(IdxMark));
            memset(idx_mark, 0, m_n_64 * sizeof(IdxMark));
            u32 mark = 0;

            size_t curr_m      = 0;
            m_neighborhoods[0] = 0;
            for (vertex_t old_u = 0; old_u < g.get_n(); ++old_u) {
                vertex_t old_v = matching.get_partner(old_u);

                if (old_u > old_v) { continue; }
                vertex_t new_u = matching.get_n(old_u);

                mark += 1;

                for (size_t i = 0; i < g.size(old_u); ++i) {
                    vertex_t vv         = g.neighbor(old_u, i);
                    weight_t ww         = g.weight(old_u, i);
                    vertex_t vv_partner = matching.get_partner(vv);
                    if (vv == old_v) { continue; }

                    // if the vv vertex is matched, then make an edge to the neighbor vertex
                    vv = std::min(vv, vv_partner);

                    // map to the new node range
                    vv = matching.get_n(vv);

                    if (idx_mark[vv].mark == mark) {
                        size_t idx = idx_mark[vv].idx;
                        m_edges_w[idx] += ww;
                    } else {
                        idx_mark[vv].idx  = curr_m;
                        idx_mark[vv].mark = mark;
                        m_edges_v[curr_m] = vv;
                        m_edges_w[curr_m] = ww;
                        curr_m += 1;
                    }
                }

                if (old_u < old_v) {
                    for (size_t i = 0; i < g.size(old_v); ++i) {
                        vertex_t vv         = g.neighbor(old_v, i);
                        vertex_t vv_partner = matching.get_partner(vv);
                        weight_t ww         = g.weight(old_v, i);
                        // do not add edge to matched vertex
                        if (vv == old_u) { continue; }

                        // if the vv vertex is matched, then make an edge to the neighbor vertex
                        vv = std::min(vv, vv_partner);

                        // map to the new node range
                        vv = matching.get_n(vv);

                        if (idx_mark[vv].mark == mark) {
                            size_t idx = idx_mark[vv].idx;
                            m_edges_w[idx] += ww;
                        } else {
                            idx_mark[vv].idx  = curr_m;
                            idx_mark[vv].mark = mark;
                            m_edges_v[curr_m] = vv;
                            m_edges_w[curr_m] = ww;
                            curr_m += 1;
                        }
                    }
                }
                m_neighborhoods[new_u + 1] = curr_m;
            }

            for (vertex_t old_u = 0; old_u < g.get_n(); ++old_u) {
                if (old_u == matching.get_partner(old_u)) {
                    vertex_t new_u     = matching.get_n(old_u);
                    m_v_weights[new_u] = g.weight(old_u);
                } else if (old_u < matching.get_partner(old_u)) {
                    vertex_t new_u     = matching.get_n(old_u);
                    vertex_t old_v     = matching.get_partner(old_u);
                    m_v_weights[new_u] = g.weight(old_u) + g.weight(old_v);
                }
            }

            free(idx_mark);

            m_m = curr_m;
        }

        // Move constructor
        GraphSplit(GraphSplit&& other) noexcept {
            m_n = other.m_n;
            m_m = other.m_m;

            m_graph_weight = other.m_graph_weight;

            m_v_weights     = other.m_v_weights;
            m_neighborhoods = other.m_neighborhoods;
            m_edges_v       = other.m_edges_v;
            m_edges_w       = other.m_edges_w;

            other.m_v_weights     = nullptr;
            other.m_neighborhoods = nullptr;
            other.m_edges_v       = nullptr;
            other.m_edges_w       = nullptr;
        }

        // Optionally disable copying.
        GraphSplit(const GraphSplit&) = delete;

        GraphSplit& operator=(const GraphSplit&) = delete;

        ~GraphSplit() override {
            free(m_v_weights);
            free(m_neighborhoods);
            free(m_edges_v);
            free(m_edges_w);
        }

        vertex_t get_n() const override { return m_n; }

        vertex_t get_m() const override { return m_m; }

        weight_t weight() const override { return m_graph_weight; }

        weight_t weight(const vertex_t u) const override { return m_v_weights[u]; }

        size_t size(const vertex_t u) const override { return m_neighborhoods[u + 1] - m_neighborhoods[u]; }

        vertex_t neighbor(const vertex_t u, const size_t idx) const override { return m_edges_v[m_neighborhoods[u] + idx]; }
        weight_t weight(const vertex_t u, const size_t idx) const override { return m_edges_w[m_neighborhoods[u] + idx]; }
    };
}

#endif //HEIPROMAP_GRAPH_SPLIT_H
