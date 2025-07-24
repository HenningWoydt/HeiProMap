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

#ifndef HEIPROMAP_DEEP_CSR_GRAPH_H
#define HEIPROMAP_DEEP_CSR_GRAPH_H

#include <omp.h>
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include "../../commons/aligned_array.h"
#include "../../commons/definitions.h"
#include "../../commons/matching.h"
#include "../../commons/utils.h"

namespace HeiProMap {
    /**
    * Standard undirected Graph that can hold vertex and edge weights.
    */
    class DeepCSRGraph {
        vertex_t m_n = 0;
        vertex_t m_m = 0;

        weight_t m_graph_weight = 0;

        AlignedArray<weight_t> m_v_weights;
        AlignedArray<size_t> m_neighborhoods;
        AlignedArray<vertex_t> m_edges_v;
        AlignedArray<weight_t> m_edges_w;

    public:
        DeepCSRGraph() = default;

        explicit DeepCSRGraph(const std::string &graph_in) {
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

            m_neighborhoods.initialize(m_n + 1);
            m_neighborhoods[0] = 0;
            m_edges_v.initialize(m_m + 1);
            m_edges_w.initialize(m_m + 1);

            size_t curr_m = 0;
            vertex_t u = 0;
            if (fmt_1 == '0' && fmt_2 == '0') {
                m_v_weights.initialize(m_n + 1, 1);

                m_graph_weight = (weight_t) m_n;
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
                        for (; i < file_size && file_arr[i] != ' ' && file_arr[i] != '\n'; ++i) {
                            v = v * 10 + (file_arr[i] - '0');
                        }
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
                m_v_weights.initialize(m_n + 1, 0);

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
                        while (i < file_size && file_arr[i] != ' ' && file_arr[i] != '\n') {
                            u_w = u_w * 10 + (file_arr[i++] - '0');
                        }
                        move_while(file_arr, i, ' ', file_size); // move to the next number
                    }
                    m_v_weights[u] = u_w;
                    m_graph_weight += u_w;

                    while (i < file_size && file_arr[i] != '\n') {
                        // read in the edges
                        vertex_t v = 0;
                        while (i < file_size && file_arr[i] != ' ' && file_arr[i] != '\n') {
                            v = v * 10 + (file_arr[i++] - '0');
                        }
                        move_while(file_arr, i, ' ', file_size); // move to the next number

                        weight_t w = 1;
                        if (fmt_2 == '1') {
                            w = 0;
                            while (i < file_size && file_arr[i] != ' ' && file_arr[i] != '\n') {
                                w = w * 10 + (file_arr[i++] - '0');
                            }
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

        void initialize(const DeepCSRGraph &g,
                        Matching &matching) {
            matching.set_translation();

            m_n = matching.get_n_coarse_nodes();
            const vertex_t temp_m = g.get_m() + 1;

            m_graph_weight = g.m_graph_weight;
            m_v_weights.initialize(m_n + 1);

            m_neighborhoods.initialize(m_n + 1);
            m_edges_v.initialize(temp_m);
            m_edges_w.initialize(temp_m);

            struct IdxMark {
                vertex_t idx;
                u32 mark;
            };

            AlignedArray<IdxMark> idx_mark;
            idx_mark.initialize(m_n, {0, 0});
            u32 mark = 0;

            size_t curr_m = 0;
            m_neighborhoods[0] = 0;
            for (vertex_t old_u = 0; old_u < g.get_n(); ++old_u) {
                vertex_t old_v = matching.get_partner(old_u);

                if (old_u > old_v) { continue; }
                vertex_t new_u = matching.get_n(old_u);

                mark += 1;

                for (size_t i = 0; i < g.size(old_u); ++i) {
                    vertex_t vv = g.neighbor(old_u, i);
                    weight_t ww = g.weight(old_u, i);
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
                        idx_mark[vv].idx = curr_m;
                        idx_mark[vv].mark = mark;
                        m_edges_v[curr_m] = vv;
                        m_edges_w[curr_m] = ww;
                        curr_m += 1;
                    }
                }

                if (old_u < old_v) {
                    for (size_t i = 0; i < g.size(old_v); ++i) {
                        vertex_t vv = g.neighbor(old_v, i);
                        vertex_t vv_partner = matching.get_partner(vv);
                        weight_t ww = g.weight(old_v, i);
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
                            idx_mark[vv].idx = curr_m;
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
                    vertex_t new_u = matching.get_n(old_u);
                    m_v_weights[new_u] = g.weight(old_u);
                } else if (old_u < matching.get_partner(old_u)) {
                    vertex_t new_u = matching.get_n(old_u);
                    vertex_t old_v = matching.get_partner(old_u);
                    m_v_weights[new_u] = g.weight(old_u) + g.weight(old_v);
                }
            }

            m_m = curr_m;
        }

        void parallel_initialize(const DeepCSRGraph &g,
                                 Matching &matching,
                                 u64 threads) {
            matching.set_translation();

            m_n = matching.get_n_coarse_nodes();
            const vertex_t temp_m = g.get_m() + 1;

            m_graph_weight = g.m_graph_weight;
            m_v_weights.initialize(m_n + 1);

            m_neighborhoods.initialize(m_n + 1);
            m_edges_v.initialize(temp_m);
            m_edges_w.initialize(temp_m);

            struct thread_info {
                vertex_t start_vertex = 0;
                vertex_t n_assigned_vertices = 0;
                vertex_t n_actual_vertices = 0;
                vertex_t curr_m = 0;
                std::vector<size_t> neighborhood;
                std::vector<vertex_t> edges_v;
                std::vector<weight_t> edges_w;
            };
            std::vector<thread_info> thread_infos(threads);

            double alpha = 2.0; // example scaling parameter > 1.0

            // Step 1: Compute weights
            std::vector<double> weights(threads);
            double total_weight = 0.0;
            for (size_t i = 0; i < threads; ++i) {
                weights[i] = std::pow(static_cast<double>(i + 1), alpha);
                total_weight += weights[i];
            }

            vertex_t n_total_vertices = g.get_n();
            vertex_t base = n_total_vertices / threads;
            vertex_t rem = n_total_vertices % threads;

            // Step 2: Assign number of vertices per thread
            vertex_t current_start = 0;
            vertex_t assigned_total = 0;
            for (size_t i = 0; i < threads; ++i) {
                double fraction = weights[i] / total_weight;
                vertex_t n_assign = static_cast<vertex_t>(std::round(fraction * n_total_vertices));

                // Ensure last thread takes any rounding residual
                if (i == threads - 1) {
                    n_assign = n_total_vertices - assigned_total;
                }

                thread_infos[i].start_vertex = current_start;
                thread_infos[i].n_assigned_vertices = n_assign;
                current_start += n_assign;
                assigned_total += n_assign;
            }

            // each thread determines the neighborhood one their assigned vertices
#pragma omp parallel num_threads(threads)
            {
                u64 t_id = omp_get_thread_num();
                thread_infos[t_id].neighborhood.push_back(0);
                for (vertex_t old_u = thread_infos[t_id].start_vertex;
                     old_u < thread_infos[t_id].start_vertex + thread_infos[t_id].n_assigned_vertices; ++old_u) {
                    vertex_t old_v = matching.get_partner(old_u);

                    if (old_u > old_v) { continue; }

                    for (size_t i = 0; i < g.size(old_u); ++i) {
                        vertex_t vv = g.neighbor(old_u, i);
                        weight_t ww = g.weight(old_u, i);
                        vertex_t vv_partner = matching.get_partner(vv);
                        if (vv == old_v) { continue; }

                        // if the vv vertex is matched, then make an edge to the neighbor vertex
                        vv = std::min(vv, vv_partner);

                        // map to the new node range
                        vv = matching.get_n(vv);

                        bool found = false;
                        for (size_t j = thread_infos[t_id].neighborhood.back(); j < thread_infos[t_id].curr_m; ++j) {
                            if (thread_infos[t_id].edges_v[j] == vv) {
                                thread_infos[t_id].edges_w[j] += ww;
                                found = true;
                                break;
                            }
                        }
                        if (!found) {
                            thread_infos[t_id].edges_v.push_back(vv);
                            thread_infos[t_id].edges_w.push_back(ww);
                            thread_infos[t_id].curr_m += 1;
                        }
                    }

                    if (old_u < old_v) {
                        for (size_t i = 0; i < g.size(old_v); ++i) {
                            vertex_t vv = g.neighbor(old_v, i);
                            vertex_t vv_partner = matching.get_partner(vv);
                            weight_t ww = g.weight(old_v, i);
                            // do not add edge to matched vertex
                            if (vv == old_u) { continue; }

                            // if the vv vertex is matched, then make an edge to the neighbor vertex
                            vv = std::min(vv, vv_partner);

                            // map to the new node range
                            vv = matching.get_n(vv);

                            bool found = false;
                            for (size_t j = thread_infos[t_id].neighborhood.back();
                                 j < thread_infos[t_id].curr_m; ++j) {
                                if (thread_infos[t_id].edges_v[j] == vv) {
                                    thread_infos[t_id].edges_w[j] += ww;
                                    found = true;
                                    break;
                                }
                            }
                            if (!found) {
                                thread_infos[t_id].edges_v.push_back(vv);
                                thread_infos[t_id].edges_w.push_back(ww);
                                thread_infos[t_id].curr_m += 1;
                            }
                        }
                    }
                    thread_infos[t_id].neighborhood.push_back(thread_infos[t_id].curr_m);
                }
                thread_infos[t_id].n_actual_vertices = thread_infos[t_id].neighborhood.size() - 1;
            }

            // each thread copies its data to the correct place in the real neighborhood
            m_neighborhoods[0] = 0;
#pragma omp parallel num_threads(threads)
            {
                u64 t_id = omp_get_thread_num();

                // determine how many vertices and edges come before
                vertex_t previous_m = 0;
                vertex_t previous_n = 0;
                for (size_t i = 0; i < t_id; ++i) {
                    previous_m += thread_infos[i].curr_m;
                    previous_n += thread_infos[i].n_actual_vertices;
                }

                // copy neighborhood sizes
                for (size_t i = 0; i < thread_infos[t_id].n_actual_vertices; ++i) {
                    m_neighborhoods[previous_n + i + 1] = thread_infos[t_id].neighborhood[i + 1] + previous_m;
                }

                // copy all edges
                for (size_t i = 0; i < thread_infos[t_id].curr_m; ++i) {
                    m_edges_v[previous_m + i] = thread_infos[t_id].edges_v[i];
                    m_edges_w[previous_m + i] = thread_infos[t_id].edges_w[i];
                }

                if (t_id == threads - 1) {
                    m_m = thread_infos[t_id].curr_m + previous_m;
                }

            }

// #pragma omp parallel for num_threads(threads)
            for (vertex_t old_u = 0; old_u < g.get_n(); ++old_u) {
                vertex_t old_v = matching.get_partner(old_u);
                if (old_u == old_v) {
                    vertex_t new_u = matching.get_n(old_u);
                    m_v_weights[new_u] = g.weight(old_u);
                } else if (old_u < old_v) {
                    vertex_t new_u = matching.get_n(old_u);
                    m_v_weights[new_u] = g.weight(old_u) + g.weight(old_v);
                }
            }
        }

        // Move constructor
        DeepCSRGraph(DeepCSRGraph &&other) noexcept {
            m_n = other.m_n;
            m_m = other.m_m;

            m_graph_weight = other.m_graph_weight;

            std::swap(m_v_weights, other.m_v_weights);
            std::swap(m_neighborhoods, other.m_neighborhoods);
            std::swap(m_edges_v, other.m_edges_v);
            std::swap(m_edges_w, other.m_edges_w);
        }

        // Optionally disable copying.
        DeepCSRGraph(const DeepCSRGraph &) = delete;

        DeepCSRGraph &operator=(const DeepCSRGraph &) = delete;

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

#endif //HEIPROMAP_DEEP_CSR_GRAPH_H
