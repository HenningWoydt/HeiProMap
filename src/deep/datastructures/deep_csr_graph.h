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

        void parallel_initialize_slower(const DeepCSRGraph &g,
                                        Matching &matching,
                                        u64 threads) {
            matching.set_translation();

            m_n = matching.get_n_coarse_nodes();
            m_m = g.get_m();

            m_graph_weight = g.m_graph_weight;
            m_v_weights.initialize(m_n);

            m_neighborhoods.initialize(m_n + 1);

            // determine maximum neighborhood size
            std::vector<size_t> max_neighborhood_size(m_n);
#pragma omp parallel num_threads(threads) default(none) shared(g, matching, max_neighborhood_size)
            for (vertex_t u = 0; u < g.m_n; ++u) {
                vertex_t u_partner = matching.get_partner(u);

                if (u > u_partner) { continue; }

                vertex_t new_u = matching.get_n(u);

                u64 size_u = g.size(u);
                u64 size_v = (u == u_partner) ? 0 : g.size(u_partner);
                max_neighborhood_size[new_u] = size_u + size_v;

                weight_t u_weight = g.weight(u);
                weight_t v_weight = (u == u_partner) ? 0 : g.weight(u_partner);
                m_v_weights[new_u] = u_weight + v_weight;
            }

            // determine the offsets into the hashing array
            std::vector<size_t> hash_offset(m_n + 1);
            hash_offset[0] = 0;
            for (size_t i = 0; i < m_n; ++i) {
                hash_offset[i + 1] = hash_offset[i] + max_neighborhood_size[i];
            }

            std::vector<u64> real_neighborhood_size(m_n);
            std::vector<u64> hash_keys(m_m, m_n);
            std::vector<weight_t> hash_vals(m_m);

#pragma omp parallel for num_threads(threads) default(none) shared(g, matching, hash_offset, hash_keys, hash_vals, real_neighborhood_size) schedule(dynamic)
            forall_gu(g, u)
                {
                    vertex_t u_partner = matching.get_partner(u);
                    if (u > u_partner) { continue; } // another iteration will handle it
                    if (u == u_partner) {
                        // vertex was not matched, only consider its neighbors
                        vertex_t u_new = matching.get_n(u);

                        // we need to insert v_new into the hash table of u_new
                        u64 offset = hash_offset[u_new];
                        u64 hash_table_size = hash_offset[u_new + 1] - hash_offset[u_new];
                        u64 n_neighbors = 0;
                        forall_guivw(g, u, i, v, w)
                            {
                                vertex_t v_new = matching.get_n(v);
                                u64 key = (v_new * 11400714819323198485ull) % hash_table_size;

                                // linear probing to find empty position or the key
                                for (u64 j = 0; j < hash_table_size; ++j) {
                                    u64 idx = offset + ((key + j) % hash_table_size);

                                    if (hash_keys[idx] == m_n) {
                                        hash_keys[idx] = v_new;
                                        hash_vals[idx] = w;
                                        n_neighbors += 1;
                                        break;
                                    } else if (hash_keys[idx] == v_new) {
                                        hash_vals[idx] += w;
                                        break;
                                    }
                                }

                            }
                        endfor
                        real_neighborhood_size[u_new] = n_neighbors;
                    } else {
                        // vertex was matched to u_partner, consider both neighbors
                        vertex_t u_new = matching.get_n(u);

                        // we need to insert v_new into the hash table of u_new
                        u64 offset = hash_offset[u_new];
                        u64 hash_table_size = hash_offset[u_new + 1] - hash_offset[u_new];
                        u64 n_neighbors = 0;
                        forall_guivw(g, u, i, v, w)
                            {
                                if (u_partner == v) { continue; } // this edge vanishes since it is matched

                                vertex_t v_new = matching.get_n(v);
                                u64 key = (v_new * 11400714819323198485ull) % hash_table_size;

                                // linear probing to find empty position or the key
                                for (u64 j = 0; j < hash_table_size; ++j) {
                                    u64 idx = offset + ((key + j) % hash_table_size);

                                    if (hash_keys[idx] == m_n) {
                                        hash_keys[idx] = v_new;
                                        hash_vals[idx] = w;
                                        n_neighbors += 1;
                                        break;
                                    } else if (hash_keys[idx] == v_new) {
                                        hash_vals[idx] += w;
                                        break;
                                    }
                                }

                            }
                        endfor

                        forall_guivw(g, u_partner, i, v, w)
                            {
                                if (u == v) { continue; } // this edge vanishes since it is matched

                                vertex_t v_new = matching.get_n(v);
                                u64 key = (v_new * 11400714819323198485ull) % hash_table_size;

                                // linear probing to find empty position or the key
                                for (u64 j = 0; j < hash_table_size; ++j) {
                                    u64 idx = offset + ((key + j) % hash_table_size);

                                    if (hash_keys[idx] == m_n) {
                                        hash_keys[idx] = v_new;
                                        hash_vals[idx] = w;
                                        n_neighbors += 1;
                                        break;
                                    } else if (hash_keys[idx] == v_new) {
                                        hash_vals[idx] += w;
                                        break;
                                    }
                                }

                            }
                        endfor
                        real_neighborhood_size[u_new] = n_neighbors;
                    }
                }
            endfor

            // prefix sum to create new neighborhoods
            m_neighborhoods[0] = 0;
            for (size_t i = 0; i < m_n; ++i) {
                m_neighborhoods[i + 1] = m_neighborhoods[i] + real_neighborhood_size[i];
            }

            // set true m
            m_m = m_neighborhoods[m_n];

            m_edges_v.initialize(m_m);
            m_edges_w.initialize(m_m);

#pragma omp parallel for default(none) shared(hash_offset, hash_keys, hash_vals)
            for (vertex_t u = 0; u < m_n; ++u) {
                size_t idx = m_neighborhoods[u];
                for (u64 i = hash_offset[u]; i < hash_offset[u + 1]; ++i) {
                    if (hash_keys[i] != m_n) {
                        m_edges_v[idx] = hash_keys[i];
                        m_edges_w[idx] = hash_vals[i];
                        idx += 1;
                    }
                }
            }
        }

        void parallel_initialize(const DeepCSRGraph &g,
                                 Matching &matching,
                                 u64 threads) {
            matching.set_translation();

            m_n = matching.get_n_coarse_nodes();
            m_graph_weight = g.m_graph_weight;
            m_v_weights.initialize(m_n);
            m_neighborhoods.initialize(m_n + 1);

            struct thread_info {
                vertex_t s_vertex = 0;
                vertex_t n_assigned_vertices = 0;
                vertex_t n_actual_vertices = 0;
                vertex_t curr_m = 0;
                std::vector<size_t> neighborhood;
                std::vector<vertex_t> edges_v;
                std::vector<weight_t> edges_w;
            };
            std::vector<thread_info> t_infos(threads);

            f64 alpha = 2.0; // example scaling parameter > 1.0

            // Step 1: Compute weights
            std::vector<f64> weights(threads);
            f64 total_weight = 0.0;
            for (size_t i = 0; i < threads; ++i) {
                weights[i] = std::pow((f64) (i + 1), alpha);
                total_weight += weights[i];
            }

            vertex_t n_total_vertices = g.get_n();

            // Step 2: Assign number of vertices per thread
            vertex_t current_start = 0;
            vertex_t assigned_total = 0;
            for (size_t i = 0; i < threads; ++i) {
                f64 fraction = weights[i] / total_weight;
                vertex_t n_assign = (vertex_t) (std::round(fraction * (f64) n_total_vertices));

                // Ensure last thread takes any rounding residual
                if (i == threads - 1) { n_assign = n_total_vertices - assigned_total; }

                t_infos[i].s_vertex = current_start;
                t_infos[i].n_assigned_vertices = n_assign;
                current_start += n_assign;
                assigned_total += n_assign;
            }

            // each thread determines the neighborhood of their assigned vertices
#pragma omp parallel num_threads(threads) default(none) shared(g, t_infos, matching)
            {
                u64 t_id = omp_get_thread_num();
                t_infos[t_id].neighborhood.push_back(0);
                for (vertex_t old_u = t_infos[t_id].s_vertex; old_u < t_infos[t_id].s_vertex + t_infos[t_id].n_assigned_vertices; ++old_u) {
                    vertex_t old_v = matching.get_partner(old_u);

                    if (old_u > old_v) { continue; }

                    weight_t old_u_w = g.weight(old_u);
                    weight_t old_v_w = old_u == old_v ? 0 : g.weight(old_v);

                    vertex_t new_u = matching.get_n(old_u);
                    m_v_weights[new_u] = old_u_w + old_v_w;

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
                        for (size_t j = t_infos[t_id].neighborhood.back(); j < t_infos[t_id].curr_m; ++j) {
                            if (t_infos[t_id].edges_v[j] == vv) {
                                t_infos[t_id].edges_w[j] += ww;
                                found = true;
                                break;
                            }
                        }
                        if (!found) {
                            t_infos[t_id].edges_v.push_back(vv);
                            t_infos[t_id].edges_w.push_back(ww);
                            t_infos[t_id].curr_m += 1;
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
                            for (size_t j = t_infos[t_id].neighborhood.back();
                                 j < t_infos[t_id].curr_m; ++j) {
                                if (t_infos[t_id].edges_v[j] == vv) {
                                    t_infos[t_id].edges_w[j] += ww;
                                    found = true;
                                    break;
                                }
                            }
                            if (!found) {
                                t_infos[t_id].edges_v.push_back(vv);
                                t_infos[t_id].edges_w.push_back(ww);
                                t_infos[t_id].curr_m += 1;
                            }
                        }
                    }
                    t_infos[t_id].neighborhood.push_back(t_infos[t_id].curr_m);
                }
                t_infos[t_id].n_actual_vertices = t_infos[t_id].neighborhood.size() - 1;
            }

            // determine the number of edges
            m_m = 0;
            for (size_t i = 0; i < threads; ++i) {
                m_m += t_infos[i].curr_m;
            }
            m_edges_v.initialize(m_m);
            m_edges_w.initialize(m_m);

            // each thread copies its data to the correct place in the real neighborhood
            m_neighborhoods[0] = 0;
#pragma omp parallel num_threads(threads) default(none) shared(t_infos) firstprivate(threads)
            {
                u64 t_id = omp_get_thread_num();

                // determine how many vertices and edges come before
                vertex_t previous_m = 0;
                vertex_t previous_n = 0;
                for (size_t i = 0; i < t_id; ++i) {
                    previous_m += t_infos[i].curr_m;
                    previous_n += t_infos[i].n_actual_vertices;
                }

                // copy neighborhood sizes
                for (size_t i = 0; i < t_infos[t_id].n_actual_vertices; ++i) {
                    m_neighborhoods[previous_n + i + 1] = t_infos[t_id].neighborhood[i + 1] + previous_m;
                }

                // copy all edges
                for (size_t i = 0; i < t_infos[t_id].curr_m; ++i) {
                    m_edges_v[previous_m + i] = t_infos[t_id].edges_v[i];
                    m_edges_w[previous_m + i] = t_infos[t_id].edges_w[i];
                }

                if (t_id == threads - 1) {
                    m_m = t_infos[t_id].curr_m + previous_m;
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
