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

#ifndef HEIPROMAP_DYN_GRAPH_H
#define HEIPROMAP_DYN_GRAPH_H

#include <vector>
#include <fstream>
#include <iostream>
#include <functional>

#include "../definitions.h"
#include "../utility/profiler.h"
#include "../utility/utils.h"

namespace HeiProMap {
    struct Neighbor {
        vertex_t u;
        weight_t w;

        Neighbor(vertex_t t_u, weight_t t_w) {
            u = t_u;
            w = t_w;
        }
    };

    class DynGraph {
    public:
        vertex_t n = 0;
        vertex_t m = 0;
        weight_t g_weight = 0;

        std::vector<weight_t> v_weights;
        std::vector<vertex_t> dirty_list;
        std::vector<s64> dirty_indices;
        std::vector<vertex_t> dirty_list_partition;
        std::vector<s64> dirty_indices_partition;
        std::vector<vertex_t> new_vertices;
        std::vector<std::vector<Neighbor> > neighbors;

        std::function<void(vertex_t)> dirty_callback = nullptr;

        void mark_dirty(vertex_t v) {
            if (v < dirty_indices.size() && dirty_indices[v] == -1) {
                dirty_indices[v] = (s64) dirty_list.size();
                dirty_list.push_back(v);
            }
            if (v < dirty_indices_partition.size() && dirty_indices_partition[v] == -1) {
                dirty_indices_partition[v] = (s64) dirty_list_partition.size();
                dirty_list_partition.push_back(v);
            }
            if (dirty_callback) dirty_callback(v);
        }

        void add_vertex(vertex_t v, weight_t w = 1) {
            if (v >= n) {
                v_weights.resize(v + 1, 0);
                dirty_indices.resize(v + 1, -1);
                dirty_indices_partition.resize(v + 1, -1);
                neighbors.resize(v + 1);
                n = v + 1;
                new_vertices.push_back(v);
            }
            g_weight += w;
            v_weights[v] += w;
            mark_dirty(v);
        }

        void add_edge(vertex_t u, vertex_t v, weight_t w = 1) {
            if (!vertex_exists(u)) { add_vertex(u); }
            if (!vertex_exists(v)) { add_vertex(v); }

            bool found = false;
            for (auto &edge: neighbors[u]) {
                if (edge.u == v) {
                    edge.w += w;
                    found = true;
                    break;
                }
            }
            if (!found) {
                neighbors[u].emplace_back(v, w);
                m += 1;
            }

            found = false;
            for (auto &edge: neighbors[v]) {
                if (edge.u == u) {
                    edge.w += w;
                    found = true;
                    break;
                }
            }
            if (!found) {
                neighbors[v].emplace_back(u, w);
                m += 1;
            }
            mark_dirty(u);
            mark_dirty(v);
        }

        void remove_vertex(vertex_t v) {
            if (v >= n) { return; }
            for (auto [nv, nw]: neighbors[v]) {
                mark_dirty(nv);
                for (size_t i = 0; i < neighbors[nv].size(); i++) {
                    if (neighbors[nv][i].u == v) {
                        std::swap(neighbors[nv][i], neighbors[nv].back());
                        neighbors[nv].pop_back();
                        m -= 1;
                        break;
                    }
                }
            }
            g_weight -= v_weights[v];
            v_weights[v] = 0;
            m -= neighbors[v].size();
            neighbors[v].clear();
        }

        void remove_vertex(vertex_t v, weight_t w) {
            if (v >= n) { return; }
            for (auto [nv, nw]: neighbors[v]) {
                mark_dirty(nv);
                for (size_t i = 0; i < neighbors[nv].size(); i++) {
                    if (neighbors[nv][i].u == v) {
                        neighbors[nv][i].w -= w;
                        if (neighbors[nv][i].w == 0) {
                            std::swap(neighbors[nv][i], neighbors[nv].back());
                            neighbors[nv].pop_back();
                            m -= 1;
                        }
                        break;
                    }
                }
            }
            g_weight -= v_weights[v];
            v_weights[v] = 0;
            m -= neighbors[v].size();
            neighbors[v].clear();
        }

        void remove_edge(vertex_t u, vertex_t v) {
            if (u >= n || v >= n) { return; }

            for (size_t i = 0; i < neighbors[u].size(); i++) {
                if (neighbors[u][i].u == v) {
                    std::swap(neighbors[u][i], neighbors[u].back());
                    neighbors[u].pop_back();
                    m -= 1;
                    break;
                }
            }
            for (size_t i = 0; i < neighbors[v].size(); i++) {
                if (neighbors[v][i].u == u) {
                    std::swap(neighbors[v][i], neighbors[v].back());
                    neighbors[v].pop_back();
                    m -= 1;
                    break;
                }
            }
            mark_dirty(u);
            mark_dirty(v);
        }

        void remove_edge(vertex_t u, vertex_t v, weight_t w) {
            if (u >= n || v >= n) { return; }

            for (size_t i = 0; i < neighbors[u].size(); i++) {
                if (neighbors[u][i].u == v) {
                    neighbors[u][i].w -= w;
                    if (neighbors[u][i].w == 0) {
                        std::swap(neighbors[u][i], neighbors[u].back());
                        neighbors[u].pop_back();
                        m -= 1;
                    }
                    break;
                }
            }
            for (size_t i = 0; i < neighbors[v].size(); i++) {
                if (neighbors[v][i].u == u) {
                    neighbors[v][i].w -= w;
                    if (neighbors[v][i].w == 0) {
                        std::swap(neighbors[v][i], neighbors[v].back());
                        neighbors[v].pop_back();
                        m -= 1;
                    }
                    break;
                }
            }
            mark_dirty(u);
            mark_dirty(v);
        }

        void set_vertex_weight(vertex_t v, weight_t w) {
            if (!vertex_exists(v)) {
                add_vertex(v, w);
                return;
            }

            g_weight -= v_weights[v];
            g_weight += w;
            v_weights[v] = w;
            mark_dirty(v);
        }

        void set_edge_weight(vertex_t u, vertex_t v, weight_t w) {
            if (!edge_exists(u, v)) {
                add_edge(u, v, w);
                return;
            }

            for (auto &edge: neighbors[u]) {
                if (edge.u == v) {
                    edge.w = w;
                    break;
                }
            }
            for (auto &edge: neighbors[v]) {
                if (edge.u == u) {
                    edge.w = w;
                    break;
                }
            }
            mark_dirty(u);
            mark_dirty(v);
        }

        void clear_dirty_status() {
            for (vertex_t v: dirty_list) {
                dirty_indices[v] = -1;
            }
            dirty_list.clear();
        }

        void clear_dirty_status_partition() {
            for (vertex_t v: dirty_list_partition) {
                dirty_indices_partition[v] = -1;
            }
            dirty_list_partition.clear();
        }

        void clear_new_vertices() {
            new_vertices.clear();
        }

        void save_to_metis(const std::string &filename) const {
            std::ofstream file(filename);
            if (!file.is_open()) {
                std::cerr << "Could not open file: " << filename << std::endl;
                return;
            }
            file << n << " " << (m / 2) << std::endl;
            for (vertex_t v = 0; v < n; ++v) {
                for (auto edge: neighbors[v]) {
                    file << (edge.u + 1) << " ";
                }
                file << std::endl;
            }
            file.close();
        }

        weight_t get_vertex_weight(vertex_t v) const {
            if (v >= n) return 0;
            return v_weights[v];
        }

        weight_t get_edge_weight(vertex_t u, vertex_t v) const {
            for (auto [nv, nw]: neighbors[u]) {
                if (nv == v) { return nw; }
            }
            return 0;
        }

        bool vertex_exists(vertex_t v) const {
            if (v >= n) { return false; }
            if (v_weights[v] == 0) { return false; }
            return true;
        }

        bool edge_exists(vertex_t u, vertex_t v) const {
            for (auto [nv, nw]: neighbors[u]) {
                if (nv == v) { return true; }
            }
            return false;
        }
    };
}

#endif //HEIPROMAP_DYN_GRAPH_H
