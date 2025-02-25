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

#ifndef HEIPROMAP_SIMPLE_GRAPH_H
#define HEIPROMAP_SIMPLE_GRAPH_H

#include <iostream>
#include <regex>
#include <vector>

#include "../../definitions.h"
#include "../../macros.h"
#include "../interfaces/ISerialGraph.h"
#include "../utility/utils.h"

namespace HeiProMap {
    /**
        * Standard undirected Graph that can hold vertex and edge weights.
        */
    class SimpleGraph final : public ISerialGraph {
        vertex_t n = 0;
        vertex_t m = 0;

        weight_t vertex_weights = 0;
        weight_t edge_weights   = 0;

        std::vector<weight_t>            v_weights;
        std::vector<std::vector<EdgeVW>> adj;
        size_t                           curr_m = 0;

    public:
        explicit SimpleGraph(const std::string &graph_in) {
            if (!file_exists(graph_in)) {
                std::cerr << "File " << graph_in << " does not exist!" << std::endl;
                exit(EXIT_FAILURE);
            }

            std::ifstream file(graph_in);
            std::string   line(64, ' ');
            if (file.is_open()) {
                bool has_v_weights  = false;
                bool has_e_weights  = false;
                u64  expected_edges = 0;

                // read in header
                while (std::getline(file, line)) {
                    if (line[0] == '%') { continue; }

                    // remove leading and trailing whitespaces, replace double whitespaces
                    line.erase(0, line.find_first_not_of(' ')).erase(line.find_last_not_of(' ') + 1);
                    line = std::regex_replace(line, std::regex("\\s{2,}"), " ");

                    // read in header
                    std::vector<std::string> header = split(line, ' ');
                    n              = std::stoi(header[0]);
                    m              = 0;
                    expected_edges = std::stoi(header[1]) * 2;

                    // allocate space
                    v_weights.resize(n, 1);
                    vertex_weights = (weight_t) n;
                    adj.resize(n);

                    // read in header
                    std::string fmt = "000";
                    if (header.size() == 3 && header[2].size() == 3) {
                        fmt = header[2];
                    }
                    has_v_weights   = fmt[1] == '1';
                    has_e_weights   = fmt[2] == '1';

                    break;
                }

                // read in edges
                vertex_t         u = 0;
                std::vector<u64> ints;

                while (std::getline(file, line)) {
                    if (line[0] == '%') { continue; }
                    // remove leading and trailing whitespaces, replace double whitespaces
                    str_to_ints(line, ints);

                    u64 i = 0;

                    // check if vertex weights
                    if (has_v_weights) {
                        weight_t w = (weight_t) ints[i++];
                        vertex_weights = vertex_weights - v_weights[u] + w;
                        v_weights[u] = w;
                    }

                    while (i < ints.size()) {
                        vertex_t v = (vertex_t) ints[i++] - 1;

                        weight_t w = 1;

                        // check if edge weights
                        if (has_e_weights) { w = (weight_t) ints[i++]; }

                        adj[u].emplace_back(v, w);
                        m += 1;
                    }

                    u += 1;
                }

                if (expected_edges != m) {
                    std::cerr << "Number of expected edges " << expected_edges << " not equal to number edges " << m << " found!" << std::endl;
                    exit(EXIT_FAILURE);
                }
            } else {
                std::cerr << "Could not open file " << graph_in << "!" << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        void set_vertex_weight(const vertex_t u, const weight_t weight = 1) {
            ASSERT(u < n);
            vertex_weights = vertex_weights - v_weights[u] + weight;
            v_weights[u] = weight;
        }

        SimpleGraph(const SimpleGraph &g,
                    EdgeUV *matches,
                    size_t &matches_size) {
            matches = ASSUME_ALIGNED(EdgeUV*, matches, 64);

            n = g.get_n();
            m = 0;
            v_weights.resize(n);
            adj.resize(n + 1);
            vertex_weights = g.vertex_weights;

            // define the state of each vertex
            constexpr u8          NOT_MATCHED    = 0;
            constexpr u8          FIRST_MATCHED  = 1;
            constexpr u8          SECOND_MATCHED = 2;
            std::vector<u8>       vertex_state(n, NOT_MATCHED);
            std::vector<vertex_t> vertex_neighbor(n);

            // check the matching
            for (size_t i = 0; i < matches_size; ++i) {
                const auto [u, v] = matches[i];

                vertex_state[u]    = FIRST_MATCHED;
                vertex_state[v]    = SECOND_MATCHED;
                vertex_neighbor[u] = v;
                vertex_neighbor[v] = u;

                v_weights[v] = 0;
                v_weights[u] = g.get_weight(u) + g.get_weight(v);
            }

            for (vertex_t u = 0; u < n; ++u) {
                if (vertex_state[u] == NOT_MATCHED) {
                    // copy it to the next graph
                    v_weights[u] = g.get_weight(u);
                    for (size_t i = 0; i < g.size(u); ++i) {
                        vertex_t vv = g.neighbor(u, i);
                        weight_t ww = g.get_weight(u, i);

                        // if the vv vertex is matched, then make an edge to the neighbor vertex
                        if (vertex_state[vv] == SECOND_MATCHED) { vv = vertex_neighbor[vv]; }

                        // if exists, add weight, otherwise append
                        bool found = false;
                        for (auto &e: adj[u]) {
                            if (e.v == vv) {
                                found = true;
                                e.w += ww;
                            }
                        }
                        if (!found) { adj[u].emplace_back(vv, ww); }
                    }
                } else if (vertex_state[u] == FIRST_MATCHED) {
                    // the vertex gets all neighbors of v
                    vertex_t v = vertex_neighbor[u];

                    // v_weights[u] = g.get_weight(u) + g.get_weight(v);
                    for (size_t i = 0; i < g.size(u); ++i) {
                        vertex_t vv = g.neighbor(u, i);
                        weight_t ww = g.get_weight(u, i);

                        // do not add edge to matched vertex
                        if (vv == v) { continue; }

                        // if the vv vertex is matched, then make an edge to the neighbor vertex
                        if (vertex_state[vv] == SECOND_MATCHED) { vv = vertex_neighbor[vv]; }

                        // if exists, add weight, otherwise append
                        bool found = false;
                        for (auto &e: adj[u]) {
                            if (e.v == vv) {
                                found = true;
                                e.w += ww;
                            }
                        }
                        if (!found) { adj[u].emplace_back(vv, ww); }
                    }
                    for (size_t i = 0; i < g.size(v); ++i) {
                        vertex_t vv = g.neighbor(v, i);
                        weight_t ww = g.get_weight(v, i);

                        // do not add edge to matched vertex
                        if (vv == u) { continue; }

                        // if the vv vertex is matched, then make an edge to the neighbor vertex
                        if (vertex_state[vv] == SECOND_MATCHED) { vv = vertex_neighbor[vv]; }

                        // if exists, add weight, otherwise append
                        bool found = false;
                        for (auto &e: adj[u]) {
                            if (e.v == vv) {
                                found = true;
                                e.w += ww;
                            }
                        }
                        if (!found) { adj[u].emplace_back(vv, ww); }
                    }
                }
            }

            for (vertex_t u = 0; u < n; ++u) {
                m += adj[u].size();
            }
        }

        vertex_t get_n() const override { return n; }

        vertex_t get_m() const override { return m; }

        weight_t get_weight() const override { return vertex_weights; }

        weight_t get_weight(const vertex_t u) const override { return v_weights[u]; }

        size_t size(const vertex_t u) const override { return adj[u].size(); }

        vertex_t neighbor(const vertex_t u, const size_t idx) const override { return adj[u][idx].v; }

        weight_t get_weight(const vertex_t u, const size_t idx) const override { return adj[u][idx].w; }

        // edge manipulation
        bool edge_exists(const vertex_t u, const vertex_t v) const override {
            for (auto e: adj[u]) {
                if (e.v == v) {
                    return true;
                }
            }
            return false;
        }

        class NeighborhoodIterator {
            std::vector<EdgeVW> &edges;

        public:
            explicit NeighborhoodIterator(std::vector<EdgeVW> &t_edges) : edges(t_edges) {}

            class Iterator {
                std::vector<EdgeVW> &edges;
                size_t idx;

            public:
                // Constructor
                Iterator(std::vector<EdgeVW> &edges, size_t idx) : edges(edges), idx(idx) {}

                // Dereference operator
                EdgeVW operator*() const { return edges[idx]; }

                // Pre-increment operator
                Iterator &operator++() {
                    idx++;
                    return *this;
                }

                bool operator!=(const Iterator &other) const { return idx != other.idx; }
            };

            Iterator begin() const { return {edges, 0}; }

            Iterator end() const { return {edges, edges.size()}; }
        };

        NeighborhoodIterator operator[](const vertex_t u) {
            return NeighborhoodIterator(adj[u]);
        }
    };
}

#endif //HEIPROMAP_SIMPLE_GRAPH_H
