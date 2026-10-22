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
#include <algorithm>
#include <unordered_map>
#include <utility>
#include <limits>

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
        std::vector<std::vector<Neighbor> > neighbors;

        std::vector<vertex_t> dirty_list;
        std::vector<s64> dirty_indices;
        std::vector<vertex_t> dirty_list_partition;
        std::vector<s64> dirty_indices_partition;
        
        // Filled by apply_pending(): vertices that came into existence (with their weight) and
        // vertices that disappeared (with the weight they had), accumulated until cleared.
        std::vector<std::pair<vertex_t, weight_t> > new_vertices;
        std::vector<std::pair<vertex_t, weight_t> > removed_vertices;

        std::function<void(vertex_t)> dirty_callback = nullptr;

        // The commands use external vertex ids. Internally the graph works on dense slots
        // 0..n, and the slot of a removed vertex is reused by the next vertex that appears.
        static constexpr vertex_t INVALID_ID = std::numeric_limits<vertex_t>::max();
        std::vector<vertex_t> ext_to_int;
        std::vector<vertex_t> int_to_ext;
        std::vector<vertex_t> free_slots;

        vertex_t find_slot(vertex_t ext) const {
            return ext < ext_to_int.size() ? ext_to_int[ext] : INVALID_ID;
        }

        vertex_t get_or_create_slot(vertex_t ext) {
            if (ext >= ext_to_int.size()) ext_to_int.resize(ext + 1, INVALID_ID);
            if (ext_to_int[ext] != INVALID_ID) return ext_to_int[ext];

            vertex_t slot;
            if (!free_slots.empty()) {
                slot = free_slots.back();
                free_slots.pop_back();
            } else {
                slot = n;
                ensure_capacity(n + 1);
            }
            ext_to_int[ext] = slot;
            int_to_ext[slot] = ext;
            return slot;
        }

        // Frees the slot of a vertex whose weight dropped to 0, so it can be reused.
        void release_if_empty(vertex_t slot) {
            if (slot >= n || v_weights[slot] != 0 || int_to_ext[slot] == INVALID_ID) return;
            if (!neighbors[slot].empty()) remove_vertex(slot);
            ext_to_int[int_to_ext[slot]] = INVALID_ID;
            int_to_ext[slot] = INVALID_ID;
            free_slots.push_back(slot);
        }

        void ensure_capacity(vertex_t cap) {
            if (cap <= n) return;
            v_weights.resize(cap, 0);
            dirty_indices.resize(cap, -1);
            dirty_indices_partition.resize(cap, -1);
            neighbors.resize(cap);
            int_to_ext.resize(cap, INVALID_ID);
            n = cap;
        }

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

        // A graph modification that is buffered until the next command that needs the graph.
        struct GraphOp {
            enum Type { ADD_V, ADD_E, REMOVE_V, REMOVE_E, SET_VW, SET_EW } type;
            vertex_t u = 0;
            vertex_t v = 0;
            weight_t w = 1;
            bool has_w = false;
        };

        std::vector<GraphOp> pending_ops;

        // Net weight change of one vertex caused by a batch of modifications.
        struct VertexChange {
            vertex_t v;
            weight_t before;
            weight_t after;
            bool is_new;  // the vertex came into existence (again) during the batch
        };

        // Applies all buffered modifications to the graph in the order they were issued and
        // records which vertices ended up new or removed (see new_vertices / removed_vertices).
        // Returns the vertices whose weight changed, so dependent structures can follow.
        std::vector<VertexChange> apply_pending() {
            std::vector<VertexChange> changes;
            if (pending_ops.empty()) return changes;
            HEIPROMAP_PROFILE_SCOPE("solver", "apply_pending", "apply_pending");

            struct Track {
                weight_t before = 0; // weight before the batch
                bool zeroed = false; // weight dropped to 0 at some point in the batch
            };
            std::unordered_map<vertex_t, Track> tracked;

            auto weight_of = [this](vertex_t x) -> weight_t { return x < n ? v_weights[x] : 0; };
            auto begin_touch = [&](vertex_t x) {
                tracked.try_emplace(x, Track{weight_of(x), false});
                return weight_of(x);
            };
            auto end_touch = [&](vertex_t x, weight_t prev) {
                if (prev > 0 && weight_of(x) == 0) tracked[x].zeroed = true;
            };

            for (const GraphOp &op: pending_ops) {
                bool two = (op.type == GraphOp::ADD_E || op.type == GraphOp::REMOVE_E || op.type == GraphOp::SET_EW);
                bool create = (op.type != GraphOp::REMOVE_V && op.type != GraphOp::REMOVE_E);

                // Translate the external ids; removals of unknown vertices are ignored.
                vertex_t su = create ? get_or_create_slot(op.u) : find_slot(op.u);
                vertex_t sv = !two ? INVALID_ID : (create ? get_or_create_slot(op.v) : find_slot(op.v));
                if (su == INVALID_ID || (two && sv == INVALID_ID)) continue;

                weight_t prev_u = begin_touch(su);
                weight_t prev_v = two ? begin_touch(sv) : 0;

                switch (op.type) {
                    case GraphOp::ADD_V: add_vertex(su, op.w);
                        break;
                    case GraphOp::ADD_E: add_edge(su, sv, op.w);
                        break;
                    case GraphOp::REMOVE_V:
                        if (op.has_w) remove_vertex(su, op.w);
                        else remove_vertex(su);
                        break;
                    case GraphOp::REMOVE_E:
                        if (op.has_w) remove_edge(su, sv, op.w);
                        else remove_edge(su, sv);
                        break;
                    case GraphOp::SET_VW: set_vertex_weight(su, op.w);
                        break;
                    case GraphOp::SET_EW: set_edge_weight(su, sv, op.w);
                        break;
                }

                end_touch(su, prev_u);
                if (two) end_touch(sv, prev_v);

                release_if_empty(su);
                if (two) release_if_empty(sv);
            }
            pending_ops.clear();

            // Deterministic order: by vertex id.
            std::vector<vertex_t> ids;
            ids.reserve(tracked.size());
            for (const auto &entry: tracked) ids.push_back(entry.first);
            std::sort(ids.begin(), ids.end());

            for (vertex_t x: ids) {
                const Track &t = tracked[x];
                weight_t after = weight_of(x);
                bool removed = t.before > 0 && (t.zeroed || after == 0);
                bool is_new = after > 0 && (t.before == 0 || t.zeroed);
                if (removed) removed_vertices.emplace_back(x, t.before);
                if (is_new) new_vertices.emplace_back(x, after);
                if (after != t.before || removed || is_new) changes.push_back({x, t.before, after, is_new});
            }
            return changes;
        }

        // Marks every vertex in new_vertices.
        std::vector<char> new_vertex_mask() const {
            std::vector<char> mask(n, 0);
            for (const auto &entry: new_vertices) mask[entry.first] = 1;
            return mask;
        }

        void add_vertex(vertex_t v, weight_t w = 1) {
            ensure_capacity(v + 1);
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
            weight_t removed = std::min(w, v_weights[v]);
            v_weights[v] -= removed;
            g_weight -= removed;
            mark_dirty(v);
            if (v_weights[v] == 0) {
                remove_vertex(v);
            }
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

        void clear_removed_vertices() {
            removed_vertices.clear();
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
