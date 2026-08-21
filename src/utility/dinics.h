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

#ifndef HEIPROMAP_DINICS_H
#define HEIPROMAP_DINICS_H

#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <limits>

#include "macros.h"
#include "memory_stack.h"
#include "push_relabel.h"

namespace HeiProMap {

template<class Cap>
class Dinics {
    struct Arc {
        int to;
        int rev;
        Cap cap;
    };

    int n, s, t;
    int total_arcs;
    Array<Arc> arcs;
    Array<int> first;       // n+1
    Array<int> level;       // n
    Array<int> ptr;         // n
    Array<int> bfs_queue;   // n

    struct TmpEdge { int u, v; Cap cap, rev_cap; };
    Array<TmpEdge> tmp_edges;
    int tmp_count;

    MemoryStack* mem;

    bool bfs() {
        for (int i = 0; i < n; i++) {
            level[i] = -1;
        }
        level[s] = 0;
        int qh = 0, qt = 0;
        bfs_queue[qt++] = s;

        while (qh < qt) {
            int u = bfs_queue[qh++];
            for (int i = first[u]; i < first[u + 1]; i++) {
                const Arc& a = arcs[i];
                if (a.cap > 0 && level[a.to] == -1) {
                    level[a.to] = level[u] + 1;
                    bfs_queue[qt++] = a.to;
                }
            }
        }
        return level[t] != -1;
    }

    Cap dfs(int u, Cap pushed) {
        if (pushed == 0) return 0;
        if (u == t) return pushed;

        for (int& i = ptr[u]; i < first[u + 1]; i++) {
            Arc& a = arcs[i];
            int trg = a.to;
            if (level[u] + 1 != level[trg] || a.cap == 0) {
                continue;
            }
            Cap tr = dfs(trg, std::min(pushed, a.cap));
            if (tr == 0) {
                continue;
            }
            a.cap -= tr;
            arcs[a.rev].cap += tr;
            return tr;
        }
        return 0;
    }

public:
    Dinics() : n(0), s(0), t(0), total_arcs(0), tmp_count(0), mem(nullptr) {}

    void init(int n_, int s_, int t_, int m, MemoryStack& stack) {
        n = n_; s = s_; t = t_;
        mem = &stack;
        tmp_count = 0;
        total_arcs = 2 * m;

        size_t tmp_or_bfs = std::max(align64(m * sizeof(TmpEdge)), align64(n * sizeof(int)));

        size_t need = align64(total_arcs * sizeof(Arc))
                    + align64((n + 1) * sizeof(int))   // first
                    + align64(n * sizeof(int))          // level
                    + align64(n * sizeof(int))          // ptr
                    + tmp_or_bfs;                       // tmp_edges / bfs_queue

        mem->ensure(need);
        mem->clear();

        arcs.init(mem->get_memory(total_arcs * sizeof(Arc)));
        first.init(mem->get_memory((n + 1) * sizeof(int)));
        level.init(mem->get_memory(n * sizeof(int)));
        ptr.init(mem->get_memory(n * sizeof(int)));

        void* shared = mem->get_memory(tmp_or_bfs);
        tmp_edges.init(shared);
        bfs_queue.init(shared);
    }

    void add_edge(int u, int v, Cap cap, Cap rev_cap = 0) {
        tmp_edges[tmp_count++] = {u, v, cap, rev_cap};
    }

    void build() {
        int m = tmp_count;
        std::memset(first.data(), 0, (n + 1) * sizeof(int));

        for (int i = 0; i < m; i++) {
            first[tmp_edges[i].u + 1]++;
            first[tmp_edges[i].v + 1]++;
        }
        for (int i = 1; i <= n; i++)
            first[i] += first[i - 1];

        // Temporary array for write positions using level array's memory (since level is not used during build)
        Array<int> pos;
        pos.init(level.data());
        std::memcpy(pos.data(), first.data(), (n + 1) * sizeof(int));

        for (int i = 0; i < m; i++) {
            auto& e = tmp_edges[i];
            int fwd = pos[e.u]++;
            int bwd = pos[e.v]++;
            arcs[fwd] = {e.v, bwd, e.cap};
            arcs[bwd] = {e.u, fwd, e.rev_cap};
        }
    }

    Cap maxflow() {
        build();

        Cap flow = 0;
        while (bfs()) {
            for (int i = 0; i < n; i++) {
                ptr[i] = first[i];
            }
            while (Cap pushed = dfs(s, std::numeric_limits<Cap>::max())) {
                flow += pushed;
            }
        }
        ASSERT(validate());
        return flow;
    }

    bool validate() const {
        // 1. Source target reachability: Target must be unreachable in the level graph (min-cut property)
        if (level[t] != -1) return false;
        if (level[s] == -1) return false;

        // 2. Residual capacity consistency across min-cut & valid level bounds
        for (int u = 0; u < n; u++) {
            if (level[u] < -1 || level[u] >= n) return false;
            for (int i = first[u]; i < first[u + 1]; i++) {
                const Arc& a = arcs[i];
                // Reverse arc index valid check
                if (a.rev < 0 || a.rev >= total_arcs) return false;
                // Residual capacity non-negative
                if (a.cap < 0) return false;
                // Cut check: No residual arc from SOURCE to SINK
                if (level[u] != -1 && level[a.to] == -1 && a.cap > 0) return false;
            }
        }
        return true;
    }

    Label what_label(int v) const {
        return level[v] != -1 ? SOURCE : SINK;
    }

    int get_source() const { return s; }
    int get_target() const { return t; }

    // Residual graph accessors for Sanders & Schulz SCC construction
    int arc_first(int u) const { return first[u]; }
    int arc_last(int u) const { return first[u + 1]; }
    int arc_to(int arc_idx) const { return arcs[arc_idx].to; }
    Cap arc_cap(int arc_idx) const { return arcs[arc_idx].cap; }
};

} // namespace HeiProMap

#endif // HEIPROMAP_DINICS_H
