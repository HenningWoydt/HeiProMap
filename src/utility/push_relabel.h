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

#ifndef HEIPROMAP_PUSH_RELABEL_H
#define HEIPROMAP_PUSH_RELABEL_H

#include <cstdlib>
#include <cstring>
#include <algorithm>

#include "memory_stack.h"

namespace HeiProMap {

enum Label : int { SOURCE = 0, SINK = 1 };

template<class T>
class Array {
    T* p = nullptr;
public:
    Array() = default;
    void init(void* mem) { p = static_cast<T*>(mem); }
    T& operator[](int i) { return p[i]; }
    const T& operator[](int i) const { return p[i]; }
    T* data() { return p; }
    const T* data() const { return p; }
};

template<class Cap>
class PushRelabel {
    struct Arc {
        int to;
        int rev;
        Cap cap;
    };

    int n, s, t;
    int total_arcs;
    Array<Arc> arcs;
    Array<int> first;       // n+1
    Array<Cap> excess;      // n
    Array<int> height;      // n
    Array<int> cur;         // n
    Array<int> cnt;         // 2n+1

    // Highest-label buckets
    Array<int> bucket_head; // 2n+1
    Array<int> bucket_next; // n
    int max_active;

    // FIFO queue for phase 1
    Array<int> fifo;        // n
    int fifo_head, fifo_tail;
    Array<int> in_fifo;     // n — flag to avoid duplicate enqueues

    // BFS queue (shares space with tmp_edges after build)
    Array<int> bfs_queue;   // n

    int work;
    int work_limit;
    bool phase2; // true = highest-label phase

    struct TmpEdge { int u, v; Cap cap, rev_cap; };
    Array<TmpEdge> tmp_edges;
    int tmp_count;

    Array<int> init_side;  // n: 0=sink, 1=source (warm-start)

    MemoryStack* mem;

    void hl_enqueue(int u) {
        bucket_next[u] = bucket_head[height[u]];
        bucket_head[height[u]] = u;
        if (height[u] > max_active) max_active = height[u];
    }

    void fifo_enqueue(int u) {
        if (in_fifo[u]) return;
        in_fifo[u] = 1;
        fifo[fifo_tail++] = u;
        if (fifo_tail == n) fifo_tail = 0;
    }

    void activate(int u) {
        if (phase2)
            hl_enqueue(u);
        else
            fifo_enqueue(u);
    }

    void push(int u, Arc& a) {
        Cap d = std::min(excess[u], a.cap);
        if (d == 0) return;
        a.cap -= d;
        arcs[a.rev].cap += d;
        excess[u] -= d;
        excess[a.to] += d;
        if (excess[a.to] == d && a.to != s && a.to != t)
            activate(a.to);
    }

    void gap(int h) {
        for (int i = 0; i < n; i++) {
            if (height[i] > h && height[i] < n && i != s && i != t) {
                cnt[height[i]]--;
                height[i] = n + 1;
                cnt[height[i]]++;
            }
        }
    }

    void relabel(int u) {
        work++;
        cnt[height[u]]--;
        int old_h = height[u];
        height[u] = 2 * n;
        for (int i = first[u]; i < first[u + 1]; i++)
            if (arcs[i].cap > 0 && height[arcs[i].to] + 1 < height[u])
                height[u] = height[arcs[i].to] + 1;
        cnt[height[u]]++;
        cur[u] = first[u];
        if (cnt[old_h] == 0)
            gap(old_h);
    }

    void discharge(int u) {
        while (excess[u] > 0) {
            int end = first[u + 1];
            for (int& i = cur[u]; i < end; i++) {
                Arc& a = arcs[i];
                if (a.cap > 0 && height[u] == height[a.to] + 1) {
                    push(u, a);
                    if (excess[u] == 0) return;
                }
            }
            relabel(u);
            if (height[u] >= n) return;
        }
    }

    void global_relabel() {
        int bkt_size = 2 * n + 1;
        std::memset(bucket_head.data(), -1, bkt_size * sizeof(int));
        std::memset(cnt.data(), 0, bkt_size * sizeof(int));
        max_active = 0;

        for (int i = 0; i < n; i++) height[i] = n;
        height[t] = 0;
        height[s] = n;

        int qh = 0, qt = 0;
        bfs_queue[qt++] = t;

        while (qh < qt) {
            int v = bfs_queue[qh++];
            int new_h = height[v] + 1;
            for (int i = first[v]; i < first[v + 1]; i++) {
                int w = arcs[i].to;
                if (height[w] == n && w != s) {
                    if (arcs[arcs[i].rev].cap > 0) {
                        height[w] = new_h;
                        cnt[new_h]++;
                        cur[w] = first[w];
                        if (excess[w] > 0 && w != t) {
                            if (phase2)
                                hl_enqueue(w);
                            else {
                                fifo_enqueue(w);
                            }
                        }
                        bfs_queue[qt++] = w;
                    }
                }
            }
        }
        cnt[0] = 1;
        work = 0;
    }

public:
    PushRelabel() : n(0), s(0), t(0), total_arcs(0), max_active(0),
        fifo_head(0), fifo_tail(0),
        work(0), work_limit(0), phase2(false), tmp_count(0),
         mem(nullptr) {}

    void init(int n_, int s_, int t_, int m, MemoryStack& stack) {
        n = n_; s = s_; t = t_;
        mem = &stack;
        tmp_count = 0;
        total_arcs = 2 * m;

        int bkt_size = 2 * n + 1;
        size_t tmp_or_bfs = std::max(align64(m * sizeof(TmpEdge)), align64(n * sizeof(int)));

        size_t need = align64(total_arcs * sizeof(Arc))
                    + align64((n + 1) * sizeof(int))   // first
                    + align64(n * sizeof(Cap))          // excess
                    + align64(n * sizeof(int))          // height
                    + align64(n * sizeof(int))          // cur
                    + align64(bkt_size * sizeof(int))   // cnt
                    + align64(bkt_size * sizeof(int))   // bucket_head
                    + align64(n * sizeof(int))          // bucket_next
                    + align64(n * sizeof(int))          // fifo
                    + align64(n * sizeof(int))          // in_fifo
                    + align64(n * sizeof(int))          // init_side
                    + tmp_or_bfs;                       // tmp_edges / bfs_queue

        mem->ensure(need);
        mem->clear();

        arcs.init(mem->get_memory(total_arcs * sizeof(Arc)));
        first.init(mem->get_memory((n + 1) * sizeof(int)));
        excess.init(mem->get_memory(n * sizeof(Cap)));
        height.init(mem->get_memory(n * sizeof(int)));
        cur.init(mem->get_memory(n * sizeof(int)));
        cnt.init(mem->get_memory(bkt_size * sizeof(int)));
        bucket_head.init(mem->get_memory(bkt_size * sizeof(int)));
        bucket_next.init(mem->get_memory(n * sizeof(int)));
        fifo.init(mem->get_memory(n * sizeof(int)));
        in_fifo.init(mem->get_memory(n * sizeof(int)));
        init_side.init(mem->get_memory(n * sizeof(int)));
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

        int* pos = bucket_head.data();
        std::memcpy(pos, first.data(), (n + 1) * sizeof(int));

        for (int i = 0; i < m; i++) {
            auto& e = tmp_edges[i];
            int fwd = pos[e.u]++;
            int bwd = pos[e.v]++;
            arcs[fwd] = {e.v, bwd, e.cap};
            arcs[bwd] = {e.u, fwd, e.rev_cap};
        }

        std::memset(excess.data(), 0, n * sizeof(Cap));
        std::memset(height.data(), 0, n * sizeof(int));
        int bkt_size = 2 * n + 1;
        std::memset(cnt.data(), 0, bkt_size * sizeof(int));
        std::memset(bucket_head.data(), -1, bkt_size * sizeof(int));
        std::memset(in_fifo.data(), 0, n * sizeof(int));
        for (int i = 0; i < n; i++) cur[i] = first[i];
        max_active = 0;
        fifo_head = fifo_tail = 0;
        phase2 = false;
        work = 0;
        work_limit = (int)((long long)n * 6 / 10);
        if (work_limit < 1) work_limit = 1;
    }

    Cap maxflow() {
        build();

        height[s] = n;
        cnt[0] = n - 1;
        cnt[n] = 1;

        // Saturate all source arcs (required for preflow correctness)
        for (int i = first[s]; i < first[s + 1]; i++) {
            Arc& a = arcs[i];
            if (a.cap > 0) {
                excess[a.to] += a.cap;
                arcs[a.rev].cap += a.cap;
                a.cap = 0;
                if (a.to != t)
                    fifo_enqueue(a.to);
            }
        }

        global_relabel();

        // Phase 1: FIFO
        while (fifo_head != fifo_tail) {
            int u = fifo[fifo_head++];
            if (fifo_head == n) fifo_head = 0;
            in_fifo[u] = 0;
            if (excess[u] > 0)
                discharge(u);
            if (work >= work_limit)
                global_relabel();
        }

        // Phase 2: Highest-label for remaining excess
        phase2 = true;
        // Move all active nodes into HL buckets
        for (int i = 0; i < n; i++) {
            if (excess[i] > 0 && i != s && i != t && height[i] < n)
                hl_enqueue(i);
        }

        while (max_active >= 0) {
            while (max_active >= 0 && bucket_head[max_active] == -1)
                max_active--;
            if (max_active < 0) break;
            int u = bucket_head[max_active];
            bucket_head[max_active] = bucket_next[u];
            if (excess[u] > 0)
                discharge(u);
            if (work >= work_limit)
                global_relabel();
        }

        return excess[t];
    }

    Label what_label(int v) const {
        return height[v] >= n ? SOURCE : SINK;
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

#endif // HEIPROMAP_PUSH_RELABEL_H
