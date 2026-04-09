#ifndef FLOW_MY_EIBFS_I_H
#define FLOW_MY_EIBFS_I_H

#include <cstdint>
#include <algorithm>
#include <cstring>

#define MY_LIKELY(x)   __builtin_expect(!!(x), 1)
#define MY_UNLIKELY(x) __builtin_expect(!!(x), 0)

namespace my_reimpls {
    static constexpr uint32_t INVALID_NODE = UINT32_MAX;
    static constexpr uint32_t INVALID_ARC = UINT32_MAX;

    inline uint64_t align64(uint64_t n) { return (n + 63) & ~uint64_t(63); }

    // O(1) operations this class can algorithmically not be sped up
    class MemoryStack {
        static constexpr size_t ALIGNMENT = 64;
        char *buffer = nullptr;
        size_t capacity = 0;
        size_t offset = 0;

    public:
        ~MemoryStack() { std::free(buffer); }

        MemoryStack() = default;

        MemoryStack(const MemoryStack &) = delete;

        MemoryStack &operator=(const MemoryStack &) = delete;

        MemoryStack(MemoryStack &&other) noexcept : buffer(other.buffer), capacity(other.capacity), offset(other.offset) {
            other.buffer = nullptr;
            other.capacity = 0;
            other.offset = 0;
        }

        MemoryStack &operator=(MemoryStack &&other) noexcept {
            if (this != &other) {
                std::free(buffer);
                buffer = other.buffer;
                capacity = other.capacity;
                offset = other.offset;
                other.buffer = nullptr;
                other.capacity = 0;
                other.offset = 0;
            }
            return *this;
        }

        void ensure(size_t total_bytes) {
            total_bytes = align64(total_bytes);
            if (total_bytes > capacity) {
                std::free(buffer);
                buffer = static_cast<char *>(std::aligned_alloc(ALIGNMENT, total_bytes));
                capacity = total_bytes;
            }
        }

        void *get_memory(size_t size) {
            size = align64(size);
            void *ptr = buffer + offset;
            offset += size;
            return ptr;
        }

        void clear() { offset = 0; }
    };

    // O(1) operations this class can algorithmically not be sped up
    template<class T>
    class Array {
        T *__restrict__ data_;

        static T *assume_aligned(T *p) { return static_cast<T *>(__builtin_assume_aligned(p, 64)); }

    public:
        Array() { data_ = nullptr; }
        Array(T *mem) { data_ = assume_aligned(mem); }
        void init(T *mem) { data_ = assume_aligned(mem); }
        T &operator[](size_t i) { return data_[i]; }
        const T &operator[](size_t i) const { return data_[i]; }
        T *data() { return data_; }
    };

    // O(1) operations this class can algorithmically not be sped up
    template<class T>
    class FixedVector {
        T *__restrict__ data_;
        size_t size_;

        static T *assume_aligned(T *p) { return static_cast<T *>(__builtin_assume_aligned(p, 64)); }

    public:
        FixedVector() {
            data_ = nullptr;
            size_ = 0;
        }

        FixedVector(T *mem) {
            data_ = assume_aligned(mem);
            size_ = 0;
        }

        void init(T *mem) {
            data_ = assume_aligned(mem);
            size_ = 0;
        }

        void push_back(T v) { data_[size_++] = v; }
        void clear() { size_ = 0; }
        bool empty() const { return size_ == 0; }
        T *begin() { return data_; }
        T *end() { return data_ + size_; }

        friend void swap(FixedVector &a, FixedVector &b) noexcept {
            std::swap(a.data_, b.data_);
            std::swap(a.size_, b.size_);
        }
    };

    // O(1) operations this class can algorithmically not be sped up
    template<class Cap>
    struct ArcData {
        static constexpr uint32_t REV_RESIDUAL_BIT = uint32_t(1) << (sizeof(uint32_t) * 8 - 1);
        static constexpr uint32_t REV_MASK = ~REV_RESIDUAL_BIT;
        uint32_t head;
        uint32_t rev; // high bit = isRevResidual
        Cap rCap;

        uint32_t getrev() const { return rev & REV_MASK; }
        bool isRevResidual() const { return rev & REV_RESIDUAL_BIT; }

        void setRevResidual(bool v) {
            if (v) rev |= REV_RESIDUAL_BIT;
            else rev &= REV_MASK;
        }
    };

    // O(1) operations this class can algorithmically not be sped up
    template<class Term>
    struct VertexData {
        int32_t label; // label > 0: distance from s, label < 0: -distance from t
        Term excess; // excess > 0: capacity from s, excess < 0: -capacity to t
        uint32_t parent;
        uint32_t firstSon;
        uint32_t nextNode;
        uint32_t prevSibling;
        int32_t lastAugTimestamp;
        bool isParentCurr;
    };

    // O(1) operations this class can algorithmically not be sped up
    template<class Term>
    class BucketsOneSided {
    public:
        Array<uint32_t> buckets;
        int64_t maxBucket = 0;

        BucketsOneSided() = default;

        BucketsOneSided(uint32_t *t_buckets,
                        int64_t n) {
            buckets = t_buckets;
            maxBucket = 0;
            std::fill_n(t_buckets, n + 1, INVALID_NODE);
        }

        void add_source(VertexData<Term> *nodes, uint32_t u) {
            VertexData<Term> &data = nodes[u];
            data.nextNode = buckets[data.label];
            buckets[data.label] = u;
            if (data.label > maxBucket) { maxBucket = data.label; }
        }

        void add_sink(VertexData<Term> *nodes, uint32_t u) {
            VertexData<Term> &data = nodes[u];
            int64_t bucket = -data.label;
            data.nextNode = buckets[bucket];
            buckets[bucket] = u;
            if (bucket > maxBucket) { maxBucket = bucket; }
        }

        uint32_t popFront(VertexData<Term> *nodes, int64_t bucket) {
            uint32_t u = buckets[bucket];
            if (u != INVALID_NODE) { buckets[bucket] = nodes[u].nextNode; }
            return u;
        }
    };

    // O(1) operations this class can algorithmically not be sped up
    template<class Term>
    class Buckets3Pass {
    public:
        Array<uint32_t> buckets;
        int64_t maxBucket = 0;

        Buckets3Pass() = default;

        Buckets3Pass(uint32_t *t_buckets, int64_t n) {
            buckets = t_buckets;
            maxBucket = 0;
            std::fill_n(t_buckets, n + 1, INVALID_NODE);
        }

        void add_source(VertexData<Term> *nodes, uint32_t u) {
            VertexData<Term> &data = nodes[u];
            data.nextNode = buckets[data.label];
            if (data.nextNode != INVALID_NODE) { nodes[data.nextNode].firstSon = u; }
            buckets[data.label] = u;
            if (data.label > maxBucket) { maxBucket = data.label; }
        }

        void add_sink(VertexData<Term> *nodes, uint32_t u) {
            VertexData<Term> &data = nodes[u];
            int64_t bucket = -data.label;
            data.nextNode = buckets[bucket];
            if (data.nextNode != INVALID_NODE) { nodes[data.nextNode].firstSon = u; }
            buckets[bucket] = u;
            if (bucket > maxBucket) { maxBucket = bucket; }
        }

        uint32_t popFront(VertexData<Term> *nodes, int64_t bucket) {
            uint32_t u = buckets[bucket];
            if (u != INVALID_NODE) {
                buckets[bucket] = nodes[u].nextNode;
                nodes[u].firstSon = INVALID_NODE;
            }
            return u;
        }

        void remove_source(VertexData<Term> *nodes, uint32_t u) {
            VertexData<Term> &data = nodes[u];
            int64_t bucket = data.label;
            if (buckets[bucket] == u) {
                buckets[bucket] = data.nextNode;
            } else {
                nodes[nodes[u].firstSon].nextNode = data.nextNode;
                if (data.nextNode != INVALID_NODE) { nodes[data.nextNode].firstSon = nodes[u].firstSon; }
            }
            nodes[u].firstSon = INVALID_NODE;
        }

        void remove_sink(VertexData<Term> *nodes, uint32_t u) {
            VertexData<Term> &data = nodes[u];
            int64_t bucket = -data.label;
            if (buckets[bucket] == u) {
                buckets[bucket] = data.nextNode;
            } else {
                nodes[nodes[u].firstSon].nextNode = data.nextNode;
                if (data.nextNode != INVALID_NODE) { nodes[data.nextNode].firstSon = nodes[u].firstSon; }
            }
            nodes[u].firstSon = INVALID_NODE;
        }

        bool isEmpty(int64_t bucket) { return buckets[bucket] == INVALID_NODE; }
    };

    // O(1) operations this class can algorithmically not be sped up
    template<class Term>
    class ExcessBuckets {
    public:
        Array<uint32_t> buckets;
        Array<uint32_t> ptrs;
        int64_t maxBucket = 0;
        int64_t minBucket = 0;

        ExcessBuckets() = default;

        ExcessBuckets(uint32_t *t_buckets,
                      uint32_t *t_ptrs,
                      int64_t n) {
            buckets = t_buckets;
            ptrs = t_ptrs;
            std::fill_n(t_buckets, n + 1, INVALID_NODE);
            reset();
        }

        void add_source(VertexData<Term> *nodes, uint32_t u) {
            VertexData<Term> &data = nodes[u];
            ptrs[u * 2] = buckets[data.label];
            if (buckets[data.label] != INVALID_NODE) { ptrs[buckets[data.label] * 2 + 1] = u; }
            buckets[data.label] = u;
            if (data.label > maxBucket) { maxBucket = data.label; }
            if (data.label != 0 && data.label < minBucket) { minBucket = data.label; }
        }

        void add_sink(VertexData<Term> *nodes, uint32_t u) {
            VertexData<Term> &data = nodes[u];
            int64_t bucket = -data.label;
            ptrs[u * 2] = buckets[bucket];
            if (buckets[bucket] != INVALID_NODE) { ptrs[buckets[bucket] * 2 + 1] = u; }
            buckets[bucket] = u;
            if (bucket > maxBucket) { maxBucket = bucket; }
            if (bucket != 0 && bucket < minBucket) { minBucket = bucket; }
        }

        uint32_t popFront(int64_t bucket) {
            uint32_t u = buckets[bucket];
            if (u != INVALID_NODE) { buckets[bucket] = ptrs[u * 2]; }
            return u;
        }

        void remove_source(VertexData<Term> *nodes, uint32_t u) {
            VertexData<Term> &data = nodes[u];
            int64_t bucket = data.label;
            if (buckets[bucket] == u) {
                buckets[bucket] = ptrs[u * 2];
            } else {
                ptrs[ptrs[u * 2 + 1] * 2] = ptrs[u * 2];
                if (ptrs[u * 2] != INVALID_NODE) { ptrs[ptrs[u * 2] * 2 + 1] = ptrs[u * 2 + 1]; }
            }
        }

        void remove_sink(VertexData<Term> *nodes, uint32_t u) {
            VertexData<Term> &data = nodes[u];
            int64_t bucket = -data.label;
            if (buckets[bucket] == u) {
                buckets[bucket] = ptrs[u * 2];
            } else {
                ptrs[ptrs[u * 2 + 1] * 2] = ptrs[u * 2];
                if (ptrs[u * 2] != INVALID_NODE) { ptrs[ptrs[u * 2] * 2 + 1] = ptrs[u * 2 + 1]; }
            }
        }

        void incMaxBucket(int64_t bucket) { if (maxBucket < bucket) { maxBucket = bucket; } }
        bool empty() const { return maxBucket < minBucket; }

        void reset() {
            maxBucket = 0;
            minBucket = INT64_MAX;
        }
    };

    template<class Cap, class Term, class Flow>
    class IBFSGraph {
        int64_t n = 0;
        int64_t m = 0;
        Array<VertexData<Term> > vertices;
        Array<uint32_t> rows;
        Array<ArcData<Cap> > arcs;

        Flow flow = 0;
        int64_t augTimestamp = 0;
        int64_t topLevel[2] = {0, 0}; // [0]=sink, [1]=source
        FixedVector<uint32_t> active, tree_active[2]; // [0]=sink, [1]=source
        Buckets3Pass<Term> orphan3PassBuckets;
        BucketsOneSided<Term> orphanBuckets;
        ExcessBuckets<Term> excessBuckets;
        int64_t n_orphans[2] = {0, 0}; // [0]=sink, [1]=source

    public:
        IBFSGraph() = default;

        IBFSGraph(int64_t t_n, int64_t t_n_edges, MemoryStack &mem_stack) {
            n = t_n;
            m = t_n_edges * 2;

            uint64_t arcMemsize = align64(sizeof(ArcData<Cap>) * t_n_edges * 2);
            uint64_t nodeMemsize = align64(sizeof(VertexData<Term>) * (t_n + 1));
            uint64_t rowMemsize = align64(sizeof(uint32_t) * (t_n + 1));
            uint64_t ptrMemsize = align64(sizeof(uint32_t) * t_n * 2);
            uint64_t bucketMemsize = align64(sizeof(uint32_t) * (t_n + 1)) * 3;
            uint64_t activeMemsize = align64(sizeof(uint32_t) * t_n) * 3;
            uint64_t totalMemsize = arcMemsize + ptrMemsize + nodeMemsize + rowMemsize + bucketMemsize + activeMemsize;
            mem_stack.ensure(totalMemsize);
            mem_stack.clear();

            vertices.init((VertexData<Term> *) mem_stack.get_memory(nodeMemsize));
            std::memset((void *) &vertices[0], 0, nodeMemsize);
            rows.init((uint32_t *) mem_stack.get_memory(rowMemsize));
            arcs.init((ArcData<Cap> *) mem_stack.get_memory(arcMemsize));

            active.init((uint32_t *) mem_stack.get_memory(sizeof(uint32_t) * t_n));
            tree_active[0].init((uint32_t *) mem_stack.get_memory(sizeof(uint32_t) * t_n));
            tree_active[1].init((uint32_t *) mem_stack.get_memory(sizeof(uint32_t) * t_n));

            excessBuckets = ExcessBuckets<Term>((uint32_t *) mem_stack.get_memory(sizeof(uint32_t) * (t_n + 1)), (uint32_t *) mem_stack.get_memory(sizeof(uint32_t) * t_n * 2), t_n);
            orphan3PassBuckets = Buckets3Pass<Term>((uint32_t *) mem_stack.get_memory(sizeof(uint32_t) * (t_n + 1)), t_n);
            orphanBuckets = BucketsOneSided<Term>((uint32_t *) mem_stack.get_memory(sizeof(uint32_t) * (t_n + 1)), t_n);

            flow = 0;
            n_orphans[0] = 0;
            n_orphans[1] = 0;
            augTimestamp = 0;
        }

        void addNode(uint32_t u, Term cap_source, Term cap_sink, uint32_t deg) {
            if (u == 0) { rows[0] = 0; }
            vertices[u].parent = rows[u]; // arc cursor (temporary)
            rows[u + 1] = rows[u] + deg;
            flow += std::min(cap_source, cap_sink);
            vertices[u].excess = cap_source - cap_sink;

            // initNodes work
            vertices[u].firstSon = INVALID_NODE;
            vertices[u].nextNode = INVALID_NODE;
            vertices[u].prevSibling = INVALID_NODE;
            vertices[u].lastAugTimestamp = 0;
            vertices[u].isParentCurr = false;
            if (vertices[u].excess == 0) {
                vertices[u].label = 0;
            } else if (vertices[u].excess > 0) {
                vertices[u].label = 1;
                tree_active[1].push_back(u);
            } else {
                vertices[u].label = -1;
                tree_active[0].push_back(u);
            }
        }

        void addEdge(uint32_t u, uint32_t v, Cap cap, Cap rev_cap) {
            ArcData<Cap> &arc_u = arcs[vertices[u].parent];
            arc_u.rev = vertices[v].parent;
            arc_u.head = v;
            arc_u.rCap = cap;
            arc_u.setRevResidual(rev_cap != 0);

            ArcData<Cap> &arc_v = arcs[vertices[v].parent];
            arc_v.rev = vertices[u].parent;
            arc_v.head = u;
            arc_v.rCap = rev_cap;
            arc_v.setRevResidual(cap != 0);

            vertices[u].parent += 1;
            vertices[v].parent += 1;
        }

        // --- Direct two-pass construction API ---
        // Pass 1: call initNodes() then incDeg() for each arc endpoint
        // Pass 2: call finalizeRows() then addEdgeDirect/setTerminal, then finalize()

        void initNodes(int64_t t_n, int64_t t_n_edges, MemoryStack &t_mem_stack) {
            n = t_n;
            m = t_n_edges * 2;

            uint64_t arcMemsize = align64(sizeof(ArcData<Cap>) * t_n_edges * 2);
            uint64_t nodeMemsize = align64(sizeof(VertexData<Term>) * (t_n + 1));
            uint64_t rowMemsize = align64(sizeof(uint32_t) * (t_n + 1));
            uint64_t ptrMemsize = align64(sizeof(uint32_t) * t_n * 2);
            uint64_t bucketMemsize = align64(sizeof(uint32_t) * (t_n + 1)) * 3;
            uint64_t activeMemsize = align64(sizeof(uint32_t) * t_n) * 3;
            uint64_t totalMemsize = arcMemsize + ptrMemsize + nodeMemsize + rowMemsize + bucketMemsize + activeMemsize;
            t_mem_stack.ensure(totalMemsize);
            t_mem_stack.clear();

            vertices.init((VertexData<Term> *) t_mem_stack.get_memory(nodeMemsize));
            rows.init((uint32_t *) t_mem_stack.get_memory(rowMemsize));
            std::memset((void *) &rows[0], 0, rowMemsize);
            arcs.init((ArcData<Cap> *) t_mem_stack.get_memory(arcMemsize));

            active.init((uint32_t *) t_mem_stack.get_memory(sizeof(uint32_t) * t_n));
            tree_active[0].init((uint32_t *) t_mem_stack.get_memory(sizeof(uint32_t) * t_n));
            tree_active[1].init((uint32_t *) t_mem_stack.get_memory(sizeof(uint32_t) * t_n));

            excessBuckets = ExcessBuckets<Term>((uint32_t *) t_mem_stack.get_memory(sizeof(uint32_t) * (t_n + 1)), (uint32_t *) t_mem_stack.get_memory(sizeof(uint32_t) * t_n * 2), t_n);
            orphan3PassBuckets = Buckets3Pass<Term>((uint32_t *) t_mem_stack.get_memory(sizeof(uint32_t) * (t_n + 1)), t_n);
            orphanBuckets = BucketsOneSided<Term>((uint32_t *) t_mem_stack.get_memory(sizeof(uint32_t) * (t_n + 1)), t_n);

            flow = 0;
            n_orphans[0] = 0;
            n_orphans[1] = 0;
            augTimestamp = 0;
        }

        // Call during pass 1 for each undirected edge (u,v)
        void incDeg(uint32_t u, uint32_t v) {
            rows[u + 1]++;
            rows[v + 1]++;
        }

        // Set degree for node u directly (alternative to incDeg)
        void setDeg(uint32_t u, uint32_t deg) {
            rows[u + 1] = deg;
        }

        // Call after pass 1, before pass 2
        void finalizeRows() {
            for (uint32_t u = 0; u < (uint32_t) n; ++u) {
                rows[u + 1] += rows[u];
            }
            // Set parent as arc cursor for addEdge
            for (uint32_t u = 0; u < (uint32_t) n; ++u) {
                vertices[u].parent = rows[u];
            }
        }

        // Call during pass 2 for each undirected edge — same as addEdge
        void addEdgeDirect(uint32_t u, uint32_t v, Cap cap, Cap rev_cap) {
            ArcData<Cap> &arc_u = arcs[vertices[u].parent];
            arc_u.rev = vertices[v].parent;
            arc_u.head = v;
            arc_u.rCap = cap;
            arc_u.setRevResidual(rev_cap != 0);

            ArcData<Cap> &arc_v = arcs[vertices[v].parent];
            arc_v.rev = vertices[u].parent;
            arc_v.head = u;
            arc_v.rCap = rev_cap;
            arc_v.setRevResidual(cap != 0);

            vertices[u].parent += 1;
            vertices[v].parent += 1;
        }

        // Call during pass 2 for each node's terminal caps, after all edges for that node
        void setTerminal(uint32_t u, Term cap_source, Term cap_sink) {
            flow += std::min(cap_source, cap_sink);
            vertices[u].excess = cap_source - cap_sink;
            vertices[u].firstSon = INVALID_NODE;
            vertices[u].nextNode = INVALID_NODE;
            vertices[u].prevSibling = INVALID_NODE;
            vertices[u].lastAugTimestamp = 0;
            vertices[u].isParentCurr = false;
            if (vertices[u].excess == 0) {
                vertices[u].label = 0;
            } else if (vertices[u].excess > 0) {
                vertices[u].label = 1;
                tree_active[1].push_back(u);
            } else {
                vertices[u].label = -1;
                tree_active[0].push_back(u);
            }
        }

        Flow computeMaxFlow() {
            vertices[n].label = m;
            // Reset parent (was used as arc cursor during addEdge)
            for (uint32_t u = 0; u <= n; ++u) {
                vertices[u].parent = INVALID_ARC;
            }
            topLevel[0] = 1;
            topLevel[1] = 1;

            bool from_source = true;
            while (true) {
                std::swap(active, tree_active[from_source]);
                topLevel[from_source] += 1;
                if (from_source) {
                    growth<true>();
                } else {
                    growth<false>();
                }

                if (tree_active[1].empty() || tree_active[0].empty()) { break; }
                from_source = (n_orphans[1] < n_orphans[0]) || (n_orphans[1] == n_orphans[0] && !from_source);
            }
            return flow;
        }

        Flow getFlow() const noexcept { return flow; }
        size_t getNumNodes() const noexcept { return n; }
        size_t getNumArcs() const noexcept { return m; }
        int isNodeOnSrcSide(uint32_t node) const { return vertices[node].label > 0 ? 1 : 0; }
        Term getExcess(uint32_t node) const { return vertices[node].excess; }
        uint32_t arcBegin(uint32_t u) const { return rows[u]; }
        uint32_t arcEnd(uint32_t u) const { return rows[u + 1]; }
        uint32_t arcHead(uint32_t i) const { return arcs[i].head; }
        Cap arcResCap(uint32_t i) const { return arcs[i].rCap; }

    private:
        void augment(uint32_t bridge_idx) {
            ArcData<Cap> &bridge = arcs[bridge_idx];
            ArcData<Cap> &rev_bridge = rev_arc(bridge_idx);
            Cap bottleneck, bottleneckT, bottleneckS;
            bool forceBottleneck = false;
            bottleneck = bridge.rCap;
            bottleneckS = bridge.rCap;
            if (bottleneck != 1) {
                uint32_t u;
                for (u = rev_bridge.head; !vertices[u].excess; u = arcs[vertices[u].parent].head) {
                    Cap rev_cap = rev_arc(vertices[u].parent).rCap;
                    if (bottleneckS > rev_cap) { bottleneckS = rev_cap; }
                }

                VertexData<Term> &data = vertices[u];
                if (bottleneckS > data.excess) { bottleneckS = data.excess; }
                if (data.label != 1) { forceBottleneck = true; }
                if (u == rev_bridge.head) { bottleneck = bottleneckS; }
            }

            if (bottleneck != 1) {
                bottleneckT = bridge.rCap;
                uint32_t u;
                for (u = bridge.head; !vertices[u].excess; u = arcs[vertices[u].parent].head) {
                    const Cap cap = arcs[vertices[u].parent].rCap;
                    if (bottleneckT > cap) { bottleneckT = cap; }
                }
                VertexData<Term> &data = vertices[u];
                if (bottleneckT > -data.excess) { bottleneckT = -data.excess; }
                if (data.label != -1) { forceBottleneck = true; }
                if (u == bridge.head && bottleneck > bottleneckT) { bottleneck = bottleneckT; }
                if (forceBottleneck) { bottleneck = bottleneckS < bottleneckT ? bottleneckS : bottleneckT; }
            }

            rev_bridge.rCap += bottleneck;
            bridge.setRevResidual(true);
            bridge.rCap -= bottleneck;
            if (bridge.rCap == 0) { rev_bridge.setRevResidual(false); }

            flow -= bottleneck;
            if (bottleneck == 1 || forceBottleneck) {
                int64_t minOrphanLevel = augmentPath<false>(bridge.head, bottleneck);
                adoption<false>(minOrphanLevel, true);
            } else {
                int64_t minOrphanLevel = augmentExcess<false>(bridge.head, bottleneck);
                adoption<false>(minOrphanLevel, false);
                augmentExcesses<false>();
            }

            if (bottleneck == 1 || forceBottleneck) {
                int64_t minOrphanLevel = augmentPath<true>(rev_bridge.head, bottleneck);
                adoption<true>(minOrphanLevel, true);
            } else {
                int64_t minOrphanLevel = augmentExcess<true>(rev_bridge.head, bottleneck);
                adoption<true>(minOrphanLevel, false);
                augmentExcesses<true>();
            }
        }

        // @ret: minimum orphan level
        template<bool from_source>
        int64_t augmentPath(uint32_t u, Cap push) {
            int64_t orphanMinLevel = topLevel[from_source] + 1;
            augTimestamp += 1;
            while (vertices[u].excess == 0) {
                ArcData<Cap> &arc = arcs[vertices[u].parent];
                ArcData<Cap> &rev_arc = arcs[arc.getrev()];
                if (from_source) {
                    arc.rCap += push;
                    rev_arc.setRevResidual(true);
                    rev_arc.rCap -= push;
                } else {
                    rev_arc.rCap += push;
                    arc.setRevResidual(true);
                    arc.rCap -= push;
                }
                if ((from_source ? rev_arc.rCap : arc.rCap) == 0) {
                    if (from_source) {
                        arc.setRevResidual(false);
                    } else {
                        rev_arc.setRevResidual(false);
                    }
                    remove_sibling(u);
                    orphanMinLevel = from_source ? vertices[u].label : -vertices[u].label;

                    if (from_source) {
                        orphanBuckets.add_source(vertices.data(), u);
                    } else {
                        orphanBuckets.add_sink(vertices.data(), u);
                    }
                }
                u = arc.head;
            }

            VertexData<Term> &data = vertices[u];
            data.excess += from_source ? -push : push;
            if (data.excess == 0) {
                orphanMinLevel = from_source ? data.label : -data.label;
                if (from_source) {
                    orphanBuckets.add_source(vertices.data(), u);
                } else {
                    orphanBuckets.add_sink(vertices.data(), u);
                }
            }

            flow += push;
            return orphanMinLevel;
        }

        // @ret: minimum level in which created an orphan
        template<bool from_source>
        int64_t augmentExcess(uint32_t u, Cap push) {
            int64_t orphanMinLevel = topLevel[from_source] + 1;
            augTimestamp += 1;
            while (from_source ? vertices[u].excess <= 0 : vertices[u].excess >= 0) {
                VertexData<Term> &data = vertices[u];
                ArcData<Cap> &arc = arcs[data.parent];
                ArcData<Cap> &rev_arc = arcs[arc.getrev()];

                if (from_source ? rev_arc.rCap < push - data.excess : arc.rCap < data.excess + push) {
                    data.excess += from_source ? rev_arc.rCap - push : push - arc.rCap;
                    push = from_source ? rev_arc.rCap : arc.rCap;
                } else {
                    push += from_source ? -data.excess : data.excess;
                    data.excess = 0;
                }

                if (from_source) {
                    arc.rCap += push;
                    rev_arc.setRevResidual(true);
                    rev_arc.rCap -= push;
                } else {
                    rev_arc.rCap += push;
                    arc.setRevResidual(true);
                    arc.rCap -= push;
                }

                if ((from_source ? rev_arc.rCap : arc.rCap) == 0) {
                    if (from_source) {
                        arc.setRevResidual(false);
                    } else {
                        rev_arc.setRevResidual(false);
                    }

                    remove_sibling(u);
                    orphanMinLevel = from_source ? data.label : -data.label;
                    if (from_source) {
                        orphanBuckets.add_source(vertices.data(), u);
                    } else {
                        orphanBuckets.add_sink(vertices.data(), u);
                    }
                    if (data.excess) { excessBuckets.incMaxBucket(from_source ? data.label : -data.label); }
                }

                u = arc.head;
                if (from_source ? vertices[u].excess < 0 : vertices[u].excess > 0) {
                    if (from_source) {
                        excessBuckets.remove_source(vertices.data(), u);
                    } else {
                        excessBuckets.remove_sink(vertices.data(), u);
                    }
                }
            }

            VertexData<Term> &data = vertices[u];
            if (push <= (from_source ? data.excess : -data.excess)) {
                flow += push;
            } else {
                flow += from_source ? data.excess : -data.excess;
            }
            data.excess += from_source ? -push : push;

            if (from_source ? data.excess <= 0 : data.excess >= 0) {
                orphanMinLevel = from_source ? data.label : -data.label;
                if (from_source) {
                    orphanBuckets.add_source(vertices.data(), u);
                } else {
                    orphanBuckets.add_sink(vertices.data(), u);
                }

                if (data.excess) { excessBuckets.incMaxBucket(from_source ? data.label : -data.label); }
            }
            return orphanMinLevel;
        }

        template<bool from_source>
        void augmentExcesses() {
            int64_t adoptedUpToLevel = excessBuckets.maxBucket;
            if (!excessBuckets.empty()) {
                for (; excessBuckets.maxBucket != excessBuckets.minBucket - 1; --excessBuckets.maxBucket) {
                    uint32_t u;
                    while ((u = excessBuckets.popFront(excessBuckets.maxBucket)) != INVALID_NODE) {
                        int64_t minOrphanLevel = augmentExcess<from_source>(u, 0);
                        if (adoptedUpToLevel < minOrphanLevel) { minOrphanLevel = adoptedUpToLevel; }

                        adoption<from_source>(minOrphanLevel, false);
                        adoptedUpToLevel = excessBuckets.maxBucket;
                    }
                }
            }
            excessBuckets.reset();
            if (orphanBuckets.maxBucket != 0) {
                adoption<from_source>(adoptedUpToLevel + 1, true);
            }

            uint32_t u;
            while ((u = excessBuckets.popFront(0)) != INVALID_NODE) {
                orphanFree<from_source>(u);
            }
        }

        template<bool from_source>
        void adoption(int64_t fromLevel, bool toTop) {
            uint32_t i;
            int64_t minLabel;
            int64_t level;
            int64_t threePassLevel = 0;
            int64_t numOrphans = 0;
            int64_t numOrphansUniq = 0;

            for (level = fromLevel; level <= orphanBuckets.maxBucket && (toTop || threePassLevel || level <= excessBuckets.maxBucket); level++) {
                uint32_t u;
                while ((u = orphanBuckets.popFront(vertices.data(), level)) != INVALID_NODE) {
                    VertexData<Term> &data = vertices[u];
                    numOrphans += 1;
                    if (data.lastAugTimestamp != augTimestamp) {
                        data.lastAugTimestamp = augTimestamp;
                        n_orphans[from_source] += 1;

                        numOrphansUniq += 1;
                    }
                    if (threePassLevel == 0 && numOrphans >= 3 * numOrphansUniq) { threePassLevel = 1; }

                    if (data.isParentCurr) {
                        i = data.parent;
                    } else {
                        i = rows[u];
                        data.isParentCurr = true;
                    }
                    data.parent = INVALID_ARC;

                    if (data.label != (from_source ? 1 : -1)) {
                        minLabel = data.label - (from_source ? 1 : -1);
                        for (; i < rows[u + 1]; ++i) {
                            ArcData<Cap> &arc = arcs[i];
                            uint32_t v = arc.head;
                            VertexData<Term> &data_v = vertices[v];
                            if ((from_source ? arc.isRevResidual() : arc.rCap) != 0 && data_v.label == minLabel) {
                                data.parent = i;
                                add_sibling(u, v);
                                break;
                            }
                        }
                    }

                    if (data.parent != INVALID_ARC) {
                        if (data.excess) {
                            if (from_source) excessBuckets.add_source(vertices.data(), u);
                            else excessBuckets.add_sink(vertices.data(), u);
                        }
                        continue;
                    }

                    if (data.label == (from_source ? topLevel[1] : -topLevel[0])) {
                        orphanFree<from_source>(u);
                        continue;
                    }

                    uint32_t next;
                    for (uint32_t v = data.firstSon; v != INVALID_NODE; v = next) {
                        next = vertices[v].nextNode;
                        if (vertices[v].excess) {
                            if (from_source) {
                                excessBuckets.remove_source(vertices.data(), v);
                            } else {
                                excessBuckets.remove_sink(vertices.data(), v);
                            }
                        }
                        if (from_source) {
                            orphanBuckets.add_source(vertices.data(), v);
                        } else {
                            orphanBuckets.add_sink(vertices.data(), v);
                        }
                    }
                    data.firstSon = INVALID_NODE;

                    if (threePassLevel) {
                        data.label += from_source ? 1 : -1;
                        if (from_source) {
                            orphan3PassBuckets.add_source(vertices.data(), u);
                        } else {
                            orphan3PassBuckets.add_sink(vertices.data(), u);
                        }
                        if (threePassLevel == 1) { threePassLevel = level + 1; }
                        continue;
                    }

                    minLabel = from_source ? topLevel[1] : -topLevel[0];
                    if (data.label != minLabel) {
                        for (i = rows[u]; i < rows[u + 1]; ++i) {
                            const ArcData<Cap> &arc = arcs[i];
                            uint32_t v = arc.head;
                            const VertexData<Term> &data_v = vertices[v];
                            if ((from_source ? arc.isRevResidual() : arc.rCap) && (from_source ? (data_v.label > 0) : (data_v.label < 0)) && (from_source ? (data_v.label < minLabel) : (data_v.label > minLabel))) {
                                minLabel = data_v.label;
                                data.parent = i;
                                if (minLabel == data.label) { break; }
                            }
                        }
                    }

                    if (data.parent != INVALID_ARC) {
                        data.label = minLabel + (from_source ? 1 : -1);
                        add_sibling(u, arcs[data.parent].head);
                        if (from_source) {
                            if (data.label == topLevel[1]) { tree_active[1].push_back(u); }
                        } else {
                            if (data.label == -topLevel[0]) { tree_active[0].push_back(u); }
                        }
                        if (data.excess) {
                            if (from_source) {
                                excessBuckets.add_source(vertices.data(), u);
                            } else {
                                excessBuckets.add_sink(vertices.data(), u);
                            }
                        }
                    } else {
                        orphanFree<from_source>(u);
                    }
                }
            }
            if (level > orphanBuckets.maxBucket) { orphanBuckets.maxBucket = 0; }
            if (threePassLevel) { adoption3Pass<from_source>(threePassLevel); }
        }

        template<bool from_source>
        void adoption3Pass(int64_t minBucket) {
            int64_t minLabel;
            for (int64_t level = minBucket; level <= orphan3PassBuckets.maxBucket; level++) {
                uint32_t u;
                while ((u = orphan3PassBuckets.popFront(vertices.data(), level)) != INVALID_NODE) {
                    VertexData<Term> &data = vertices[u];
                    if (data.parent == INVALID_ARC) {
                        minLabel = from_source ? topLevel[1] : -topLevel[0];
                        int64_t destLabel = data.label - (from_source ? 1 : -1);
                        for (uint32_t i = rows[u]; i != rows[u + 1]; ++i) {
                            ArcData<Cap> &arc = arcs[i];
                            uint32_t v = arc.head;
                            VertexData<Term> &data_v = vertices[v];

                            if ((from_source ? arc.isRevResidual() : arc.rCap) && ((from_source ? data_v.excess > 0 : data_v.excess < 0) || data_v.parent != INVALID_ARC) && (from_source ? (data_v.label > 0) : (data_v.label < 0)) && (from_source ? (data_v.label < minLabel) : (data_v.label > minLabel))) {
                                data.parent = i;
                                minLabel = data_v.label;
                                if (minLabel == destLabel) { break; }
                            }
                        }

                        if (data.parent == INVALID_ARC) {
                            data.label = 0;
                            if (data.excess) {
                                if (from_source) {
                                    excessBuckets.add_source(vertices.data(), u);
                                } else {
                                    excessBuckets.add_sink(vertices.data(), u);
                                }
                            }
                            continue;
                        }

                        data.label = minLabel + (from_source ? 1 : -1);
                        if (data.label != (from_source ? level : -level)) {
                            if (from_source) {
                                orphan3PassBuckets.add_source(vertices.data(), u);
                            } else {
                                orphan3PassBuckets.add_sink(vertices.data(), u);
                            }
                            continue;
                        }
                    }

                    if (data.label != (from_source ? topLevel[1] : -topLevel[0])) {
                        minLabel = data.label + (from_source ? 1 : -1);
                        for (uint32_t i = rows[u]; i != rows[u + 1]; ++i) {
                            ArcData<Cap> &arc = arcs[i];
                            uint32_t v = arc.head;
                            VertexData<Term> &data_v = vertices[v];

                            if ((from_source ? arc.rCap : arc.isRevResidual()) && (data_v.label == 0 || (from_source ? (minLabel < data_v.label) : (minLabel > data_v.label)))) {
                                if (data_v.label != 0) {
                                    if (from_source) {
                                        orphan3PassBuckets.remove_source(vertices.data(), v);
                                    } else {
                                        orphan3PassBuckets.remove_sink(vertices.data(), v);
                                    }
                                } else if (data_v.excess) {
                                    if (from_source) {
                                        excessBuckets.remove_source(vertices.data(), v);
                                    } else {
                                        excessBuckets.remove_sink(vertices.data(), v);
                                    }
                                }
                                data_v.label = minLabel;
                                data_v.parent = arc.getrev();
                                if (from_source) {
                                    orphan3PassBuckets.add_source(vertices.data(), v);
                                } else {
                                    orphan3PassBuckets.add_sink(vertices.data(), v);
                                }
                            }
                        }
                    }

                    add_sibling(u, arcs[data.parent].head);
                    data.isParentCurr = false;
                    if (data.excess) {
                        if (from_source) {
                            excessBuckets.add_source(vertices.data(), u);
                        } else {
                            excessBuckets.add_sink(vertices.data(), u);
                        }
                    }
                    if (from_source) {
                        if (data.label == topLevel[1]) { tree_active[1].push_back(u); }
                    } else {
                        if (data.label == -topLevel[0]) { tree_active[0].push_back(u); }
                    }
                }
            }
            orphan3PassBuckets.maxBucket = 0;
        }

        template<bool from_source>
        void growth() {
            for (auto u: active) {
                VertexData<Term> &data = vertices[u];
                if (data.label != (from_source ? topLevel[1] - 1 : -(topLevel[0] - 1))) { continue; }

                for (uint32_t i = rows[u]; i < rows[u + 1]; ++i) {
                    ArcData<Cap> &arc = arcs[i];
                    if (!(from_source ? arc.rCap : arc.isRevResidual())) { continue; }

                    uint32_t v = arc.head;
                    VertexData<Term> &data_v = vertices[v];
                    if (MY_UNLIKELY(data_v.label == 0)) {
                        data_v.isParentCurr = false;
                        data_v.label = data.label + (from_source ? 1 : -1);
                        data_v.parent = arc.getrev();
                        add_sibling(v, u);
                        tree_active[from_source].push_back(v);
                    } else if (MY_UNLIKELY(from_source ? data_v.label < 0 : data_v.label > 0)) {
                        augment(from_source ? i : arc.getrev());
                        if (data.label != (from_source ? topLevel[1] - 1 : -(topLevel[0] - 1))) { break; }
                        if (from_source ? arc.rCap : arc.isRevResidual()) { i -= 1; }
                    }
                }
            }
            active.clear();
        }

        void remove_sibling(uint32_t u) {
            if (vertices[u].prevSibling != INVALID_NODE) {
                vertices[vertices[u].prevSibling].nextNode = vertices[u].nextNode;
            } else {
                parent_node(u).firstSon = vertices[u].nextNode;
            }

            if (vertices[u].nextNode != INVALID_NODE) {
                vertices[vertices[u].nextNode].prevSibling = vertices[u].prevSibling;
            }
        }

        void add_sibling(uint32_t u, uint32_t parent) {
            vertices[u].prevSibling = INVALID_NODE;
            vertices[u].nextNode = vertices[parent].firstSon;
            if (vertices[parent].firstSon != INVALID_NODE) {
                vertices[vertices[parent].firstSon].prevSibling = u;
            }
            vertices[parent].firstSon = u;
        }

        uint32_t parent_idx(uint32_t u) const { return arcs[vertices[u].parent].head; }
        VertexData<Term> &parent_node(uint32_t u) { return vertices[parent_idx(u)]; }

        ArcData<Cap> &rev_arc(uint32_t arc) { return arcs[arcs[arc].getrev()]; }
        const ArcData<Cap> &rev_arc(uint32_t arc) const { return arcs[arcs[arc].getrev()]; }

        template<bool from_source>
        void orphanFree(uint32_t u) {
            VertexData<Term> &data = vertices[u];
            if (data.excess) {
                data.label = from_source ? -topLevel[0] : topLevel[1];
                tree_active[!from_source].push_back(u);
                data.isParentCurr = false;
            } else {
                data.label = 0;
            }
        }
    };
} // namespace my_reimpls

#endif //FLOW_MY_EIBFS_I_H
