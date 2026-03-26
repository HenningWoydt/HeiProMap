#ifndef HEIPROMAP_FLOW_ALGORITHM_ADAPTERS_H
#define HEIPROMAP_FLOW_ALGORITHM_ADAPTERS_H

#include <vector>
#include <string>
#include <memory>

#include "flow_interface.h"

#include "../../extern/maxflow_algorithms/reimpls/mbk.h"
#include "../../extern/maxflow_algorithms/ibfs/ibfs.h"

namespace HeiProMap {
    template<typename captype = int, typename tcaptype = int, typename flowtype = int>
    class BKAdapter : public IFlowAlgorithm<captype, tcaptype, flowtype> {
    private:
        reimpls::Graph<captype, tcaptype, flowtype> *g;
        vertex_t n;
        vertex_t source;
        vertex_t target;

    public:
        BKAdapter() : g(nullptr), n(0), source(0), target(0) {
        }

        ~BKAdapter() override {
            if (g) delete g;
        }

        void initialize(size_t t_n) override {
            n = t_n;
            if (g) delete g;
            g = new reimpls::Graph<captype, tcaptype, flowtype>(n + 2, (n + 2) * 4 + 1024);
            g->add_node(n + 2);
            source = n;
            target = n + 1;
        }

        void add(vertex_t u, vertex_t v, weight_t w) override {
            ASSERT(u < n);
            ASSERT(v < n);
            ASSERT(w >= 0);
            g->add_edge(u, v, w, w);
        }

        void add_s_edge(vertex_t v, weight_t w) override {
            ASSERT(v < n);
            ASSERT(w >= 0);
            g->add_edge(source, v, w, 0);
        }

        void add_t_edge(vertex_t v, weight_t w) override {
            ASSERT(v < n);
            ASSERT(w >= 0);
            g->add_edge(v, target, w, 0);
        }

        void solve() override {
            const int INF = std::numeric_limits<int>::max() / 2;
            g->add_tweights(source, INF, 0);
            g->add_tweights(target, 0, INF);
            g->maxflow();
        }

        void get_cut(std::vector<u8> &is_left) override {
            is_left.resize(n);
            for (vertex_t u = 0; u < n; ++u) {
                is_left[u] = g->what_segment(u) == reimpls::Graph<captype, tcaptype, flowtype>::SOURCE;
            }
        }

        void build_residual_network(ResidualFlowNetwork &residual_g) override {
            residual_g.initialize(n);

            int n_edges = g->get_arc_num();
            typename reimpls::Graph<captype, tcaptype, flowtype>::arc_id arc = g->get_first_arc();

            int u, v;
            for (int i = 0; i < n_edges; ++i) {
                g->get_arc_ends(arc, u, v);
                weight_t w = g->get_rcap(arc);

                if (w > 0) {
                    if (u == (int) source) {
                        residual_g.add_edge_from_source(v, w);
                    } else if (v == (int) source) {
                        residual_g.add_edge_to_source(u, w);
                    } else if (u == (int) target) {
                        residual_g.add_edge_from_target(v, w);
                    } else if (v == (int) target) {
                        residual_g.add_edge_to_target(u, w);
                    } else {
                        residual_g.add_directed_edge(u, v, w);
                    }
                }

                arc = g->get_next_arc(arc);
            }
        }

        void print() const override {
        }
    };

    template<typename captype = int64_t, typename tcaptype = int64_t, typename flowtype = int64_t>
    class IBFSAdapter : public IFlowAlgorithm<captype, tcaptype, flowtype> {
    private:
        struct Edge {
            vertex_t u;
            vertex_t v;
            captype cap_uv;
            captype cap_vu;
        };

        using GraphT = ibfs::IBFSGraph<captype, tcaptype, flowtype>;

        GraphT *g = nullptr;
        vertex_t n = 0;

        // Buffered problem instance
        std::vector<Edge> edges;
        std::vector<tcaptype> source_caps; // SOURCE -> v
        std::vector<tcaptype> sink_caps; // v -> SINK

        // Reuse bookkeeping
        bool capacity_initialized = false;
        size_t reserved_n = 0;
        size_t reserved_m = 0;

        bool built = false;
        bool solved = false;

    public:
        IBFSAdapter()
            : g(new GraphT(GraphT::IB_INIT_FAST)) {
        }

        ~IBFSAdapter() override {
            delete g;
        }

        void initialize(size_t t_n) override {
            n = static_cast<vertex_t>(t_n);

            edges.clear();
            source_caps.assign(n, 0);
            sink_caps.assign(n, 0);

            built = false;
            solved = false;
        }

        void add(vertex_t u, vertex_t v, weight_t w) override {
            ASSERT(u < n);
            ASSERT(v < n);
            ASSERT(u != v);
            ASSERT(w >= 0);

            edges.push_back({
                u,
                v,
                static_cast<captype>(w),
                static_cast<captype>(w)
            });
        }

        void add_s_edge(vertex_t v, weight_t w) override {
            ASSERT(v < n);
            ASSERT(w >= 0);
            source_caps[v] += static_cast<tcaptype>(w);
        }

        void add_t_edge(vertex_t v, weight_t w) override {
            ASSERT(v < n);
            ASSERT(w >= 0);
            sink_caps[v] += static_cast<tcaptype>(w);
        }

        void solve() override {
            build_if_needed();
            g->computeMaxFlow();
            solved = true;
        }

        void get_cut(std::vector<u8> &is_left) override {
            ASSERT(solved);

            is_left.resize(n);
            for (vertex_t u = 0; u < n; ++u) {
                is_left[u] = static_cast<u8>(g->isNodeOnSrcSide(static_cast<int64_t>(u), 0) == 1);
            }
        }

        void build_residual_network(ResidualFlowNetwork &residual_g) override {
            ASSERT(solved);

            residual_g.initialize(n);

            // Terminal residuals from node excess:
            //   excess > 0  => residual SOURCE -> node
            //   excess < 0  => residual node -> SINK
            for (vertex_t u = 0; u < n; ++u) {
                const auto &node = g->nodes[u];
                if (node.excess > 0) {
                    residual_g.add_edge_from_source(u, static_cast<weight_t>(node.excess));
                } else if (node.excess < 0) {
                    residual_g.add_edge_to_target(u, static_cast<weight_t>(-node.excess));
                }
            }

            // Residual ordinary arcs
            for (vertex_t u = 0; u < n; ++u) {
                auto *a = g->nodes[u].firstArc;
                auto *a_end = (g->nodes + (u + 1))->firstArc;
                for (; a != a_end; ++a) {
                    if (a->rCap > 0) {
                        vertex_t v = static_cast<vertex_t>(a->head - g->nodes);
                        residual_g.add_directed_edge(u, v, static_cast<weight_t>(a->rCap));
                    }
                }
            }
        }

        void print() const override {
        }

    private:
        void ensure_capacity() {
            const size_t m = edges.size();

            if (!capacity_initialized) {
                g->initSize(static_cast<int64_t>(n), static_cast<int64_t>(m));
                reserved_n = n;
                reserved_m = m;
                capacity_initialized = true;
                return;
            }

            // IBFS reset() only reuses already allocated storage for the original size.
            // If the new instance exceeds that size, we must recreate the graph.
            if (n > reserved_n || m > reserved_m) {
                delete g;
                g = new GraphT(GraphT::IB_INIT_FAST);
                g->initSize(static_cast<int64_t>(n), static_cast<int64_t>(m));
                reserved_n = n;
                reserved_m = m;
                capacity_initialized = true;
                return;
            }

            g->reset();
        }

        void build_if_needed() {
            if (built) return;

            ensure_capacity();

            for (vertex_t u = 0; u < n; ++u) {
                g->addNode(
                    static_cast<int64_t>(u),
                    source_caps[u],
                    sink_caps[u]
                );
            }

            for (const auto &e: edges) {
                g->addEdge(
                    static_cast<int64_t>(e.u),
                    static_cast<int64_t>(e.v),
                    e.cap_uv,
                    e.cap_vu
                );
            }

            g->initGraph();

            built = true;
        }
    };
} // namespace HeiProMap

#endif // HEIPROMAP_FLOW_ALGORITHM_ADAPTERS_H
