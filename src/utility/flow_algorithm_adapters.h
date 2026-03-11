#ifndef HEIPROMAP_FLOW_ALGORITHM_ADAPTERS_H
#define HEIPROMAP_FLOW_ALGORITHM_ADAPTERS_H

#include <vector>
#include <string>
#include <memory>

#include "flow_interface.h"

#include "../../extern/maxflow_algorithms/bk/graph.h"
#include "../../extern/maxflow_algorithms/hi_pr/hi_pr.h"

namespace HeiProMap {

template <typename captype = int, typename tcaptype = int, typename flowtype = int>
class BKAdapter : public IFlowAlgorithm<captype, tcaptype, flowtype> {
private:
    bk::Graph<captype, tcaptype, flowtype> *g;
    vertex_t n;
    vertex_t source;
    vertex_t target;

public:
    BKAdapter() : g(nullptr), n(0), source(0), target(0) {}

    ~BKAdapter() override {
        if (g) delete g;
    }

    void initialize(size_t t_n) override {
        n = t_n;
        if (g) delete g;
        g = new bk::Graph<captype, tcaptype, flowtype>(n + 2, (n + 2) * 4 + 1024);
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
            is_left[u] = g->what_segment(u) == bk::SOURCE;
        }
    }

    void build_residual_network(ResidualFlowNetwork &residual_g) override {
        residual_g.initialize(n);

        int n_edges = g->get_arc_num();
        typename bk::Graph<captype, tcaptype, flowtype>::arc_id arc = g->get_first_arc();

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

    void print() const override {}
};

template <typename captype = int, typename tcaptype = int, typename flowtype = int>
class HiPrAdapter : public IFlowAlgorithm<captype, tcaptype, flowtype> {
private:
    hi_pr::HiPr solver;
    vertex_t n;
    std::vector<int> edges;
    std::vector<captype> caps;
    std::vector<tcaptype> excess;

public:
    HiPrAdapter() : n(0) {}

    void initialize(size_t t_n) override {
        n = t_n;
        edges.clear();
        caps.clear();
        excess.assign(n + 2, 0);
    }

    void add(vertex_t u, vertex_t v, weight_t w) override {
        ASSERT(u < n);
        ASSERT(v < n);
        ASSERT(w >= 0);

        // one edge record, two capacities
        edges.push_back(u);
        edges.push_back(v);
        caps.push_back(w); // u -> v
        caps.push_back(w); // v -> u
    }

    void add_s_edge(vertex_t v, weight_t w) override {
        ASSERT(v < n);
        ASSERT(w >= 0);
        excess[v] += w;
    }

    void add_t_edge(vertex_t v, weight_t w) override {
        ASSERT(v < n);
        ASSERT(w >= 0);
        excess[v] -= w;
    }

    void solve() override {
        solver.construct(
            static_cast<unsigned int>(n),
            static_cast<unsigned int>(edges.size() / 2),
            edges.data(),
            caps.data(),
            excess.data()
        );
        solver.stageOne();
    }

    void get_cut(std::vector<u8> &is_left) override {
        is_left.resize(n);
        for (vertex_t u = 0; u < n; ++u) {
            is_left[u] = solver.is_weak_source(&solver.nodes[u]);
        }
    }

    void build_residual_network(ResidualFlowNetwork &residual_g) override {
        residual_g.initialize(n);
        for (vertex_t u = 0; u < n; ++u) {
            node *node_u = &solver.nodes[u];
            arc *stopA = (u + 1 < n)
                ? solver.nodes[u + 1].first
                : solver.arcs + 2 * solver.m;

            for (arc *a = node_u->first; a != stopA; ++a) {
                if (a->resCap > 0) {
                    vertex_t v = static_cast<vertex_t>(a->head - solver.nodes);
                    if (v < n) {
                        residual_g.add_directed_edge(u, v, a->resCap);
                    }
                }
            }
        }
    }

    void print() const override {}
};


} // namespace HeiProMap

#endif // HEIPROMAP_FLOW_ALGORITHM_ADAPTERS_H
