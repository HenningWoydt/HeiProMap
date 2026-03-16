#ifndef HEIPROMAP_FLOW_ALGORITHM_ADAPTERS_H
#define HEIPROMAP_FLOW_ALGORITHM_ADAPTERS_H

#include <vector>
#include <string>
#include <memory>

#include "flow_interface.h"

#include "../../extern/maxflow_algorithms/nbk/graph.h"

namespace HeiProMap {

template <typename captype = int, typename tcaptype = int, typename flowtype = int>
class BKAdapter : public IFlowAlgorithm<captype, tcaptype, flowtype> {
private:
    nbk::Graph<captype, tcaptype, flowtype> *g;
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
        g = new nbk::Graph<captype, tcaptype, flowtype>(n + 2, (n + 2) * 4 + 1024);
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
            is_left[u] = g->what_segment(u) == nbk::SOURCE;
        }
    }

    void build_residual_network(ResidualFlowNetwork &residual_g) override {
        residual_g.initialize(n);

        int n_edges = g->get_arc_num();
        typename nbk::Graph<captype, tcaptype, flowtype>::arc_id arc = g->get_first_arc();

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

} // namespace HeiProMap

#endif // HEIPROMAP_FLOW_ALGORITHM_ADAPTERS_H
