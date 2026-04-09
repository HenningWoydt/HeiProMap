#ifndef HEIPROMAP_FLOW_ALGORITHM_ADAPTERS_H
#define HEIPROMAP_FLOW_ALGORITHM_ADAPTERS_H

#include <vector>
#include <string>
#include <memory>
#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <cstdint>

#include "flow_interface.h"
#include "my_eibfs_i.h"

namespace HeiProMap {
    template<typename captype = int, typename tcaptype = int, typename flowtype = int>
    class EIBFSAdapter : public IFlowAlgorithm<captype, tcaptype, flowtype> {
    private:
        struct Edge {
            vertex_t u, v;
            weight_t w;
        };

        struct TermEdge {
            vertex_t v;
            weight_t w;
        };

        vertex_t n = 0;
        std::vector<Edge> edges;
        std::vector<TermEdge> s_edges, t_edges;
        std::vector<tcaptype> s_cap, t_cap; // per-node terminal caps
        std::vector<uint32_t> deg;

        my_reimpls::MemoryStack &mem_stack;
        std::unique_ptr<my_reimpls::IBFSGraph<captype, tcaptype, flowtype> > g;

    public:
        explicit EIBFSAdapter(my_reimpls::MemoryStack &t_mem_stack) : mem_stack(t_mem_stack) {}

        ~EIBFSAdapter() override = default;

        void initialize(size_t t_n) override {
            n = static_cast<vertex_t>(t_n);
            edges.clear();
            s_edges.clear();
            t_edges.clear();
            s_cap.assign(n, 0);
            t_cap.assign(n, 0);
            deg.assign(n, 0);
            g.reset();
        }

        void add(vertex_t u, vertex_t v, weight_t w) override {
            ASSERT(u < n && v < n && w >= 0);
            edges.push_back({u, v, w});
            deg[u]++;
            deg[v]++;
        }

        void add_s_edge(vertex_t v, weight_t w) override {
            ASSERT(v < n && w >= 0);
            s_cap[v] += static_cast<tcaptype>(w);
            s_edges.push_back({v, w});
        }

        void add_t_edge(vertex_t v, weight_t w) override {
            ASSERT(v < n && w >= 0);
            t_cap[v] += static_cast<tcaptype>(w);
            t_edges.push_back({v, w});
        }

        void solve() override {
            g = std::make_unique<my_reimpls::IBFSGraph<captype, tcaptype, flowtype> >(static_cast<int64_t>(n), static_cast<int64_t>(edges.size()), mem_stack);

            for (uint32_t u = 0; u < n; ++u) {
                g->addNode(u, s_cap[u], t_cap[u], deg[u]);
            }
            for (auto &e: edges) {
                g->addEdge(static_cast<uint32_t>(e.u), static_cast<uint32_t>(e.v), static_cast<captype>(e.w), static_cast<captype>(e.w));
            }
            g->computeMaxFlow();
        }

        flowtype get_flow_value() const { return static_cast<flowtype>(g->getFlow()); }

        void get_cut(std::vector<u8> &is_left) override {
            is_left.resize(n);
            for (vertex_t u = 0; u < n; ++u) {
                is_left[u] = g->isNodeOnSrcSide(static_cast<uint32_t>(u)) ? 1 : 0;
            }
        }

        void build_residual_network(ResidualFlowNetwork &residual_g) override {
            residual_g.initialize(n);

            // Internal arcs: iterate CSR and emit arcs with positive residual capacity
            for (uint32_t u = 0; u < n; ++u) {
                for (uint32_t i = g->arcBegin(u); i < g->arcEnd(u); ++i) {
                    uint32_t v = g->arcHead(i);
                    captype rc = g->arcResCap(i);
                    if (rc > 0 && v < n) {
                        residual_g.add_directed_edge(u, v, static_cast<weight_t>(rc));
                    }
                }
            }

            // Terminal arcs: reconstruct from original caps and final excess.
            // After max-flow, excess[u] = (s_cap[u] - t_cap[u]) - (flow_pushed_from_s[u] - flow_pushed_to_t[u])
            // Residual s->u = s_cap[u] - flow_from_s[u], residual u->s = flow_from_s[u]
            // Residual u->t = t_cap[u] - flow_to_t[u],   residual t->u = flow_to_t[u]
            //
            // Let orig_excess = s_cap - t_cap, final_excess = excess after solve.
            // net_flow_into_u_from_terminals = orig_excess - final_excess
            // For a source-side node: net_flow = flow_from_s (t_cap was 0 or fully cancelled)
            // We use: flow_from_s = max(0, min(s_cap, s_cap - (final_excess where final > 0 means leftover)))
            for (uint32_t u = 0; u < n; ++u) {
                tcaptype sc = s_cap[u];
                tcaptype tc = t_cap[u];
                tcaptype excess = g->getExcess(u);
                // Original excess = sc - tc. After solve, excess = (sc - tc) - pushed_net.
                // pushed_net = (sc - tc) - excess
                // The min(sc,tc) was already cancelled in addNode, so:
                // effective_s = sc - min(sc,tc), effective_t = tc - min(sc,tc)
                // But the residual on terminal arcs:
                tcaptype mn = std::min(sc, tc);
                tcaptype eff_s = sc - mn; // effective source cap
                tcaptype eff_t = tc - mn; // effective sink cap
                // After solve, excess tells us how much is left.
                // If eff_s > 0: flow_from_s = eff_s - max(0, excess)
                // If eff_t > 0: flow_to_t = eff_t - max(0, -excess)
                tcaptype flow_from_s = eff_s > 0 ? eff_s - std::max<tcaptype>(0, excess) : 0;
                tcaptype flow_to_t = eff_t > 0 ? eff_t - std::max<tcaptype>(0, -excess) : 0;

                // Residual s->u (remaining source capacity)
                tcaptype res_s_to_u = eff_s - flow_from_s;
                if (res_s_to_u > 0) residual_g.add_edge_from_source(u, static_cast<weight_t>(res_s_to_u));
                // Residual u->s (flow that can be returned)
                if (flow_from_s > 0) residual_g.add_edge_to_source(u, static_cast<weight_t>(flow_from_s));
                // Residual u->t (remaining sink capacity)
                tcaptype res_u_to_t = eff_t - flow_to_t;
                if (res_u_to_t > 0) residual_g.add_edge_to_target(u, static_cast<weight_t>(res_u_to_t));
                // Residual t->u (flow that can be returned)
                if (flow_to_t > 0) residual_g.add_edge_from_target(u, static_cast<weight_t>(flow_to_t));
            }
        }
    };
} // namespace HeiProMap

#endif // HEIPROMAP_FLOW_ALGORITHM_ADAPTERS_H
