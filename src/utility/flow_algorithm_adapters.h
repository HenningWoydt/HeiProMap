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
#include "my_hpf_hl.h"
#include "profiler.h"

namespace HeiProMap {
    template<typename captype = int, typename tcaptype = int, typename flowtype = int>
    class HPF_HLAdapter : public IFlowAlgorithm<captype, tcaptype, flowtype> {
        struct Edge {
            vertex_t u, v;
            weight_t w;
        };

        vertex_t n = 0;
        std::vector<Edge> edges;
        std::vector<tcaptype> s_cap, t_cap;

        MemoryStack &mem_stack;
        my_reimpls_hpf::Hpf<captype> g;

    public:
        explicit HPF_HLAdapter(MemoryStack &t_mem_stack) : mem_stack(t_mem_stack), g(0, 0, t_mem_stack) {
        }

        ~HPF_HLAdapter() override = default;

        void initialize(size_t t_n) override {
            n = t_n;
            edges.clear();
            s_cap.assign(n, 0);
            t_cap.assign(n, 0);
        }

        void add(vertex_t u, vertex_t v, weight_t w) override {
            ASSERT(u < n && v < n && w >= 0);
            edges.push_back({u, v, w});
        }

        void add_s_edge(vertex_t v, weight_t w) override {
            ASSERT(v < n && w >= 0);
            s_cap[v] += static_cast<tcaptype>(w);
        }

        void add_t_edge(vertex_t v, weight_t w) override {
            ASSERT(v < n && w >= 0);
            t_cap[v] += static_cast<tcaptype>(w);
        }

        void solve() override {
            size_t s_count = 0;
            size_t t_count = 0;
            for (vertex_t u = 0; u < n; ++u) {
                if (s_cap[u] > 0) s_count++;
                if (t_cap[u] > 0) t_count++;
            }

            {
                ScopedTimer _t("refinement", "FlowBasedRefinement", "solve_construct");
                g.reset(n + 2, edges.size() * 2 + s_count + t_count, mem_stack);
                g.add_node(n + 2);
                g.set_source(n);
                g.set_sink(n + 1);
            }

            {
                ScopedTimer _t("refinement", "FlowBasedRefinement", "solve_addEdge");
                for (const auto &e: edges) {
                    g.add_edge(e.u, e.v, static_cast<captype>(e.w));
                    g.add_edge(e.v, e.u, static_cast<captype>(e.w));
                }
                for (vertex_t u = 0; u < n; ++u) {
                    if (s_cap[u] > 0) g.add_edge(n, u, static_cast<captype>(s_cap[u]));
                    if (t_cap[u] > 0) g.add_edge(u, n + 1, static_cast<captype>(t_cap[u]));
                }
            }

            {
                ScopedTimer _t("refinement", "FlowBasedRefinement", "solve_maxflow");
                g.mincut();
            }
        }

        flowtype get_flow_value() const override {
            flowtype cut_value = 0;

            for (const auto &e: edges) {
                const bool u_on_source_side = g.what_label(e.u) == my_reimpls_hpf::Hpf<captype>::source_side;
                const bool v_on_source_side = g.what_label(e.v) == my_reimpls_hpf::Hpf<captype>::source_side;
                if (u_on_source_side != v_on_source_side) {
                    cut_value += static_cast<flowtype>(e.w);
                }
            }

            for (vertex_t u = 0; u < n; ++u) {
                if (g.what_label(u) == my_reimpls_hpf::Hpf<captype>::source_side) {
                    cut_value += static_cast<flowtype>(t_cap[u]);
                } else {
                    cut_value += static_cast<flowtype>(s_cap[u]);
                }
            }

            return cut_value;
        }

        void get_cut(std::vector<u8> &is_left) override {
            is_left.resize(n);
            for (vertex_t u = 0; u < n; ++u) {
                is_left[u] = (g.what_label(u) == my_reimpls_hpf::Hpf<captype>::source_side) ? 1 : 0;
            }
        }

        void build_residual_network(ResidualFlowNetwork &residual_g) override {
            residual_g.initialize(n);

            // HPF is a pseudoflow solver — arc flow values are NOT valid network flow.
            // Use cut labels to determine the residual network.
            // Cut edges (source→target side): saturated, residual only backward.
            // Same-side edges: conservatively assume residual in both directions.
            for (const auto &e : edges) {
                bool u_src = (g.what_label(e.u) == my_reimpls_hpf::Hpf<captype>::source_side);
                bool v_src = (g.what_label(e.v) == my_reimpls_hpf::Hpf<captype>::source_side);
                if (u_src == v_src) {
                    residual_g.add_directed_edge(e.u, e.v, 1);
                    residual_g.add_directed_edge(e.v, e.u, 1);
                } else if (u_src) {
                    // u on source side, v on target: saturated u→v, residual only v→u
                    residual_g.add_directed_edge(e.v, e.u, 1);
                } else {
                    // v on source side, u on target: saturated v→u, residual only u→v
                    residual_g.add_directed_edge(e.u, e.v, 1);
                }
            }

            // Terminal edges: s→u crosses cut if u on target side, u→t crosses cut if u on source side
            for (vertex_t u = 0; u < n; ++u) {
                bool u_src = (g.what_label(u) == my_reimpls_hpf::Hpf<captype>::source_side);
                if (s_cap[u] > 0) {
                    if (u_src) {
                        // s→u same side (both source), residual in both directions
                        residual_g.add_edge_from_source(u, 1);
                        residual_g.add_edge_to_source(u, 1);
                    } else {
                        // s→u crosses cut (saturated), residual only u→s
                        residual_g.add_edge_to_source(u, 1);
                    }
                }
                if (t_cap[u] > 0) {
                    if (!u_src) {
                        // u→t same side (both target), residual in both directions
                        residual_g.add_edge_to_target(u, 1);
                        residual_g.add_edge_from_target(u, 1);
                    } else {
                        // u→t crosses cut (saturated), residual only t→u
                        residual_g.add_edge_from_target(u, 1);
                    }
                }
            }
        }

        flowtype compute_cut_value(const std::vector<u8> &is_left) const override {
            flowtype val = 0;
            for (const auto &e : edges) {
                if (is_left[e.u] != is_left[e.v]) val += static_cast<flowtype>(e.w);
            }
            for (vertex_t u = 0; u < n; ++u) {
                if (is_left[u]) val += static_cast<flowtype>(t_cap[u]);
                else val += static_cast<flowtype>(s_cap[u]);
            }
            return val;
        }
    };

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

        MemoryStack &mem_stack;
        my_reimpls::IBFSGraph<captype, tcaptype, flowtype> g;

    public:
        explicit EIBFSAdapter(MemoryStack &t_mem_stack) : mem_stack(t_mem_stack) {
        }

        ~EIBFSAdapter() override = default;

        void initialize(size_t t_n) override {
            n = static_cast<vertex_t>(t_n);
            edges.clear();
            s_edges.clear();
            t_edges.clear();
            s_cap.assign(n, 0);
            t_cap.assign(n, 0);
            deg.assign(n, 0);
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
            {
                ScopedTimer _t("refinement", "FlowBasedRefinement", "solve_construct");
                g = my_reimpls::IBFSGraph<captype, tcaptype, flowtype>(static_cast<int64_t>(n), static_cast<int64_t>(edges.size()), mem_stack);
            }
            {
                ScopedTimer _t("refinement", "FlowBasedRefinement", "solve_addNode");
                for (uint32_t u = 0; u < n; ++u) {
                    g.addNode(u, s_cap[u], t_cap[u], deg[u]);
                }
            }
            {
                ScopedTimer _t("refinement", "FlowBasedRefinement", "solve_addEdge");
                for (auto &e: edges) {
                    g.addEdge(static_cast<uint32_t>(e.u), static_cast<uint32_t>(e.v), static_cast<captype>(e.w), static_cast<captype>(e.w));
                }
            }
            {
                ScopedTimer _t("refinement", "FlowBasedRefinement", "solve_maxflow");
                g.computeMaxFlow();
            }
        }

        flowtype get_flow_value() const { return static_cast<flowtype>(g.getFlow()); }

        void get_cut(std::vector<u8> &is_left) override {
            is_left.resize(n);
            for (vertex_t u = 0; u < n; ++u) {
                is_left[u] = g.isNodeOnSrcSide(static_cast<uint32_t>(u)) ? 1 : 0;
            }
        }

        void build_residual_network(ResidualFlowNetwork &residual_g) override {
            residual_g.initialize(n);

            // Internal arcs: iterate CSR and emit arcs with positive residual capacity
            for (uint32_t u = 0; u < n; ++u) {
                for (uint32_t i = g.arcBegin(u); i < g.arcEnd(u); ++i) {
                    uint32_t v = g.arcHead(i);
                    captype rc = g.arcResCap(i);
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
                tcaptype excess = g.getExcess(u);
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

        flowtype compute_cut_value(const std::vector<u8> &is_left) const override {
            flowtype val = 0;
            for (const auto &e : edges) {
                if (is_left[e.u] != is_left[e.v]) val += static_cast<flowtype>(e.w);
            }
            for (vertex_t u = 0; u < n; ++u) {
                if (is_left[u]) val += static_cast<flowtype>(t_cap[u]);
                else val += static_cast<flowtype>(s_cap[u]);
            }
            return val;
        }
    };
} // namespace HeiProMap

#endif // HEIPROMAP_FLOW_ALGORITHM_ADAPTERS_H
