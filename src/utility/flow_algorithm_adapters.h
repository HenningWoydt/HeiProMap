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

#include "../../extern/maxflow_algorithms/nbk/graph.h"
#include "../../extern/maxflow_algorithms/ibfs/ibfs.h"

namespace HeiProMap {
    struct CanonicalEdge {
        vertex_t u;
        vertex_t v;
        weight_t w;

        bool operator<(const CanonicalEdge &other) const {
            if (u != other.u) return u < other.u;
            if (v != other.v) return v < other.v;
            return w < other.w;
        }
    };

    static uint64_t fnv1a_64(const std::string &data) {
        uint64_t hash = 14695981039346656037ull;
        for (unsigned char c: data) {
            hash ^= static_cast<uint64_t>(c);
            hash *= 1099511628211ull;
        }
        return hash;
    }

    static std::string to_hex(uint64_t x) {
        std::ostringstream oss;
        oss << std::hex << std::setw(16) << std::setfill('0') << x;
        return oss.str();
    }

    template<typename captype = int, typename tcaptype = int, typename flowtype = int>
    class BKAdapter : public IFlowAlgorithm<captype, tcaptype, flowtype> {
    private:
        nbk::Graph<captype, tcaptype, flowtype> *g;
        vertex_t n;
        vertex_t source;
        vertex_t target;

        struct Edge {
            vertex_t u;
            vertex_t v;
            weight_t w;
        };

        struct TerminalEdge {
            vertex_t v;
            weight_t w;
        };

        // Original graph edges between normal vertices
        std::vector<Edge> edges;

        // Optional: keep terminal connections too
        std::vector<TerminalEdge> s_edges;
        std::vector<TerminalEdge> t_edges;

    public:
        BKAdapter() : g(nullptr), n(0), source(0), target(0) {
        }

        ~BKAdapter() override {
            if (g) delete g;
        }

        void initialize(size_t t_n) override {
            n = static_cast<vertex_t>(t_n);

            if (g) delete g;
            g = new nbk::Graph<captype, tcaptype, flowtype>(n + 2, (n + 2) * 4 + 1024);
            g->add_node(n + 2);

            source = n;
            target = n + 1;

            edges.clear();
            s_edges.clear();
            t_edges.clear();
        }

        void add(vertex_t u, vertex_t v, weight_t w) override {
            ASSERT(u < n);
            ASSERT(v < n);
            ASSERT(w >= 0);

            g->add_edge(u, v, w, w);
            edges.push_back({u, v, w});
        }

        void add_s_edge(vertex_t v, weight_t w) override {
            ASSERT(v < n);
            ASSERT(w >= 0);

            g->add_edge(source, v, w, 0);
            s_edges.push_back({v, w});
        }

        void add_t_edge(vertex_t v, weight_t w) override {
            ASSERT(v < n);
            ASSERT(w >= 0);

            g->add_edge(v, target, w, 0);
            t_edges.push_back({v, w});
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
                is_left[u] = (g->what_segment(u) == nbk::SOURCE) ? 1 : 0;
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

        void save_graph(const std::string &filename) const {
            std::ofstream out(filename);
            if (!out) {
                throw std::runtime_error("Failed to open METIS output file: " + filename);
            }

            const vertex_t total_n = n + 2;
            const vertex_t source_id = n;
            const vertex_t target_id = n + 1;

            // adjacency list
            std::vector<std::vector<std::pair<vertex_t, weight_t> > > adj(total_n);

            // normal edges (undirected)
            for (const auto &e: edges) {
                adj[e.u].push_back({e.v, e.w});
                adj[e.v].push_back({e.u, e.w});
            }

            // source edges
            for (const auto &e: s_edges) {
                adj[source_id].push_back({e.v, e.w});
                adj[e.v].push_back({source_id, e.w});
            }

            // target edges
            for (const auto &e: t_edges) {
                adj[e.v].push_back({target_id, e.w});
                adj[target_id].push_back({e.v, e.w});
            }

            // count undirected edges
            size_t m = edges.size() + s_edges.size() + t_edges.size();

            // header: edge-weighted graph
            out << total_n << " " << m << " 001\n";

            // write adjacency (1-based indexing)
            for (vertex_t u = 0; u < total_n; ++u) {
                for (size_t i = 0; i < adj[u].size(); ++i) {
                    out << (adj[u][i].first + 1) << " " << adj[u][i].second;
                    if (i + 1 < adj[u].size()) out << " ";
                }
                out << "\n";
            }
        }

        void save_cut(const std::vector<u8> &is_left, const std::string &filename) const {
            std::ofstream out(filename);
            if (!out) {
                throw std::runtime_error("Failed to open cut output file: " + filename);
            }

            for (vertex_t u = 0; u < n; ++u) {
                out << (is_left[u] ? 1 : 0) << "\n";
            }
        }

        std::string canonical_graph_string() const {
            std::vector<CanonicalEdge> all_edges;
            all_edges.reserve(edges.size() + s_edges.size() + t_edges.size());

            // Normal undirected edges: normalize as (min(u,v), max(u,v))
            for (const auto &e: edges) {
                vertex_t a = std::min(e.u, e.v);
                vertex_t b = std::max(e.u, e.v);
                all_edges.push_back({a, b, e.w});
            }

            // Artificial source = n, target = n + 1
            for (const auto &e: s_edges) {
                all_edges.push_back({source, e.v, e.w});
            }

            for (const auto &e: t_edges) {
                all_edges.push_back({e.v, target, e.w});
            }

            std::sort(all_edges.begin(), all_edges.end());

            std::ostringstream oss;
            oss << "n=" << n << ";";
            oss << "source=" << source << ";";
            oss << "target=" << target << ";";

            for (const auto &e: all_edges) {
                oss << e.u << "," << e.v << "," << e.w << ";";
            }

            return oss.str();
        }

        std::string graph_hash() const {
            return to_hex(fnv1a_64(canonical_graph_string()));
        }

        void save_graph_and_cut_with_hash(const std::vector<u8> &is_left,
                                          const std::string &directory,
                                          const std::string &prefix = "graph") const {
            namespace fs = std::filesystem;

            fs::create_directories(directory);

            const std::string hash = graph_hash();

            const fs::path graph_file = fs::path(directory) / (prefix + "_" + hash + ".metis");
            const fs::path cut_file = fs::path(directory) / (prefix + "_" + hash + ".cut");

            save_graph(graph_file.string());
            save_cut(is_left, cut_file.string());
        }
    };
} // namespace HeiProMap

#endif // HEIPROMAP_FLOW_ALGORITHM_ADAPTERS_H
