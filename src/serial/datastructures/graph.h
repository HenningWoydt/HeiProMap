#ifndef HEIDELBERGPROCESSMAPPING_GRAPH_H
#define HEIDELBERGPROCESSMAPPING_GRAPH_H

#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <sys/stat.h>

#include "../../definitions.h"
#include "../interfaces/ISerialGraph.h"
#include "../utility/utils.h"

namespace HeiProMap {
    /**
    * Standard undirected Graph that can hold vertex and edge weights.
    */
    class Graph final : public ISerialGraph {
        vertex_t n        = 0; // number of vertices
        vertex_t m        = 0; // number of edges
        weight_t g_weight = 0; // graph weight

        std::vector<weight_t> v_weights; // vertex weight

        std::vector<std::vector<EdgeVW>> m_adj; // adjacency and edge weights

    public:
        /**
         * Default constructor.
         */
        Graph() = default;

        /**
         * Initializes the graph via a file.
         *
         * @param graph_in Path to the file.
         */
        void initialize(const std::string& graph_in) override {
            std::ifstream file(graph_in);
            std::string line(64, ' ');

            if (file.is_open()) {
                bool has_v_weights = false;
                bool has_e_weights = false;
                std::vector<u64> ints;

                while (std::getline(file, line)) {
                    if (line[0] == '%') { continue; }
                    line_to_ints(line, ints);

                    n = ints[0];
                    m = ints[1];

                    v_weights.reserve(n);
                    g_weight = 0;
                    m_adj.reserve(n);

                    if (ints.size() == 3) {
                        if (ints[2] == 1) { has_e_weights = true; }
                        if (ints[2] == 10) { has_v_weights = true; }
                        if (ints[2] == 11) {
                            has_v_weights = true;
                            has_e_weights = true;
                        }
                    }

                    break;
                }

                // read in edges
                while (std::getline(file, line)) {
                    if (line[0] == '%') { continue; }
                    line_to_ints(line, ints);
                    size_t idx = 0;

                    // check vertex weight
                    weight_t u_w = 1;
                    if (has_v_weights) {
                        u_w = (weight_t)ints[idx];
                        idx += 1;
                    }
                    v_weights.push_back(u_w);
                    g_weight += u_w;

                    // reserve enough space for the edges
                    m_adj.emplace_back();
                    m_adj.back().reserve(ints.size());

                    while (idx < ints.size()) {
                        vertex_t v = (vertex_t)ints[idx] - 1;
                        idx += 1;

                        weight_t w = 1;
                        if (has_e_weights) {
                            w = (weight_t)ints[idx];
                            idx += 1;
                        }
                        m_adj.back().emplace_back(v, w);
                    }
                    std::sort(m_adj.back().begin(), m_adj.back().end());
                }
            } else {
                std::cerr << "Could not open file " << graph_in << "!" << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        /**
         * Gets the number of vertices.
         *
         * @return NUmber of vertices.
         */
        vertex_t get_n() const override { return n; }

        /**
         * Gets the number of edges.
         *
         * @return Number of edges.
         */
        vertex_t get_m() const override { return m; }

        /**
         * Gets the graph weight, e.g. the sum of all vertices.
         *
         * @return The graph weight.
         */
        weight_t get_weight() const override { return g_weight; }

        /**
         * Gets the weight of a vertex.
         *
         * @param u The vertex.
         * @return Weight of the vertex.
         */
        weight_t get_weight(const vertex_t u) const override { return v_weights[u]; }

        /**
         * Gets the number of edges of a vertex.
         *
         * @param u The vertex.
         * @return The number of adjacent edges.
         */
        size_t size(const vertex_t u) const override { return m_adj[u].size(); }

        /**
         * Gets the neighbor of a vertex.
         *
         * @param u The vertex.
         * @param idx The nth neighbor.
         * @return The neighboring vertex.
         */
        vertex_t neighbor(const vertex_t u, const size_t idx) const override { return m_adj[u][idx].v; }

        /**
         * Gets the weight of the edge.
         *
         * @param u The vertex.
         * @param idx The nth edge.
         * @return Weight of the nth edge.
         */
        weight_t get_weight(const vertex_t u, const size_t idx) const override { return m_adj[u][idx].w; }

        /**
         * Checks if an edge between vertices exists.
         *
         * @param u The first vertex.
         * @param v The second vertex.
         * @return True if an edge exists, false else.
         */
        bool edge_exists(const vertex_t u, const vertex_t v) const override { return std::any_of(m_adj[u].begin(), m_adj[u].end(), [&](const EdgeVW& e) { return e.v == v; }); }

        /**
         * Contracts the edge between two vertices.
         *
         * @param u The first vertex.
         * @param v The second vertex.
         */
        void contract(const vertex_t u, const vertex_t v) override {
            // add weight of v to u
            v_weights[u] += v_weights[v];

            // remove v from all its neighbors
            for (const EdgeVW& e : m_adj[v]) {
                remove_edge(e.v, v); // remove v from e.v
            }

            // connect neighbors of v to u, but not u
            for (const EdgeVW& e : m_adj[v]) {
                if (u != e.v) {
                    add_edge_with_weight(u, e.v, e.w);
                    add_edge_with_weight(e.v, u, e.w);
                }
            }
        }

        /**
         * Uncontracts an edge between two vertices.
         *
         * @param u The first vertex.
         * @param v The second vertex.
         */
        void uncontract(const vertex_t u, const vertex_t v) override {
            // remove neighbors of v from u
            for (const EdgeVW& e : m_adj[v]) {
                if (u != e.v) {
                    remove_edge_with_weight(u, e.v, e.w);
                    remove_edge_with_weight(e.v, u, e.w);
                }
            }

            // connect v to all its neighbors
            for (const EdgeVW& e : m_adj[v]) {
                add_edge_with_weight_guaranteed(e.v, v, e.w); // add v to e.v
            }

            // subtract weight of v from u
            v_weights[u] -= v_weights[v];
        }

    private:
        void remove_edge(const vertex_t u, const vertex_t v) {
            const size_t lower_idx = own_lower_bound_guaranteed(m_adj[u], v);
            m_adj[u].erase(m_adj[u].cbegin() + lower_idx);
        }

        void add_edge_with_weight(const vertex_t u, const vertex_t v, const weight_t weight = 1) {
            const size_t lower_idx = own_lower_bound_not_guaranteed(m_adj[u], v);
            if (lower_idx != m_adj[u].size() && m_adj[u][lower_idx].v == v) {
                m_adj[u][lower_idx].w += weight;
            } else {
                m_adj[u].insert(m_adj[u].begin() + lower_idx, {v, weight});
            }
        }

        void add_edge_with_weight_guaranteed(const vertex_t u, const vertex_t v, const weight_t weight = 1) {
            const size_t lower_idx = own_lower_bound_not_guaranteed(m_adj[u], v);
            m_adj[u].insert(m_adj[u].begin() + lower_idx, {v, weight});
        }

        void remove_edge_with_weight(const vertex_t u, const vertex_t v, const weight_t weight = 1) {
            const size_t lower_idx = own_lower_bound_guaranteed(m_adj[u], v);
            m_adj[u][lower_idx].w -= weight;
            if (m_adj[u][lower_idx].w == 0) {
                m_adj[u].erase(m_adj[u].begin() + lower_idx);
            }
        }
    };
}

#endif //HEIDELBERGPROCESSMAPPING_GRAPH_H
