#ifndef HEIDELBERGPROCESSMAPPING_GRAPH_H
#define HEIDELBERGPROCESSMAPPING_GRAPH_H

#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../utility/definitions.h"
#include "../utility/utils.h"
#include "../utility/macros.h"
#include "../interfaces/IGraph.h"

namespace HeiProMap {

    /**
    * Standard undirected Graph that can hold vertex and edge weights.
    */
    class Graph : public IGraph {

    private:
        vertex_t m_n = 0; // original number of vertices
        vertex_t m_m = 0; // original number of edges

        // vertex weight
        std::vector<weight_t> m_v_weights;

        weight_t m_g_weight = 0;

        // adjacency and edge weights
        std::vector<std::vector<EdgeW>> m_adj;

    public:
        Graph() = default;

        // initialization
        void initialize(const std::string &graph_in, [[maybe_unused]] u64 n_threads) final {
            if (!file_exists(graph_in)) {
                std::cerr << "File " << graph_in << " does not exist!" << std::endl;
                exit(EXIT_FAILURE);
            }

            bool has_v_weights = false;
            bool has_e_weights = false;
            vertex_t expected_edges = 0;

            std::ifstream file(graph_in);
            std::string line(64, ' ');
            if (file.is_open()) {
                while (std::getline(file, line)) {
                    if (line[0] == '%') { continue; }

                    // remove leading and trailing whitespaces, replace double whitespaces
                    line.erase(0, line.find_first_not_of(' ')).erase(line.find_last_not_of(' ') + 1);
                    line = std::regex_replace(line, std::regex("\\s{2,}"), " ");

                    // read in header
                    std::vector<std::string> header = split(line, ' ');
                    m_n = std::stoi(header[0]);
                    m_m = 0;
                    expected_edges = std::stoi(header[1]);
                    m_m = expected_edges;

                    m_v_weights.resize(m_n, 1);
                    m_g_weight = m_n;
                    m_adj.resize(m_n);

                    std::string fmt = "000";
                    if (header.size() == 3 && header[2].size() == 3) {
                        fmt = header[2];
                    }

                    has_v_weights = fmt[1] == '1';
                    has_e_weights = fmt[2] == '1';

                    break;
                }

                // read in edges
                vertex_t u = 0;
                std::vector<u64> ints;

                while (std::getline(file, line)) {
                    if (line[0] == '%') { continue; }
                    line_to_ints(line, ints);
                    size_t idx = 0;

                    // check vertex weight
                    if (has_v_weights) {
                        set_weight(u, ints[idx]);
                        idx += 1;
                    }

                    // reserve enough space for the edges
                    reserve_space_for_edges(u, ints.size());

                    while (idx < ints.size()) {
                        vertex_t v = (vertex_t) ints[idx] - 1;
                        idx += 1;

                        weight_t w = 1;
                        if (has_e_weights) {
                            w = (weight_t) ints[idx];
                            idx += 1;
                        }
                        add_edge(u, v, w);
                    }
                    sort_neighborhood(u);
                    u += 1;
                }

                if (expected_edges != m_m) {
                    std::cerr << "Number of expected edges " << expected_edges << " not equal to number edges " << m_m << " found!" << std::endl;
                    exit(EXIT_FAILURE);
                }

            } else {
                std::cerr << "Could not open file " << graph_in << "!" << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        // graph properties
        vertex_t get_n() const final { return m_n; }

        vertex_t get_m() const final { return m_m; }

        weight_t get_weight() const final { return m_g_weight; }

        // vertex weights
        void set_weight(vertex_t u, weight_t weight) final {
            m_g_weight = m_g_weight - m_v_weights[u] + weight;
            m_v_weights[u] = weight;
        }

        weight_t get_weight(vertex_t u) const final { return m_v_weights[u]; }

        size_t n_neighbors(vertex_t u) const final { return m_adj[u].size(); }

        EdgeW &neighbor(vertex_t u, size_t idx) final { return m_adj[u][idx]; }

        const EdgeW &neighbor(vertex_t u, size_t idx) const final { return m_adj[u][idx]; }

        // edge manipulation
        void add_edge(vertex_t u, vertex_t v, weight_t w) final { m_adj[u].emplace_back(v, w); }

        bool edge_exists(vertex_t u, vertex_t v) const final { return std::any_of(m_adj[u].begin(), m_adj[u].end(), [&](const EdgeW &e) { return e.v == v; }); }

        // coarsing and uncoarsing
        void contract(vertex_t u, vertex_t v) final {
            // add weight of v to u
            m_v_weights[u] += m_v_weights[v];

            // remove v from all its neighbors
            for (const EdgeW &e: m_adj[v]) {
                remove_edge(e.v, v); // remove v from e.v
            }

            // connect neighbors of v to u, but not u
            for (const EdgeW &e: m_adj[v]) {
                if (u != e.v) {
                    add_edge_with_weight(u, e.v, e.w);
                    add_edge_with_weight(e.v, u, e.w);
                }
            }
        }

        void uncontract(vertex_t u, vertex_t v) final {
            // remove neighbors of v from u
            for (const EdgeW &e: m_adj[v]) {
                if (u != e.v) {
                    remove_edge_with_weight(u, e.v, e.w);
                    remove_edge_with_weight(e.v, u, e.w);
                }
            }

            // connect v to all its neighbors
            for (const EdgeW &e: m_adj[v]) {
                add_edge_with_weight_guaranteed(e.v, v, e.w); // add v to e.v
            }

            // subtract weight of v from u
            m_v_weights[u] -= m_v_weights[v];
        }

    private:
        // initialization
        void reserve_space_for_edges(vertex_t u, size_t size) { m_adj[u].reserve(size); }

        void sort_neighborhood(vertex_t u) { std::sort(m_adj[u].begin(), m_adj[u].end()); }

        // coarsing
        void remove_edge(vertex_t u, vertex_t v) {
            size_t lower_idx = own_lower_bound_guaranteed(m_adj[u], v);
            m_adj[u].erase(m_adj[u].cbegin() + lower_idx);
        }

        void add_edge_with_weight(vertex_t u, vertex_t v, weight_t weight = 1) {
            size_t lower_idx = own_lower_bound_not_guaranteed(m_adj[u], v);
            if (lower_idx != m_adj[u].size() && m_adj[u][lower_idx].v == v) {
                m_adj[u][lower_idx].w += weight;
            } else {
                m_adj[u].insert(m_adj[u].begin() + lower_idx, {v, weight});
            }
        }

        // uncoarsing
        void add_edge_with_weight_guaranteed(vertex_t u, vertex_t v, weight_t weight = 1) {
            size_t lower_idx = own_lower_bound_not_guaranteed(m_adj[u], v);
            m_adj[u].insert(m_adj[u].begin() + lower_idx, {v, weight});
        }

        void remove_edge_with_weight(vertex_t u, vertex_t v, weight_t weight = 1) {
            size_t lower_idx = own_lower_bound_guaranteed(m_adj[u], v);
            m_adj[u][lower_idx].w -= weight;
            if (m_adj[u][lower_idx].w == 0) {
                m_adj[u].erase(m_adj[u].begin() + lower_idx);
            }
        }

    };

}

#endif //HEIDELBERGPROCESSMAPPING_GRAPH_H
