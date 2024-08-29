#ifndef HEIDELBERGPROCESSMAPPING_GRAPH_H
#define HEIDELBERGPROCESSMAPPING_GRAPH_H

#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>
#include <iostream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "../../definitions.h"
#include "../utility/utils.h"
#include "../../macros.h"
#include "../../interfaces/IGraph.h"
#include "../interfaces/ISerialGraph.h"

namespace HeiProMap {

    /**
    * Standard undirected Graph that can hold vertex and edge weights.
    */
    class Graph : public ISerialGraph {

    private:
        vertex_t m_n = 0; // original number of vertices
        vertex_t m_m = 0; // original number of edges

        // vertex weight
        std::vector<weight_t> m_v_weights;

        weight_t m_g_weight = 0;

        // adjacency and edge weights
        std::vector<std::vector<EdgeVW>> m_adj;

    public:
        Graph() = default;

        // initialization
        void initialize(const std::string &graph_in) final {
            bool has_v_weights = false;
            bool has_e_weights = false;

            std::ifstream file(graph_in);
            std::string line(64, ' ');
            std::vector<u64> ints;
            if (file.is_open()) {
                while (std::getline(file, line)) {
                    if (line[0] == '%') { continue; }
                    line_to_ints(line, ints);

                    m_n = ints[0];
                    m_m = ints[1];

                    m_v_weights.reserve(m_n);
                    m_g_weight = 0;
                    m_adj.reserve(m_n);

                    if(ints.size() == 3){
                        if(ints[2] == 1){has_e_weights = true; }
                        if(ints[2] == 10){has_v_weights = true; }
                        if(ints[2] == 11){has_v_weights = true; has_e_weights = true; }
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
                        u_w = (weight_t) ints[idx];
                        idx += 1;
                    }
                    m_v_weights.push_back(u_w);
                    m_g_weight += u_w;

                    // reserve enough space for the edges
                    m_adj.emplace_back();
                    m_adj.back().reserve(ints.size());

                    while (idx < ints.size()) {
                        vertex_t v = (vertex_t) ints[idx] - 1;
                        idx += 1;

                        weight_t w = 1;
                        if (has_e_weights) {
                            w = (weight_t) ints[idx];
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

        // graph properties
        vertex_t get_n() const final { return m_n; }

        vertex_t get_m() const final { return m_m; }

        weight_t get_weight() const final { return m_g_weight; }

        // vertex weights
        weight_t get_weight(vertex_t u) const final { return m_v_weights[u]; }

        size_t size(vertex_t u) const final { return m_adj[u].size(); }

        vertex_t neighbor(vertex_t u, size_t idx) const final { return m_adj[u][idx].v; }
        weight_t get_weight(vertex_t u, size_t idx) const final { return m_adj[u][idx].w; }

        // edge manipulation
        bool edge_exists(vertex_t u, vertex_t v) const final { return std::any_of(m_adj[u].begin(), m_adj[u].end(), [&](const EdgeVW &e) { return e.v == v; }); }

        // coarsing and uncoarsing
        void contract(vertex_t u, vertex_t v) final {
            // add weight of v to u
            m_v_weights[u] += m_v_weights[v];

            // remove v from all its neighbors
            for (const EdgeVW &e: m_adj[v]) {
                remove_edge(e.v, v); // remove v from e.v
            }

            // connect neighbors of v to u, but not u
            for (const EdgeVW &e: m_adj[v]) {
                if (u != e.v) {
                    add_edge_with_weight(u, e.v, e.w);
                    add_edge_with_weight(e.v, u, e.w);
                }
            }
        }

        void uncontract(vertex_t u, vertex_t v) final {
            // remove neighbors of v from u
            for (const EdgeVW &e: m_adj[v]) {
                if (u != e.v) {
                    remove_edge_with_weight(u, e.v, e.w);
                    remove_edge_with_weight(e.v, u, e.w);
                }
            }

            // connect v to all its neighbors
            for (const EdgeVW &e: m_adj[v]) {
                add_edge_with_weight_guaranteed(e.v, v, e.w); // add v to e.v
            }

            // subtract weight of v from u
            m_v_weights[u] -= m_v_weights[v];
        }

    private:
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
