#ifndef MT_RECPROMAP_GRAPH_H
#define MT_RECPROMAP_GRAPH_H

#include <vector>
#include <fstream>
#include <regex>
#include <numeric>
#include <random>

#include "../utility/definitions.h"
#include "../utility/utils.h"
#include "../utility/macros.h"

namespace SPM {

    /**
    * Standard undirected Graph that can hold vertex and edge weights.
    */
    class Graph {

    private:
        vertex_t m_n = 0; // number of vertices
        vertex_t m_m = 0; // number of edges

        // edges per vertex
        std::vector<vertex_t> m_v_edges;

        // vertex weight
        std::vector<weight_t> m_v_weights;
        weight_t m_g_weight = 0;

        // adjacency and edge weights
        std::vector<std::vector<EdgeW>> m_adj;

        // vertex active
        std::vector<u8> m_v_active;
        std::vector<vertex_t> m_active_vertices;

    public:
        /******************************
         GRAPH INITIALIZATION
         ******************************/

        /**
         * Reads in a graph, that is in METIS format.
         *
         * @param file_path Path to the file.
         */
        explicit Graph(const std::string &file_path) {
            if (!file_exists(file_path)) {
                std::cerr << "File " << file_path << " does not exist!" << std::endl;
                exit(EXIT_FAILURE);
            }

            bool has_v_weights = false;
            bool has_e_weights = false;
            vertex_t expected_edges = 0;

            std::ifstream file(file_path);
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

                    m_v_edges.resize(m_n, 0);
                    m_v_weights.resize(m_n, 1);
                    m_v_active.resize(m_n, 1);
                    m_active_vertices.resize(m_n);
                    std::iota(m_active_vertices.begin(), m_active_vertices.end(), 0);
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
                    if(has_v_weights){
                        set_vertex_weight(u, ints[idx]);
                        idx += 1;
                    }

                    // reserve enough space for the edges
                    reserve_space_for_edges(u, ints.size());

                    while(idx < ints.size()){
                        vertex_t v = (vertex_t) ints[idx] - 1;
                        idx += 1;

                        weight_t w = 1;
                        if(has_e_weights){
                            w = (weight_t) ints[idx];
                            idx += 1;
                        }
                        add_edge_1way(u,v,w);
                    }
                    sort_neighborhood(u);
                    u += 1;
                }

                if (expected_edges != m_m) {
                    std::cerr << "Number of expected edges " << expected_edges << " not equal to number edges " << m_m << " found!" << std::endl;
                    exit(EXIT_FAILURE);
                }

            } else {
                std::cerr << "Could not open file " << file_path << "!" << std::endl;
                exit(EXIT_FAILURE);
            }

            assert_graph();
        }

        /**
         * Constructor for the graph.
         *
         * @param n Number of vertices.
         * @param has_vertex_w Whether the graph has vertex weights.
         * @param has_edge_w Whether the graph has edge weights.
         */

        Graph() = default;

        explicit Graph(vertex_t n) {
            m_n = n;
            m_m = 0;

            m_v_edges.resize(m_n, 0);
            m_v_weights.resize(m_n, 1);
            m_v_active.resize(m_n, 1);
            m_g_weight = m_n;
            m_adj.resize(m_n);
        }

        /**
         * Sets the vertex weight of vertex u.
         *
         * @param u The vertex.
         * @param weight The weight.
         */
        void set_vertex_weight(vertex_t u, weight_t weight = 1) {
            ASSERT(u < m_n);
            m_g_weight = m_g_weight - m_v_weights[u] + weight;
            m_v_weights[u] = weight;
        }

        /**
         * Adds an edge.
         *
         * @param u Vertex u.
         * @param v Vertex v.
         * @param weight The weight of the edge.
         */
        void add_edge_2way(vertex_t u, vertex_t v, weight_t weight = 1) {
            ASSERT(u < m_n);
            ASSERT(v < m_n);
            ASSERT(weight > 0);
            ASSERT(!edge_exists(u, v));

            m_adj[u].emplace_back(v, weight);
            m_adj[v].emplace_back(u, weight);
            m_v_edges[u] += 1;
            m_v_edges[v] += 1;
            m_m += 1;

            ASSERT(m_adj[u].size() < m_n);
            ASSERT(m_adj[v].size() < m_n);
        }

        void add_edge_1way(vertex_t u, vertex_t v, weight_t weight = 1){
            ASSERT(u < m_n);
            ASSERT(v < m_n);
            ASSERT(weight > 0);
            ASSERT(!edge_exists(u, v));

            m_adj[u].emplace_back(v, weight);
            m_v_edges[u] += 1;

            ASSERT(m_adj[u].size() < m_n);
        }

        /**
         * Adds an edge.
         *
         * @param u Vertex u.
         * @param v Vertex v.
         * @param weight The weight of the edge.
         */
        void add_edge_2way_if_not_exist(vertex_t u, vertex_t v, weight_t weight = 1) {
            ASSERT(u < m_n);
            ASSERT(v < m_n);
            ASSERT(weight > 0);

            bool exists;
            if (m_adj[u].size() > m_adj[v].size()) {
                exists = edge_exists(v, u);
            } else {
                exists = edge_exists(u, v);
            }
            if (!exists) {
                add_edge_2way(u, v, weight);
            }
        }

        void reserve_space_for_edges(vertex_t u, size_t size){
            m_adj[u].reserve(size);
        }

        void sort_neighborhood(vertex_t u) {
            std::sort(m_adj[u].begin(), m_adj[u].end());
        }

        /******************************
         GRAPH UTILITY
         ******************************/

        /**
         * Returns the number of vertices.
         *
         * @return Number of vertices.
         */
        vertex_t get_n() const {
            return m_n;
        }

        /**
         * Returns the number of vertices.
         *
         * @return Number of vertices.
         */
        vertex_t get_n_active() const {
            vertex_t n = 0;
            for (u8 a: m_v_active) {
                n += a;
            }
            return n;
        }

        /**
         * Returns the number of edges.
         *
         * @return Number of edges.
         */
        vertex_t get_m() const {
            return m_m;
        }

        /**
         * Returns the number of active edges.
         *
         * @return Number of edges.
         */
        vertex_t get_m_active() const {
            vertex_t m = 0;
            for (vertex_t u = 0; u < m_n; ++u) {
                if (m_v_active[u] == 1) {
                    m += m_adj[u].size();
                }
            }
            return m;
        }

        /**
         * Determines if an edge exists between u and v in the graph.
         *
         * @param u Vertex u.
         * @param v Vertex v.
         * @return True if the edge exists, false else.
         */
        bool edge_exists(vertex_t u, vertex_t v) const {
            ASSERT(u < m_n);
            ASSERT(v < m_n);

            for (const EdgeW &e: m_adj[u]) {
                if (e.v == v) {
                    return true;
                }
            }
            return false;
        }

        /**
         * Determines if an edge exists between u and v in the graph.
         *
         * @param u Vertex u.
         * @param v Vertex v.
         * @return True if the edge exists, false else.
         */
        bool edge_exists_2way(vertex_t u, vertex_t v) const {
            ASSERT(u < m_n);
            ASSERT(v < m_n);

            return edge_exists(u, v) && edge_exists(v, u);
        }

        /**
         * Returns the weight of vertex u.
         *
         * @param u The vertex.
         * @return The weight.
         */
        vertex_t get_vertex_n_edge(vertex_t u) const {
            ASSERT(m_n == m_v_edges.size());
            ASSERT(u < m_n);
            return m_v_edges[u];
        }

        /**
         * Returns the weight of vertex u.
         *
         * @param u The vertex.
         * @return The weight.
         */
        vertex_t get_vertex_state(vertex_t u) const {
            ASSERT(m_n == m_v_edges.size());
            ASSERT(u < m_n);
            return m_v_active[u];
        }

        /**
         * Returns the weight of vertex u.
         *
         * @param u The vertex.
         * @return The weight.
         */
        void set_vertex_state(vertex_t u, u8 state) {
            ASSERT(m_n == m_v_edges.size());
            ASSERT(u < m_n);
            m_v_active[u] = state;
        }

        /**
         * Returns the weight of vertex u.
         *
         * @param u The vertex.
         * @return The weight.
         */
        weight_t get_vertex_weight(vertex_t u) const {
            ASSERT(m_n == m_v_weights.size());
            ASSERT(u < m_n);
            return m_v_weights[u];
        }

        /**
         * Returns the sum of vertex weights.
         *
         * @return The sum of all vertex weights.
         */
        weight_t get_sum_vertex_weights() const {
            weight_t s = 0;
            for (vertex_t u = 0; u < m_n; ++u) {
                if (m_v_active[u] == 1) {
                    s += m_v_weights[u];
                }
            }
            ASSERT(s == m_g_weight);
            return m_g_weight;
        }

        /**
         * Returns the weight of edge between vertices u and v. Undefined
         * behaviour if edge does not exist.
         *
         * @param u The vertex u.
         * @param v The vertex v.
         * @return The weight.
         */
        weight_t get_edge_weight(vertex_t u, vertex_t v) const {
            ASSERT(u < m_n);
            ASSERT(v < m_n);
            ASSERT(edge_exists(u, v));

            vertex_t min = std::min(u, v);
            vertex_t max = std::max(u, v);

            for (auto &e: m_adj[min]) {
                if (e.v == max) {
                    return e.w;
                }
            }
            // unreachable
            abort();
        }

        /**
         * Returns the sum of vertex weights.
         *
         * @return The sum of all vertex weights.
         */
        weight_t get_sum_edge_weights() const {
            weight_t edge_weights = 0;
            for (vertex_t u = 0; u < m_n; ++u) {
                for (const EdgeW &e: m_adj[u]) {
                    edge_weights += e.w;
                }
            }
            return edge_weights;
        }

        /**
         * Get the adjacency of vertex u.
         *
         * @param u The vertex.
         * @return Reference to the adjacency.
         */
        std::vector<EdgeW> &operator[](vertex_t u) {
            ASSERT(u < m_n);

            return m_adj[u];
        }

        /**
         * Get the adjacency of vertex u.
         *
         * @param u The vertex.
         * @return Reference to the adjacency.
         */
        const std::vector<EdgeW> &operator[](vertex_t u) const {
            ASSERT(u < m_n);

            return m_adj[u];
        }

        std::vector<u8> &get_active_states(){
            return m_v_active;
        }

        std::vector<vertex_t> &get_active_vertices(){
            return m_active_vertices;
        }

        /******************************
         GRAPH EDGE CONTRACTION
         ******************************/

        /**
         * Contracts the edge from u to v.
         *
         * @param u Vertex to keep.
         * @param v Vertex that will be removed.
         */
        void contract_edge(vertex_t u, vertex_t v) {
            // u has to be valid
            ASSERT(m_v_active[u] == 1);
            ASSERT(has_only_active_edges(u));
            ASSERT(no_duplicate_edges(u));
            ASSERT(no_self_loop(u));
            ASSERT(is_sorted(u));

            // v has to be valid
            ASSERT(m_v_active[v] == 1);
            ASSERT(has_only_active_edges(v));
            ASSERT(no_duplicate_edges(v));
            ASSERT(no_self_loop(v));
            ASSERT(is_sorted(v));

            ASSERT(edge_exists_2way(u, v));

            // disable second vertex
            m_v_active[v] = 0;

            // add weight of v to u
            m_v_weights[u] += m_v_weights[v];

            // remove v from all its neighbors
            for (const EdgeW &e: m_adj[v]) {
                remove_edge(e.v, v); // remove v from e.v

                ASSERT(no_duplicate_edges(e.v));
                ASSERT(has_only_active_edges(e.v));
                ASSERT(no_self_loop(e.v));
                // ASSERT(is_sorted(e.v));
            }

            // connect neighbors of v to u, but not u
            for (const EdgeW &e: m_adj[v]) {
                if (u != e.v) {
                    add_edge_with_weight(u, e.v, e.w);
                    add_edge_with_weight(e.v, u, e.w);

                    ASSERT(has_only_active_edges(u));
                    ASSERT(no_duplicate_edges(u));
                    ASSERT(no_self_loop(u));
                    // ASSERT(is_sorted(u));

                    ASSERT(has_only_active_edges(e.v));
                    ASSERT(no_duplicate_edges(e.v));
                    ASSERT(no_self_loop(e.v));
                    // ASSERT(is_sorted(e.v));
                }
            }

            // u has to be valid
            ASSERT(m_v_active[u] == 1);
            ASSERT(has_only_active_edges(u));
            ASSERT(no_duplicate_edges(u));
            ASSERT(no_self_loop(u));
            ASSERT(is_sorted(u));

            // v is inactive
            ASSERT(m_v_active[v] == 0);
        }

        /**
         * Removes v from u' adjacency.
         * @param u
         * @param v
         */
        void remove_edge(vertex_t u, vertex_t v) {
            ASSERT(u < m_n);
            ASSERT(v < m_n);
            ASSERT(edge_exists(u, v));

            size_t lower_idx = own_lower_bound_guaranteed(m_adj[u], v);
            ASSERT(m_adj[u][lower_idx].v == v);
            m_adj[u].erase(m_adj[u].cbegin() + lower_idx);
            return;

            size_t size = m_adj[u].size();
            size_t idx = size - 1;
            for (u64 i = 0; i < size - 1; ++i) {
                if(m_adj[u][i].v == v){
                    idx = i;
                    break;
                }
            }
            m_adj[u][idx] = m_adj[u][size - 1];
            m_adj[u].pop_back();
        }

        /**
         * Adds an edge.
         *
         * @param u Vertex u.
         * @param v Vertex v.
         * @param weight The weight of the edge.
         */
        void add_edge_with_weight(vertex_t u, vertex_t v, weight_t weight = 1) {
            ASSERT(u != v);
            ASSERT(u < m_n);
            ASSERT(v < m_n);
            ASSERT(weight > 0);

            size_t lower_idx = own_lower_bound_not_guaranteed(m_adj[u], v);
            if(lower_idx != m_adj[u].size() && m_adj[u][lower_idx].v == v){
                m_adj[u][lower_idx].w += weight;
            } else{
                m_adj[u].insert(m_adj[u].begin() + lower_idx, {v, weight});
            }
            return;

            for (EdgeW &e: m_adj[u]) {
                if (e.v == v) {
                    e.w += weight;
                    return;
                }
            }
            m_adj[u].emplace_back(v, weight);

            ASSERT(m_adj[u].size() < m_n);
            ASSERT(m_adj[v].size() < m_n);
        }

        /**
         * Uncontracts the edge from u to v.
         *
         * @param u Vertex that was kept.
         * @param v Vertex that was removed.
         */
        void uncontract_edge(vertex_t u, vertex_t v) {
            ASSERT(m_v_active[u] == 1);
            ASSERT(m_v_active[v] == 0);

            ASSERT(has_only_active_edges(u));
            ASSERT(no_duplicate_edges(u));

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

            // activate second vertex
            m_v_active[v] = 1;
            m_active_vertices.push_back(v);

            // subtract weight of v from u
            m_v_weights[u] -= m_v_weights[v];

            ASSERT(m_v_active[u] == 1);
            ASSERT(m_v_active[v] == 1);
            ASSERT(edge_exists_2way(u, v));
            ASSERT(has_only_active_edges(u));
            ASSERT(has_only_active_edges(v));
            ASSERT(no_duplicate_edges(u));
            ASSERT(no_duplicate_edges(v));
            for (const EdgeW &e: m_adj[v]) {
                ASSERT(has_only_active_edges(e.v));
                ASSERT(no_duplicate_edges(e.v));
            }
        }

        void add_edge_with_weight_guaranteed(vertex_t u, vertex_t v, weight_t weight = 1) {
            ASSERT(u < m_n);
            ASSERT(v < m_n);
            ASSERT(weight > 0);

            size_t lower_idx = own_lower_bound_not_guaranteed(m_adj[u], v);
            m_adj[u].insert(m_adj[u].begin() + lower_idx, {v, weight});
            return;

            m_adj[u].emplace_back(v, weight);

            ASSERT(m_adj[u].size() < m_n);
            ASSERT(m_adj[v].size() < m_n);
        }

        /**
         * Removes an edge.
         *
         * @param u Vertex u.
         * @param v Vertex v.
         * @param weight The weight of the edge.
         */
        void remove_edge_with_weight(vertex_t u, vertex_t v, weight_t weight = 1) {
            ASSERT(u < m_n);
            ASSERT(v < m_n);
            ASSERT(weight > 0);
            ASSERT(edge_exists(u, v));

            size_t lower_idx = own_lower_bound_guaranteed(m_adj[u], v);
            ASSERT(m_adj[u][lower_idx].v == v);
            m_adj[u][lower_idx].w -= weight;
            if(m_adj[u][lower_idx].w == 0){
                m_adj[u].erase(m_adj[u].begin() + lower_idx);
            }
            return;

            // remove v from u
            size_t size = m_adj[u].size();
            size_t idx = size - 1;
            for (u64 i = 0; i < size - 1; ++i) {
                if (m_adj[u][i].v == v) {
                    idx = i;
                    break;
                }
            }
            m_adj[u][idx].w -= weight;
            if (m_adj[u][idx].w == 0) {
                std::swap(m_adj[u][idx], m_adj[u][size - 1]);
                m_adj[u].pop_back();
            }

            ASSERT(m_adj[u].size() < m_n);
            ASSERT(m_adj[v].size() < m_n);
        }

        /******************************
         GRAPH MISC
         ******************************/

        /**
         * Used for asserting that the graph is correctly set.
         */
        void write_metis_graph(const std::string &file_path) const {
            std::stringstream ss;
            ss << m_n << " " << m_m << " 011" << std::endl;
            for (vertex_t u = 0; u < m_n; ++u) {
                ss << m_v_weights[u] << " ";
                for (auto &e: m_adj[u]) {
                    ss << e.v + 1 << " " << e.w << " ";
                }
                ss << std::endl;
            }

            std::ofstream file(file_path);
            file << ss.rdbuf();
            file.close();
        };

        Graph copy() {
            Graph g(m_n);

            for (vertex_t u = 0; u < m_n; ++u) {
                g.set_vertex_weight(u, get_vertex_weight(u));
                g.set_vertex_state(u, get_vertex_state(u));

                for (const EdgeW &e: m_adj[u]) {
                    g.add_edge_2way_if_not_exist(u, e.v, e.w);
                }
            }

            return g;
        }

        /******************************
         GRAPH ASSERTION
         ******************************/

        bool has_only_active_edges(vertex_t u) {
            for (const EdgeW &e: m_adj[u]) {
                if (m_v_active[e.v] == 0) {
                    return false;
                }
            }
            return true;
        }

        bool no_duplicate_edges(vertex_t u) {
            for (u64 i = 0; i < m_adj[u].size(); ++i) {
                for (u64 j = i + 1; j < m_adj[u].size(); ++j) {
                    if (m_adj[u][i].v == m_adj[u][j].v) {
                        return false;
                    }
                }
            }
            return true;
        }

        bool no_self_loop(vertex_t u) {
            for (const EdgeW &e: m_adj[u]) {
                if (e.v == u) {
                    return false;
                }
            }
            return true;
        }

        bool is_sorted(vertex_t u){
            return std::is_sorted(m_adj[u].begin(), m_adj[u].end());
        }

        bool assert_graph() {
            // for each active vertex, check that it is only connected to active vertices
            for (vertex_t u = 0; u < m_n; ++u) {
                if (m_v_active[u] == 1) {
                    // only active edges
                    if(!has_only_active_edges(u)){
                        return false;
                    }

                    // no duplicates
                    if(!no_duplicate_edges(u)){
                        return false;
                    }

                    // no self loop
                    if(!no_self_loop(u)){
                        return false;
                    }
                }
            }
            return true;
        }
    };

}

#endif //MT_RECPROMAP_GRAPH_H
