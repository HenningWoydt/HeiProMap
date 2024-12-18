#ifndef HEIDELBERGPROCESSMAPPING_GRAPH_CSR_H
#define HEIDELBERGPROCESSMAPPING_GRAPH_CSR_H

#include <fcntl.h>
#include <iostream>
#include <unistd.h>
#include <vector>
#include <sys/mman.h>
#include <sys/stat.h>

#include "arena_allocator.h"
#include "small_vector.h"
#include "../../definitions.h"
#include "../interfaces/ISerialGraph.h"
#include "../utility/utils.h"

namespace HeiProMap {
    /**
    * Standard undirected Graph that can hold vertex and edge weights.
    */
    class GraphCSR final : public ISerialGraph {
    private:
        vertex_t n = 0;
        vertex_t m = 0;

        weight_t vertex_weights = 0;
        weight_t edge_weights   = 0;

        std::vector<weight_t> v_weights;
        std::vector<size_t>   neighborhoods;
        std::vector<EdgeVW>   edges;
        size_t                curr_m = 0;


    public:
        GraphCSR() = default;

        void allocate(const vertex_t new_n, const vertex_t new_m) {
            n = new_n;
            m = new_m;
            v_weights.resize(n);
            neighborhoods.resize(n + 1);
            neighborhoods[0] = 0;
            edges.resize(m);
        }

        void set_weight(const vertex_t u, const weight_t w) {
            v_weights[u] = w;
        }

        void add_edge_on_stack(const vertex_t v, const weight_t w) {
            edges[curr_m++] = {v, w};
        }

        void add_weight_or_edge_on_stack(const vertex_t u, const vertex_t v, const weight_t w) {
            for (size_t i = neighborhoods[u]; i < neighborhoods[u + 1]; ++i) {
                if (edges[i].v == v) {
                    edges[i].w += w;
                    return;
                }
            }
            edges[curr_m++] = {v, w};
        }

        void set_size(const vertex_t u) {
            neighborhoods[u + 1] = curr_m;
        }


        // initialization
        void initialize(const std::string &graph_in) override {
            // Open the file
            int fd = open(graph_in.c_str(), O_RDONLY);
            if (fd == -1) {
                std::cerr << "File " << graph_in << " does not exist!" << std::endl;
                exit(EXIT_FAILURE);
            }

            // Get the file size
            struct stat fileInfo;
            if (fstat(fd, &fileInfo) == -1) {
                std::cerr << "File " << graph_in << " Could not get file size!" << std::endl;
                close(fd);
                exit(EXIT_FAILURE);
            }
            size_t file_size = fileInfo.st_size;

            // Memory-map the file
            char *file_arr = static_cast<char *>(mmap(nullptr, file_size, PROT_READ, MAP_PRIVATE, fd, 0));
            if (file_arr == MAP_FAILED) {
                std::cerr << "File " << graph_in << " Could not map the file!" << std::endl;
                close(fd);
                exit(EXIT_FAILURE);
            }

            size_t i = 0;
            while (file_arr[i] == '%') {
                // skip line, since it is a comment
                move_while_not(file_arr, i, '\n', file_size);
                ++i;
            }

            // skip leading white spaces
            move_while(file_arr, i, ' ', file_size);

            // read n
            n = 0;
            while (file_arr[i] != ' ') { n = n * 10 + (file_arr[i++] - '0'); }

            // skip whitespaces
            move_while(file_arr, i, ' ', file_size);

            // read m
            m = 0;
            while (file_arr[i] != ' ' && file_arr[i] != '\n') { m = m * 10 + (file_arr[i++] - '0'); }
            m *= 2;

            // skip whitespaces
            move_while(file_arr, i, ' ', file_size);

            // read fmt, since its optional special code
            // char fmt_0 = '0';
            char fmt_1 = '0';
            char fmt_2 = '0';
            while (file_arr[i] != '\n') {
                if (file_arr[i] != ' ') {
                    // found one fmt number
                    fmt_2 = file_arr[i++];

                    if (i < file_size && file_arr[i] != ' ' && file_arr[i] != '\n') {
                        fmt_1 = fmt_2;
                        fmt_2 = file_arr[i++];

                        if (i < file_size && file_arr[i] != ' ' && file_arr[i] != '\n') {
                            // fmt_0 = fmt_1;
                            fmt_1 = fmt_2;
                            fmt_2 = file_arr[i++];
                        }
                    }
                    break;
                }
                i += 1;
            }
            // now only ' ' expected until '\n'
            move_while(file_arr, i, ' ', file_size);
            ++i; // now on the next line

            // v_weights.resize(n, 1);
            neighborhoods.resize(n + 1);
            neighborhoods[0] = 0;
            edges.resize(m);

            vertex_t u = 0;
            if (fmt_1 == '0' && fmt_2 == '0') {
                v_weights.resize(n, 1);
                vertex_weights = (weight_t) n;
                while (true) {
                    if (file_arr[i] == '%') {
                        // this line is a comment, ignore it
                        while (file_arr[i] != '\n') { ++i; }
                        ++i;
                        continue;
                    }
                    // this line contains vertex information

                    while (file_arr[i] == ' ') { ++i; }

                    while (file_arr[i] != '\n') {
                        // read in the edges
                        vertex_t v = 0;
                        while (file_arr[i] != ' ' && file_arr[i] != '\n') { v = v * 10 + (file_arr[i++] - '0'); }
                        while (file_arr[i] == ' ') { ++i; }

                        edges[curr_m++] = {v, 1};
                    }

                    ++i;
                    neighborhoods[u + 1] = curr_m;
                    u += 1;

                    if (u + 32 >= n) {
                        break;
                    }
                }
                while (true) {
                    if (file_arr[i] == '%') {
                        // this line is a comment, ignore it
                        while (i < file_size && file_arr[i] != '\n') { ++i; }
                        // move_while_not(file_arr, i, '\n', file_size);
                        ++i;
                        continue;
                    }
                    // this line contains vertex information

                    while (i < file_size && file_arr[i] == ' ') { ++i; }

                    while (i < file_size && file_arr[i] != '\n') {
                        // read in the edges
                        vertex_t v = 0;
                        for (; i < file_size && file_arr[i] != ' ' && file_arr[i] != '\n'; ++i) { v = v * 10 + (file_arr[i] - '0'); }
                        while (i < file_size && file_arr[i] == ' ') { ++i; }

                        edges[curr_m++] = {v, 1};
                    }

                    ++i;
                    neighborhoods[u + 1] = curr_m;
                    u += 1;

                    if (u == n) {
                        break;
                    }
                }
            } else {
                v_weights.resize(n, 0);
                while (true) {
                    if (file_arr[i] == '%') {
                        // this line is a comment, ignore it
                        move_while_not(file_arr, i, '\n', file_size);
                        ++i;
                        continue;
                    }
                    // this line contains vertex information

                    move_while(file_arr, i, ' ', file_size); // skip leading whitespaces

                    weight_t u_w = 1;
                    if (fmt_1 == '1') {
                        // read in vertex weight
                        u_w = 0;
                        while (i < file_size && file_arr[i] != ' ' && file_arr[i] != '\n') { u_w = u_w * 10 + (file_arr[i++] - '0'); }
                        move_while(file_arr, i, ' ', file_size); // move to next number
                    }
                    v_weights[u] = u_w;
                    vertex_weights += u_w;

                    while (i < file_size && file_arr[i] != '\n') {
                        // read in the edges
                        vertex_t v = 0;
                        while (i < file_size && file_arr[i] != ' ' && file_arr[i] != '\n') { v = v * 10 + (file_arr[i++] - '0'); }
                        move_while(file_arr, i, ' ', file_size); // move to next number

                        weight_t w = 1;
                        if (fmt_2 == '1') {
                            w = 0;
                            while (i < file_size && file_arr[i] != ' ' && file_arr[i] != '\n') { w = w * 10 + (file_arr[i++] - '0'); }
                            move_while(file_arr, i, ' ', file_size); // move to next number
                        }

                        edges[curr_m++] = {v, w};
                    }

                    ++i;
                    neighborhoods[u + 1] = curr_m;
                    u += 1;

                    if (u == n) {
                        break;
                    }
                }
            }

            // Clean up
            munmap(file_arr, file_size);
            close(fd);
        }

        // graph properties
        vertex_t get_n() const override { return n; }

        vertex_t get_m() const override { return m; }
        vertex_t get_true_m() const { return curr_m; }

        weight_t get_weight() const override { return vertex_weights; }

        // vertex weights
        weight_t get_weight(const vertex_t u) const override { return v_weights[u]; }

        size_t size(const vertex_t u) const override { return neighborhoods[u + 1] - neighborhoods[u]; }

        vertex_t neighbor(const vertex_t u, const size_t idx) const override { return edges[neighborhoods[u] + idx].v; }

        weight_t get_weight(const vertex_t u, const size_t idx) const override { return edges[neighborhoods[u] + idx].w; }

        // edge manipulation
        bool edge_exists(const vertex_t u, const vertex_t v) const override {
            for (size_t i = neighborhoods[u]; i < neighborhoods[u + 1]; ++i) {
                if (edges[i].v == v) {
                    return true;
                }
            }
            return false;
        }

        // coarsening and uncoarsening
        void contract(const vertex_t u, const vertex_t v) override {
            return;
            /*
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
             */
        }

        void uncontract(const vertex_t u, const vertex_t v) override {
            return;
            /*
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
             */
        }

        GraphCSR contract_matching(const std::vector<EdgeUV> &matching) const {
            // allocate space for the new graph
            GraphCSR g;
            g.allocate(n, curr_m);

            // define the state of each vertex
            const u8              NOT_MATCHED    = 0;
            const u8              FIRST_MATCHED  = 1;
            const u8              SECOND_MATCHED = 2;
            std::vector<u8>       vertex_state(n, NOT_MATCHED);
            std::vector<vertex_t> vertex_neighbor(n);

            // check the matching
            for (const auto &[u, v]: matching) {
                vertex_state[u]    = FIRST_MATCHED;
                vertex_state[v]    = SECOND_MATCHED;
                vertex_neighbor[u] = v;
                vertex_neighbor[v] = u;
            }

            for (vertex_t u = 0; u < n; ++u) {
                if (vertex_state[u] == NOT_MATCHED) {
                    // copy it to the next graph
                    g.set_weight(u, v_weights[u]);
                    for (size_t i = neighborhoods[u]; i < neighborhoods[u + 1]; ++i) {
                        g.add_edge_on_stack(edges[i].v, edges[i].w);
                    }
                    g.set_size(u);
                } else if (vertex_state[u] == SECOND_MATCHED) {
                    // the vertex will not be present in the next graph
                    g.set_weight(u, 0); // no weight
                    g.set_size(u); // no edges
                } else {
                    // the vertex gets all neighbors of v
                    vertex_t v = vertex_neighbor[u];

                    g.set_weight(u, v_weights[u] + v_weights[v]); // add both weights
                    for (size_t i = neighborhoods[u]; i < neighborhoods[u + 1]; ++i) {
                        g.add_edge_on_stack(edges[i].v, edges[i].w); // add all edges of u
                    }
                    g.set_size(u);
                    for (size_t i = neighborhoods[v]; i < neighborhoods[v + 1]; ++i) {
                        g.add_weight_or_edge_on_stack(u, edges[i].v, edges[i].w); // add all edges of v
                    }
                    g.set_size(u);
                }
            }

            for (vertex_t u = 0; u < n; ++u) {
                if (vertex_state[u] == NOT_MATCHED) {
                    if(v_weights[u] != g.get_weight(u)){
                        std::cout << "a" << std::endl;
                    }
                    for (size_t i = neighborhoods[u]; i < neighborhoods[u + 1]; ++i) {
                        if(!g.edge_exists(u, edges[i].v)){
                            std::cout << "b" << std::endl;
                        }
                    }
                    if(g.size(u) != size(u)){
                        std::cout << "c" << std::endl;
                    }
                } else if (vertex_state[u] == SECOND_MATCHED) {
                    if(0 != g.get_weight(u)){
                        std::cout << "d" << std::endl;
                    }
                    if(g.size(u) != 0){
                        std::cout << "f" << std::endl;
                    }
                } else {
                    vertex_t v = vertex_neighbor[u];

                    if(v_weights[u] + v_weights[v] != g.get_weight(u)){
                        std::cout << "g" << std::endl;
                    }
                }
            }

            return g;
        }

    };
}

#endif //HEIDELBERGPROCESSMAPPING_GRAPH_CSR_H
