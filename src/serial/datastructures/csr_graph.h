#ifndef HEIDELBERGPROCESSMAPPING_CSR_GRAPH_H
#define HEIDELBERGPROCESSMAPPING_CSR_GRAPH_H

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
#include "arena_allocator.h"
#include "small_vector.h"
#include "../interfaces/ISerialGraph.h"

namespace HeiProMap {

    /**
    * Standard undirected Graph that can hold vertex and edge weights.
    */
    class CSRGraph : public ISerialGraph {

    private:
        vertex_t m_n = 0; // original number of vertices
        vertex_t m_m = 0; // original number of edges

        // vertex weight
        std::vector<weight_t> m_v_weights;

        weight_t m_g_weight = 0;

        // Arena memory allocator
        ArenaAllocator<EdgeVW> arena_allocator;

        // adjacency and edge weights
        std::vector<SmallVector<EdgeVW>> m_adj;

    public:
        CSRGraph() = default;

        // initialization
        void initialize(const std::string &graph_in) final {
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
            m_n = 0;
            while (file_arr[i] != ' ') { m_n = m_n * 10 + (file_arr[i++] - '0'); }

            // skip whitespaces
            move_while(file_arr, i, ' ', file_size);

            // read m
            m_m = 0;
            while (file_arr[i] != ' ' && file_arr[i] != '\n') { m_m = m_m * 10 + (file_arr[i++] - '0'); }

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

            arena_allocator.initialize(m_m * 2);
            EdgeVW *curr_ptr = arena_allocator.request(m_m * 2);
            m_v_weights.reserve(m_n);
            m_adj.reserve(m_n);

            vertex_t u = 0;
            if(fmt_1 == '0' && fmt_2 == '0'){
                m_v_weights.resize(m_n, 1);
                m_g_weight = (weight_t) m_n;
                while (true) {
                    if (file_arr[i] == '%') {
                        // this line is a comment, ignore it
                        while(file_arr[i] != '\n'){ ++i; }
                        ++i;
                        continue;
                    }
                    // this line contains vertex information

                    while(file_arr[i] == ' '){ ++i; }

                    EdgeVW *start_ptr = curr_ptr;
                    while (file_arr[i] != '\n') {
                        // read in the edges
                        vertex_t v = 0;
                        while(file_arr[i] != ' ' && file_arr[i] != '\n') { v = v * 10 + (file_arr[i++] - '0'); }
                        while(file_arr[i] == ' '){ ++i; }

                        curr_ptr->v = v - 1;
                        curr_ptr->w = 1;
                        curr_ptr += 1;
                    }
                    m_adj.emplace_back(start_ptr, curr_ptr);
                    std::sort(start_ptr, curr_ptr);

                    ++i;
                    u += 1;

                    if (u+32 >= m_n) {
                        break;
                    }
                }
                while (true) {
                    if (file_arr[i] == '%') {
                        // this line is a comment, ignore it
                        while(i < file_size && file_arr[i] != '\n'){ ++i; }
                        // move_while_not(file_arr, i, '\n', file_size);
                        ++i;
                        continue;
                    }
                    // this line contains vertex information

                    while(i < file_size && file_arr[i] == ' '){ ++i; }

                    EdgeVW *start_ptr = curr_ptr;
                    while (i < file_size && file_arr[i] != '\n') {
                        // read in the edges
                        vertex_t v = 0;
                        for(;i < file_size && file_arr[i] != ' ' && file_arr[i] != '\n';++i){v = v * 10 + (file_arr[i] - '0');}
                        while(i < file_size && file_arr[i] == ' '){ ++i; }

                        curr_ptr->v = v - 1;
                        curr_ptr->w = 1;
                        curr_ptr += 1;
                    }
                    m_adj.emplace_back(start_ptr, curr_ptr);
                    std::sort(start_ptr, curr_ptr);

                    ++i;
                    u += 1;

                    if (u == m_n) {
                        break;
                    }
                }
            } else {
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
                    m_v_weights.push_back(u_w);
                    m_g_weight += u_w;

                    EdgeVW *start_ptr = curr_ptr;
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

                        curr_ptr->v = v - 1;
                        curr_ptr->w = 1;
                        curr_ptr += 1;
                    }
                    m_adj.emplace_back(start_ptr, curr_ptr);
                    std::sort(start_ptr, curr_ptr);

                    ++i;
                    u += 1;

                    if (u == m_n) {
                        break;
                    }
                }
            }

            // Clean up
            munmap(file_arr, file_size);
            close(fd);
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
        bool edge_exists(vertex_t u, vertex_t v) const final {
            for(size_t i = 0; i < m_adj[u].size(); ++i){
                if(m_adj[u][i].v == v){ return true; }
            }
            return false;
        }

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
            m_adj[u].erase(lower_idx);
        }

        void add_edge_with_weight(vertex_t u, vertex_t v, weight_t weight = 1) {
            size_t lower_idx = own_lower_bound_not_guaranteed(m_adj[u], v);
            if (lower_idx != m_adj[u].size() && m_adj[u][lower_idx].v == v) {
                m_adj[u][lower_idx].w += weight;
            } else {
                m_adj[u].insert(lower_idx, {v, weight}, arena_allocator);
            }
        }

        // uncoarsing
        void add_edge_with_weight_guaranteed(vertex_t u, vertex_t v, weight_t weight = 1) {
            size_t lower_idx = own_lower_bound_not_guaranteed(m_adj[u], v);
            m_adj[u].insert(lower_idx, {v, weight}, arena_allocator);
        }

        void remove_edge_with_weight(vertex_t u, vertex_t v, weight_t weight = 1) {
            size_t lower_idx = own_lower_bound_guaranteed(m_adj[u], v);
            m_adj[u][lower_idx].w -= weight;
            if (m_adj[u][lower_idx].w == 0) {
                m_adj[u].erase(lower_idx);
            }
        }

    };

}

#endif //HEIDELBERGPROCESSMAPPING_CSR_GRAPH_H
