#ifndef HEIDELBERGPROCESSMAPPING_GRAPH_CSR_H
#define HEIDELBERGPROCESSMAPPING_GRAPH_CSR_H

#include <fcntl.h>
#include <iostream>
#include <unistd.h>
#include <vector>
#include <sys/mman.h>
#include <sys/stat.h>

#include "../../definitions.h"
#include "../interfaces/ISerialGraph.h"
#include "../utility/utils.h"

namespace HeiProMap {
    /**
    * Standard undirected Graph that can hold vertex and edge weights.
    */
    class GraphCSR final : public ISerialGraph {
        vertex_t n = 0;
        vertex_t m = 0;

        weight_t vertex_weights = 0;
        weight_t edge_weights   = 0;

        std::vector<weight_t> v_weights;
        std::vector<size_t> neighborhoods;
        std::vector<EdgeVW> edges;
        size_t curr_m = 0;

    public:
        explicit GraphCSR(const std::string& graph_in) {
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
            char* file_arr = static_cast<char*>(mmap(nullptr, file_size, PROT_READ, MAP_PRIVATE, fd, 0));
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
                vertex_weights = (weight_t)n;
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

                        edges[curr_m++] = {v - 1, 1};
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

                        edges[curr_m++] = {v - 1, 1};
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

                        edges[curr_m++] = {v - 1, w};
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

        GraphCSR(const GraphCSR& g, const std::vector<EdgeUV>& matching) {
            n      = g.get_n();
            m      = g.get_true_m();
            curr_m = 0;
            v_weights.resize(n);
            neighborhoods.resize(n + 1);
            neighborhoods[0] = 0;
            edges.resize(m + 1);
            vertex_weights = g.vertex_weights;

            // define the state of each vertex
            constexpr u8 NOT_MATCHED    = 0;
            constexpr u8 FIRST_MATCHED  = 1;
            constexpr u8 SECOND_MATCHED = 2;
            std::vector<u8> vertex_state(n, NOT_MATCHED);
            std::vector<vertex_t> vertex_neighbor(n);

            // check the matching
            for (const auto& [u, v] : matching) {
                vertex_state[u]    = FIRST_MATCHED;
                vertex_state[v]    = SECOND_MATCHED;
                vertex_neighbor[u] = v;
                vertex_neighbor[v] = u;

                v_weights[v] = 0;
                v_weights[u] = g.get_weight(u) + g.get_weight(v);
            }

            for (vertex_t u = 0; u < n; ++u) {
                if (vertex_state[u] == NOT_MATCHED) {
                    // copy it to the next graph
                    v_weights[u] = g.get_weight(u);
                    for (size_t i = 0; i < g.size(u); ++i) {
                        vertex_t vv = g.neighbor(u, i);
                        weight_t ww = g.get_weight(u, i);

                        // if the vv vertex is matched, then make an edge to the neighbor vertex
                        if (vertex_state[vv] == SECOND_MATCHED) { vv = vertex_neighbor[vv]; }

                        edges[curr_m].v   = vv;
                        edges[curr_m++].w = ww;

                        // if the edge is present then add the weight, else expand it
                        for (size_t j = neighborhoods[u]; j < curr_m - 1; ++j) {
                            if (edges[j].v == vv) {
                                edges[j].w += ww;
                                curr_m -= 1;
                                break;
                            }
                        }
                    }
                } else if (vertex_state[u] == FIRST_MATCHED) {
                    // the vertex gets all neighbors of v
                    vertex_t v = vertex_neighbor[u];

                    // v_weights[u] = g.get_weight(u) + g.get_weight(v);
                    for (size_t i = 0; i < g.size(u); ++i) {
                        vertex_t vv = g.neighbor(u, i);
                        weight_t ww = g.get_weight(u, i);

                        // do not add edge to matched vertex
                        if (vv == v) { continue; }

                        // if the vv vertex is matched, then make an edge to the neighbor vertex
                        if (vertex_state[vv] == SECOND_MATCHED) { vv = vertex_neighbor[vv]; }

                        edges[curr_m].v   = vv;
                        edges[curr_m++].w = ww;

                        // if the edge is present then add the weight, else expand it
                        for (size_t j = neighborhoods[u]; j < curr_m - 1; ++j) {
                            if (edges[j].v == vv) {
                                edges[j].w += ww;
                                curr_m -= 1;
                                break;
                            }
                        }
                    }
                    for (size_t i = 0; i < g.size(v); ++i) {
                        vertex_t vv = g.neighbor(v, i);
                        weight_t ww = g.get_weight(v, i);

                        // do not add edge to matched vertex
                        if (vv == u) { continue; }

                        // if the vv vertex is matched, then make an edge to the neighbor vertex
                        if (vertex_state[vv] == SECOND_MATCHED) { vv = vertex_neighbor[vv]; }

                        edges[curr_m].v   = vv;
                        edges[curr_m++].w = ww;

                        // if the edge is present then add the weight, else expand it
                        for (size_t j = neighborhoods[u]; j < curr_m - 1; ++j) {
                            if (edges[j].v == vv) {
                                edges[j].w += ww;
                                curr_m -= 1;
                                break;
                            }
                        }
                    }
                }

                neighborhoods[u + 1] = curr_m;
            }
        }

        vertex_t get_n() const override { return n; }
        vertex_t get_m() const override { return m; }
        vertex_t get_true_m() const { return curr_m; }
        weight_t get_weight() const override { return vertex_weights; }
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
    };
}

#endif //HEIDELBERGPROCESSMAPPING_GRAPH_CSR_H
