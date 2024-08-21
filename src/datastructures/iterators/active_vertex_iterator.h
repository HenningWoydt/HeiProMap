#ifndef SERIALPROCESSMAPPING_ACTIVE_VERTEX_ITERATOR_H
#define SERIALPROCESSMAPPING_ACTIVE_VERTEX_ITERATOR_H

#include <vector>
#include <fstream>
#include <regex>
#include <numeric>

#include "../../utility/definitions.h"
#include "../../utility/utils.h"
#include "../../utility/macros.h"
#include "../graph.h"

namespace SPM {

    class ActiveVertexIterator {
    private:
        std::vector<u8> &states;
        std::vector<vertex_t> &active_vertices;
        size_t idx;

        // random
        std::mt19937 rng;
        bool do_shuffle;
        size_t shuffle_size = 1024 * 32;

    public:
        explicit ActiveVertexIterator(Graph &g, bool do_shuffle = false) : states(g.get_active_states()), active_vertices(g.get_active_vertices()), do_shuffle(do_shuffle) {
            idx = 0;

            rng.seed(0);
        }

        vertex_t get(){
            return active_vertices[idx];
        }

        void next(){
            idx += 1;
        }

        bool not_end(){
            if(do_shuffle && idx % (shuffle_size / 2) == 0){
                shuffle();
            }

            while(idx < active_vertices.size() && states[active_vertices[idx]] == 0){
                active_vertices[idx] = active_vertices.back();
                active_vertices.pop_back();
            }

            if(idx >= active_vertices.size()){
                return false;
            }
            return true;
        }

    private:
        void shuffle() {
            size_t start = idx;
            size_t end = std::min(start + shuffle_size, active_vertices.size());
            std::shuffle(active_vertices.begin() + start, active_vertices.begin() + end, rng);
        }
    };

}

#endif //SERIALPROCESSMAPPING_ACTIVE_VERTEX_ITERATOR_H
